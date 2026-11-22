"""Tests of the rerun scorer's matching and orientation (PROTOCOL.md of groupbench, addendum item 16)."""
import unittest

from .score import ABSOLUTE, match, raw_result


def run_of(case, control, clock='Horvathv1', direction=1.0, p=0.01, effect=0.8):
    return {'arms': {('A vs B', 'A'): {'case': list(case), 'control': list(control)}},
            'rows': [{'comparison': 'A vs B', 'condition': 'A', 'clock': clock, 'p': p, 'effect': effect,
                      'direction': direction, 'test': 't-test'}]}


def record_of(case, control, clock='Horvathv1'):
    return {'clock_clockbase': clock, 'groups': {'case_gsms': list(case), 'control_gsms': list(control)}}


class MatchTest(unittest.TestCase):
    case, control = ['G1', 'G2', 'G3', 'G4'], ['G5', 'G6', 'G7', 'G8']

    def test_same_orientation(self):
        row, swapped, overlap = match(record_of(self.case, self.control), run_of(self.case, self.control))
        self.assertFalse(swapped)
        self.assertEqual(overlap, 1.0)
        self.assertEqual(raw_result(row, swapped, overlap)['direction'], 'higher')

    def test_swapped_arms_reverse_the_reported_direction(self):
        row, swapped, overlap = match(record_of(self.case, self.control), run_of(self.control, self.case))
        self.assertTrue(swapped)
        self.assertEqual(raw_result(row, swapped, overlap)['direction'], 'lower')

    def test_partial_overlap_at_the_threshold(self):
        found = match(record_of(self.case, self.control), run_of(['G1', 'G2', 'G3', 'X1', 'X2'], self.control))
        self.assertIsNotNone(found)                      # Jaccard 3 / 6 = 0.5
        self.assertIsNone(match(record_of(self.case, self.control), run_of(['G1', 'G2', 'X1', 'X2', 'X3'], self.control)))

    def test_other_clock_does_not_match(self):
        self.assertIsNone(match(record_of(self.case, self.control, 'PhenoAge'), run_of(self.case, self.control)))

    def test_direction_falls_back_to_the_effect_sign(self):
        row, swapped, overlap = match(record_of(self.case, self.control), run_of(self.case, self.control, direction=None, effect=-0.4))
        self.assertEqual(raw_result(row, swapped, overlap)['direction'], 'lower')
        row['effect'] = None
        self.assertFalse(raw_result(row, swapped, overlap)['covered'])

    def test_absolute_paths_in_commands(self):
        self.assertEqual(ABSOLUTE.findall('python /Entropy/kying/x/run.py data/in.csv'), ['/Entropy/kying/x/run.py'])
        self.assertEqual(ABSOLUTE.findall('open(base + "/clock_output/a.csv")'), ['/clock_output/a.csv'])  # filtered later
        self.assertEqual(ABSOLUTE.findall('pd.read_csv("data/clock_output/a.csv")'), [])


if __name__ == '__main__':
    unittest.main()


class ClaudeAuditTest(unittest.TestCase):
    """PROTOCOL.md item 22: the Claude Code audit sees paths in shell commands and in file-tool arguments."""

    def _run(self, blocks, result=None):
        import json
        import tempfile
        from pathlib import Path
        tmp = Path(tempfile.mkdtemp(dir='/Entropy/kying/tmp'))
        run = tmp / 'claude-claude-opus-5-5' / 'GSE1'
        run.mkdir(parents=True)
        events = [{'type': 'assistant', 'message': {'content': blocks}}]
        if result:
            events.append(result)
        (run / 'events.jsonl').write_text('\n'.join(json.dumps(e) for e in events))
        (run / 'status').write_text('0')
        return run

    def test_paths_inside_the_run_are_allowed(self):
        from .score import audit_claude
        run = self._run([{'type': 'tool_use', 'name': 'Read', 'input': {'file_path': 'results/x.csv'}}])
        run2 = self._run([{'type': 'tool_use', 'name': 'Bash', 'input': {'command': 'ls data'}}])
        self.assertEqual(audit_claude(run, '/usr/bin/python3')['outside_paths'], [])
        self.assertEqual(audit_claude(run2, '/usr/bin/python3')['commands'], 1)

    def test_file_tool_outside_the_run_is_caught(self):
        from .score import audit_claude
        run = self._run([{'type': 'tool_use', 'name': 'Read',
                          'input': {'file_path': '/Entropy/kying/agent-worktrees/x/README.md'}}])
        self.assertEqual(audit_claude(run, '/usr/bin/python3')['outside_paths'],
                         ['/Entropy/kying/agent-worktrees/x/README.md'])

    def test_command_outside_the_run_is_caught(self):
        from .score import audit_claude
        run = self._run([{'type': 'tool_use', 'name': 'Bash', 'input': {'command': 'cat /home/kying/notes.txt'}}])
        self.assertEqual(audit_claude(run, '/usr/bin/python3')['outside_paths'], ['/home/kying/notes.txt'])

    def test_web_and_subagent_tools_are_counted(self):
        from .score import audit_claude
        run = self._run([{'type': 'tool_use', 'name': 'WebSearch', 'input': {'query': 'x'}},
                         {'type': 'tool_use', 'name': 'Agent', 'input': {'prompt': 'x'}}])
        out = audit_claude(run, '/usr/bin/python3')
        self.assertEqual((out['web_searches'], out['spawned_agents']), (1, 1))

    def test_usage_from_the_result_event(self):
        from .score import audit_claude
        run = self._run([], {'type': 'result', 'total_cost_usd': 1.5,
                             'usage': {'input_tokens': 10, 'cache_read_input_tokens': 90,
                                       'cache_creation_input_tokens': 5, 'output_tokens': 7}})
        usage = audit_claude(run, '/usr/bin/python3')['usage']
        self.assertEqual((usage['input_tokens'], usage['cached_input_tokens'], usage['output_tokens']), (105, 90, 7))


class BiomniAuditTest(unittest.TestCase):
    """PROTOCOL.md item 24: the Biomni audit reads the executed blocks of A1's messages."""

    def _run(self, contents, usage=()):
        import json
        import tempfile
        from pathlib import Path
        tmp = Path(tempfile.mkdtemp(dir='/Entropy/kying/tmp'))
        run = tmp / 'biomni-gpt-6-astra' / 'GSE1'
        run.mkdir(parents=True)
        events = [{'type': 'message', 'role': 'ai', 'content': c} for c in contents] + list(usage)
        (run / 'events.jsonl').write_text('\n'.join(json.dumps(e) for e in events))
        (run / 'status').write_text('0')
        return run

    def test_executed_block_rules(self):
        from .score import executed_block
        self.assertEqual(executed_block('plan <execute>print(1)</execute>'), 'print(1)')
        self.assertEqual(executed_block('<execute>print(2)'), 'print(2)')                 # A1 closes the tag
        self.assertIsNone(executed_block('<execute>x</execute><solution>done</solution>'))  # the solution ends the run
        self.assertIsNone(executed_block('<observation>/home/x</observation>'))
        self.assertIsNone(executed_block('no tags'))

    def test_paths_and_usage(self):
        from .score import audit_biomni
        run = self._run(['<execute>import pandas as pd\npd.read_csv("data/x.csv")</execute>',
                         '<execute>open("/home/kying/notes.txt")</execute>',
                         '<observation>Error: /home/kying/other.txt</observation>'],
                        [{'type': 'usage', 'input_tokens': 10, 'cached_input_tokens': 4, 'output_tokens': 2},
                         {'type': 'usage', 'input_tokens': 5, 'cached_input_tokens': 0, 'output_tokens': 1}])
        out = audit_biomni(run, '/usr/bin/python3')
        self.assertEqual(out['commands'], 2)
        self.assertEqual(out['outside_paths'], ['/home/kying/notes.txt'])
        self.assertEqual(out['usage'], {'input_tokens': 15, 'cached_input_tokens': 4, 'output_tokens': 3})

    def test_network_and_tool_modules_count_as_web_access(self):
        from .score import audit_biomni
        run = self._run(['<execute>from biomni.tool.database import query_geo\nquery_geo("x")</execute>',
                         '<execute>#!BASH\ncurl -s example.org</execute>',
                         '<execute>from biomni.tool.support_tools import run_python_repl</execute>'])
        out = audit_biomni(run, '/usr/bin/python3')
        self.assertEqual((out['web_searches'], out['network_attempts'], out['tool_module_calls']), (2, 1, 1))


class BiomniStopTest(unittest.TestCase):
    """PROTOCOL.md item 25: replies are cut at A1's first stop sequence."""

    def test_truncate_at_stop(self):
        from .biomni_runner import truncate_at_stop
        self.assertEqual(truncate_at_stop('plan <execute>x = 1</execute> <solution>done</solution>'), 'plan <execute>x = 1')
        self.assertEqual(truncate_at_stop('<solution>done</solution><execute>x</execute>'), '<solution>done')
        self.assertEqual(truncate_at_stop('<execute>x = 1'), '<execute>x = 1')
