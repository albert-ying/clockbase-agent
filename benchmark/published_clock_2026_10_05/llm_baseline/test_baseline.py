"""S4 checks with synthetic text fixtures and a fake transport; no model calls."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from .answers import parse_answer
from .cli import main
from .prompts import (CLOCKS, StudyDescriptions, build_prompt, check_blinding,
                      clock_description)
from .runner import MAX_TOKENS, request_hash, run
from .scorer import read_jsonl


def synthetic_record(**changes):
    record = {'id': 'synthetic-GSE123-horvath-1', 'tier': 'T1', 'gse': 'GSE123',
              'modality': 'human DNAm', 'clock_reported': 'Horvath DNAm age acceleration',
              'clock_clockbase': 'Horvath', 'measure': 'acceleration_residual',
              'case_label': 'condition', 'control_label': 'control',
              'case_definition': 'Participants with condition A.',
              'control_definition': 'Participants without condition A.', 'n_case': 12, 'n_control': 14,
              'direction': 'higher', 'significant': True, 'threshold': 'P < 0.01', 'p_value': 0.00731,
              'effect': {'value': 2.314, 'unit': 'years', 'type': 'mean difference'},
              'source': {'title': 'Synthetic paper result title', 'citation': 'Synthetic author 2024',
                         'quote': 'The clock result was higher by 2.314 years (P = 0.00731).'}}
    record.update(changes)
    return record


def synthetic_study():
    return {'title': 'DNA methylation in condition A',
            'summary': 'Blood samples were collected from participants with and without condition A.',
            'overall_design': 'The study compared two groups.', 'n_samples': 26, 'sources': []}


def valid_answer(**changes):
    answer = {'prediction': 'lower', 'significance': False, 'confidence': 0.7,
              'recognizes_publication': False, 'publication': None}
    answer.update(changes)
    return json.dumps(answer)


class FakeResponse:
    def __init__(self, answer):
        self.payload = {'model': 'fake-returned-model',
                        'usage': {'prompt_tokens': 420, 'completion_tokens': 60, 'total_tokens': 480},
                        'choices': [{'message': {'content': answer}}]}

    def read(self):
        return json.dumps(self.payload).encode()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class BaselineTest(unittest.TestCase):
    def setUp(self):
        # All temporary test artifacts stay within the package this worker owns.
        workspace = tempfile.TemporaryDirectory(dir=Path(__file__).parent)
        self.addCleanup(workspace.cleanup)
        self.root = Path(workspace.name)
        self.geo = self.root / 'geo'
        self.geo.mkdir()
        self.ground_truth = self.root / 'ground_truth.jsonl'
        self.write_records([synthetic_record()])
        self.write_index()
        env = patch.dict(os.environ, {'AUDIT_CACHE_ONLY': '0', 'TOKENSWITCH_API_KEY': 'synthetic-test-key',
                                      'TOKENSWITCH_BASE_URL': 'https://example.invalid',
                                      'BENCH_S4_MAX_REQUESTS': '3000'}, clear=False)
        env.start()
        self.addCleanup(env.stop)
        sleeper = patch('audit_system.label_adjudication_2026_09_27.time.sleep')
        sleeper.start()
        self.addCleanup(sleeper.stop)

    def write_records(self, records):
        self.ground_truth.write_text(''.join(json.dumps(row) + '\n' for row in records))

    def write_index(self, filename='gse_me_hs.csv', rows=None):
        if rows is None:
            study = synthetic_study()
            rows = [{'accession': 'GSE123', 'title': study['title'], 'description': study['summary'], 'samples': 26}]
        with (self.geo / filename).open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['accession', 'title', 'description', 'samples'])
            writer.writeheader()
            writer.writerows(rows)

    def execute(self, **kwargs):
        return run(self.ground_truth, self.root / 'out', descriptions=StudyDescriptions(self.geo),
                   models=kwargs.pop('models', ['fake-model']), **kwargs)

    def test_prompt_excludes_source_and_result_fields(self):
        record = synthetic_record()
        prepared = build_prompt(record, synthetic_study())
        self.assertEqual(check_blinding(prepared['stimulus'], record), [])
        for forbidden in (record['source']['quote'], record['source']['title'], record['source']['citation'],
                          record['id'], str(record['p_value']), str(record['effect']['value']), record['threshold']):
            self.assertNotIn(forbidden, prepared['prompt'])
        self.assertNotIn('accelerat', prepared['stimulus'].lower())
        self.assertIn('Residual of clock output', prepared['prompt'])
        self.assertIn('Series samples: 26', prepared['prompt'])
        self.assertEqual(prepared['flags'], [])

    def test_prompt_is_invariant_to_published_direction_and_significance(self):
        first = build_prompt(synthetic_record(), synthetic_study())
        second = build_prompt(synthetic_record(direction='lower', significant=False), synthetic_study())
        self.assertEqual(first['prompt'], second['prompt'])

    def test_geo_result_sentences_and_paper_title_are_removed(self):
        record, study = synthetic_record(), synthetic_study()
        study['title'] = record['source']['title']
        study['summary'] += (' Age was accelerated in cases. ' + record['source']['quote'] +
                             ' Epigenetic age did not differ between groups.')
        prepared = build_prompt(record, study)
        self.assertIn('geo_title_redacted', prepared['flags'])
        self.assertIn('geo_summary_redacted', prepared['flags'])
        self.assertEqual(check_blinding(prepared['stimulus'], record), [])
        self.assertNotIn('did not differ', prepared['stimulus'])

    def test_revealing_group_definitions_are_kept_and_flagged(self):
        record = synthetic_record(case_definition='Participants selected for accelerated epigenetic aging.')
        prepared = build_prompt(record, synthetic_study())
        self.assertIn(record['case_definition'], prepared['stimulus'])
        self.assertIn('case_control_answer_leakage', prepared['flags'])
        self.assertEqual(prepared['blinding_issues'], ['result_wording'])

    def test_exact_source_quote_in_group_definition_stops_request(self):
        record = synthetic_record()
        record['case_definition'] = record['source']['quote']
        self.write_records([record])
        with patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen') as transport:
            summary = self.execute()
        transport.assert_not_called()
        self.assertEqual(summary['statuses'], {'PROMPT_FAILED': 1})
        result = next(read_jsonl(self.root / 'out' / 'predictions.jsonl'))
        self.assertIn('case_control_answer_leakage', result['prompt_flags'])

    def test_numeric_result_collision_omits_sample_count(self):
        record = synthetic_record(effect={'value': 26})
        prepared = build_prompt(record, synthetic_study())
        self.assertIn('Series samples: [withheld]', prepared['stimulus'])
        self.assertIn('sample_count_result_collision', prepared['flags'])
        self.assertEqual(check_blinding(prepared['stimulus'], record), [])

    def test_clocks_have_neutral_descriptions(self):
        for name in CLOCKS:
            with self.subTest(name=name):
                record = synthetic_record(clock_clockbase=name, clock_reported=name)
                chosen, description = clock_description(record)
                self.assertEqual(chosen, name)
                self.assertTrue(description.startswith(('Estimates ', 'Combines ')))
        self.assertEqual(clock_description(synthetic_record(clock_clockbase='Horvath2'))[0], 'Skin and blood')
        self.assertEqual(clock_description(synthetic_record(clock_clockbase='Horvathv1',
                         clock_reported='Horvath 2013 (353-CpG epigenetic clock)'))[0], 'Horvath pan-tissue age')
        self.assertEqual(clock_description(synthetic_record(clock_clockbase='GrimAge',
                         clock_reported='GrimAge 2019'))[0], 'GrimAge')
        with self.assertRaisesRegex(ValueError, 'No neutral clock description'):
            clock_description(synthetic_record(clock_clockbase=None, clock_reported='Unregistered clock'))

    def test_modality_indices_and_provenance(self):
        for modality, filename in [('human DNAm', 'gse_me_hs.csv'), ('human RNA-seq', 'gse_rna_hs.csv'),
                                   ('mouse RNA-seq', 'gse_rna_mm.csv')]:
            self.write_index(filename)
            study = StudyDescriptions(self.geo).get(synthetic_record(modality=modality))
            self.assertEqual(study['n_samples'], 26)
            source = study['sources'][0]
            self.assertEqual(source['sha256'], hashlib.sha256((self.geo / filename).read_bytes()).hexdigest())
            self.assertIsInstance(source['mtime_ns'], int)

    def test_series_header_fallback_stops_before_matrix_data(self):
        (self.geo / 'gse_me_hs.csv').unlink()
        folder = self.geo / 'GEO_hs_dnam' / 'series_matrix'
        folder.mkdir(parents=True)
        with gzip.open(folder / 'GSE123-GPL1_series_matrix.txt.gz', 'wt') as handle:
            handle.write('!Series_title\t"Synthetic header title"\n'
                         '!Series_summary\t"A description."\n!Series_summary\t"Second paragraph."\n'
                         '!Series_overall_design\t"Two groups."\n'
                         '!Sample_geo_accession\t"GSM1"\t"GSM2"\n'
                         '!series_matrix_table_begin\n!Series_summary\t"do not read matrix data"\n')
        study = StudyDescriptions(self.geo).get(synthetic_record())
        self.assertEqual(study['title'], 'Synthetic header title')
        self.assertEqual(study['summary'], 'A description. Second paragraph.')
        self.assertEqual(study['overall_design'], 'Two groups.')
        self.assertEqual(study['n_samples'], 2)

    def test_missing_geo_source_stops_without_model_request(self):
        (self.geo / 'gse_me_hs.csv').unlink()
        with patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen') as transport:
            with self.assertRaisesRegex(RuntimeError, 'Cannot proceed - required data file not found:'):
                self.execute()
        transport.assert_not_called()

    def test_json_fenced_json_and_aliases(self):
        for text in (valid_answer(), '```json\n' + valid_answer() + '\n```',
                     json.dumps({'direction': 'not different', 'significant': 'no', 'confidence': '0.6',
                                 'recognizes_publication': 'yes', 'publication_reference': 'Synthetic author 2024'})):
            with self.subTest(text=text):
                self.assertIsNotNone(parse_answer(text))
        self.assertEqual(parse_answer(text)['prediction'], 'none')

    def test_malformed_ambiguous_and_out_of_range_answers_are_rejected(self):
        invalid = ['lower, not significant', '{"prediction": "lower",}', '[]', 'null',
                   'Here is my answer: ' + valid_answer(), valid_answer(confidence=1.1),
                   valid_answer(confidence=float('nan')), valid_answer(confidence=True),
                   valid_answer(significance='maybe'), valid_answer(prediction='increased'),
                   valid_answer(recognizes_publication=True), valid_answer(publication='Synthetic paper'),
                   valid_answer(prediction='none', significance=True),
                   valid_answer().replace('"prediction": "lower"', '"prediction": "higher", "prediction": "lower"'),
                   '{"prediction":"lower","confidence":0.8}']
        for text in invalid:
            with self.subTest(text=text):
                self.assertIsNone(parse_answer(text))

    def test_ledger_caching_hashing_costs_and_engine_export(self):
        with patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen',
                   return_value=FakeResponse(valid_answer())) as transport:
            first = self.execute()
            second = self.execute()
        self.assertEqual(transport.call_count, 1)
        body = json.loads(transport.call_args.args[0].data)
        self.assertEqual(body['temperature'], 0)
        self.assertEqual(body['max_tokens'], MAX_TOKENS)
        digest = request_hash(body['model'], body['messages'][0]['content'])
        ledger = list(read_jsonl(self.root / 'out' / 'ledger.jsonl'))
        self.assertEqual(ledger[0]['request_sha256'], digest)
        self.assertEqual(ledger[0]['model_requested'], 'fake-model')
        self.assertEqual(ledger[0]['model_returned'], 'fake-returned-model')
        self.assertEqual(ledger[0]['response_text'], valid_answer())
        self.assertEqual(first['requests_sent_this_run'], 1)
        self.assertEqual(second['requests_sent_this_run'], 0)
        self.assertEqual(second['cached_responses'], 1)
        self.assertEqual(second['usage_by_model']['fake-model']['total_tokens'], 480)
        exported = list(read_jsonl(self.root / 'out' / 's4_predictions.jsonl'))
        self.assertEqual(len(exported), 1)
        self.assertEqual(exported[0]['direction'], 'lower')
        self.assertIs(exported[0]['significant'], False)
        self.assertEqual(exported[0]['request_hash'], digest)
        # Check the engine's reader directly; no GEO/sample analysis or model calls run here.
        from benchmark.published_clock_2026_10_05.engine.cli import load_s4
        engine = load_s4(self.root / 'out' / 's4_predictions.jsonl')
        self.assertTrue(engine[synthetic_record()['id'], 'fake-model']['covered'])
        self.assertEqual(engine[synthetic_record()['id'], 'fake-model']['confidence'], 0.7)

    def test_malformed_response_is_cached_without_repurchasing(self):
        with patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen',
                   return_value=FakeResponse('not JSON')) as transport:
            self.execute()
            summary = self.execute()
        self.assertEqual(transport.call_count, 1)
        self.assertEqual(summary['statuses'], {'UNPARSEABLE': 1})
        self.assertEqual(summary['cached_responses'], 1)
        result = next(read_jsonl(self.root / 'out' / 's4_predictions.jsonl'))
        self.assertIsNone(result['direction'])
        self.assertIn('request_sha256=', result['raw_response_reference'])
        from benchmark.published_clock_2026_10_05.engine.cli import load_s4
        engine = load_s4(self.root / 'out' / 's4_predictions.jsonl')
        self.assertFalse(engine[synthetic_record()['id'], 'fake-model']['covered'])

    def test_budget_stop_including_restart_and_cached_response(self):
        with patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen',
                   return_value=FakeResponse(valid_answer())) as transport:
            first = self.execute(models=['fake-a', 'fake-b'], max_requests=1)
            second = self.execute(models=['fake-a', 'fake-b'], max_requests=1)
        self.assertEqual(transport.call_count, 1)
        self.assertEqual(first['statuses'], {'OK': 1, 'BUDGET_STOP': 1})
        self.assertEqual(second['cached_responses'], 1)
        stopped = list(read_jsonl(self.root / 'out' / 'predictions.jsonl'))[-1]
        self.assertEqual(stopped['status'], 'BUDGET_STOP')
        self.assertIsNotNone(stopped['request_hash'])
        self.assertIsNone(stopped['raw_response_reference'])

    def test_cache_only_miss_and_hit_never_send(self):
        with patch.dict(os.environ, {'AUDIT_CACHE_ONLY': '1'}), \
                patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen') as transport:
            missing = self.execute()
        transport.assert_not_called()
        self.assertEqual(missing['statuses'], {'NOT_CACHED': 1})
        self.assertFalse((self.root / 'out' / 'ledger.jsonl').exists())
        with patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen',
                   return_value=FakeResponse(valid_answer())):
            self.execute()
        with patch.dict(os.environ, {'AUDIT_CACHE_ONLY': '1'}), \
                patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen') as transport:
            replay = self.execute()
        transport.assert_not_called()
        self.assertEqual(replay['statuses'], {'OK': 1})
        self.assertEqual(replay['cached_responses'], 1)

    def test_environment_models_budget_limit_and_zero_limit(self):
        with patch.dict(os.environ, {'BENCH_S4_MODELS': 'fake-a,fake-b', 'BENCH_S4_MAX_REQUESTS': '0'}), \
                patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen') as transport:
            summary = run(self.ground_truth, self.root / 'out', descriptions=StudyDescriptions(self.geo))
            empty = self.execute(limit=0)
        transport.assert_not_called()
        self.assertEqual(summary['statuses'], {'BUDGET_STOP': 2})
        self.assertEqual(empty['records'], 0)

    def test_transport_failures_are_ledgered_and_count_toward_budget(self):
        with patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen',
                   side_effect=OSError('synthetic transport failure')) as transport:
            summary = self.execute(models=['fake-a', 'fake-b'], max_requests=1)
        self.assertEqual(transport.call_count, 3)
        self.assertEqual(summary['statuses'], {'TRANSPORT_ERROR': 1, 'BUDGET_STOP': 1})
        entry = next(read_jsonl(self.root / 'out' / 'ledger.jsonl'))
        self.assertEqual(entry['attempts'], 3)
        self.assertIn('synthetic transport failure', entry['error'])
        self.assertEqual(summary['usage_by_model']['fake-a']['usage_missing'], 1)

    def test_malformed_transport_payload_is_saved_as_client_error(self):
        response = FakeResponse(valid_answer())
        response.payload['choices'] = []
        with patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen',
                   return_value=response):
            summary = self.execute()
        self.assertEqual(summary['statuses'], {'CLIENT_ERROR': 1})
        self.assertEqual(next(read_jsonl(self.root / 'out' / 'ledger.jsonl'))['status'], 'CLIENT_ERROR')

    def test_zero_budget_never_sends_and_direct_openai_names_are_recorded(self):
        with patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen') as transport:
            self.assertEqual(self.execute(max_requests=0)['statuses'], {'BUDGET_STOP': 1})
            self.assertEqual(self.execute(models=['openai:direct'], max_requests=0)['statuses'], {'BUDGET_STOP': 1})
        transport.assert_not_called()
        manifest = json.loads((self.root / 'out' / 'manifest.json').read_text())
        self.assertIn('no temperature sent', manifest['sampling']['openai:'])
        self.assertIn('temperature 0', manifest['sampling']['gateway'])

    def test_unowned_outputs_and_changed_ground_truth_are_not_overwritten(self):
        out = self.root / 'out'
        out.mkdir()
        prediction = out / 'predictions.jsonl'
        prediction.write_text('someone else owns this\n')
        with self.assertRaisesRegex(RuntimeError, 'Refusing to replace unowned'):
            self.execute(limit=0)
        self.assertEqual(prediction.read_text(), 'someone else owns this\n')
        prediction.unlink()
        self.execute(limit=0)
        self.write_records([synthetic_record(id='another-id')])
        with self.assertRaisesRegex(ValueError, 'output identity mismatch'):
            self.execute(limit=0)

    def test_cli_missing_source_reports_failure(self):
        self.assertEqual(main(['run', '--ground-truth', str(self.root / 'missing.jsonl'),
                               '--out', str(self.root / 'out')]), 1)

    def test_real_curated_prompt_without_model_calls(self):
        path = Path('/Entropy/kying/clockbase-benchmark-2026-10-05/curated/records.jsonl')
        if not path.exists():
            self.skipTest('Curated source file is not available')
        raw = next(read_jsonl(path))
        if 'id' in raw:
            record = raw
        else:
            # The extraction table uses earlier field names; its note states the residual measure.
            self.assertTrue(any('residual' in note for note in raw['notes']))
            record = dict(raw, id=raw['record_id'], gse=raw['gse_id'], clock_reported=raw['clock_source'],
                          measure='acceleration_residual', source={'quote': raw['quote']},
                          p_value=raw['statistic_fields']['p_value']['value'])
        with patch('audit_system.label_adjudication_2026_09_27.urllib.request.urlopen') as transport:
            prepared = build_prompt(record, StudyDescriptions().get(record))
        transport.assert_not_called()
        self.assertEqual(check_blinding(prepared['stimulus'], record), [])
        self.assertNotIn(record['source']['quote'], prepared['prompt'])
        self.assertNotIn(str(record['p_value']), prepared['prompt'])
        self.assertTrue(prepared['sources'])
        self.assertGreater(prepared['estimated_prompt_tokens'], 0)


if __name__ == '__main__':
    unittest.main()
