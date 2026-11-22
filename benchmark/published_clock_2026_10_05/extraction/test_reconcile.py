import json
import tempfile
import unittest
from pathlib import Path

from .reconcile import assemble, compare, problems


TEXT = {'paper_id': 'MED:1', 'fulltext_status': 'ok', 'blocks': [
    {'text': 'Patients had higher Horvath age acceleration than controls (P = 0.003; n = 20 and 25).'},
    {'text': 'Table 2 GSE1 DunedinPACE 1.05 0.98 0.21'}]}


def record(**changes):
    row = {'tier': 'T1', 'gse': 'GSE1', 'modality': 'human DNAm', 'tissue': 'blood', 'clock_reported': 'Horvath',
           'clock_clockbase': 'Horvathv1', 'measure': 'acceleration_residual', 'case_label': 'patients',
           'control_label': 'controls', 'case_definition': 'patients', 'control_definition': 'controls',
           'n_case': 20, 'n_control': 25, 'direction': 'higher', 'significant': True, 'threshold': None,
           'p_value': 0.003, 'effect': None, 'test': None, 'covariates': [],
           'source': {'quote': 'Patients had higher Horvath age acceleration than controls (P = 0.003; n = 20 and 25).',
                      'location': 'Results'}}
    row.update(changes)
    return row


class ReconcileTest(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        text = self.root / 'text.json'
        text.write_text(json.dumps(TEXT))
        for name in ('inputs', 'first', 'second'):
            (self.root / name).mkdir()
        paper = {'paper_id': 'MED:1', 'ids': {'id': '1', 'source': 'MED', 'pmid': '1', 'pmcid': None, 'doi': None},
                 'title': 'A paper', 'year': '2020', 'text_path': str(text)}
        (self.root / 'inputs/MED_1.json').write_text(json.dumps({'paper': paper, 'series': [{'gse': 'GSE1', 'modality': 'human DNAm'}]}))
        self.texts = [b['text'] for b in TEXT['blocks']]

    def write(self, label, records):
        (self.root / label / 'MED_1.json').write_text(json.dumps({'paper_id': 'MED:1', 'records': records}))

    def test_checks_quotes_and_statistics(self):
        self.assertEqual(problems(record(), self.texts), [])
        self.assertTrue(any('statistic 0.03' in p for p in problems(record(p_value=0.03), self.texts)))
        quoted = record(source={'quote': 'Patients had lower clocks', 'location': 'x'}, p_value=None, n_case=None,
                        n_control=None)
        self.assertTrue(any('absent' in p for p in problems(quoted, self.texts)))
        across = record(source={'quote': 'than controls (P = 0.003; n = 20 and 25). Table 2 GSE1', 'location': 'x'})
        self.assertEqual(problems(across, self.texts), [])
        joined = record(source={'quote': 'Patients had higher Horvath age acceleration than controls (P = 0.003; n = 20 '
                                         'and 25). Table 2 GSE1 DunedinPACE 1.05', 'location': 'x'})
        self.assertEqual(problems(joined, list(reversed(self.texts))), [])
        stitched = record(source={'quote': 'Patients had higher Horvath age acceleration GSE1 DunedinPACE 1.05 0.98 0.21',
                                  'location': 'x'}, p_value=None, n_case=None, n_control=None)
        self.assertEqual(problems(stitched, self.texts), [])
        invented = record(source={'quote': 'Patients had a much lower GrimAge than controls', 'location': 'x'},
                          p_value=None, n_case=None, n_control=None)
        self.assertTrue(problems(invented, self.texts))
        curly = record(source={'quote': 'Patients had higher Horvath age acceleration than controls (P = 0.003; '
                                        'n = 20 and 25)', 'location': 'x'})
        self.assertEqual(problems(curly, [t.replace('-', '\u2013') for t in self.texts]), [])

    def test_agreement_requires_identical_results_and_review_gets_the_rest(self):
        self.write('first', [record(), record(clock_clockbase='DunedinPACE', p_value=None, n_case=None,
                                              n_control=None, effect={'value': 0.21, 'unit': None, 'type': None},
                                              source={'quote': 'Table 2 GSE1 DunedinPACE 1.05 0.98 0.21',
                                                      'location': 'Table 2'})])
        self.write('second', [record(case_label='cases'), record(direction='lower')])
        summary = compare(self.root, 'first', 'second')
        self.assertEqual(summary['agreed_records'], 1)
        review = json.loads((self.root / 'review/MED_1.json').read_text())
        self.assertEqual(len(review['unmatched']['first']), 1)
        self.assertEqual(review['unmatched']['second'][0]['record']['direction'], 'lower')

    def test_collect_keeps_only_reviewed_records_that_pass_checks(self):
        from .reconcile import collect
        self.write('first', [record()])
        self.write('second', [record(direction='lower')])
        compare(self.root, 'first', 'second')
        (self.root / 'reviewer').mkdir()
        good = dict(record(), extractors=['first', 'codex-reviewer'], agreement='single', notes='')
        bad = dict(record(p_value=0.5), extractors=['second'], agreement='single', notes='')
        (self.root / 'reviewer/MED_1.json').write_text(json.dumps({'records': [good, bad], 'rejected': []}))
        summary = collect(self.root, 'reviewer')
        self.assertEqual((summary['records_single'], summary['records_failing_checks']), (1, 1))
        kept = json.loads((self.root / 'reviewed/MED_1.json').read_text())['records']
        self.assertEqual([r['p_value'] for r in kept], [0.003])

    def test_clock_names_map_to_clockbase(self):
        from .reconcile import clockbase_clock
        self.assertEqual(clockbase_clock('Horvathv2 (skin and blood)'), 'Horvathv2')
        self.assertIsNone(clockbase_clock('DeepMAge'))
        self.assertIsNone(clockbase_clock(None))

    def test_assemble_assigns_ids_and_keeps_the_base(self):
        self.write('first', [record()])
        self.write('second', [record()])
        compare(self.root, 'first', 'second')
        base = self.root / 'base.jsonl'
        base.write_text(json.dumps({'id': 'old-0'}) + '\n')
        result = assemble(self.root, base, self.root / 'out.jsonl')
        rows = [json.loads(line) for line in (self.root / 'out.jsonl').read_text().splitlines()]
        self.assertEqual(result['added'], {'agreed': 1})
        self.assertEqual([r['id'] for r in rows], ['old-0', 'MED-1-GSE1-Horvathv1-0'])
        self.assertEqual(rows[1]['extraction']['agreement'], 'agreed')
        self.assertEqual(rows[1]['source']['pmid'], '1')
        base.write_text(json.dumps({'id': 'old-0', 'gse': 'GSE1', 'source': {'paper_id': 'MED:1'}}) + '\n')
        result = assemble(self.root, base, self.root / 'out2.jsonl')
        self.assertEqual((result['added'], result['skipped']), ({}, {'paper_series_already_in_base': 1}))


if __name__ == '__main__':
    unittest.main()
