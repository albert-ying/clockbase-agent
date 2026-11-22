"""Check cached source identities, extraction boundaries and atlas linkage."""
import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import openpyxl

from benchmark.published_clock_2026_10_05.curated import assemble_2026_10_05 as assembly
from benchmark.published_clock_2026_10_05.curated.retrieve_2026_10_05 import ROOT, retrieve


class CuratedSourcesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cache = ROOT / 'downloads'
        cls.records = (assembly.classic_records(cls.cache) + assembly.phenoage_records(cls.cache)
                       + assembly.computage_records(cls.cache) + assembly.signature_records(cls.cache)
                       + assembly.biolearn_records(cls.cache))
        assembly.bind_atlas(cls.records, assembly.ATLAS)
        assembly.bind_evidence(cls.records, cls.cache)

    def test_every_quote_has_verified_download_bytes(self):
        for row in self.records:
            with self.subTest(record=row['record_id']):
                self.assertTrue(row['quote_evidence'])
                for quote in row['quote_evidence']:
                    payload = (self.cache / quote['file']).read_bytes()
                    self.assertEqual(hashlib.sha256(payload).hexdigest(), quote['sha256'])
                    self.assertTrue(quote['url'].startswith('https://'))
                    self.assertTrue(quote['retrieved_at'])
                    if quote['file'].endswith('.txt') or quote['file'].endswith('-soft'):
                        self.assertIn(assembly.normalized(quote['text']), assembly.normalized(payload.decode()))

    def test_expected_directions_do_not_acquire_reported_statistics(self):
        rows = [r for r in self.records if r['source_id'] == 'computage']
        self.assertTrue(rows)
        self.assertTrue(any(not r['control_available_within_series'] for r in rows))
        for row in rows:
            self.assertEqual(row['ground_truth_type'], 'literature_expected_direction')
            self.assertEqual(row['expected_direction'], 'case higher')
            self.assertIsNone(row['reported_direction'])
            self.assertEqual(row['reported_statistics'], [])
            self.assertIsNone(row['clock_clockbase'])
            self.assertIsNone(row['n_case'])
            self.assertIsNone(row['n_control'])
        with (self.cache / 'hf-meta-pinned.tsv').open() as handle:
            metadata = list(csv.DictReader(handle, delimiter='\t'))
        for row in rows:
            expected = sum(m['DatasetID'] == row['gse_id'] and m['Tissue'] == row['tissue']
                           and m['Condition'] == row['condition_code'] for m in metadata)
            self.assertEqual(expected, row['metadata_n_case'])

    def test_signature_scores_are_copied_and_not_mapped_to_clocks(self):
        workbook = openpyxl.load_workbook(self.cache / 'tysh2023-mmc6.xlsx', read_only=True, data_only=True)
        sheet = workbook['A (Public datasets)']
        row = next(r for r in self.records if r['source_id'] == 'tysh2023'
                   and r['gse_id'] == 'GSE15891' and r['signature_source'] == 'Interventions')
        stats = {s['name']: s for s in row['reported_statistics']}
        self.assertEqual(stats['normalized_enrichment_score']['value'], sheet['D3'].value)
        self.assertEqual(stats['p_value']['value'], sheet['E3'].value)
        self.assertEqual(stats['adjusted_p_value']['value'], sheet['F3'].value)
        workbook.close()
        for row in self.records:
            if row['ground_truth_type'].startswith('reported_signature'):
                self.assertIsNone(row['clock_source'])
                self.assertIsNone(row['clock_clockbase'])
                self.assertIsNone(row['reported_direction'])
                self.assertFalse(row['direct_case_control_statistic'])

    def test_hiv_consensus_and_nonsignificant_horvath_remain_distinct(self):
        rows = [r for r in self.records if r['source_id'] == 'hivgross']
        for row in rows:
            if row['clock_source'].startswith('Consensus'):
                self.assertIsNone(row['clock_clockbase'])
        row = next(r for r in rows if r['clock_source'] == 'Horvath model')
        self.assertEqual(row['reported_direction'], 'no significant difference')
        self.assertEqual(row['clock_clockbase'], 'Horvathv1')
        pvalue = next(s for s in row['reported_statistics'] if s['name'] == 'p_value')
        self.assertEqual(pvalue, assembly.statistic('p_value', 0.05, 'probability', '>'))

    def test_pooled_parkinson_results_cannot_become_per_series_effects(self):
        rows = [r for r in self.records if r['source_id'] == 'parkinson']
        self.assertEqual(len({r['pooling_id'] for r in rows}), 2)
        for row in rows:
            self.assertEqual(row['ground_truth_type'], 'reported_pooled_case_control_statistic')
            self.assertFalse(row['direct_case_control_statistic'])
            self.assertEqual(row['pooled_gse_ids'], ['GSE72774', 'GSE72776'])

    def test_atlas_matches_require_both_clock_and_modality(self):
        for row in self.records:
            matched = [r for r in row['atlas']['rows'] if r['clock'] == row['clock_clockbase']
                       and r['modality'] == row['modality']]
            self.assertEqual(row['atlas']['series_clock_match'], bool(matched))
            self.assertIsNone(row['atlas']['comparison_match'])
        absent = next(r for r in self.records if r['source_id'] == 'huntington')
        self.assertFalse(absent['atlas']['series_present'])

    def test_source_values_and_unreported_directions_are_preserved(self):
        rows = [r for r in self.records if r['source_id'] == 'phenoage']
        pd = next(r for r in rows if r['gse_id'] == 'GSE72774')
        self.assertIsNone(pd['reported_direction'])
        self.assertEqual((pd['n_case'], pd['n_control']), (289, 219))
        self.assertEqual(pd['statistic_fields']['p_value']['value'], 0.028)
        liver = next(r for r in rows if r['gse_id'] == 'GSE48325' and r['statistic_fields']['beta'] is not None
                     and r['clock_clockbase'] == 'PhenoAge')
        self.assertEqual(liver['n_total'], 85)
        self.assertEqual(liver['statistic_fields']['beta']['value'], 0.13)
        self.assertEqual(liver['ground_truth_type'], 'reported_continuous_association')

    def test_get_cache_refuses_identity_changes_without_network(self):
        validation = ROOT / 'validation'
        validation.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=validation) as directory:
            root = Path(directory)
            cache = root / 'downloads'
            cache.mkdir()
            payload = (self.cache / 'hf-readme.md').read_bytes()
            (cache / 'source.md').write_bytes(payload)
            (cache / 'source.md.receipt.json').write_text(json.dumps({
                'url': 'https://example.invalid/source', 'sha256': hashlib.sha256(payload).hexdigest()}))
            self.assertEqual(retrieve('https://example.invalid/source', 'source.md', root), cache / 'source.md')
            with self.assertRaises(ValueError):
                retrieve('https://example.invalid/changed', 'source.md', root)
            with self.assertRaises(ValueError):
                retrieve('https://example.invalid/source', '../outside.md', root)


if __name__ == '__main__':
    unittest.main()
