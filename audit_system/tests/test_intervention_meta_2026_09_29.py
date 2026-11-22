"""Unit tests for the intervention meta-analysis. All fixtures are synthetic."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from audit_system import intervention_meta_2026_09_29 as im


class InterventionMetaTests(unittest.TestCase):
    def test_hedges_g_correction(self):
        g, var = im.hedges_g(1.0, 10, 10)
        j = 1 - 3 / (4 * 18 - 1)
        self.assertAlmostEqual(g, j)
        self.assertAlmostEqual(var, j ** 2 * (0.2 + 1 / 40))

    def test_dersimonian_laird_homogeneous_equals_fixed(self):
        out = im.dersimonian_laird([0.5, 0.5, 0.5], [0.1, 0.2, 0.4])
        self.assertAlmostEqual(out['estimate'], 0.5)
        self.assertEqual(out['tau2'], 0.0)
        self.assertEqual(out['k'], 3)

    def test_dersimonian_laird_reference_values(self):
        # Hand-computed: y = [0, 1], v = [0.1, 0.1]; fixed = 0.5, Q = 5, c = 10, tau2 = 0.4
        out = im.dersimonian_laird([0.0, 1.0], [0.1, 0.1])
        self.assertAlmostEqual(out['tau2'], 0.4)
        self.assertAlmostEqual(out['estimate'], 0.5)
        self.assertAlmostEqual(out['se'], np.sqrt(1 / (2 / 0.5)))
        self.assertAlmostEqual(out['i2'], 0.8)

    def test_one_effect_per_study(self):
        rows = pd.DataFrame({'gse_id': ['A', 'A', 'B'], 'g': [1.0, 3.0, 2.0], 'var_g': [0.2, 0.4, 0.3],
                             'comparison_name': ['x', 'x', 'y'], 'condition_name': ['x', 'x', 'y'],
                             'v3_n_case': [5, 5, 4], 'v3_n_control': [5, 5, 4]})
        studies = im.study_effects(rows).set_index('gse_id')
        self.assertEqual(len(studies), 2)
        self.assertAlmostEqual(studies.loc['A', 'g'], 2.0)
        self.assertAlmostEqual(studies.loc['A', 'var'], 0.3)

    def test_one_effect_per_project(self):
        rows = pd.DataFrame({'gse_id': ['GSE1', 'GSE2', 'GSE3'], 'project': ['GSE1', 'GSE1', 'GSE3'],
                             'g': [1.0, 3.0, 2.0], 'var_g': [0.2, 0.4, 0.3], 'comparison_name': ['x', 'x', 'y'],
                             'condition_name': ['x', 'x', 'y'], 'v3_n_case': [5, 5, 4], 'v3_n_control': [5, 5, 4]})
        units = im.study_effects(rows, 'project').set_index('project')
        self.assertEqual(len(units), 2)
        self.assertAlmostEqual(units.loc['GSE1', 'g'], 2.0)
        self.assertEqual(units.loc['GSE1', 'series'], 'GSE1;GSE2')
        self.assertEqual(units.loc['GSE1', 'comparisons'], 2)

    def test_hartung_knapp_reference_values(self):
        # y = [0, 1], v = [0.1, 0.1]: tau2 = 0.4, w* = 2 each, estimate 0.5,
        # se_hk = sqrt((2 * 0.25 * 2) / (1 * 4)) = 0.5, half width = t(0.975, 1) * 0.5
        out = im.dersimonian_laird([0.0, 1.0], [0.1, 0.1])
        half = stats.t.ppf(0.975, 1) * 0.5
        self.assertAlmostEqual(out['ci_low_hk'], 0.5 - half)
        self.assertAlmostEqual(out['ci_high_hk'], 0.5 + half)
        self.assertTrue(np.isnan(im.dersimonian_laird([0.3], [0.1])['ci_low_hk']))

    def test_parse_geo_record(self):
        text = ('!Series_title = A study\n!Series_geo_accession = GSE10\n!Series_pubmed_id = 123\n'
                '!Series_contact_name = Jane,J.,Doe\n!Series_contact_institute = Some University\n'
                '!Series_relation = SubSeries of: GSE99\n!Series_relation = BioProject: https://x/PRJNA1\n')
        record = im.parse_geo_record(text)
        self.assertEqual(record['related'], ['GSE99'])
        self.assertEqual(record['pubmed'], ['123'])
        self.assertEqual(record['contact'], 'jane doe | some university')

    def test_projects_link_superseries_pubmed_contact_transitively(self):
        def record(related=(), pubmed=(), contact=''):
            return {'related': list(related), 'pubmed': list(pubmed), 'contact': contact}
        records = {'GSE1': record(related=['GSE100']), 'GSE2': record(related=['GSE100'], pubmed=['7']),
                   'GSE3': record(pubmed=['7']), 'GSE4': record(contact='a b | x'), 'GSE5': record(contact='a b | x'),
                   'GSE6': record(), 'GSE7': None}
        project, kinds = im.projects_from_records(records)
        self.assertEqual({project[s] for s in ('GSE1', 'GSE2', 'GSE3')}, {'GSE1'})
        self.assertEqual(project['GSE4'], project['GSE5'])
        self.assertNotEqual(project['GSE4'], project['GSE1'])
        self.assertEqual(project['GSE6'], 'GSE6')
        self.assertEqual(project['GSE7'], 'GSE7')  # a record that was not read stays its own project
        self.assertEqual(kinds['GSE1'], ['pubmed', 'series'])
        without, _ = im.projects_from_records(records, link_contact=False)
        self.assertNotEqual(without['GSE4'], without['GSE5'])
        self.assertEqual(without['GSE3'], 'GSE1')

    def test_lifespan_expectations(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            genage, drugage = Path(tmp) / 'g.csv', Path(tmp) / 'd.csv'
            pd.DataFrame({'symbol': ['Igf1r', 'Sirt6', 'Foo'], 'organism': ['Mus musculus'] * 3,
                          'longevity influence': ['Anti-Longevity', 'Pro-Longevity', 'Unclear']}).to_csv(genage, index=False)
            pd.DataFrame({'compound_name': ['Rapamycin', 'Rapamycin', 'Mixed', 'Mixed', 'Toxin', 'Weak'],
                          'species': ['Mus musculus'] * 6,
                          'avg_lifespan_change_percent': [20, -5, 10, -10, -30, 5],
                          'avg_lifespan_significance': ['S', 'NS', 'S', 'S', 'S', 'NS']}).to_csv(drugage, index=False)
            e = im.lifespan_expectations(genage, drugage)
        self.assertEqual(e[('GENE', 'igf1r', 'DECREASE')], -1)   # loss of an anti-longevity gene: lower outputs
        self.assertEqual(e[('GENE', 'sirt6', 'DECREASE')], 1)
        self.assertEqual(e[('GENE', 'sirt6', 'INCREASE')], -1)
        self.assertNotIn(('GENE', 'foo', 'DECREASE'), e)
        self.assertEqual(e[('EXPOSURE', 'rapamycin', 'EXPOSURE')], -1)  # extending: lower outputs
        self.assertEqual(e[('EXPOSURE', 'toxin', 'EXPOSURE')], 1)
        self.assertNotIn(('EXPOSURE', 'mixed', 'EXPOSURE'), e)
        self.assertNotIn(('EXPOSURE', 'weak', 'EXPOSURE'), e)

    def test_superseries_itself_joins_its_subseries(self):
        records = {'GSE20': {'related': ['GSE22'], 'pubmed': [], 'contact': ''},
                   'GSE22': {'related': ['GSE20', 'GSE21'], 'pubmed': [], 'contact': ''},
                   'GSE21': {'related': ['GSE22'], 'pubmed': [], 'contact': ''}}
        project, _ = im.projects_from_records(records)
        self.assertEqual(set(project.values()), {'GSE20'})

    def test_group_keys_separate_the_atlas_categories(self):
        """Correction of 2026-10-15: one intervention's comparisons are pooled only within one atlas category."""
        self.assertEqual(im.GROUP_KEYS + ['condition_category'],
                         ['species', 'canonical_target', 'intervention_class', 'manipulation', 'condition_category'])
        rows = pd.DataFrame({'species': ['mouse'] * 3, 'canonical_target': ['plx5622'] * 3,
                             'intervention_class': ['EXPOSURE'] * 3, 'manipulation': ['EXPOSURE'] * 3,
                             'condition_category': ['Drug', 'Drug', 'Environment'],
                             'project': ['GSE1', 'GSE2', 'GSE3'], 'gse_id': ['GSE1', 'GSE2', 'GSE3'],
                             'g': [-1.0, -1.0, 2.0], 'var_g': [0.2, 0.2, 0.2],
                             'comparison_name': ['drug vs vehicle', 'drug vs vehicle', 'repopulation vs control'],
                             'condition_name': ['drug', 'drug', 'repopulation'],
                             'v3_n_case': [4, 4, 4], 'v3_n_control': [4, 4, 4]})
        keys = im.GROUP_KEYS + ['condition_category']
        groups = {k[-1]: len(im.study_effects(g, 'project')) for k, g in rows.groupby(keys, dropna=False)}
        self.assertEqual(groups, {'Drug': 2, 'Environment': 1})
        self.assertEqual(sum(len(g) for _, g in rows.groupby(keys, dropna=False)), len(rows))

    def test_design_flags_without_an_audit_mark_nothing(self):
        usable = pd.DataFrame({'modality': ['mouse RNA-seq'], 'gse_id': ['GSE1'], 'comparison_name': ['a vs b'],
                               'condition_name': ['a']})
        flagged = im.design_flags(usable, ())
        self.assertFalse(bool(flagged.design_audited.iloc[0]))
        self.assertFalse(bool(flagged.design_tumor.iloc[0]))

    def test_design_flags_mark_tumour_and_cultured_samples(self):
        """Correction of 2026-10-17: both models must agree before a comparison is dropped."""
        audit = pd.DataFrame({'modality': ['mouse RNA-seq'] * 4, 'gse_id': ['GSE1', 'GSE2', 'GSE3', 'GSE4'],
                              'comparison_name': ['a', 'b', 'c', 'd'], 'condition_name': ['x'] * 4,
                              'tumor': ['YES', 'NO', None, 'NO'],
                              'specimen': ['ORGANISM_DIRECT', 'CELL_LINE', 'ORGANISM_DIRECT', 'ORGANISM_DIRECT']})
        usable = audit[im.KEY_COLUMNS].copy() if hasattr(im, 'KEY_COLUMNS') else audit[
            ['modality', 'gse_id', 'comparison_name', 'condition_name']].copy()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'run'
            path.mkdir()
            audit.to_json(path / 'comparisons.jsonl', orient='records', lines=True)
            flagged = im.design_flags(usable, [path]).set_index('gse_id')
        self.assertTrue(bool(flagged.loc['GSE1', 'design_tumor']))              # both models: tumour
        self.assertTrue(bool(flagged.loc['GSE2', 'design_not_from_organism']))  # both models: cell line
        self.assertFalse(bool(flagged.loc['GSE3', 'design_audited']))           # disagreed: stays, unflagged
        self.assertFalse(bool(flagged.loc['GSE4', 'design_tumor']))
        self.assertFalse(bool(flagged.loc['GSE4', 'design_not_from_organism']))

    def test_grouping_keeps_effects_whose_keys_are_missing(self):
        """dropna=False: an effect with no category is its own group instead of vanishing from the analysis."""
        rows = pd.DataFrame({'species': ['mouse', 'mouse'], 'canonical_target': ['x', 'x'],
                             'intervention_class': ['EXPOSURE'] * 2, 'manipulation': ['EXPOSURE'] * 2,
                             'condition_category': ['Drug', np.nan]})
        keys = im.GROUP_KEYS + ['condition_category']
        self.assertEqual(sum(len(g) for _, g in rows.groupby(keys, dropna=False)), 2)
        self.assertEqual(sum(len(g) for _, g in rows.groupby(keys)), 1)  # the old default dropped it


if __name__ == '__main__':
    unittest.main()
