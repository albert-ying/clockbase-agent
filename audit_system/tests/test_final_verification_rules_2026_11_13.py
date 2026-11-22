"""Synthetic software fixtures for final verification rules; no empirical results."""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from audit_system import final_verification_rules_2026_11_13 as rules
from audit_system import design_estimand_2026_10_10 as de
from audit_system.row_repairs_2026_09_28 import recompute, assign_fdr
from audit_system.error_census_2026_09_29 import classify, MATERIAL


def fixture(a=None, b=None, clock='Horvathv1'):
    # All values and identifiers are synthetic software fixtures.
    a = [20., 22., 27., 32.] if a is None else a
    b = [18., 23., 25., 28.] if b is None else b
    case, control = [f'A{i}' for i in range(len(a))], [f'B{i}' for i in range(len(b))]
    values = pd.Series(a + b, index=case + control)
    record = {'row_ordinal': 0, 'clock': clock, 'comparison_name': 'post vs pre',
              'condition_name': 'post', 'gse_id': 'GSE_SYNTHETIC', 'final_class': 'VERIFIED',
              'reasons': [], 'actions': [], 'soft_flags': [],
              'inputs': {'case_gsms': case, 'control_gsms': control},
              'recomputed': recompute(values, case, control)}
    series = {'subjects': {}, 'ages': {}, 'files': [], 'verification_metadata':
              {'headers': {}, 'titles': {}, 'individual_fields': {}}}
    return record, values, series, case, control


class FinalRules(unittest.TestCase):
    def test_gse126803_record_bytes_across_hash_seeds(self):
        # Integration fixture uses deposited study records and GEO metadata.
        from audit_system import paths
        source = Path(os.environ.get('CLOCKBASE_FINAL', paths.VERIFIER_ROOT / 'final-2026-11-13'))
        source = source / 'public-export/studies/hs_rna/GSE126803.json'
        if not source.is_file():
            self.skipTest(f'Deposited integration fixture unavailable: {source}')
        script = '''
import json, sys
from pathlib import Path
from audit_system import paths, final_verification_rules_2026_11_13 as rules
from audit_system.v3_public_export_2026_09_30 import write_json
study = json.loads(Path(sys.argv[1]).read_text())
meta = rules.metadata(paths.GEO, study['modality'], study['gse_id'])
rebuilt = []
for row in study['rows']:
    evidence = row.get('verification_rules') or {}
    if 'title_pair_conflict' in evidence:
        conflict = rules.title_pairs(meta, row['case_gsms'], row['control_gsms'])
        assert conflict and conflict['conflicting_individual_fields']
        evidence['title_pair_conflict'] = conflict
        rebuilt.append(conflict)
assert rebuilt
# Check construction order before the exporter canonicalizes mapping keys.
Path(sys.argv[3]).write_text(json.dumps(rebuilt))
write_json(Path(sys.argv[2]), study)
'''
        with tempfile.TemporaryDirectory() as folder:
            outputs, constructions = [], []
            for seed in (de.MIN_UNITS, de.MIN_UNITS + de.MIN_UNITS):
                output, construction = (Path(folder) / f'{seed}{suffix}'
                                        for suffix in ('.json', '.construction.json'))
                subprocess.run([sys.executable, '-c', script, str(source), str(output), str(construction)],
                               env={**os.environ, 'PYTHONHASHSEED': str(seed)}, check=True,
                               capture_output=True)
                outputs.append(output.read_bytes())
                constructions.append(construction.read_bytes())
            self.assertEqual(*constructions)
            self.assertEqual(*outputs)

    def test_year_bounds_inclusive_and_rate_scores_exempt(self):
        for mean, invalid in [(rules.LOWER_YEARS, False), (rules.UPPER_YEARS, False),
                              (rules.LOWER_YEARS - 1, True), (rules.UPPER_YEARS + 1, True)]:
            r, v, s, a, b = fixture([mean] * 4)
            rules.apply(r, v, s, 'dnam')
            self.assertEqual(r['final_class'] == 'PENDING_INPUT', invalid)
        for clock in ('DunedinPACE', 'DNAmPAI1', 'Hazard2_All_All_All_scaleddiff_EN_tAge'):
            r, v, s, a, b = fixture([-300.] * 4, clock=clock)
            rules.apply(r, v, s, 'dnam')
            self.assertEqual(r['final_class'], 'VERIFIED')

    def test_floor_requires_threshold_in_either_arm_and_defined_transform(self):
        for n_floor in (8, 9):
            r, v, s, a, b = fixture([-1.] * n_floor + [20.] * (10 - n_floor))
            rules.apply(r, v, s, 'dnam')
            self.assertEqual(r['final_class'] == 'PENDING_INPUT', n_floor == 9)
        r, v, s, a, b = fixture([-1. + 1e-11] * 4, [-1.] * 4)
        rules.apply(r, v, s, 'dnam')
        self.assertEqual(r['pending_reason'], [rules.OUTPUT_REASON])
        self.assertNotIn('p', r['recomputed'])
        r, v, s, a, b = fixture([-1.] * 4, clock='Hannum')
        rules.apply(r, v, s, 'dnam')
        self.assertNotEqual(r['final_class'], 'PENDING_INPUT')

    def test_single_cells_without_individuals_have_no_test(self):
        r, v, s, a, b = fixture(clock='Chronoage_All_All_WT_scaleddiff_EN_tAge')
        s['verification_metadata']['headers']['title'] = 'Single cell RNA-seq software fixture'
        s['verification_metadata']['titles'] = {g: f'Single cell {i}' for i, g in enumerate(a+b)}
        rules.apply(r, v, s, 'hs_rna')
        self.assertEqual(r['final_class'], 'DESCRIPTIVE')
        self.assertEqual(r['pending_reason'], [rules.CELL_REASON])
        self.assertNotIn('p', r['recomputed'])
        self.assertEqual(assign_fdr([r]), [])

    def test_individual_means_and_few_individuals(self):
        for n in (2, 3):
            r, v, s, a, b = fixture([2., 4., 5., 7., 9., 12.][:2*n], [1., 2., 4., 5., 6., 8.][:2*n], clock='RNA_SYNTHETIC')
            s['verification_metadata']['headers']['design'] = 'Smart-seq software fixture'
            s['verification_metadata']['titles'] = {g: f'Single cell {i}' for i, g in enumerate(a+b)}
            mapping = {**{g: f'case{i//2}' for i,g in enumerate(a)}, **{g: f'control{i//2}' for i,g in enumerate(b)}}
            s['verification_metadata']['individual_fields']['animal id'] = mapping
            rules.apply(r, v, s, 'mm_rna')
            self.assertEqual(r['recomputed']['n_case'], n)
            self.assertEqual(r['final_class'], 'VERIFIED' if n >= de.MIN_UNITS else 'DESCRIPTIVE')
            if n >= de.MIN_UNITS:
                aa, bb = np.asarray(v[a]).reshape(n, -1).mean(axis=1), np.asarray(v[b]).reshape(n, -1).mean(axis=1)
                self.assertAlmostEqual(r['recomputed']['p'], stats.ttest_ind(aa, bb, equal_var=False).pvalue)

    def test_incomplete_individual_map_is_not_padded_with_cells(self):
        r, v, s, a, b = fixture(clock='RNA_SYNTHETIC')
        s['verification_metadata']['headers']['title'] = 'Single nucleus RNA-seq software fixture'
        s['verification_metadata']['titles'] = {g: f'Single nucleus {i}' for i, g in enumerate(a+b)}
        s['verification_metadata']['individual_fields']['donor'] = {g: g for g in (a+b)[:-1]}
        rules.apply(r, v, s, 'hs_rna')
        self.assertEqual(r['pending_reason'], [rules.CELL_REASON])

    def test_explicit_mouse_tokens_in_cell_titles_identify_individuals(self):
        r, v, s, a, b = fixture(clock='RNA_SYNTHETIC')
        meta = s['verification_metadata']
        meta['headers']['title'] = 'Single cell RNA-seq software fixture'
        meta['titles'] = {**{g: f'single cell_KO_m{i//2}_A{i+1}' for i,g in enumerate(a)},
                          **{g: f'single cell_WT_m{i//2}_A{i+1}' for i,g in enumerate(b)}}
        rules.apply(r, v, s, 'mm_rna')
        self.assertEqual(r['recomputed']['n_case'], len(a)//2)
        self.assertEqual(r['recomputed']['n_control'], len(b)//2)
        self.assertEqual(r['design']['method'], 'UNIT_MEANS')

    def test_numbers_in_bulk_titles_do_not_establish_single_cells(self):
        r, v, s, a, b = fixture(clock='RNA_SYNTHETIC')
        meta = s['verification_metadata']
        meta['headers'] = {'title': 'Bulk RNA sequencing software fixture', 'summary': 'Accompanying single-cell RNA-seq is described elsewhere'}
        meta['titles'] = {g: f'bulk replicate {i}' for i,g in enumerate(a+b)}
        found = rules.cell_evidence(meta, a, b)
        self.assertFalse(found['detected'])
        self.assertTrue(found['ambiguous'])
        meta['headers'] = {}
        meta['titles'] = {g: f'sample_{i}' for i,g in enumerate(a+b)}
        self.assertFalse(rules.cell_evidence(meta, a, b)['detected'])

    def test_summary_and_droplet_libraries_are_not_single_cells(self):
        r, v, s, a, b = fixture(clock='RNA_SYNTHETIC')
        meta = s['verification_metadata']
        for protocol in ('10x single nucleus RNA-seq', 'InDrop single cell RNA-seq',
                         'Spatial transcriptomics', 'Smart-seq whole embryos', 'Bulk organoids'):
            meta['headers'] = {'title': protocol, 'summary': 'We performed single-cell RNA sequencing.'}
            meta['titles'] = {g: f'Donor library {i}' for i, g in enumerate(a+b)}
            self.assertFalse(rules.cell_evidence(meta, a, b)['detected'], protocol)
        meta['headers'] = {'title': 'Single-cell RNA-seq'}
        meta['titles'] = {g: f'Single cell {i}' for i, g in enumerate(a)}
        self.assertFalse(rules.cell_evidence(meta, a, b)['detected'])

    def test_cell_counts_clones_and_donor_library_replicates_are_not_cells(self):
        r, v, s, a, b = fixture(clock='RNA_SYNTHETIC')
        meta = s['verification_metadata']
        meta['headers'] = {'summary': 'We performed single-cell RNA sequencing.'}
        meta['characteristics'] = {g: {'cell number': str(len(a+b))} for g in a+b}
        for template in ('Single cell RNAseq replicate {}', 'single cell clone {}', 'ES cell {}'):
            meta['titles'] = {g: template.format(i) for i, g in enumerate(a+b)}
            self.assertFalse(rules.cell_evidence(meta, a, b)['detected'], template)
        meta['characteristics'] = {}
        meta['titles'] = {g: f'Single cell RNAseq_{i}' for i, g in enumerate(a+b)}
        meta['individual_fields'] = {'donor': {g: g for g in a+b}}
        self.assertFalse(rules.cell_evidence(meta, a, b)['detected'])

    def test_whole_embryo_stages_and_unnumbered_assay_libraries_are_not_cells(self):
        r, v, s, a, b = fixture(clock='RNA_SYNTHETIC')
        meta = s['verification_metadata']
        meta['headers'] = {'summary': 'We performed single-cell RNA sequencing.'}
        meta['titles'] = {g: f'control {len(a)}-cell {i} embryo' for i, g in enumerate(a+b)}
        self.assertFalse(rules.cell_evidence(meta, a, b)['detected'])
        meta['titles'] = {g: 'Primary single cell RNA-seq' for g in a+b}
        self.assertFalse(rules.cell_evidence(meta, a, b)['detected'])

    def test_output_screen_uses_arm_age_coverage_and_inclusive_deviation(self):
        r, v, s, a, b = fixture(clock='PhenoAge')
        s['ages'] = {g: float(v[g]) + rules.AGE_ERROR_YEARS for g in a}
        evidence = rules.output_evidence(r['clock'], v, a, b, s)
        self.assertTrue(evidence['arms']['case']['uses_recorded_age'])
        self.assertFalse(evidence['arms']['case']['outside'])
        s['ages'] = {g: age + 1 for g, age in s['ages'].items()}
        self.assertTrue(rules.output_evidence(r['clock'], v, a, b, s)['detected'])
        s['ages'] = {}
        s['verification_metadata']['age_fields_years'] = {'age (years)': {g: v[g] for g in a[:len(a)//2]}}
        self.assertTrue(rules.output_evidence(r['clock'], v, a, b, s)['arms']['case']['uses_recorded_age'])
        s['verification_metadata']['age_fields_years']['age (years)'].pop(a[0])
        self.assertFalse(rules.output_evidence(r['clock'], v, a, b, s)['arms']['case']['uses_recorded_age'])

    def test_non_assay_single_cell_phrases_and_small_nuclear_rna_are_not_cell_sequencing(self):
        for text in ('RNA was extracted from single-cell suspensions generated by FACS.',
                     'We used single cell phenotyping to identify cells.',
                     'We studied gene function at single cell resolution.',
                     'The small nuclear RNA (snRNA) spliceosome was altered.'):
            r,v,s,a,b=fixture(clock='RNA_SYNTHETIC')
            s['verification_metadata']['headers']['summary']=text
            self.assertFalse(rules.cell_evidence(s['verification_metadata'],a,b)['detected'],text)

    def test_background_assay_and_bulk_library_protocol_are_ambiguous(self):
        for headers in ({'summary':'Single-cell RNA sequencing recently identified a candidate gene.',
                         'design':'RNA was extracted from sorted cells using TRIzol.'},
                        {'summary':'We performed single nucleus transcriptomics.',
                         'design':'Stranded mRNA-seq profiles were obtained from muscle tissue.'},
                        {'title':'Whole colon using Takara SMART-Seq Stranded kit.'},
                        {'summary':'Using single cell RNA-seq technology on whole embryos.',
                         'design':'Global mRNA profiles of whole preimplantation embryos.'}):
            r,v,s,a,b=fixture(clock='RNA_SYNTHETIC')
            s['verification_metadata']['headers']=headers
            found=rules.cell_evidence(s['verification_metadata'],a,b)
            self.assertFalse(found['detected'],headers)
            self.assertTrue(found['ambiguous'],headers)

    def test_enumerated_animal_libraries_are_individuals(self):
        # Synthetic design enumerates one animal per library in each group.
        titles={f'G{i}':f'{sex}_{treatment}_{j}' for i,(sex,treatment,j) in enumerate(
            (sex,treatment,j) for sex in ('Female','Male') for treatment in ('Saline','LPS') for j in range(1,4))}
        meta={'titles':titles,'headers':{'design':'10x libraries. 3 saline female, 3 saline male, 3 LPS female, 3 LPS male animals'}}
        mapping=rules.animal_library_individuals(meta)
        self.assertEqual(set(mapping),set(titles))
        self.assertEqual(len(set(mapping.values())),len(titles))
        meta['headers']['design']='10x libraries from mice; sample titles have replicate numbers.'
        self.assertFalse(rules.animal_library_individuals(meta))

    def test_well_coordinates_require_a_cell_or_plate_signal(self):
        r, v, s, a, b = fixture(clock='RNA_SYNTHETIC')
        meta = s['verification_metadata']
        meta['titles'] = {g: f'SS2_software_A{i+1}' for i,g in enumerate(a+b)}
        self.assertTrue(rules.cell_evidence(meta,a,b)['detected'])
        meta['titles'] = {g: f'bulk_A{i+1}' for i,g in enumerate(a+b)}
        found = rules.cell_evidence(meta,a,b)
        self.assertFalse(found['detected'])
        self.assertTrue(found['ambiguous'])

    def test_title_pair_bijection_and_paired_statistic(self):
        r, v, s, a, b = fixture(clock='RNA_SYNTHETIC')
        s['verification_metadata']['titles'] = {**{g: f'{i}B Young HIIT PostTraining' for i,g in enumerate(a)},
                                                **{g: f'{i}A Young HIIT PreTraining' for i,g in enumerate(b)}}
        rules.apply(r, v, s, 'hs_rna')
        self.assertEqual(r['design']['method'], 'PAIRED')
        self.assertAlmostEqual(r['recomputed']['p'], stats.ttest_rel(v[a], v[b]).pvalue)
        self.assertEqual(len(r['verification_rules']['title_pairs']['tokens']),len(a))
        s['verification_metadata']['titles'][b[0]] = '99A Young HIIT PreTraining'
        self.assertIsNone(rules.title_pairs(s['verification_metadata'],a,b))
        self.assertIsNone(rules.title_pairs(s['verification_metadata'],a,b[:-1]))

    def test_title_numbers_without_temporal_labels_are_not_pairs(self):
        r, v, s, a, b = fixture(clock='RNA_SYNTHETIC')
        s['verification_metadata']['titles'] = {**{g:f'{i}B KO' for i,g in enumerate(a)}, **{g:f'{i}A WT' for i,g in enumerate(b)}}
        self.assertIsNone(rules.title_pairs(s['verification_metadata'],a,b))

    def test_title_pairs_cannot_override_a_contradictory_GEO_subject_field(self):
        r,v,s,a,b=fixture(clock='RNA_SYNTHETIC')
        meta=s['verification_metadata']
        meta['titles']={**{g:f'{i}B PostTraining' for i,g in enumerate(a)},
                        **{g:f'{i}A PreTraining' for i,g in enumerate(b)}}
        meta['individual_fields']['subject']={g:g for g in a+b}
        original=dict(r['recomputed'])
        rules.apply(r,v,s,'hs_rna')
        self.assertEqual(r['recomputed'],original)
        self.assertIn('title_pair_conflict',r['verification_rules'])
        self.assertNotIn('title_pairs',r['verification_rules'])

    def test_error_census_includes_both_defects(self):
        for code in (rules.OUTPUT_CODE, rules.CELL_CODE):
            r, _, _, _, _ = fixture()
            r['reasons'] = ['MEMBERSHIP_RECOVERED',code]
            defects, reproduced = classify(r,{})
            self.assertIn(code,defects)
            self.assertIn(code,MATERIAL)
            self.assertTrue(reproduced)

    def test_one_recorded_animal_per_library_is_not_a_pseudoreplication_defect(self):
        r,v,s,a,b=fixture(clock='RNA_SYNTHETIC')
        s['verification_metadata']['headers']['title']='Single cell RNA-seq software fixture'
        s['verification_metadata']['individual_fields']['animal id']={g:g for g in a+b}
        rules.apply(r,v,s,'mm_rna')
        self.assertEqual(r['final_class'],'VERIFIED')
        self.assertNotIn(rules.CELL_CODE,r['reasons'])
        self.assertEqual(r['recomputed']['n_case'],len(a))

    def test_complete_donor_field_takes_precedence_over_title_tokens(self):
        r,v,s,a,b=fixture(clock='RNA_SYNTHETIC')
        meta=s['verification_metadata']
        meta['headers']['title']='Single cell RNA-seq software fixture'
        meta['titles']={g:f'single cell_patient{i}_A1' for i,g in enumerate(a+b)}
        meta['individual_fields']['donor']={**dict.fromkeys(a,'case donor'),**dict.fromkeys(b,'control donor')}
        rules.apply(r,v,s,'hs_rna')
        self.assertEqual(r['final_class'],'DESCRIPTIVE')
        self.assertEqual(r['recomputed']['n_case'],1)
        self.assertEqual(r['verification_rules']['single_cells']['individual_field'],'donor')

    def test_deposited_arm_decisions_retain_precedence_on_overlapping_rules(self):
        from audit_system.clock_output_arms_2026_11_12 import apply, comparison_key
        r, v, s, _, _ = fixture([-300.] * 4)
        rules.apply(r, v, s, 'dnam')
        key = comparison_key('dnam', r)
        decision = {'final_class': 'NOT_REPAIRABLE', 'reason': 'arms defined by a clock output',
                    'reason_code': 'ARMS_DEFINED_BY_CLOCK_OUTPUT'}
        apply([r], {key: decision}, 'dnam')
        self.assertEqual(r['final_class'], decision['final_class'])
        self.assertIn(rules.OUTPUT_CODE, r['reasons'])
        self.assertIn(decision['reason_code'], r['reasons'])


if __name__ == '__main__':
    unittest.main()
