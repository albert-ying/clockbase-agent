"""Unit tests for the streaming DNAm repair runner. All fixtures are synthetic."""
import gzip
import json
import tempfile
import threading
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

from src import dnam_stream_repair_2026_09_26 as repair


SAMPLES = {
    'GSM1': {'title': ['Patient A [S01]'], 'characteristics_ch1': ['age: 40', 'sample id: S01'],
             'supplementary_file': ['ftp://x/GSM1_2001_R01C01_Grn.idat.gz'], 'platform_id': ['GPL13534']},
    'GSM2': {'title': ['Patient B [S02]'], 'characteristics_ch1': ['age: 60', 'sample id: S02'],
             'supplementary_file': ['ftp://x/GSM2_2001_R02C01_Grn.idat.gz'], 'platform_id': ['GPL13534']},
    'GSM3': {'title': ['Patient C [S03]'], 'characteristics_ch1': ['age: 80', 'sample id: S03'],
             'supplementary_file': ['ftp://x/GSM3_2001_R03C01_Grn.idat.gz'], 'platform_id': ['GPL13534']},
}


def write_table(path, header, rows, sep='\t'):
    with gzip.open(path, 'wt') as handle:
        handle.write(sep.join(header) + '\n')
        for row in rows:
            handle.write(sep.join(str(v) for v in row) + '\n')
    return path



def approx(value, tol=1e-6):
    class _Approx:
        def __eq__(self, other):
            return abs(float(other) - value) <= tol * max(1.0, abs(value))
    return _Approx()


class StreamRepairTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_aliases_are_unique_and_drop_shared_tokens(self):
        samples = {**SAMPLES, 'GSM4': {'title': ['Patient A [S01]'], 'characteristics_ch1': []}}
        aliases = repair.sample_aliases(samples)
        assert aliases['gsm1'][0] == 'GSM1'
        assert 's02' in aliases and aliases['s02'][0] == 'GSM2'
        assert 's01' not in aliases  # owned by GSM1 and GSM4
        assert aliases['2001_r03c01'] == ('GSM3', 'idat_basename')


    def test_table_maps_titles_ignores_detection_and_annotation(self):
        path = write_table(self.tmp / 't.txt.gz',
                           ['ID_REF', 'S01.AVG_Beta', 'S01.Detection Pval', 'S02.AVG_Beta', 'S02.Detection Pval',
                            'S03.AVG_Beta', 'S03.Detection Pval', 'CHR'],
                           [['cg0000001', .1, .001, .2, .001, .3, .5, 1], ['cg0000002', .4, 0, 'NA', 0, .6, 0, 2]])
        frame, evidence, sample_map = repair.read_beta_table(path, SAMPLES, kind='beta')
        assert list(frame.columns) == ['GSM1', 'GSM2', 'GSM3']
        assert np.isnan(frame.loc['cg0000002', 'GSM2'])
        assert evidence['detection_p_columns_ignored'] == 3 and evidence['annotation_columns'] == ['CHR']
        assert evidence['value_scale'] == 'beta' and set(sample_map.gsm) == {'GSM1', 'GSM2', 'GSM3'}


    def test_unmapped_numeric_column_fails_route(self):
        path = write_table(self.tmp / 't.txt.gz', ['ID_REF', 'GSM1', 'GSM2', 'mystery'],
                           [['cg1', .1, .2, .3]])
        with self.assertRaises(repair.RouteFailure) as failure:
            repair.read_beta_table(path, SAMPLES, kind='beta')
        assert failure.exception.code == 'UNMAPPED_COLUMNS'


    def test_mvalues_convert_and_percent_scale_rejected(self):
        path = write_table(self.tmp / 'm.txt.gz', ['ID_REF', 'GSM1', 'GSM2'], [['cg1', -2, 0], ['cg2', 3, 1]])
        frame, evidence, _ = repair.read_beta_table(path, SAMPLES, kind='mvalue')
        assert evidence['value_scale'] == 'mvalue_converted'
        assert frame.loc['cg1', 'GSM1'] == approx(0.2)
        assert frame.loc['cg1', 'GSM2'] == approx(0.5)
        path = write_table(self.tmp / 'p.txt.gz', ['ID_REF', 'GSM1', 'GSM2'], [['cg1', 10, 90]])
        with self.assertRaises(repair.RouteFailure) as failure:
            repair.read_beta_table(path, SAMPLES, kind='beta')
        assert failure.exception.code == 'UNRECOGNISED_VALUE_SCALE'


    def test_unlabelled_feature_column_and_geo_preamble(self):
        path = self.tmp / 'g.txt.gz'
        with gzip.open(path, 'wt') as handle:
            handle.write('!Series_title\t"x"\n!series_matrix_table_begin\n')
            handle.write('GSM1\tGSM2\ncg1\t0.1\t0.2\ncg2\t0.3\t0.4\n!series_matrix_table_end\n')
        frame, evidence, _ = repair.read_beta_table(path, SAMPLES, kind='beta')
        assert frame.shape == (2, 2) and evidence['features'] == 2


    def test_epicv2_replicates_average(self):
        frame = pd.DataFrame({'GSM1': [.2, .4, .5]}, index=['cg00000029_TC21', 'cg00000029_TC22', 'cg00000100_BC11'])
        collapsed, removed = repair.collapse_feature_ids(frame)
        assert removed == 1
        assert collapsed.loc['cg00000029', 'GSM1'] == approx(.3)


    def test_duplicate_sample_columns_fail(self):
        aliases = repair.sample_aliases(SAMPLES)
        with self.assertRaises(repair.RouteFailure) as failure:
            repair.classify_columns(['ID_REF', 'GSM1', 'S01'], aliases, list(SAMPLES))
        assert failure.exception.code == 'DUPLICATE_SAMPLE_COLUMNS'


    def test_routes_prefer_beta_then_mvalue_then_idat(self):
        row = {'class': 'IDAT_AVAILABLE',
               'files': 'a_M.txt.gz[50,mvalue] | a_beta.txt.gz[100,beta] | a_RAW.tar[900] | sig.txt.gz[70,signals]'}
        routes = repair.routes_for(row)
        assert [r['file'] for r in routes] == ['a_beta.txt.gz', 'a_M.txt.gz', 'a_RAW.tar']
        assert repair.routes_for({**row, 'class': 'BETA_TABLE_AVAILABLE'})[-1]['type'] == 'table'


    def test_verdict_rules_coverage_and_degeneracy(self):
        beta = pd.DataFrame(np.nan, index=[f'cg{i}' for i in range(10)], columns=['GSM1', 'GSM2', 'GSM3'])
        beta.iloc[:9] = 0.5          # 90% of sites observed
        sites = {'Good': [f'cg{i}' for i in range(10)], 'Half': ['cg0', 'cg1', 'cg9', 'cgX'],
                 'Flat': ['cg0'], 'SexEstimation': ['cg0']}
        predictions = pd.DataFrame({'Good': [1., 2., 3.], 'Half': [1., 2., 4.], 'Flat': [5., 5., 5.]},
                                   index=beta.columns)
        rows, coverage = repair.verdicts(beta, predictions, sites, pd.Series(dtype=float))
        by = {r['clock']: r for r in rows}
        assert by['Good']['verdict'] == 'VERIFIED'
        assert by['Half']['verdict'] == 'WITHHELD' and by['Half']['coverage_min'] == approx(.5)
        assert by['Flat']['verdict'] == 'WITHHELD' and 'DEGENERATE_OUTPUT' in by['Flat']['reasons']
        assert 'SexEstimation' not in by and 'SexEstimation' in coverage
        sites['Mid'] = ['cg0', 'cg1', 'cg9', 'cg2', 'cg3', 'cg4', 'cg5', 'cg6', 'cg7', 'cgX']
        predictions['Mid'] = [1., 2., 3.]
        rows, _ = repair.verdicts(beta, predictions, sites, pd.Series(dtype=float))
        assert {r['clock']: r['verdict'] for r in rows}['Mid'] == 'VERIFIED'
        sites['Mid'] = ['cg0', 'cg1', 'cg2', 'cgX', 'cgY', 'cg3', 'cg4']  # 5/7 = 0.71
        rows, _ = repair.verdicts(beta, predictions, sites, pd.Series(dtype=float))
        assert {r['clock']: r['verdict'] for r in rows}['Mid'] == 'PROVISIONAL'


    def test_recorded_ages_parse_only_age_keys(self):
        samples = {'GSM1': {'characteristics_ch1': ['age: 40', 'gestational age: 12']},
                   'GSM2': {'characteristics_ch1': ['Age (years): 61.5']},
                   'GSM3': {'characteristics_ch1': ['age: unknown']}}
        ages = repair.recorded_ages(samples)
        assert ages.to_dict() == {'GSM1': 40.0, 'GSM2': 61.5}


    def test_disk_budget_blocks_until_release(self):
        patcher = mock.patch.object(repair.shutil, 'disk_usage', lambda _: type('U', (), {'free': 10 ** 15})())
        patcher.start()
        self.addCleanup(patcher.stop)
        budget = repair.DiskBudget(self.tmp, max_inflight=100, min_free=0)
        stop = threading.Event()
        assert budget.acquire(80, stop)
        started = threading.Event()

        def second():
            started.set()
            budget.acquire(50, stop)

        thread = threading.Thread(target=second)
        thread.start()
        started.wait()
        thread.join(timeout=0.3)
        assert thread.is_alive()          # 80 + 50 exceeds 100
        budget.release(80)
        thread.join(timeout=5)
        assert not thread.is_alive() and budget.reserved == 50
        stop.set()
        assert budget.acquire(500, stop) is False  # would exceed the budget while another study runs


    def test_idat_pairs_accept_green_red_channel_names(self):
        import io
        import tarfile
        archive = self.tmp / 'raw.tar'
        with tarfile.open(archive, 'w') as tar:
            for name in ['GSM1_x_wt_green.idat.gz', 'GSM1_x_wt_red.idat.gz', 'GSM2_y_Grn.idat.gz', 'GSM2_y_Red.idat.gz',
                         'GSM3_z_Grn.idat.gz']:
                info = tarfile.TarInfo(name); info.size = 4
                tar.addfile(info, io.BytesIO(b'IDAT'))
        pairs, problems = repair.idat_pairs(archive, self.tmp / 'idat', SAMPLES)
        assert set(pairs) == {'GSM1', 'GSM2'} and 'GSM3' in problems
        assert Path(pairs['GSM1'] + '_Grn.idat.gz').exists() and Path(pairs['GSM1'] + '_Red.idat.gz').exists()


    def test_title_parts_and_separator_normalisation_map_columns(self):
        samples = {'GSM7': {'title': ['202242420201: Primed_NS_control']},
                   'GSM8': {'title': ['202242420201: Primed DNMT3A KD1']}}
        aliases = repair.sample_aliases(samples)
        assert aliases['primed_ns_control'][0] == 'GSM7'
        assert aliases[repair.normalize_token('Primed_DNMT3A-KD1')][0] == 'GSM8'
        assert '202242420201' not in aliases  # shared chip barcode is not an identity


    def test_sample_qc_excludes_bad_samples_from_clock_verdicts(self):
        beta = pd.DataFrame(0.5, index=[f'cg{i}' for i in range(20)], columns=['GSM1', 'GSM2', 'GSM3', 'GSM4'])
        beta.iloc[:12, 3] = np.nan                 # GSM4 is 60% missing
        status = repair.sample_qc(beta)
        assert status.loc['GSM4', 'status'] == 'WITHHELD' and (status.status == 'RETAINED').sum() == 3
        retained = status.index[status.status == 'RETAINED']
        sites = {'C': [f'cg{i}' for i in range(20)], 'Broken': ['cg0']}
        predictions = pd.DataFrame({'C': [1., 2., 3., 9.], 'Broken': np.nan}, index=beta.columns)
        rows, _ = repair.verdicts(beta, predictions, sites, pd.Series(dtype=float), retained,
                                  {'Broken': 'Missing required CpG sites'})
        by = {r['clock']: r for r in rows}
        assert by['C']['verdict'] == 'VERIFIED' and by['C']['samples_assessed'] == 3
        assert by['Broken']['verdict'] == 'WITHHELD' and by['Broken']['reasons'] == ['MISSING_REQUIRED_SITES']


    def test_status_json_replaces_nan_and_inf(self):
        path = self.tmp / 's.json'
        repair.write_json(path, {'r': float('nan'), 'x': [np.float64('inf'), np.int64(3), 1.5]})
        assert path.read_text().count('null') == 2 and '"x": [\n    null,\n    3,\n    1.5' in path.read_text()


    def test_enrichment_assays_are_withheld_even_with_good_coverage(self):
        samples = {'GSM1': {'title': ['5hmCIP_450K (0015-C-HMC)']}, 'GSM2': {'title': ['noIP_450K (0015-C)']},
                   'GSM3': {'title': ['blood'], 'characteristics_ch1': ['assay: oxBS']},
                   'GSM4': {'title': ['Chip 3 position R01C01'], 'description': ['whole blood']},
                   'GSM5': {'title': ['H3K4me3 ChIP input']},
                   'GSM6': {'title': ['monocyte'], 'characteristics_ch1': ['ChIP: 6264509100']}}
        beta = pd.DataFrame(0.5, index=['cg1', 'cg2'], columns=list(samples))
        status = repair.sample_qc(beta, samples)
        assert status.status.to_dict() == {'GSM1': 'WITHHELD', 'GSM2': 'RETAINED', 'GSM3': 'WITHHELD',
                                           'GSM4': 'RETAINED', 'GSM5': 'WITHHELD', 'GSM6': 'RETAINED'}
        assert status.loc['GSM1', 'reason'].startswith('NON_STANDARD_ASSAY')


    def test_mesa_and_horvath_style_columns_and_duplicate_headers(self):
        samples = {'GSM10': {'title': ['100001_peripheral_CD14 [methylation]']},
                   'GSM11': {'title': ['100002_peripheral_CD14 [methylation]']}}
        path = write_table(self.tmp / 'mesa.txt.gz', ['ID_REF', '100001.Mvalue', '100002.Mvalue', 'Detection Pval',
                                                      'Detection Pval'], [['cg1', -1, 1, 0, 0]])
        frame, evidence, _ = repair.read_beta_table(path, samples, kind='mvalue')
        assert list(frame.columns) == ['GSM10', 'GSM11'] and evidence['detection_p_columns_ignored'] == 2
        samples = {'GSM20': {'title': ['OldTsimane: genomic DNA from whole blood of subject Sample1'],
                             'characteristics_ch1': ['dnamage: 40']},
                   'GSM21': {'title': ['OldTsimane: genomic DNA from whole blood of subject Sample2'],
                             'characteristics_ch1': ['dnamage: 50']}}
        path = write_table(self.tmp / 'h.csv.gz', ['', 'Sample1', 'Sample2'], [['cg1', .2, .3]], sep=',')
        frame, _, _ = repair.read_beta_table(path, samples, kind='beta')
        assert list(frame.columns) == ['GSM20', 'GSM21']
        assert repair.author_dnam_ages(samples).to_dict() == {'GSM20': 40.0, 'GSM21': 50.0}


    def test_same_named_beta_and_detection_columns_are_separated_by_values(self):
        rng = np.random.default_rng(0)
        rows = [[f'cg{i}', b1, p1, b2, p2] for i, (b1, p1, b2, p2) in
                enumerate(zip(rng.uniform(0, 1, 200), rng.uniform(0, 0.001, 200), rng.uniform(0, 1, 200),
                              rng.uniform(0, 0.001, 200)))]
        path = write_table(self.tmp / 'd.csv.gz', ['ID_REF', 'S01', 'S01', 'S02', 'S02'], rows, sep=',')
        frame, evidence, _ = repair.read_beta_table(path, SAMPLES, kind='beta')
        assert list(frame.columns) == ['GSM1', 'GSM2'] and evidence['unlabelled_detection_p_columns_inferred'] == 2
        assert frame['GSM1'].median() > 0.2
        rows = [[f'cg{i}', b, 'p<0.01', b2, 'p<0.01'] for i, (b, b2) in enumerate(zip(rng.uniform(0, 1, 50), rng.uniform(0, 1, 50)))]
        path = write_table(self.tmp / 't.csv.gz', ['ID_REF', 'S01', 'S01', 'S02', 'S02'], rows, sep=',')
        frame, evidence, _ = repair.read_beta_table(path, SAMPLES, kind='beta')
        assert list(frame.columns) == ['GSM1', 'GSM2'] and evidence['unlabelled_detection_p_columns_inferred'] == 2
        rows = [[f'cg{i}', b, b2] for i, (b, b2) in enumerate(zip(rng.uniform(0, 1, 50), rng.uniform(0, 1, 50)))]
        path = write_table(self.tmp / 'e.csv.gz', ['ID_REF', 'S01', 'S01'], rows, sep=',')
        with self.assertRaises(repair.RouteFailure) as failure:
            repair.read_beta_table(path, SAMPLES, kind='beta')
        assert failure.exception.code == 'DUPLICATE_SAMPLE_COLUMNS'


    def test_meth_intensity_and_beta_pair_keeps_beta_and_blank_cells_are_missing(self):
        rows = [[f'cg{i}', 2000 + i, 0.2, 3000 + i, ' '] for i in range(30)]
        path = write_table(self.tmp / 'mb.txt.gz', ['ID_REF', 'S01_Meth', 'S01_Beta', 'S02_Meth', 'S02_Beta'], rows)
        rows[0][4] = 0.4
        path = write_table(self.tmp / 'mb.txt.gz', ['ID_REF', 'S01_Meth', 'S01_Beta', 'S02_Meth', 'S02_Beta'], rows)
        frame, evidence, _ = repair.read_beta_table(path, SAMPLES, kind='beta')
        assert list(frame.columns) == ['GSM1', 'GSM2'] and frame['GSM1'].max() <= 1
        assert frame['GSM2'].notna().sum() == 1

    def test_unreadable_idat_is_isolated_not_fatal(self):
        import subprocess
        (self.tmp / 'idat').mkdir()
        calls = []

        def fake_run(command, **kwargs):
            config = json.loads(Path(command[-1]).read_text())
            calls.append(config['sample_ids'])
            if 'GSM2' in config['sample_ids']:
                return subprocess.CompletedProcess(command, 1, '', 'invalid UTF-8 input in readChar()')
            pd.DataFrame({g: [0.5] for g in config['sample_ids']}, index=['cg1']).rename_axis('ID_REF') \
                .to_csv(config['beta_output'], sep='\t', compression='gzip')
            pd.DataFrame({'sample_id': config['sample_ids']}).to_csv(config['qc_output'], sep='\t', index=False)
            Path(config['receipt_output']).write_text('{"status": "PASS"}')
            return subprocess.CompletedProcess(command, 0, '', '')

        pairs = {'GSM1': str(self.tmp / 'idat/GSM1_a'), 'GSM2': str(self.tmp / 'idat/GSM2_b')}
        with mock.patch.object(repair.subprocess, 'run', fake_run), open(self.tmp / 'log.txt', 'w') as log:
            frame, evidence, _ = repair.process_idats(pairs, SAMPLES, self.tmp / 'idat', log)
        assert list(frame.columns) == ['GSM1'] and list(evidence['unreadable_idat_gsms']) == ['GSM2']
        assert calls[0] == ['GSM1', 'GSM2'] and ['GSM1'] in calls and ['GSM2'] in calls


if __name__ == '__main__':
    unittest.main()
