"""Software regression fixtures for merging predicted samples; no clock runs."""
import importlib.util
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch,Mock
import pandas as pd

path=Path(__file__).resolve().parents[2]/'src_rna/2_run_rna_clock.py'
spec=importlib.util.spec_from_file_location('rna_identity_guard_test',path)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class RnaIdentityGuards(unittest.TestCase):
    def test_strict_reader_rejects_width_identity_and_feature_defects(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'input.csv'
            good=',gene_a,gene_b\nGSM1,1,2\nGSM2,3,4\n'
            path.write_text(good)
            parsed=module.read_strict_rna_csv(path)
            self.assertEqual(parsed.index.tolist(),['GSM1','GSM2'])
            self.assertEqual(parsed.columns.tolist(),['gene_a','gene_b'])
            for text in [good.replace('GSM2,3,4','GSM2,3'),
                         good.replace('GSM2,3,4','GSM1,3,4'),
                         good.replace('GSM2,3,4','sample2,3,4'),
                         good.replace('GSM2,3,4','GSM２,3,4'),
                         good.replace('gene_a,gene_b','gene_a,gene_a'),
                         good.replace('gene_a,gene_b','gene_a,   '),
                         ',gene_a,gene_b\n']:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    module.read_strict_rna_csv(path)

    def test_same_width_source_mutation_between_scan_and_parse_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'input.csv'
            path.write_text(',gene_a\nGSM1,1.0\n')
            original_read_csv=module.pd.read_csv
            def mutate_then_read(*args,**kwargs):
                path.write_text(path.read_text().replace('1.0','2.0'))
                return original_read_csv(*args,**kwargs)
            with patch.object(module.pd,'read_csv',side_effect=mutate_then_read):
                with self.assertRaisesRegex(ValueError,'bytes changed during parsing'):
                    module.read_strict_rna_csv(path)

    def test_observed_truncated_yugene_diff_sources_are_rejected_unchanged(self):
        cases=[
            (Path('/Entropy/kying/GEO_data/GEO_mm_rna/processed_data/GSE107493/GSE107493_yugene_diff_data.csv'),
             'c7071e4ac4f71b48630d4c517981f1108e826912493cd8b0e895872497439f73',
             r'record 53 has 307 fields; expected 18697'),
            (Path('/Entropy/kying/GEO_data/GEO_mm_rna/processed_data/GSE152298/GSE152298_yugene_diff_data.csv'),
             '4378d78299566083aba133e3513083346d1e0984002648e2590326548f615232',
             r'record 9 has 2173 fields; expected 18697'),
        ]
        for path,expected_hash,message in cases:
            with self.subTest(path=path):
                before=hashlib.sha256(path.read_bytes()).hexdigest()
                self.assertEqual(before,expected_hash)
                with self.assertRaisesRegex(ValueError,message):
                    module.read_strict_rna_csv(path)
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),before)

    def test_existing_cache_is_structurally_checked_without_rewrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'input/GSE_FIXTURE').mkdir(parents=True);(root/'output').mkdir()
            output=root/'output/GSE_FIXTURE.csv'
            output.write_text('sample_id,clock_a\nGSM1,1.0\n')
            before=output.read_bytes()
            module.process_file_folder('GSE_FIXTURE',str(root/'input'),str(root/'output'),'mouse')
            self.assertEqual(output.read_bytes(),before)

            output.write_text('sample_id,clock_a\nNA,1.0\n')
            before=output.read_bytes()
            with self.assertRaisesRegex(ValueError,'canonical GSM'):
                module.process_file_folder('GSE_FIXTURE',str(root/'input'),str(root/'output'),'mouse')
            self.assertEqual(output.read_bytes(),before)

    def test_safe_wrapper_retains_input_and_cache_validation_receipts(self):
        for phase,make_output in [('input_validation',False),('cache_validation',True)]:
            with self.subTest(phase=phase),tempfile.TemporaryDirectory() as directory:
                root=Path(directory);(root/'input/GSE_FIXTURE').mkdir(parents=True);(root/'output').mkdir()
                if make_output:
                    target=root/'output/GSE_FIXTURE.csv'
                    target.write_text('sample_id,clock_a\nNA,1.0\n')
                else:
                    target=root/'input/GSE_FIXTURE/GSE_FIXTURE_scaled_data.csv'
                    target.write_text(',gene\nGSM1,1.0\nGSM2\n')
                before=target.read_bytes()
                module.safe_process_file_folder(
                    'GSE_FIXTURE',str(root/'input'),str(root/'output'),'mouse'
                )
                self.assertEqual(target.read_bytes(),before)
                receipts=list((root/'output').glob('*_audit_failure_*.json'))
                self.assertEqual(len(receipts),1)
                record=json.loads(receipts[0].read_text())
                self.assertEqual(record['phase'],phase)
                self.assertEqual(record['status'],'CANNOT VERIFY')
                self.assertEqual(record['existing_completed_export_preserved'],make_output)
                self.assertFalse(record['new_completed_export_written'])

    def test_all_missing_clock_features_never_reach_model_imputation(self):
        clock=Mock();clock.feature_names=['gene_a','gene_b']
        data=pd.DataFrame({'gene_a':[float('nan')],'gene_b':[float('nan')],'unrelated_gene':[1.]})
        with self.assertRaises(ValueError):module.predict_tAge(data,clock,'fixture_scaledlog_EN')
        clock.predict.assert_not_called()

    def test_reordered_samples_are_joined_by_identity(self):
        first=pd.DataFrame({'sample_id':['GSM_A','GSM_B'],'clock_a':[1.,2.]})
        second=pd.DataFrame({'sample_id':['GSM_B','GSM_A'],'clock_b':[20.,10.]})
        merged=module.merge_clock_outputs([first,second])
        self.assertEqual(merged['clock_b'].tolist(),[10.,20.])
        self.assertEqual(first['sample_id'].tolist(),['GSM_A','GSM_B'])
        self.assertEqual(second['sample_id'].tolist(),['GSM_B','GSM_A'])

    def test_membership_conflicts_and_duplicate_clock_labels_are_rejected(self):
        first=pd.DataFrame({'sample_id':['GSM_A','GSM_B'],'clock_a':[1.,2.]})
        variants=[pd.DataFrame({'sample_id':['GSM_A','GSM_C'],'clock_b':[10.,30.]}),
                  pd.DataFrame({'sample_id':['GSM_A','GSM_A'],'clock_b':[10.,30.]}),
                  first.copy()]
        for second in variants:
            with self.assertRaises(ValueError):module.merge_clock_outputs([first,second])

    def test_missing_predictions_and_identity_only_exports_are_rejected(self):
        for frames in [[],[pd.DataFrame({'sample_id':['GSM_A']})],
                       [pd.DataFrame({'sample_id':['GSM_A'],'clock_a':[float('nan')]})]]:
            with self.assertRaises(ValueError):module.merge_clock_outputs(frames)

    def test_clock_exception_is_retained_without_writing_a_completed_export(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'input/GSE_FIXTURE').mkdir(parents=True);(root/'output').mkdir()
            pd.DataFrame({'gene':[1.,2.]},index=['GSM1','GSM2']).to_csv(root/'input/GSE_FIXTURE/GSE_FIXTURE_scaled_data.csv')
            with patch.object(module,'clocks',[object()]),patch.object(module,'clock_names',['fixture_scaledlog_EN']),patch.object(module,'apply_clocks',side_effect=ValueError('fixture model failure')):
                with self.assertRaises(ValueError):module.process_file_folder('GSE_FIXTURE',str(root/'input'),str(root/'output'),'mouse')
            self.assertFalse((root/'output/GSE_FIXTURE.csv').exists())
            receipts=list((root/'output').glob('*.json'))
            self.assertEqual(len(receipts),1)
            record=json.loads(receipts[0].read_text())
            self.assertEqual(record['status'],'CANNOT VERIFY')
            self.assertEqual(record['phase'],'clock_execution')
            self.assertIn('fixture model failure',record['failures'][0]['error'])


if __name__=='__main__':unittest.main()
