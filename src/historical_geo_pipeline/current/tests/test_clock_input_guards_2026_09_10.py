"""Input/format regression tests; no model loading, clock execution or API calls.

Load the two functions through AST to avoid the historical module-level
Biolearn reference download. Fixtures below test parser behavior, not biology.
"""
import ast
import csv
import gzip
import hashlib
from pathlib import Path
import re
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

SOURCE=Path(__file__).resolve().parents[1]/'4_compute_age.py'
tree=ast.parse(SOURCE.read_text())
namespace={
    'np':np,'pd':pd,'gzip':gzip,'csv':csv,'re':re,'hashlib':hashlib,
    'GEO_TABLE_BEGIN':'!series_matrix_table_begin',
    'GEO_TABLE_END':'!series_matrix_table_end',
    'GEO_SAMPLE_ACCESSION':'!Sample_geo_accession',
}
functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in
           ('require_observed_beta_matrix','_sha256_file','_open_geo_matrix',
            '_scan_geo_matrix','read_gse_matrix')]
exec(compile(ast.Module(body=functions,type_ignores=[]),str(SOURCE),'exec'),namespace)
guard=namespace['require_observed_beta_matrix']
read=namespace['read_gse_matrix']


class InputGuards(unittest.TestCase):
    def test_invalid_inputs_are_rejected_before_population_imputation(self):
        bad=[pd.DataFrame(columns=['GSM_A']),pd.DataFrame([[np.nan]],index=['cg1']),
             pd.DataFrame([[.5,.6]],columns=['GSM_A','GSM_A']),
             pd.DataFrame([[.2],[.3]],index=['cg1','cg1']),
             pd.DataFrame([[1.01]],index=['cg1']),pd.DataFrame([[np.inf]],index=['cg1'])]
        for frame in bad:
            with self.subTest(shape=frame.shape),self.assertRaises(ValueError):guard(frame)

    def test_observed_input_and_permitted_missing_cells_are_preserved(self):
        frame=pd.DataFrame([[.2,np.nan],[.3,.4]],index=['cg1','cg2'],columns=['GSM_A','GSM_B'])
        original=frame.copy(deep=True)
        self.assertIs(guard(frame),frame)
        pd.testing.assert_frame_equal(frame,original)

    def test_source_header_and_terminators_cannot_be_silently_repaired(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'matrix.txt'
            good='!Series_title = software fixture\n!series_matrix_table_begin\n"ID_REF"\t"GSM1"\n"cg1"\t0.2\n!series_matrix_table_end\n'
            path.write_text(good)
            parsed=guard(read(str(path)))
            self.assertEqual(parsed.shape,(1,1))
            self.assertEqual(float(parsed.iloc[0,0]),.2)
            for text in [good.replace('!series_matrix_table_begin\n',''),
                         good.replace('!series_matrix_table_end\n',''),
                         good.replace('"ID_REF"\t"GSM1"','"ID_REF"\t"GSM1"\t"GSM1"'),
                         good.replace('"ID_REF"','"feature"'),
                         good.replace('"GSM1"','"sample-one"'),
                         good.replace('"GSM1"','"GSM１"'),
                         good.replace('"cg1"\t0.2','"cg1"'),
                         good.replace('"cg1"','"   "'),
                         good.replace('"cg1"','"!series_matrix_table_begin"'),
                         good.replace('!series_matrix_table_end\n','!series_matrix_table_end\textra\n'),
                         good.replace('!series_matrix_table_end\n','!series_matrix_table_end\nextra\n')]:
                path.write_text(text)
                with self.assertRaises(ValueError):read(str(path))

    def test_logical_record_scan_preserves_multiline_metadata_and_round_trip_float(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'matrix.txt'
            token='0.18292417962450002'
            path.write_text(
                '!Series_title\t"first physical line\nsecond physical line"\n'
                '!Sample_geo_accession\t"GSM1"\t"GSM2"\n'
                '!series_matrix_table_begin\n'
                '"ID_REF"\t"GSM1"\t"GSM2"\n'
                f'"cg1"\t{token}\t0.4\n'
                '!series_matrix_table_end\n'
            )
            parsed=read(path)
            self.assertEqual(parsed.columns.tolist(),['GSM1','GSM2'])
            self.assertEqual(parsed.index.tolist(),['cg1'])
            self.assertEqual(parsed.loc['cg1','GSM1'],float(token))

            path.write_text(path.read_text().replace('"GSM1"\t"GSM2"\n!series',
                                                     '"GSM2"\t"GSM1"\n!series',1))
            with self.assertRaisesRegex(ValueError,'metadata differs'):
                read(path)

    def test_same_width_source_mutation_between_scan_and_parse_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'matrix.txt'
            path.write_text(
                '!series_matrix_table_begin\n'
                '"ID_REF"\t"GSM1"\n'
                '"cg1"\t0.2\n'
                '!series_matrix_table_end\n'
            )
            original_read_csv=pd.read_csv
            def mutate_then_read(*args,**kwargs):
                path.write_text(path.read_text().replace('0.2','0.3'))
                return original_read_csv(*args,**kwargs)
            with patch.object(pd,'read_csv',side_effect=mutate_then_read):
                with self.assertRaisesRegex(ValueError,'bytes changed during parsing'):
                    read(path)

    def test_empty_geo_table_remains_empty_and_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'empty.txt.gz'
            with gzip.open(path,'wt') as handle:
                handle.write('!series_matrix_table_begin\n"ID_REF"\t"GSM1"\n!series_matrix_table_end\n')
            frame=read(str(path))
            self.assertEqual(frame.shape,(0,1))
            with self.assertRaises(ValueError):guard(frame)


if __name__=='__main__':unittest.main()
