"""Software regression fixtures, not empirical audit-accuracy observations."""
import json
import copy
import csv
import gzip
import os
import shutil
import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

from langchain_core.messages import AIMessage, HumanMessage
from src.coding_agent import CodingAgent
from src.pipeline_runtime_2026_09_10 import BoundedOpenAI, IsolatedExecutor, StageStore, PipelineIntegrityError
from src.pipeline_inputs_2026_09_10 import read_matrix, parse_samples, acquire_study_context, acquire_metadata, fetch
from src.pipeline_inputs_2026_09_10 import propose_dnam_recovery, reconstruct_dnam_input, acquire_dnam, preprocess_dnam
from audit_system.source_packet import read_study_context, verify_source_design
from audit_system.adapters.measurement_pilot import verify_clock_measurement
from audit_system.agent_workflow_2026_09_10 import sha, write_new
from src.pipeline_2026_09_10 import (numerical_audit, recompute_declared_contrast, report, parser, build_client,
                                   read_primary_outputs, PRIMARY_SCHEMA, adjust_dependent_family,
                                   CodingActionValidationClient, CodingOutputCompletionGuard, primary_analysis)
from audit_system.agent_workflow_2026_09_10 import DESIGN_SET_SCHEMA, parse_design
from audit_system.report_2026_09_10 import check_measurement_profiles, extract_clock_sites
from audit_system.report_2026_09_10 import recover_dnam_input, verify_dnam_recovery
from audit_system.evaluation_2026_09_10 import (compare_development_case, run_comparison_plan,
                                             replay_development_comparison, build_families, COMPARISON_METHODS)
import pandas as pd
from audit_system.tests import require_development_artifact


class PipelineTests(unittest.TestCase):
    def test_geo_metadata_view_selects_url_and_preserves_sample_fields(self):
        series_payload = ('^SERIES = GSE1\n!Series_status = Public on Jan 01 2020\n'
                          '!Series_sample_id = GSM1\n!Series_sample_id = GSM2\n')
        sample_payload = ''.join(
            '^SAMPLE = ' + sid + '\n!Sample_title = Fixture ' + sid +
            '\n!Sample_organism_ch1 = Homo sapiens\n!Sample_molecule_ch1 = genomic DNA\n'
            '!Sample_characteristics_ch1 = tissue: blood\n'
            for sid in ('GSM1', 'GSM2'))
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            calls = []

            def fake_fetch(url, destination, limit=20_000_000):
                calls.append(url)
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text(series_payload if 'targ=self' in url else sample_payload)
                return destination

            with patch('src.pipeline_inputs_2026_09_10.fetch', side_effect=fake_fetch), \
                    patch('src.pipeline_inputs_2026_09_10.time.sleep'):
                record = acquire_metadata('GSE1', root / 'brief', publication_methods=False,
                                          metadata_view='brief')
            self.assertIn('targ=gsm', calls[1]); self.assertIn('view=brief', calls[1])
            self.assertEqual(record['sample_metadata_view'], 'brief')
            metadata = (root / 'brief' / 'metadata.csv').read_text()
            self.assertIn('GSM1', metadata); self.assertIn('Homo sapiens', metadata)
            self.assertEqual(set(record['sample_ids']), {'GSM1', 'GSM2'})

            calls.clear()
            with patch('src.pipeline_inputs_2026_09_10.fetch', side_effect=fake_fetch), \
                    patch('src.pipeline_inputs_2026_09_10.time.sleep'):
                legacy = acquire_metadata('GSE1', root / 'full', publication_methods=False)
            self.assertIn('view=full', calls[1])
            self.assertEqual(legacy['sample_metadata_view'], 'full')

    def test_geo_metadata_view_rejects_unknown_value_before_network(self):
        with tempfile.TemporaryDirectory() as folder, patch('src.pipeline_inputs_2026_09_10.fetch') as fetch_call:
            with self.assertRaisesRegex(ValueError, 'brief, full'):
                acquire_metadata('GSE1', Path(folder), metadata_view='quick')
            fetch_call.assert_not_called()

    def test_geo_metadata_view_is_in_cli_and_binds_resume_config(self):
        required = ['--study', 'GSE1', '--modality', 'human_dnam', '--output', 'out']
        self.assertEqual(parser().parse_args(required).geo_metadata_view, 'full')
        self.assertEqual(parser().parse_args(required + ['--geo-metadata-view', 'brief']).geo_metadata_view, 'brief')
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            StageStore(root, {'geo_metadata_view': 'brief'}).run(
                'metadata', lambda directory: {'sample_metadata_view': 'brief'})
            with self.assertRaisesRegex(PipelineIntegrityError, 'Configuration'):
                StageStore(root, {'geo_metadata_view': 'full'})

    def test_prompt_cap_has_stable_default_explicit_value_and_frozen_resume_binding(self):
        required = ['--study', 'GSE1', '--modality', 'human_dnam', '--output', 'out']
        self.assertEqual(parser().parse_args(required).max_prompt_bytes, 48000)
        self.assertEqual(parser().parse_args(required + ['--max-prompt-bytes', '96000']).max_prompt_bytes, 96000)
        for value in ('0', '-1', '1.5', 'not-an-integer'):
            with self.subTest(value=value), self.assertRaises(SystemExit):
                parser().parse_args(required + ['--max-prompt-bytes', value])
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            StageStore(root, {'max_prompt_bytes': 96000}).run(
                'fixture', lambda directory: {'status': 'PASS', 'scope': 'software fixture'})
            with self.assertRaisesRegex(PipelineIntegrityError, 'Configuration'):
                StageStore(root, {'max_prompt_bytes': 48000})

    def test_prompt_cap_propagates_to_both_actor_clients(self):
        required = ['--study', 'GSE1', '--modality', 'human_dnam', '--output', 'out',
                    '--audit-model', 'gpt-5.6-luna', '--audit-input-usd-per-million', '0.5',
                    '--audit-output-usd-per-million', '2.4', '--input-usd-per-million', '5',
                    '--output-usd-per-million', '20', '--max-prompt-bytes', '96000']
        args = parser().parse_args(required)
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            design = build_client(args, output, 'design')
            primary = build_client(args, output, 'primary')
            self.assertEqual(design.prompt_bytes, 96000)
            self.assertEqual(primary.prompt_bytes, 96000)
            self.assertEqual(design.model, 'gpt-5.6-luna')
            self.assertEqual(primary.model, 'gpt-5.6-sol')

    def test_prompt_cap_allows_bounded_feedback_retry_without_truncation(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            body = {'model': 'gpt-5.6-sol',
                    'usage': {'prompt_tokens': 300, 'completion_tokens': 100},
                    'choices': [{'finish_reason': 'stop',
                                 'message': {'content': 'software fixture response'}}]}
            transport_calls = []

            def transport(*args, **kwargs):
                transport_calls.append((args, kwargs))
                return SimpleNamespace(returncode=0,
                    stdout=(json.dumps(body) + '\n200').encode())

            first = [HumanMessage(content='initial software fixture ' + ('a' * 20000))]
            retry = first + [AIMessage(content='model fixture response ' + ('b' * 20000)),
                             HumanMessage(content='formatting feedback ' + ('c' * 10000))]
            with patch.dict(os.environ, {'OPENAI_API_KEY': 'software-fixture-key'}), \
                    patch('src.pipeline_runtime_2026_09_10.subprocess.run', side_effect=transport):
                client = BoundedOpenAI(root / 'api', 'fixture', model='gpt-5.6-sol',
                                       max_calls=2, max_cost=1.0, output_tokens=100,
                                       prompt_bytes=96000)
                client.invoke(first)
                client.invoke(retry)
                with self.assertRaisesRegex(RuntimeError, 'call/cost limit'):
                    client.invoke(retry)
            self.assertEqual(len(transport_calls), 2)
            saved = json.loads((root / 'api' / 'fixture-001' / 'request.json').read_text())
            retry_bytes = len(json.dumps(saved, ensure_ascii=False).encode())
            self.assertGreater(retry_bytes, 48000)
            self.assertLessEqual(retry_bytes, 96000)
            self.assertEqual(saved['messages'][-2]['content'], retry[-2].content)
            self.assertEqual(saved['messages'][-1]['content'], retry[-1].content)
            self.assertEqual(saved['messages'], [
                {'role': 'user', 'content': retry[0].content},
                {'role': 'assistant', 'content': retry[1].content},
                {'role': 'user', 'content': retry[2].content}])

            low_root = root / 'low-cap'
            low_calls = []

            def low_transport(*args, **kwargs):
                low_calls.append((args, kwargs))
                return SimpleNamespace(returncode=0,
                    stdout=(json.dumps(body) + '\n200').encode())

            with patch.dict(os.environ, {'OPENAI_API_KEY': 'software-fixture-key'}), \
                    patch('src.pipeline_runtime_2026_09_10.subprocess.run', side_effect=low_transport):
                low_client = BoundedOpenAI(low_root / 'api', 'fixture', model='gpt-5.6-sol',
                                           max_calls=2, max_cost=1.0, output_tokens=100,
                                           prompt_bytes=48000)
                low_client.invoke(first)
                with self.assertRaisesRegex(ValueError, r'bytes > permitted 48000'):
                    low_client.invoke(retry)
            self.assertEqual(len(low_calls), 1)
            self.assertEqual(len(list((low_root / 'api').glob('*/reservation.json'))), 1)

    def test_prompt_cap_constructor_rejects_nonpositive_or_non_exact_values_before_output(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for value in (0, -1, 1.5, True, '96000'):
                with self.subTest(value=value):
                    destination = root / ('invalid-' + str(value).replace('.', '_'))
                    with self.assertRaisesRegex(ValueError, 'positive integer'):
                        BoundedOpenAI(destination, 'fixture', prompt_bytes=value)
                    self.assertFalse(destination.exists())

    def test_raw_filter_precedes_aggregation_only_in_explicit_profile(self):
        from src.pipeline_inputs_2026_09_10 import mapped_rna_counts
        # Software fixture: individually subthreshold rows must not enter the
        # published-order library through an above-threshold aggregate.
        raw = pd.DataFrame({'GSM1': [6, 6, 10, 20], 'GSM2': [6, 6, 12, 30]},
                           index=['E1', 'E2', 'E3', 'E4'])
        table = pd.DataFrame({'Entrez': [1, 1, 2, 3]}, index=raw.index)
        legacy = mapped_rna_counts(raw, table)
        published = mapped_rna_counts(raw, table, raw_filter_first=True)
        self.assertEqual(legacy.loc['1'].tolist(), [12, 12])
        self.assertEqual(published.index.tolist(), ['2', '3'])
        self.assertEqual(published.to_dict(), {'GSM1': {'2': 10, '3': 20}, 'GSM2': {'2': 12, '3': 30}})
        table['Entrez'] = [1.5, 1, 2, 3]
        with self.assertRaisesRegex(ValueError, 'positive finite integers'):
            mapped_rna_counts(raw, table, raw_filter_first=True)

    def test_mouse_preprocessing_stratum_binds_every_source_sample(self):
        from src.pipeline_inputs_2026_09_10 import rna_preprocessing_stratum
        packet = {'samples': [{'sample_id': sid, 'organism_ch1': 'Mus musculus',
            'molecule_ch1': 'total RNA', 'characteristics_ch1': 'tissue: adipose tissue'} for sid in ('GSM1', 'GSM2')]}
        result = rna_preprocessing_stratum(packet, ['GSM2', 'GSM1'])
        self.assertEqual(result['reference_ids'], ['GSM2', 'GSM1'])
        self.assertFalse(result['matched_controls_verified'])
        self.assertEqual(result['biological_applicability'], 'CANNOT VERIFY')
        with self.assertRaisesRegex(ValueError, 'membership'):
            rna_preprocessing_stratum(packet, ['GSM1'])
        packet['samples'][1]['characteristics_ch1'] = 'tissue: liver'
        with self.assertRaisesRegex(ValueError, 'distinct'):
            rna_preprocessing_stratum(packet, ['GSM1', 'GSM2'])
        packet['samples'][1]['characteristics_ch1'] = 'genotype: wild type'
        with self.assertRaisesRegex(ValueError, 'missing or ambiguous'):
            rna_preprocessing_stratum(packet, ['GSM1', 'GSM2'])

    def test_mouse_preprocessing_rejects_other_species_before_file_operations(self):
        from src.pipeline_inputs_2026_09_10 import preprocess_rna
        with self.assertRaisesRegex(ValueError, 'mouse only'):
            preprocess_rna(Path('absent-source'), 'human', Path('absent-output'), normalization='published_mouse_rle')

    def test_rle_factor_calculation_rejects_undefined_reference_and_preserves_library_scaling(self):
        from src.pipeline_inputs_2026_09_10 import rna_rle_factors
        import numpy as np
        # Explicit software fixture; no empirical performance inference.
        counts = pd.DataFrame({'GSM1': [10, 20, 0], 'GSM2': [20, 40, 7]}, index=['g1', 'g2', 'g3'])
        factors = rna_rle_factors(counts)
        scaled = counts.mul(pd.Series({'GSM1': 3, 'GSM2': 7}), axis=1)
        other = rna_rle_factors(scaled)
        self.assertTrue(np.allclose(counts.mul(factors['count_multiplier'], axis=1),
                                    scaled.mul(other['count_multiplier'], axis=1)))
        counts.iloc[1] = [20, 0]; counts.iloc[0] = [0, 40]
        with self.assertRaisesRegex(ValueError, 'observed in every sample'):
            rna_rle_factors(counts)

    def test_normalization_sensitivity_does_not_release_held_out_execution(self):
        from audit_system.report_2026_09_10 import compare_rna_normalization
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); plan = root / 'plan.json'
            plan.write_text(json.dumps({'schema': 'clockbase-rna-normalization-sensitivity-v1',
                                        'split': 'held_out', 'study_id': 'GSE1', 'species': 'mouse'}))
            with self.assertRaisesRegex(ValueError, 'exposed development'):
                compare_rna_normalization(plan, root / 'output')
            self.assertFalse((root / 'output').exists())

    def test_geo_parser_keeps_repeated_extraction_protocols_for_source_assessment(self):
        payload = ('^SAMPLE = GSM1\n!Sample_title = software fixture\n'
                   '!Sample_extract_protocol_ch1 = Culture and harvest protocol\n'
                   '!Sample_extract_protocol_ch1 = Library preparation protocol\n'
                   '!Sample_extract_protocol_ch1 = Read alignment protocol\n'
                   '!Sample_contact_email = private@example.test\n')
        rows = parse_samples(payload, 'GSE1')
        self.assertEqual([rows[0][key] for key in ('extract_protocol_ch1', 'extract_protocol_ch1.1', 'extract_protocol_ch1.2')],
                         ['Culture and harvest protocol', 'Library preparation protocol', 'Read alignment protocol'])
        self.assertNotIn('private@example.test', json.dumps(rows))

    def test_geo_method_evidence_survives_acquisition_and_both_packet_adapters(self):
        from audit_system.agent_workflow_2026_09_10 import source_packet
        from audit_system.source_packet import metadata_packet, compact_source_packet, expand_source_packet
        # Software fixture: assay values are narrative strings, not measurements.
        expected = {'label_ch1': 'Cy3', 'label_protocol_ch1': 'Antibody label',
                    'label_protocol_ch1.1': 'Wash after labeling', 'hyb_protocol': 'Gene-array hybridization',
                    'scan_protocol': 'Scanner configuration', 'scan_protocol.1': 'Scanner gain',
                    'data_processing': 'Normalized fold fluorescence',
                    'data_processing.1': 'NARRATIVE_FINDING_VISIBLE_WITH_SCOPE_WARNING'}
        payload = '^SAMPLE = GSM1\n!Sample_title = software fixture\n'
        for field, value in expected.items():
            payload += '!Sample_' + field.split('.')[0] + ' = ' + value + '\n'
        payload += '!Sample_contact_email = PRIVATE_SENTINEL\n!Sample_p_value = OUTCOME_SENTINEL\n'
        rows = parse_samples(payload, 'GSE1')
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'metadata.csv'
            with path.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
            current = source_packet(path, 'GSE1', 'Inspect assay suitability')
            legacy = metadata_packet(path, 'GSE1')
        for field, value in expected.items():
            self.assertEqual(rows[0][field], value)
            self.assertEqual(current['samples'][0][field], value)
            self.assertEqual(legacy['samples'][0]['raw.' + field], value)
        self.assertEqual(expand_source_packet(compact_source_packet(current)), current)
        for packet in (current, legacy):
            self.assertIn('published findings', packet['method_provenance']['blinding_scope'])
            self.assertNotIn('PRIVATE_SENTINEL', json.dumps(packet))
            self.assertNotIn('OUTCOME_SENTINEL', json.dumps(packet))

    def recovery_fixture(self, root):
        """Tiny source-recovery software fixture; not empirical GEO measurements."""
        study='GSE999001';ids=['GSM1','GSM2','GSM3','GSM4'];tokens=['P1','P2','X03','3A']
        original=root/'empty.csv.gz'
        with gzip.open(original,'wt') as handle:handle.write('ID_REF,'+','.join(ids)+'\n')
        supplement=root/'software-fixture-beta.tsv.gz'
        with gzip.open(supplement,'wt') as handle:
            handle.write('ID\t'+'\t'.join(tokens)+'\ncg001\t0.1\t0.2\t0.3\t0.4\ncg002\t0.5\t0.6\t0.7\t0.8\n')
        metadata=root/'metadata.csv'
        with metadata.open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=['GSM','gse_id','title']);writer.writeheader()
            writer.writerows({'GSM':sid,'gse_id':study,'title':'Software fixture ['+token+']'} for sid,token in zip(ids,tokens))
        series=root/'series.soft'
        series.write_text('^SERIES = '+study+'\n!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE999nnn/'+study+'/suppl/'+supplement.name+'\n'
            + ''.join('^SAMPLE = '+sid+'\n!Sample_title = Software fixture ['+token+']\n' for sid,token in zip(ids,tokens)))
        return study,original,supplement,metadata,series

    def xlsx_recovery_fixture(self, root):
        """Spreadsheet software fixture, never empirical GEO measurements."""
        from openpyxl import Workbook
        study, original, old_supplement, metadata, series = self.recovery_fixture(root)
        supplement = root / 'software-fixture-beta.xlsx'
        workbook = Workbook(); sheet = workbook.active; sheet.title = 'Measured beta'
        sheet.append(['ID_REF', 'P1', 'P2', 'X03', '3A'])
        sheet.append(['cg001', .1, .2, .3, .4])
        sheet.append(['00123', .5, None, .7, .8])
        workbook.save(supplement); workbook.close()
        series.write_text(series.read_text().replace(old_supplement.name, supplement.name))
        return study, original, supplement, metadata, series

    def test_xlsx_recovery_preserves_sample_order_features_missingness_and_resumes(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); inputs = self.xlsx_recovery_fixture(root)
            result = recover_dnam_input(*inputs, root/'out')
            self.assertEqual((result['status'], result['source_rows'], result['source_samples'],
                              result['values_compared'], result['missing_values']), ('PASS', 2, 4, 8, 1))
            self.assertIn('openpyxl', result['parser'])
            self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
            self.assertFalse(result['clock_execution_performed'])
            with gzip.open(result['recovered_input'], 'rt') as source:
                rows = list(csv.reader(source, delimiter='\t'))
            self.assertEqual(rows[0], ['ID_REF', 'GSM1', 'GSM2', 'GSM3', 'GSM4'])
            self.assertEqual(rows[2][0], '00123'); self.assertEqual(rows[2][2], '')
            before = {p: sha(p) for p in (root/'out').rglob('*') if p.is_file()}
            with patch('src.pipeline_inputs_2026_09_10.reconstruct_dnam_input', side_effect=AssertionError('Repeated reconstruction')):
                self.assertEqual(recover_dnam_input(*inputs, root/'out'), result)
            self.assertEqual(before, {p: sha(p) for p in before})

    def test_xlsx_recovery_rejects_aggregate_headers_extra_sheets_formulas_and_nonbeta(self):
        from openpyxl import load_workbook
        cases = ('aggregate', 'duplicate_detection_p', 'extra_sheet', 'formula', 'boolean',
                 'error', 'outside_range', 'duplicate_feature', 'missing_feature', 'numeric_feature')
        for case in cases:
            with self.subTest(case=case), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); inputs = self.xlsx_recovery_fixture(root)
                workbook = load_workbook(inputs[2]); sheet = workbook.active
                if case == 'aggregate':
                    for col, value in enumerate(['ID_REF', 'Control', 'Treated', 'Control2', 'Treated2'], 1):
                        sheet.cell(1, col, value)
                elif case == 'duplicate_detection_p':
                    for col, value in enumerate(['ID_REF', 'Control-siRNA', 'Detection Pval', 'PHGDH-siRNA', 'Detection Pval'], 1):
                        sheet.cell(1, col, value)
                elif case == 'extra_sheet': workbook.create_sheet('Ambiguous second table')
                elif case == 'formula': sheet['B3'] = '=0.5'
                elif case == 'boolean': sheet['B3'] = True
                elif case == 'error': sheet['B3'] = '#DIV/0!'
                elif case == 'outside_range': sheet['B3'] = 1.2
                elif case == 'duplicate_feature': sheet['A3'] = 'cg001'
                elif case == 'missing_feature': sheet['A3'] = None
                elif case == 'numeric_feature': sheet['A3'] = 123
                workbook.save(inputs[2]); workbook.close()
                with self.assertRaises(ValueError):
                    recover_dnam_input(*inputs, root/'out')
                self.assertFalse((root/'out/02_independent_verification/complete.json').exists())

    def test_xlsx_merged_measurements_cannot_become_values_plus_missingness(self):
        from openpyxl import load_workbook
        for cells in ('B2:C2', 'B2:B3'):
            with self.subTest(merge=cells), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); inputs = self.xlsx_recovery_fixture(root)
                # Build a valid sparse projection first, then show that even
                # a resealed source merge cannot be certified independently.
                workbook = load_workbook(inputs[2]); sheet = workbook.active
                if cells == 'B2:C2': sheet['C2'] = None
                else: sheet['B3'] = None
                sheet.append(['cg003', .4, .3, .2, .1])
                workbook.save(inputs[2]); workbook.close()
                proposal = propose_dnam_recovery(*inputs); output = root/'out'; output.mkdir()
                result = reconstruct_dnam_input(proposal, output)
                self.assertEqual(verify_dnam_recovery(proposal, result)['status'], 'PASS')
                workbook = load_workbook(inputs[2]); workbook.active.merge_cells(cells)
                workbook.save(inputs[2]); workbook.close()
                # Same visible cell values; the layout now makes the original
                # individual-cell interpretation unsupported.
                proposal['source_files']['supplement']['sha256'] = sha(inputs[2])
                with self.assertRaisesRegex(ValueError, 'merged cells'):
                    verify_dnam_recovery(proposal, result)
                fresh = root/'fresh'; fresh.mkdir()
                with self.assertRaisesRegex(ValueError, 'merged cells'):
                    reconstruct_dnam_input(proposal, fresh)

    def test_xlsx_dimensions_do_not_truncate_rows_or_sample_columns(self):
        import zipfile
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); inputs = self.xlsx_recovery_fixture(root)
            with zipfile.ZipFile(inputs[2]) as archive:
                contents = {name: archive.read(name) for name in archive.namelist()}
            self.assertIn(b'<dimension ref="A1:E3"', contents['xl/worksheets/sheet1.xml'])
            contents['xl/worksheets/sheet1.xml'] = contents['xl/worksheets/sheet1.xml'].replace(
                b'<dimension ref="A1:E3"', b'<dimension ref="A1:A1"')
            with zipfile.ZipFile(inputs[2], 'w', zipfile.ZIP_DEFLATED) as archive:
                for name, payload in contents.items(): archive.writestr(name, payload)
            result = recover_dnam_input(*inputs, root/'out')
            self.assertEqual((result['source_rows'], result['values_compared']), (2, 8))

    def test_xlsx_coordinate_ambiguity_cannot_create_or_duplicate_measurements(self):
        import zipfile
        changes = ((b'r="B3"', b'r="C3"'), (b'r="B3"', b'r="B4"'),
                   (b'<row r="3">', b'<row r="4">'), (b'r="E3"', b'r="F3"'))
        for before, after in changes:
            with self.subTest(after=after), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); inputs = self.xlsx_recovery_fixture(root)
                with zipfile.ZipFile(inputs[2]) as archive:
                    contents = {name: archive.read(name) for name in archive.namelist()}
                xml = contents['xl/worksheets/sheet1.xml']; self.assertIn(before, xml)
                # The first change fills an otherwise missing cell coordinate;
                # make a duplicate by moving D3 there as well.
                xml = xml.replace(before, after)
                if before == b'r="B3"' and after == b'r="C3"':
                    xml = xml.replace(b'r="D3"', b'r="C3"')
                contents['xl/worksheets/sheet1.xml'] = xml
                with zipfile.ZipFile(inputs[2], 'w', zipfile.ZIP_DEFLATED) as archive:
                    for name, payload in contents.items(): archive.writestr(name, payload)
                with self.assertRaises(ValueError): recover_dnam_input(*inputs, root/'out')

    def test_xlsx_recovery_rejects_resealed_value_tampering(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); inputs = self.xlsx_recovery_fixture(root)
            proposal = propose_dnam_recovery(*inputs); output = root/'out'; output.mkdir()
            result = reconstruct_dnam_input(proposal, output)
            target = Path(result['path'])
            with gzip.open(target, 'rt') as handle: rows = list(csv.reader(handle, delimiter='\t'))
            rows[1][1] = '0.9'
            with gzip.open(target, 'wt') as handle:
                writer = csv.writer(handle, delimiter='\t', lineterminator='\n'); writer.writerows(rows)
            result['sha256'] = sha(target)
            with self.assertRaisesRegex(ValueError, 'values differ'):
                verify_dnam_recovery(proposal, result)

    def test_xlsx_numeric_date_format_is_rejected_by_independent_reader(self):
        from openpyxl import load_workbook
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); inputs = self.xlsx_recovery_fixture(root)
            workbook = load_workbook(inputs[2]); workbook.active['B3'].number_format = 'yyyy-mm-dd'
            workbook.save(inputs[2]); workbook.close()
            proposal = propose_dnam_recovery(*inputs); output = root/'out'; output.mkdir()
            result = reconstruct_dnam_input(proposal, output)
            with self.assertRaisesRegex(ValueError, 'date|not numeric'):
                verify_dnam_recovery(proposal, result)

    def test_recovery_accepts_source_bound_space_table_and_idat_coordinate_aliases(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); study = 'GSE999002'
            ids = ['GSM9001', 'GSM9002']
            original = root / 'empty.csv.gz'
            with gzip.open(original, 'wt') as handle:
                handle.write('ID_REF,' + ','.join(ids) + '\n')
            supplement = root / 'Betas.txt.gz'
            with gzip.open(supplement, 'wt') as handle:
                handle.write('"ARRAY_R05C01" "ARRAY_R06C01"\n')
                handle.write('"cg0001" 0.1 0.2\n"cg0002" 0.3 0.4\n')
            metadata = root / 'metadata.csv'
            with metadata.open('w') as handle:
                writer = csv.DictWriter(handle, fieldnames=['GSM', 'gse_id', 'title'])
                writer.writeheader()
                writer.writerows({'GSM': sid, 'gse_id': study, 'title': 'Sample ' + sid} for sid in ids)
            series = root / 'series.soft'
            series.write_text(
                '^SERIES = ' + study + '\n'
                '!Series_sample_id = GSM9001\n!Series_sample_id = GSM9002\n'
                '!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE999nnn/'
                + study + '/suppl/' + supplement.name + '\n')
            samples = root / 'samples.soft'
            samples.write_text(
                '^SAMPLE = GSM9001\n!Sample_title = Sample GSM9001\n'
                '!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9nnn/'
                'GSM9001/suppl/GSM9001_ARRAY_R05C01_Grn.idat.gz\n'
                '!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9nnn/'
                'GSM9001/suppl/GSM9001_ARRAY_R05C01_Red.idat.gz\n'
                '^SAMPLE = GSM9002\n!Sample_title = Sample GSM9002\n'
                '!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9nnn/'
                'GSM9002/suppl/GSM9002_ARRAY_R06C01_Grn.idat.gz\n'
                '!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9nnn/'
                'GSM9002/suppl/GSM9002_ARRAY_R06C01_Red.idat.gz\n')
            inputs = (study, original, supplement, metadata, series)
            proposal = propose_dnam_recovery(*inputs, sample_source=samples)
            self.assertEqual(proposal['source_layout']['delimiter'], 'space')
            self.assertFalse(proposal['source_layout']['feature_column'])
            self.assertEqual([row['gsm_id'] for row in proposal['crosswalk']], ids)
            self.assertTrue(all(row['mapping_rule'] == 'exact_sample_supplementary_idat_coordinate'
                                for row in proposal['crosswalk']))
            output = root / 'out'; output.mkdir()
            result = reconstruct_dnam_input(proposal, output)
            self.assertEqual(result['feature_count'], 2)
            self.assertEqual(result['finite_values_per_sample'], {'GSM9001': 2.0, 'GSM9002': 2.0})
            self.assertEqual(verify_dnam_recovery(proposal, result)['values_compared'], 4)
            with gzip.open(result['path'], 'rt') as handle:
                self.assertEqual(handle.read(), 'ID_REF\tGSM9001\tGSM9002\ncg0001\t0.1\t0.2\ncg0002\t0.3\t0.4\n')
            recovered = recover_dnam_input(*inputs, root / 'run', sample_source=samples)
            self.assertEqual(recovered['status'], 'PASS')
            direct = root / 'direct'; direct.mkdir()
            acquired = acquire_dnam(study, ids, original, direct, metadata=metadata, series=series,
                                    supplement=supplement, sample_source=samples)
            self.assertEqual(acquired['repair_status'], 'PASS')
            forged = copy.deepcopy(proposal)
            forged['crosswalk'][0]['quote'] = forged['crosswalk'][0]['quote'].replace(
                'ARRAY_R05C01_Grn.idat.gz', 'R05C01_Grn.idat.gz')
            with self.assertRaisesRegex(ValueError, 'exact source binding'):
                verify_dnam_recovery(forged, result)

    def test_recovery_rejects_metadata_title_swaps_even_when_local_crosswalk_is_consistent(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);inputs=self.recovery_fixture(root)
            proposal=propose_dnam_recovery(*inputs)
            inputs[3].write_text(inputs[3].read_text().replace('[P1]','[TMP]').replace('[P2]','[P1]').replace('[TMP]','[P2]'))
            with self.assertRaisesRegex(ValueError,'title differs from primary'):
                propose_dnam_recovery(*inputs)
            # A forged local proposal cannot bypass the verifier's source check.
            proposal['source_files']['metadata']['sha256']=sha(inputs[3])
            proposal['crosswalk'][0]['gsm_id']='GSM2'
            proposal['crosswalk'][1]['gsm_id']='GSM1'
            output=root/'out';output.mkdir()
            result=reconstruct_dnam_input(proposal,output)
            with self.assertRaisesRegex(ValueError,'title differs from primary'):
                verify_dnam_recovery(proposal,result)

    def test_recovery_requires_primary_sample_records_and_consistent_membership(self):
        for mutation in ('missing','duplicate','membership','conflicting'):
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as folder:
                root=Path(folder);inputs=self.recovery_fixture(root);text=inputs[4].read_text()
                if mutation=='missing':text=text.split('^SAMPLE')[0]
                if mutation=='duplicate':text=text.replace('^SAMPLE = GSM2','^SAMPLE = GSM1')
                if mutation=='membership':text+='!Series_sample_id = GSM1\n'
                if mutation=='conflicting':text+='!Sample_geo_accession = GSM9\n'
                inputs[4].write_text(text)
                with self.assertRaises(ValueError):propose_dnam_recovery(*inputs)

    def test_recovery_rejects_primary_alias_ambiguity_even_when_local_metadata_omits_it(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);inputs=self.recovery_fixture(root)
            proposal=propose_dnam_recovery(*inputs)
            text=inputs[4].read_text().replace('Software fixture [P1]\n',
                'Software fixture [P1]\n!Sample_characteristics_ch1 = sample code: P2\n')
            text=text.replace('Software fixture [P2]\n',
                'Software fixture [P2]\n!Sample_characteristics_ch1 = sample code: P1\n')
            inputs[4].write_text(text)
            with self.assertRaisesRegex(ValueError,'ambiguous in primary'):propose_dnam_recovery(*inputs)
            proposal['source_files']['series']['sha256']=sha(inputs[4])
            output=root/'out';output.mkdir();result=reconstruct_dnam_input(proposal,output)
            with self.assertRaisesRegex(ValueError,'ambiguous in primary'):verify_dnam_recovery(proposal,result)

    def test_recovery_soft_text_with_embedded_tabs_is_not_a_matrix_header(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);inputs=self.recovery_fixture(root)
            with inputs[4].open('a') as handle:handle.write('!Sample_description = Software fixture\twith tab\n')
            self.assertEqual(propose_dnam_recovery(*inputs)['source_identity_verification']['samples'],4)

    def test_recovery_accepts_separate_series_and_sample_sources_and_pins_both(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);inputs=self.recovery_fixture(root)
            series_text,samples_text=inputs[4].read_text().split('^SAMPLE',1)
            samples=root/'samples.soft';samples.write_text('^SAMPLE'+samples_text)
            inputs[4].write_text(series_text+'!Series_sample_id = GSM1 GSM2 GSM3 GSM4\n')
            with self.assertRaisesRegex(ValueError,'sample identity evidence is missing'):
                propose_dnam_recovery(*inputs)
            result=recover_dnam_input(*inputs,root/'run',sample_source=samples)
            self.assertEqual(result['status'],'PASS')
            config=json.loads((root/'run/config.json').read_text())
            self.assertEqual(config['proposal']['source_files']['sample_records']['sha256'],sha(samples))
            with patch('src.pipeline_inputs_2026_09_10.reconstruct_dnam_input',side_effect=AssertionError('Repeated repair')):
                self.assertEqual(recover_dnam_input(*inputs,root/'run',sample_source=samples),result)
            output=root/'direct';output.mkdir()
            direct=acquire_dnam(inputs[0],['GSM1','GSM2','GSM3','GSM4'],inputs[1],output,
                metadata=inputs[3],series=inputs[4],supplement=inputs[2],sample_source=samples)
            self.assertEqual(direct['repair_status'],'PASS')
            # Source changes invalidate both the pinned proposal and subsequent resume.
            samples.write_text(samples.read_text().replace('[P1]','[wrong]'))
            with self.assertRaises(ValueError):recover_dnam_input(*inputs,root/'run',sample_source=samples)

    def test_separate_sample_sources_require_exact_series_membership(self):
        for membership in ('','!Series_sample_id = GSM1\n'):
            with self.subTest(membership=membership),tempfile.TemporaryDirectory() as folder:
                root=Path(folder);inputs=self.recovery_fixture(root)
                series_text,samples_text=inputs[4].read_text().split('^SAMPLE',1)
                samples=root/'samples.soft';samples.write_text('^SAMPLE'+samples_text)
                inputs[4].write_text(series_text+membership)
                with self.assertRaises(ValueError):propose_dnam_recovery(*inputs,sample_source=samples)

    def test_recovery_accepts_matrix_header_identities_and_rejects_ragged_source_fields(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);inputs=self.recovery_fixture(root)
            text=inputs[4].read_text().split('^SAMPLE')[0]
            text+='!Series_sample_id\t"GSM1 GSM2 GSM3 GSM4"\n'
            text+='!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"\n'
            text+='!Sample_title\t"Software fixture [P1]"\t"Software fixture [P2]"\t"Software fixture [X03]"\t"Software fixture [3A]"\n'
            inputs[4].write_text(text+'!series_matrix_table_begin\n')
            self.assertEqual(propose_dnam_recovery(*inputs)['source_identity_verification']['status'],'PASS')
            inputs[4].write_text(text.replace('\t"Software fixture [3A]"','')+'!series_matrix_table_begin\n')
            with self.assertRaisesRegex(ValueError,'field widths differ'):propose_dnam_recovery(*inputs)

    def test_recovery_explicit_sample_labels_must_be_present_in_primary_characteristics(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);inputs=self.recovery_fixture(root)
            tokens=['P1','P2','X03','3A']
            with inputs[3].open('w') as handle:
                writer=csv.DictWriter(handle,fieldnames=['GSM','title','characteristics_ch1']);writer.writeheader()
                writer.writerows({'GSM':f'GSM{i}','title':'Software fixture',
                    'characteristics_ch1':'sample code: '+token} for i,token in enumerate(tokens,1))
            text=inputs[4].read_text().split('^SAMPLE')[0]
            text+=''.join(f'^SAMPLE = GSM{i}\n!Sample_title = Software fixture\n!Sample_characteristics_ch1 = sample code: {token}\n'
                for i,token in enumerate(tokens,1))
            inputs[4].write_text(text)
            result=recover_dnam_input(*inputs,root/'run')
            self.assertEqual(result['status'],'PASS')
            inputs[3].write_text(inputs[3].read_text().replace('sample code: P1','sample code: wrong'))
            with self.assertRaisesRegex(ValueError,'sample label differs from primary'):propose_dnam_recovery(*inputs)

    def test_recovery_compares_numeric_feature_identifiers_as_exact_tokens(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);inputs=self.recovery_fixture(root)
            with gzip.open(inputs[2],'rt') as handle:text=handle.read().replace('cg001','001').replace('cg002','002')
            with gzip.open(inputs[2],'wt') as handle:handle.write(text)
            proposal=propose_dnam_recovery(*inputs);output=root/'out';output.mkdir()
            result=reconstruct_dnam_input(proposal,output)
            self.assertEqual(verify_dnam_recovery(proposal,result)['status'],'PASS')
            self.assertEqual(list(read_matrix(result['path'],beta=True).index),['001','002'])
            with gzip.open(result['path'],'rt') as handle:text=handle.read().replace('001\t','1\t')
            with gzip.open(result['path'],'wt') as handle:handle.write(text)
            result['sha256']=sha(result['path'])
            with self.assertRaisesRegex(ValueError,'probe identities/order changed'):verify_dnam_recovery(proposal,result)

    def test_imputation_verification_checks_axes_even_when_all_values_are_equal(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            frame=pd.DataFrame([[.5,.5]],index=['cg001'],columns=['GSM1','GSM2'])
            means=pd.Series([.5],index=frame.index)
            with patch('src.pipeline_inputs_2026_09_10.read_matrix',return_value=frame), \
                 patch('src.pipeline_inputs_2026_09_10.reference_means',return_value=(means,{})), \
                 patch('biolearn.imputation.impute_from_standard',return_value=frame[['GSM2','GSM1']]):
                with self.assertRaisesRegex(ValueError,'sample/feature identities changed'):
                    preprocess_dnam(root/'raw',root/'ref',root/'cache',root)

    def test_source_recovery_reconstructs_every_value_and_resumes_without_rewriting(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);inputs=self.recovery_fixture(root)
            before={str(path):sha(path) for path in inputs[1:]}
            ready=recover_dnam_input(*inputs,root/'run',prepare_only=True)
            self.assertTrue(ready['prepared']);self.assertEqual(ready['source_samples'],4)
            result=recover_dnam_input(*inputs,root/'run')
            self.assertEqual(result['status'],'PASS');self.assertEqual(result['values_compared'],8)
            self.assertEqual(result['scientific_status'],'CANNOT VERIFY');self.assertEqual(result['generation_calls'],0)
            repaired=read_matrix(result['recovered_input'],beta=True)
            self.assertEqual(list(repaired.columns),['GSM1','GSM2','GSM3','GSM4'])
            self.assertEqual(repaired.loc['cg002','GSM4'],.8)
            self.assertEqual(before,{str(path):sha(path) for path in inputs[1:]})
            with patch('src.pipeline_inputs_2026_09_10.reconstruct_dnam_input',side_effect=AssertionError('Completed recovery repeated')):
                self.assertEqual(recover_dnam_input(*inputs,root/'run'),result)

    def test_recovery_is_integrated_before_model_stages_and_checks_requested_membership(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);study,original,supplement,metadata,series=self.recovery_fixture(root)
            output=root/'stage';output.mkdir()
            result=acquire_dnam(study,['GSM1','GSM2','GSM3','GSM4'],original,output,
                                metadata=metadata,series=series,supplement=supplement)
            self.assertEqual(result['repair_status'],'PASS');self.assertEqual(result['feature_count'],2)
            with self.assertRaisesRegex(ValueError,'requested GEO samples'):
                acquire_dnam(study,['GSM1'],original,output,metadata=metadata,series=series,supplement=supplement)

    def test_recovery_preserves_non_cg_source_features_without_imputing_or_renaming(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);inputs=self.recovery_fixture(root)
            with gzip.open(inputs[2],'at') as handle:
                handle.write('ch.1.101940785F\t0.2\t0.3\t0.4\t0.5\nrs123\t0.1\t0.4\t0.6\t0.9\n')
            result=recover_dnam_input(*inputs,root/'run')
            self.assertEqual(result['status'],'PASS')
            self.assertEqual(result['feature_id_categories'],{'cg_digits':2,'ch_prefix':1,'other':1})
            self.assertEqual(result['values_compared'],16)
            self.assertEqual(result['feature_annotation_status'],'CANNOT VERIFY')
            frame=read_matrix(result['recovered_input'],beta=True)
            self.assertEqual(list(frame.index),['cg001','cg002','ch.1.101940785F','rs123'])
            output=root/'direct';output.mkdir()
            direct=acquire_dnam(inputs[0],list(frame.columns),Path(result['recovered_input']),output)
            self.assertEqual(direct['feature_count'],4)
            pd.testing.assert_frame_equal(frame,read_matrix(direct['path'],beta=True))

    def test_recovered_clock_table_joins_historical_values_by_gsm_and_rejects_ambiguity(self):
        from scripts.run_recovered_dnam_clocks_2026_09_16 import summarize
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);fresh=root/'fresh.csv';old=root/'old.csv'
            fresh.write_text('sample_id,Hannum\nGSM1,10\nGSM2,20\n')
            old.write_text(',Hannum,unrelated,unrelated\nGSM2,21,1,2\nGSM1,11,3,4\n')
            receipt={'path':str(fresh),'clocks':['Hannum'],'checks':[{'clock':'Hannum',
                'required_features':71,'observed_features_min':60,'numerical_status':'PASS',
                'max_abs_independent_difference':0}]}
            recovery={'source_samples':2,'study_id':'GSE999001'}
            result=summarize(receipt,recovery,old,root)
            self.assertEqual(result['rows'][0]['max_absolute_change'],1)
            self.assertEqual(result['rows'][0]['imputed_cpgs_max'],11)
            self.assertEqual(result['scientific_status'],'CANNOT VERIFY')
            old.write_text(',Hannum,Hannum\nGSM1,10,10\nGSM2,20,20\n')
            with self.assertRaisesRegex(ValueError,'missing/ambiguous'):
                summarize(receipt,recovery,old,root)

    def test_recovery_rejects_ambiguous_or_incomplete_identity_and_unlinked_sources(self):
        for variant in ('ambiguous','missing','unlinked','nonempty'):
            with self.subTest(variant=variant),tempfile.TemporaryDirectory() as folder:
                root=Path(folder);inputs=self.recovery_fixture(root)
                if variant=='ambiguous':inputs[3].write_text(inputs[3].read_text().replace('[P2]','[P1]'))
                if variant=='missing':
                    with gzip.open(inputs[2],'wt') as handle:handle.write('ID\tP1\ncg001\t0.2\n')
                if variant=='unlinked':inputs[4].write_text(inputs[4].read_text().replace('ftp.ncbi.nlm.nih.gov','example.org'))
                if variant=='nonempty':
                    with gzip.open(inputs[1],'at') as handle:handle.write('cg001,0.1,0.2,0.3,0.4\n')
                with self.assertRaises(ValueError):propose_dnam_recovery(*inputs)

    def test_recovery_refuses_ragged_duplicate_nonbeta_and_all_missing_sources(self):
        for tail in ('cg001\t0.2\t0.3\t0.4\t0.5\n',
                     'cg003\t0.2\t0.3\n','cg003\t1.2\t0.3\t0.4\t0.5\n',
                     '\t0.2\t0.3\t0.4\t0.5\n','invalid probe\t0.2\t0.3\t0.4\t0.5\n'):
            with self.subTest(tail=tail),tempfile.TemporaryDirectory() as folder:
                root=Path(folder);inputs=self.recovery_fixture(root)
                with gzip.open(inputs[2],'at') as handle:handle.write(tail)
                proposal=propose_dnam_recovery(*inputs);output=root/'out';output.mkdir()
                with self.assertRaises(ValueError):reconstruct_dnam_input(proposal,output)
                self.assertFalse((output/'reconstruction.json').exists())
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);inputs=self.recovery_fixture(root)
            with gzip.open(inputs[2],'rt') as handle:text=handle.read()
            with gzip.open(inputs[2],'wt') as handle:handle.write(text.replace('\t0.1','\tNA').replace('\t0.5','\tNA'))
            proposal=propose_dnam_recovery(*inputs);output=root/'out';output.mkdir()
            with self.assertRaisesRegex(ValueError,'no measured features'):reconstruct_dnam_input(proposal,output)

    def test_independent_recovery_verifier_rejects_value_or_mapping_tampering(self):
        for variant in ('value','mapping','source'):
            with self.subTest(variant=variant),tempfile.TemporaryDirectory() as folder:
                root=Path(folder);inputs=self.recovery_fixture(root)
                proposal=propose_dnam_recovery(*inputs);output=root/'out';output.mkdir()
                result=reconstruct_dnam_input(proposal,output)
                if variant=='value':
                    with gzip.open(result['path'],'rt') as handle:text=handle.read()
                    with gzip.open(result['path'],'wt') as handle:handle.write(text.replace('\t0.1','\t0.9'))
                    result['sha256']=sha(result['path'])
                elif variant=='mapping':
                    proposal['crosswalk'][0]['quote']='Invented source label [P1]'
                else:inputs[3].write_text(inputs[3].read_text()+'\n')
                with self.assertRaises(ValueError):verify_dnam_recovery(proposal,result)

    def test_comparison_plan_freezes_all_cases_and_enforces_budget_before_calls(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);packet,clocks,design,primary,_=self.measured_fixture(root)
            packet['study_id']='GSE1'
            proposal={**design,'study_id':'GSE1','biological_unit':'software fixture',
                      'excluded_ids':[],'evidence':[{'sample_id':s['sample_id'],'field':'characteristics_ch1.0',
                        'quote':s['characteristics_ch1.0']} for s in packet['samples']]}
            family=root/'families.json';family.write_text(json.dumps(build_families([],{'GSE1'})))
            archive=root/'archive-identities.json';archive.write_text(json.dumps({'source_sha256':{}}))
            plan={'schema':'clockbase-development-comparison-v1','split':'development',
                  'provider':'openai','model':'gpt-5.6-sol','methods':list(COMPARISON_METHODS),
                  'review_calls_per_model_arm':1,'baseline_rule':'retain_producer_RESULT_rows_without_audit',
                  'family_manifest':str(family),'archive_identity_manifest':str(archive),
                  'cases':[{'study_id':'GSE1','archive':'software-fixture'}],
                  'max_calls':2,'max_cost':.00001,'output_tokens':3500,'prompt_bytes':48000,
                  'input_per_million':5.,'output_per_million':20.}
            plan_path=root/'plan.json';plan_path.write_text(json.dumps(plan))
            prepared={'study_id':'GSE1','packet':packet,'primary':primary,'clocks':clocks,
                      'source_sha256':{clocks['path']:sha(clocks['path'])}}
            calls=[]
            class Client:
                def invoke(self,messages):
                    calls.append(messages)
                    return AIMessage(content=json.dumps(proposal))
            with patch('audit_system.report_2026_09_10.load_archived_comparison',return_value=prepared):
                with self.assertRaisesRegex(ValueError,'budget cannot reserve'):
                    run_comparison_plan(plan_path,root/'run',client_factory=lambda *args:Client())
                self.assertFalse((root/'run').exists())
                plan['max_cost']=2.;plan_path.write_text(json.dumps(plan))
                ready=run_comparison_plan(plan_path,root/'run',prepare_only=True)
                self.assertTrue(ready['prepared']);self.assertEqual(calls,[])
                result=run_comparison_plan(plan_path,root/'run',client_factory=lambda *args:Client())
                self.assertEqual(len(calls),2);self.assertEqual(len(result['records']),4)
                self.assertFalse(result['independent_reference_available'])
                self.assertEqual(run_comparison_plan(plan_path,root/'run',client_factory=lambda *args:Client()),result)
                self.assertEqual(len(calls),2)
                plan['max_cost']=3.;plan_path.write_text(json.dumps(plan))
                with self.assertRaisesRegex(ValueError,'Configuration, source or implementation changed'):
                    run_comparison_plan(plan_path,root/'run',client_factory=lambda *args:Client())

    def attempted_plan_fixture(self, root, *, ready=False):
        """Pinned failure and optional prepared case; software inputs only."""
        upstream = root / 'failed-upstream'; attempt = upstream / '02_molecular_input/attempt-000'
        attempt.mkdir(parents=True)
        config = upstream / 'config.json'; write_new(config, {'study': 'GSE2', 'scope': 'software fixture'})
        failure = attempt / 'failure.json'
        write_new(failure, {'status': 'CANNOT VERIFY', 'error': 'ValueError: software fixture prerequisite failure',
                            'elapsed_seconds': 0.2})
        partial = attempt / 'partial.txt'; partial.write_text('Unscored software fixture partial output\n')
        expected = {str(path): sha(path) for path in (config, failure, partial)}
        cases = [{'study_id': 'GSE2', 'archive': str(upstream), 'failure_receipt': str(failure)}]
        prepared = design = None
        if ready:
            source = root / 'ready'; source.mkdir()
            packet, clocks, design, primary, _ = self.measured_fixture(source)
            packet['study_id'] = 'GSE1'
            prepared = {'study_id': 'GSE1', 'packet': packet, 'primary': primary, 'clocks': clocks,
                        'source_sha256': {clocks['path']: sha(clocks['path'])}}
            expected.update(prepared['source_sha256'])
            cases.insert(0, {'study_id': 'GSE1', 'archive': 'software-fixture-ready'})
        family = root / 'families.json'
        family.write_text(json.dumps(build_families([], {case['study_id'] for case in cases})))
        identities = root / 'archive-identities.json'; write_new(identities, {'source_sha256': expected})
        plan = {'schema': 'clockbase-development-comparison-v2', 'split': 'development',
            'provider': 'openai', 'model': 'gpt-5.6-sol', 'methods': list(COMPARISON_METHODS),
            'review_calls_per_model_arm': 1, 'baseline_rule': 'retain_producer_RESULT_rows_without_audit',
            'shared_prerequisite_failure_rule': 'retain_all_methods_without_execution',
            'family_manifest': str(family), 'archive_identity_manifest': str(identities), 'cases': cases,
            'max_calls': 2 * len(cases), 'max_cost': 2., 'output_tokens': 3500, 'prompt_bytes': 48000,
            'input_per_million': 5., 'output_per_million': 20.}
        plan_path = root / 'plan.json'; write_new(plan_path, plan)
        return plan_path, prepared, design, failure

    def test_attempted_comparison_preserves_terminal_failure_and_later_recovery(self):
        from audit_system.evaluation_2026_09_10 import (freeze_reference_unavailable, read_reference_records,
                                                       evaluate, export_reference_designs)
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); plan, _, _, failure = self.attempted_plan_fixture(root)
            no_call = lambda *args: self.fail('Unavailable prerequisite invoked a model')
            run = root / 'comparison'
            result = run_comparison_plan(plan, run, client_factory=no_call)
            self.assertEqual(result['attempted_studies'], 1)
            self.assertEqual(result['attempted_method_records'], 4)
            self.assertEqual(result['generation_attempts'], 0)
            self.assertTrue(all(not row['accepted'] and not row['flagged_material_error'] for row in result['records']))
            sources = export_reference_designs(run, root / 'source-reference')
            self.assertEqual(sources['packets'], [])
            self.assertEqual(len(sources['unavailable_cases']), 1)
            frozen = freeze_reference_unavailable(run, 'GSE2', root / 'unavailable-reference')
            records = read_reference_records(frozen['frozen_outcomes'])
            self.assertTrue(all(row['reference_adjudicator'] is None for row in records))
            metrics = evaluate(reference_freezes=[{'path': frozen['frozen_outcomes'], 'sha256': frozen['sha256']}])
            for value in metrics['methods'].values():
                self.assertEqual(value['attempted_studies'], 1)
                self.assertEqual(value['reference_resolved_studies'], 0)
                self.assertEqual(value['material_error_denominator'], 0)
                self.assertEqual(value['shared_prerequisite_unavailable_studies'], 1)
                self.assertIsNone(value['accepted_result_error_rate'])
            before = {str(p): sha(p) for p in run.rglob('*') if p.is_file()}
            # A later upstream completion must not overwrite the frozen failure.
            write_new(failure.parent.parent / 'complete.json', {'scope': 'later software fixture recovery'})
            self.assertEqual(run_comparison_plan(plan, run, client_factory=no_call), result)
            self.assertEqual(before, {str(p): sha(p) for p in run.rglob('*') if p.is_file()})
            self.assertEqual(read_reference_records(frozen['frozen_outcomes']), records)
            failure.write_text(failure.read_text() + '\n')
            with self.assertRaisesRegex(ValueError, 'source identity changed'):
                run_comparison_plan(plan, run, client_factory=no_call)

    def test_attempted_comparison_mixes_ready_and_unavailable_without_losing_roster(self):
        from audit_system.evaluation_2026_09_10 import reference_cases, freeze_reference_unavailable
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); plan, prepared, design, _ = self.attempted_plan_fixture(root, ready=True)
            packet = prepared['packet']; calls = []
            proposal = {**design, 'study_id': 'GSE1', 'biological_unit': 'software fixture', 'excluded_ids': [],
                'evidence': [{'sample_id': s['sample_id'], 'field': 'characteristics_ch1.0',
                             'quote': s['characteristics_ch1.0']} for s in packet['samples']]}
            class Client:
                def invoke(self, messages):
                    calls.append(messages); return AIMessage(content=json.dumps(proposal))
            with patch('audit_system.report_2026_09_10.load_archived_comparison', return_value=prepared):
                result = run_comparison_plan(plan, root / 'run', client_factory=lambda *args: Client())
            self.assertEqual(len(calls), 2)
            self.assertEqual(len(result['records']), 8)
            self.assertEqual(len(reference_cases(root / 'run')), 2)
            self.assertTrue(any(row['accepted'] for row in result['records'] if row['study_id'] == 'GSE1'))
            self.assertTrue(all(not row['accepted'] for row in result['records'] if row['study_id'] == 'GSE2'))
            with self.assertRaisesRegex(ValueError, 'completed analyses need independent adjudication'):
                freeze_reference_unavailable(root / 'run', 'GSE1', root / 'invalid-reference')

    def test_unavailable_attempt_rejects_identity_integrity_and_completed_upstream(self):
        from audit_system.report_2026_09_10 import load_unavailable_comparison
        for mode in ('wrong_study', 'changed', 'outside', 'completed', 'unselected_attempt'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); plan_path, _, _, failure = self.attempted_plan_fixture(root)
                plan = json.loads(plan_path.read_text()); case = {**plan['cases'][0], 'family_id': 'software family'}
                expected = json.loads(Path(plan['archive_identity_manifest']).read_text())['source_sha256']
                if mode == 'wrong_study': case['study_id'] = 'GSE999'
                elif mode == 'changed': failure.write_text(failure.read_text() + '\n')
                elif mode == 'outside': case['failure_receipt'] = str(root / 'elsewhere/failure.json')
                elif mode == 'completed': write_new(failure.parent.parent / 'complete.json', {'software': True})
                elif mode == 'unselected_attempt': (failure.parent.parent / 'attempt-001').mkdir()
                with self.assertRaises(ValueError): load_unavailable_comparison(case, expected, 'fixture enrollment')

    def test_attempt_resource_summary_retains_failed_request_reservation_and_known_usage(self):
        from audit_system.evaluation_2026_09_10 import comparison_resource_record
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); actor = root / 'api/GSE1-independent_audit-000'; actor.mkdir(parents=True)
            write_new(actor / 'reservation.json', {'reserved_usd': .2, 'scope': 'Software receipt fixture'})
            row = {'study_id': 'GSE1', 'method': 'independent_audit'}
            failed = comparison_resource_record(root, row)
            self.assertEqual(failed['generation_attempts'], 1)
            self.assertEqual(failed['unknown_billing_attempts'], 1)
            self.assertIsNone(failed['billed_cost_usd'])
            write_new(actor / 'receipt.json', {'estimated_cost_usd': .05, 'returned_model': 'fixture-model',
                'usage': {'prompt_tokens': 12, 'completion_tokens': 34}})
            complete = comparison_resource_record(root, row)
            self.assertEqual(complete['unknown_billing_attempts'], 0)
            self.assertEqual(complete['estimated_known_cost_usd'], .05)
            self.assertEqual(complete['prompt_tokens'], 12)
            self.assertEqual(complete['completion_tokens'], 34)
            self.assertIsNone(complete['billed_cost_usd'])

    def test_four_arm_comparison_shares_verification_and_resumes_without_calls(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);packet,clocks,design,primary,_=self.measured_fixture(root)
            primary['rows'][0]['mean_difference']+=100000
            primary['rows'][0]['reason']='PRECEDING_ANALYSIS_SENTINEL'
            frozen=copy.deepcopy(primary)
            proposal={**design,'study_id':packet['study_id'],'biological_unit':'software fixture',
                      'excluded_ids':[],'evidence':[{'sample_id':s['sample_id'],'field':'characteristics_ch1.0',
                        'quote':s['characteristics_ch1.0']} for s in packet['samples']]}
            captured={}
            class Client:
                def __init__(self,method):self.method=method
                def invoke(self,messages):
                    captured[self.method]=messages
                    return AIMessage(content=json.dumps(proposal))
            store=StageStore(root/'comparison',{'scope':'software fixture'})
            case={'study_id':packet['study_id'],'family_id':'software-fixture-family',
                  'packet':packet,'clocks':clocks,'primary':primary}
            rows=compare_development_case(store,case,lambda study,method:Client(method))
            by_method={row['method']:row for row in rows}
            self.assertEqual(len(rows),4)
            self.assertTrue(by_method['existing_pipeline']['accepted'])
            self.assertFalse(by_method['existing_pipeline']['flagged_material_error'])
            for method in ('deterministic','independent_audit','self_review'):
                self.assertTrue(by_method[method]['accepted'])
                self.assertTrue(by_method[method]['flagged_material_error'])
                self.assertTrue(by_method[method]['repair_proposed'])
                self.assertIsNone(by_method[method]['accepted_result_error_rate'])
            self.assertNotIn('PRECEDING_ANALYSIS_SENTINEL',captured['independent_audit'][1].content)
            self.assertIn('PRECEDING_ANALYSIS_SENTINEL',captured['self_review'][1].content)
            self.assertEqual(primary,frozen)
            def no_call(*args):raise AssertionError('Completed model stage repeated')
            resumed=compare_development_case(store,case,no_call)
            self.assertEqual([r['method_result'] for r in rows],[r['method_result'] for r in resumed])

    def test_comparator_retains_model_failure_as_an_unresolved_attempt(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);packet,clocks,_,primary,_=self.measured_fixture(root)
            class FailedClient:
                def invoke(self,messages):raise TimeoutError('software fixture transport failure')
            case={'study_id':packet['study_id'],'family_id':'software-fixture-family',
                  'packet':packet,'clocks':clocks,'primary':primary}
            rows=compare_development_case(StageStore(root/'comparison',{}),case,lambda *args:FailedClient())
            self.assertEqual(len(rows),4)
            for row in rows[2:]:
                self.assertEqual(row['decision_status'],'CANNOT VERIFY')
                self.assertIn('TimeoutError',row['execution_error'])
                self.assertFalse(row['accepted'])

    def measured_fixture(self, root):
        """Construct a software fixture, not an empirical methylation study."""
        from biolearn.model_gallery import ModelGallery
        coefficients=ModelGallery().get('Hannum',imputation_method='none').coefficients.iloc[:,0]
        sites=[site for site in coefficients.index if site!='intercept']
        ids=['GSM1','GSM2','GSM3','GSM4']
        beta=pd.DataFrame({sid:[value]*len(sites) for sid,value in zip(ids,[.2,.25,.4,.5])},index=sites)
        beta_path=root/'software-fixture-beta.csv'; beta.to_csv(beta_path)
        coefficient_path=root/'coefficients.csv'; coefficients.to_csv(coefficient_path)
        predicted=coefficients.loc[sites].to_numpy() @ beta.to_numpy()+float(coefficients.get('intercept',0))
        predictions=pd.DataFrame({'sample_id':ids,'Hannum':predicted})
        prediction_path=root/'software-fixture-predictions.csv'; predictions.to_csv(prediction_path,index=False)
        packet={'study_id':'SOFTWARE_FIXTURE','target':'drugX versus vehicle','samples':[]}
        for index,sid in enumerate(ids):
            fields=['donor id: '+sid,'treatment: '+('drugX' if index<2 else 'vehicle'),
                    'age (y): 50','sex: F','genotype: declared fixture group','tissue: whole blood','batch: A']
            packet['samples'].append({'sample_id':sid,'organism_ch1':'Homo sapiens','molecule_ch1':'genomic DNA',
                                      'platform_id':'GPL13534',**{f'characteristics_ch1.{i}':v for i,v in enumerate(fields)}})
        clocks={'path':str(prediction_path),'clocks':['Hannum'],'sample_count':len(ids),'checks':[{'clock':'Hannum','numerical_status':'PASS',
                 'required_features':len(sites),'observed_features_min':len(sites),'measurement_sources':{
                 key:{'path':str(path),'sha256':sha(path)} for key,path in [('raw_beta',beta_path),('coefficients',coefficient_path)]}}]}
        design={'proposed_action':'recompute_welch','status':'DESIGN_PROPOSED','pairing':'independent',
                'case_ids':ids[:2],'control_ids':ids[2:],'findings':[],'unresolved':[]}
        original=recompute_declared_contrast(predictions.set_index('sample_id'),ids[:2],ids[2:],'Hannum','Welch')
        original['status']='RESULT'
        return packet,clocks,design,{'rows':[original]},beta

    def test_failed_auditor_proposal_is_not_counted_as_a_producer_error_in_four_arms(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary, _ = self.measured_fixture(root)
            before = copy.deepcopy(primary)
            proposal = {**design, 'study_id': packet['study_id'], 'biological_unit': 'software fixture',
                        'case_ids': design['control_ids'], 'control_ids': design['case_ids'], 'excluded_ids': [],
                        'evidence': [{'sample_id': row['sample_id'], 'field': 'characteristics_ch1.1',
                                      'quote': row['characteristics_ch1.1']} for row in packet['samples']]}
            class Client:
                def invoke(self, messages): return AIMessage(content=json.dumps(proposal))
            case = {'study_id': packet['study_id'], 'family_id': 'software-fixture-family',
                    'packet': packet, 'clocks': clocks, 'primary': primary}
            rows = compare_development_case(StageStore(root / 'comparison', {'scope': 'software fixture'}),
                                            case, lambda *args: Client())
            by_method = {row['method']: row for row in rows}
            self.assertTrue(by_method['existing_pipeline']['accepted'])
            self.assertTrue(by_method['deterministic']['accepted'])
            for method in ('independent_audit', 'self_review'):
                row = by_method[method]; check = row['method_result']['checks'][0]
                self.assertTrue(row['accepted'])
                self.assertFalse(row['flagged_material_error'])
                self.assertEqual(row['failed_proposal_results'], 1)
                self.assertEqual(check['primary_source_design_verification']['status'], 'PASS')
                self.assertEqual(check['primary_measurement_verification']['status'], 'PASS')
                self.assertEqual(check['candidate_source_design_verification']['status'], 'FAIL')
                self.assertEqual(check['producer_error_signals'], [])
                self.assertTrue(check['proposal_errors'])
                self.assertIsNone(row['accepted_result_error_rate'])
                self.assertEqual(check['repair_decision'], 'REJECTED')
                self.assertEqual(check['original_scientific_status'], 'PASS')
                self.assertEqual(check['candidate_scientific_status'], 'FAIL')
                final = row['method_result']['final_outputs'][0]
                self.assertEqual(final['selected_output'], 'original')
                self.assertFalse(final['repair_applied'])
                for field in ('case_ids', 'control_ids', 'mean_difference', 'p_value', 'effect_size'):
                    self.assertEqual(final[field], before['rows'][0][field])
                ledger = next((root / 'comparison' / packet['study_id'] / method).rglob('decision-ledger.jsonl'))
                records = [json.loads(line) for line in ledger.read_text().splitlines()]
                dependencies = [item['parent_id'] for item in records
                                if item.get('child_id') == 'claim:scientific:0' and item['requirement'] == 'REQUIRED_CHECK']
                subjects = {item['subject_id'] for item in records if item['record_id'] in dependencies}
                self.assertEqual(subjects, {'artifact:original'})
            self.assertEqual(primary, before)

    def test_proposal_computation_failure_keeps_original_verification_and_subject(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary, _ = self.measured_fixture(root)
            # Direct callers still receive a separate failure for an invalid
            # proposal even though the model parser rejects singleton arms.
            design['case_ids'] = design['case_ids'][:1]
            result = numerical_audit(design, primary, clocks, root, packet=packet)
            check = result['checks'][0]
            self.assertEqual(check['calculation_status'], 'PASS')
            self.assertEqual(check['primary_source_design_verification']['status'], 'PASS')
            self.assertFalse(check['producer_error_allegation'])
            self.assertTrue(check['proposal_attempted'])
            self.assertEqual(check['proposal_verification_status'], 'FAIL')
            self.assertIn('cannot be recomputed', check['proposal_errors'][0])
            self.assertIsNone(check['candidate_source_design_verification'])
            self.assertEqual(result['accepted_results'], 1)
            self.assertEqual(result['final_outputs'][0]['selected_output'], 'original')

    def test_unresolved_candidate_keeps_an_original_that_passes_its_own_checks(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            original = primary['rows'][0]
            original.update(recompute_declared_contrast(pd.read_csv(clocks['path']).set_index('sample_id'),
                original['case_ids'], original['control_ids'], 'Hannum', 'Student t-test'))
            design['comparisons'][0]['target'] = 'Does drugX reverse radiation resistance?'
            before = copy.deepcopy(primary)
            output = root / 'audit'; output.mkdir()
            result = numerical_audit(design, primary, clocks, output, packet=packet)
            check, final = result['checks'][0], result['final_outputs'][0]
            self.assertEqual(check['original_scientific_status'], 'PASS')
            self.assertEqual(check['candidate_source_design_verification']['status'], 'CANNOT VERIFY')
            self.assertEqual(check['candidate_scientific_status'], 'CANNOT VERIFY')
            self.assertEqual(check['repair_decision'], 'CANNOT VERIFY')
            self.assertEqual(final['selected_output'], 'original')
            self.assertTrue(final['accepted_result'])
            self.assertFalse(final['repair_applied'])
            self.assertEqual(final['mean_difference'], before['rows'][0]['mean_difference'])
            self.assertEqual(final['target'], before['rows'][0]['target'])
            self.assertNotEqual(final['target'], final['audit_target'])
            self.assertNotEqual(result['candidates'][0]['p_value'], final['p_value'])
            self.assertEqual(primary, before)

    def test_no_passing_claim_preserves_unaccepted_original_and_rejected_candidate(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary, _ = self.measured_fixture(root)
            primary['rows'][0]['mean_difference'] += 1000
            before = copy.deepcopy(primary)
            design['case_ids'], design['control_ids'] = design['control_ids'], design['case_ids']
            result = numerical_audit(design, primary, clocks, root, packet=packet)
            check, final = result['checks'][0], result['final_outputs'][0]
            self.assertEqual(check['original_scientific_status'], 'FAIL')
            self.assertEqual(check['candidate_scientific_status'], 'FAIL')
            self.assertEqual(check['repair_decision'], 'REJECTED')
            self.assertEqual(final['selected_output'], 'original')
            self.assertFalse(final['accepted_result'])
            self.assertEqual(final['mean_difference'], before['rows'][0]['mean_difference'])
            self.assertEqual(len(result['candidates']), 1)
            self.assertEqual(primary, before)

    def test_source_review_uncertainty_is_applied_to_original_before_fallback(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary, _ = self.measured_fixture(root)
            design.update(proposed_action='obtain_source_evidence', unresolved=['Unresolved source-unit conflict'])
            result = numerical_audit(design, primary, clocks, root, packet=packet)
            check = result['checks'][0]
            self.assertEqual(check['primary_source_design_verification']['status'], 'PASS')
            self.assertEqual(check['original_reviewed_design_verification']['status'], 'CANNOT VERIFY')
            self.assertEqual(check['original_scientific_status'], 'CANNOT VERIFY')
            self.assertFalse(result['final_outputs'][0]['accepted_result'])

    def test_single_comparison_adjustment_retains_unsuccessful_clock_endpoints(self):
        from scipy.stats import false_discovery_control
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary, _ = self.measured_fixture(root)
            clocks['clocks'].append('UnavailableClock')
            primary['rows'].append({**primary['rows'][0], 'clock': 'UnavailableClock', 'status': 'CANNOT VERIFY',
                                    'p_value': '', 'mean_difference': '', 'effect_size': '', 'statistical_test': ''})
            result = numerical_audit(design, primary, clocks, root, packet=packet)
            self.assertEqual(len(result['candidates']), 1)
            candidate = result['candidates'][0]
            expected = float(false_discovery_control([candidate['p_value'], 1.], method='by')[0])
            self.assertAlmostEqual(candidate['p_value_by_attempted_family'], expected, places=14)
            self.assertNotIn('q_value_run_family', candidate)

    def test_preserved_comparison_replay_checks_sources_and_resumes_without_generation(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            paths = [primary['path'], primary['manifest_path'], clocks['path']]
            paths += [item['path'] for item in clocks['checks'][0]['measurement_sources'].values()]
            case = {'study_id': packet['study_id'], 'family_id': 'software-fixture-family',
                    'packet': packet, 'clocks': clocks, 'primary': primary,
                    'source_sha256': {path: sha(path) for path in paths}}
            archive = root/'preserved'; store = StageStore(archive, {'development_only': True})
            def freeze(directory):
                write_new(directory/'prepared-cases.json', [case])
                return {'cases': [case]}
            store.run('00_inputs', freeze)
            calls = []
            class Client:
                def invoke(self, messages):
                    calls.append(messages)
                    return AIMessage(content=json.dumps(design))
            compare_development_case(store, case, lambda *args: Client())
            before = {str(path): sha(path) for path in archive.rglob('*') if path.is_file()}
            replay = root/'replay'
            with patch.object(BoundedOpenAI, 'invoke', side_effect=AssertionError('Replay must not call a model')):
                result = replay_development_comparison(archive, replay)
                saved = {str(path): sha(path) for path in replay.rglob('*') if path.is_file()}
                self.assertEqual(replay_development_comparison(archive, replay), result)
                self.assertEqual(saved, {str(path): sha(path) for path in replay.rglob('*') if path.is_file()})
            self.assertEqual(len(calls), 2)
            self.assertEqual(result['generation_calls'], 0)
            self.assertEqual(result['study_count'], 1)
            self.assertEqual(len(result['records']), 4)
            self.assertTrue(all(row['accepted'] for row in result['records']))
            self.assertEqual(before, {str(path): sha(path) for path in archive.rglob('*') if path.is_file()})
            # A separate software archive exercises the older API-only source
            # layout. These transport fixtures are not empirical model runs.
            import shutil
            legacy = root/'legacy-preserved'; shutil.copytree(archive, legacy)
            for method in ('independent_audit', 'self_review'):
                review = legacy/packet['study_id']/method/'review'
                receipt = json.loads((review/'complete.json').read_text())
                attempt = review/'attempt-000'
                normalized = json.loads((attempt/'design.json').read_text())
                model = json.loads((attempt/'model-design.json').read_text())
                reviewed = json.loads((attempt/'review-result.json').read_text())
                shutil.rmtree(attempt); attempt.mkdir()
                write_new(attempt/'design.json', normalized); write_new(attempt/'review-result.json', reviewed)
                receipt['outputs'] = {str(path.relative_to(review)): sha(path) for path in attempt.rglob('*') if path.is_file()}
                (review/'complete.json').write_text(json.dumps(receipt))
                (review/'complete.sha256').write_text(sha(review/'complete.json')+'\n')
                slot = legacy/'api'/(packet['study_id']+'-'+method+'-000'); slot.mkdir(parents=True)
                body = packet if method == 'independent_audit' else {'source_packet': packet, 'preceding_analysis': {'rows': primary['rows']}}
                write_new(slot/'request.json', {'messages': [{'role': 'system', 'content': 'Software fixture'},
                                                            {'role': 'user', 'content': json.dumps(body)}]})
                write_new(slot/'response.json', {'choices': [{'message': {'content': json.dumps(model)}}]})
                write_new(slot/'receipt.json', {'response_sha256': sha(slot/'response.json')})
            legacy_result = replay_development_comparison(legacy, root/'legacy-replay')
            self.assertEqual(len(legacy_result['records']), 4)
            self.assertTrue(all(row['accepted'] for row in legacy_result['records']))
            (slot/'response.json').write_text('{}')
            with self.assertRaisesRegex(ValueError, 'API receipt'):
                replay_development_comparison(legacy, root/'changed-response')
            self.assertFalse((root/'changed-response').exists())
            Path(clocks['path']).write_text('modified source')
            with self.assertRaisesRegex(ValueError, 'source changed'):
                replay_development_comparison(archive, root/'changed-source')
            self.assertFalse((root/'changed-source').exists())

    def test_absent_root_pooling_condition_does_not_block_separate_comparisons(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            design['findings'] = [{'category': 'SOURCE_DESIGN_CONCERN', 'statement': 'Pooling source treatments requires justification',
                                   'affected_result': None, 'condition': {'operation': 'pools_groups', 'arm': 'case',
                                   'groups': [['GSM1', 'GSM2'], ['GSM5', 'GSM6']]},
                                   'evidence': design['comparisons'][0]['evidence']}]
            directory = root / 'audit'; directory.mkdir()
            result = numerical_audit(design, primary, clocks, directory, packet=packet)
            self.assertEqual(result['accepted_results'], 2)
            self.assertTrue(all(not row['producer_error_allegation'] for row in result['checks']))
            self.assertTrue(all(row['finding_attribution'][0]['attribution'] == 'ACTION_NOT_PRESENT' for row in result['checks']))
            self.assertEqual(len(design['findings']), 1)
            # A condition supplied directly still requires checked source quotes.
            design['findings'][0]['evidence'] = []
            rejected = root / 'unquoted'; rejected.mkdir()
            with self.assertRaisesRegex(ValueError, 'quotation'):
                numerical_audit(design, primary, clocks, rejected, packet=packet)

    def test_supported_fixture_is_accepted_without_becoming_reference_truth(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); packet,clocks,design,primary,_=self.measured_fixture(root)
            result=numerical_audit(design,primary,clocks,root,packet=packet)
            self.assertEqual(result['accepted_results'],1)
            self.assertEqual(result['study_decision_status'],'PASS')
            self.assertEqual(result['checks'][0]['source_design_verification']['status'],'PASS')
            self.assertEqual(result['checks'][0]['measurement_verification']['status'],'PASS')
            self.assertEqual(result['final_outputs'][0]['independent_reference_status'],'CANNOT VERIFY')
            self.assertIsNone(result['accepted_result_error_rate'])
            self.assertEqual(json.loads((root/'final-outputs.json').read_text()),result['final_outputs'])

    def multiple_comparison_fixture(self, root):
        """Six-sample software fixture: two treatments share two controls."""
        packet, clocks, design, _, beta = self.measured_fixture(root)
        for sid, value in [('GSM5', .32), ('GSM6', .35)]:
            beta[sid] = value
            sample = copy.deepcopy(packet['samples'][0]); sample['sample_id'] = sid
            sample['characteristics_ch1.0'] = 'donor id: ' + sid
            sample['characteristics_ch1.1'] = 'treatment: drugY'
            packet['samples'].append(sample)
        packet['target'] = 'drugX versus vehicle; drugY versus vehicle'
        sources = clocks['checks'][0]['measurement_sources']
        beta.to_csv(sources['raw_beta']['path']); sources['raw_beta']['sha256'] = sha(sources['raw_beta']['path'])
        coefficients = pd.read_csv(sources['coefficients']['path'], index_col=0).iloc[:, 0]
        values = coefficients.loc[beta.index].to_numpy() @ beta.to_numpy() + float(coefficients.get('intercept', 0))
        frame = pd.DataFrame({'sample_id': list(beta.columns), 'Hannum': values})
        frame.to_csv(clocks['path'], index=False); clocks['sample_count'] = 6
        children, declarations, rows = [], [], []
        for name, ids in [('drugX', ['GSM1', 'GSM2']), ('drugY', ['GSM5', 'GSM6'])]:
            item = {'comparison_name': name, 'target': name + ' versus vehicle',
                    'case_ids': ids, 'control_ids': ['GSM3', 'GSM4']}
            declarations.append(item)
            children.append({**design, **item, 'biological_unit': 'Explicit source donor IDs',
                             'excluded_ids': sorted(set(beta.columns) - set(ids + ['GSM3', 'GSM4'])),
                             'evidence': [{'sample_id': s['sample_id'], 'field': 'characteristics_ch1.1',
                                           'quote': s['characteristics_ch1.1']} for s in packet['samples']]})
            row = recompute_declared_contrast(frame.set_index('sample_id'), ids, ['GSM3', 'GSM4'], 'Hannum', 'Welch')
            rows.append({**row, 'comparison_name': name, 'status': 'RESULT', 'reason': 'Software fixture'})
        manifest = {'schema': PRIMARY_SCHEMA, 'comparisons': declarations, 'unassigned_ids': []}
        manifest_path = root / 'comparisons.json'; manifest_path.write_text(json.dumps(manifest))
        csv_path = root / 'aging_analysis.csv'
        fields = ['comparison_name', 'clock', 'case_ids', 'control_ids', 'statistical_test', 'p_value',
                  'mean_difference', 'effect_size', 'status', 'reason']
        with csv_path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore'); writer.writeheader()
            writer.writerows({key: json.dumps(value) if isinstance(value, list) else value
                              for key, value in row.items()} for row in rows)
        primary = read_primary_outputs(csv_path, clocks, packet, manifest_path)
        design = parse_design(json.dumps({'schema': DESIGN_SET_SCHEMA, 'study_id': packet['study_id'],
                              'status': 'DESIGN_PROPOSED', 'comparisons': children,
                              'unassigned_ids': [], 'findings': [], 'unresolved': [], 'evidence': []}), packet)
        return packet, clocks, design, primary

    def test_shared_control_comparisons_are_separately_verified_and_reported(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            before = copy.deepcopy(primary)
            primary['rows'][1]['comparison_name'] = 'Different producer name'
            primary['rows'][1]['mean_difference'] += 1000
            changed = copy.deepcopy(primary)
            directory = root / 'audit'; directory.mkdir()
            result = numerical_audit(design, primary, clocks, directory, packet=packet)
            self.assertEqual(result['comparison_count'], 2)
            self.assertEqual(result['accepted_results'], 2)
            self.assertEqual(result['verified_declared_calculations'], 1)
            self.assertEqual(set(result['shared_samples']), {'GSM3', 'GSM4'})
            self.assertEqual(len({row['result_id'] for row in result['checks']}), 2)
            self.assertEqual(primary, changed)
            self.assertEqual(json.loads((directory / 'original-analysis.json').read_text()), changed)
            for final, original in zip(result['final_outputs'], before['rows'], strict=True):
                self.assertAlmostEqual(final['mean_difference'], original['mean_difference'], delta=1e-12)
            self.assertEqual(result['multiplicity']['attempted_endpoints'], 2)
            self.assertTrue(all('q_value_run_family' not in row for row in result['final_outputs']))
            report_dir = root / 'report'; report_dir.mkdir()
            summary = report(root, design, primary, result, clocks, report_dir)
            self.assertEqual(summary['comparison_count'], 2)
            self.assertIn('drugY', (report_dir / 'report.html').read_text())

    def test_missing_producer_comparison_is_new_analysis_not_an_erroneous_result(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            primary['rows'] = primary['rows'][:1]; primary['comparisons'] = primary['comparisons'][:1]
            directory = root / 'audit'; directory.mkdir()
            result = numerical_audit(design, primary, clocks, directory, packet=packet)
            missing = next(row for row in result['checks'] if not row['producer_present'])
            self.assertEqual(missing['primary_analysis_state']['state'], 'NOT_REPORTED')
            self.assertEqual(missing['calculation_status'], 'CANNOT VERIFY')
            self.assertEqual(missing['deterministic_errors'], [])
            self.assertFalse(missing['producer_error_allegation'])
            self.assertEqual(missing['repair_kind'], 'NEW_ANALYSIS_FOR_MISSING_COMPARISON')
            self.assertEqual(result['accepted_results'], 2)
            self.assertEqual(len(primary['rows']), 1)

    def test_empty_and_partial_source_requests_do_not_invent_missing_results(self):
        for partial in (False, True):
            with self.subTest(partial=partial), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
                request = design['comparisons'][1]
                request.update(case_ids=['GSM5'] if partial else [], control_ids=[], excluded_ids=[],
                               status='CANNOT VERIFY', proposed_action='obtain_source_evidence',
                               unresolved=['The control assignment needs source adjudication'])
                request.pop('comparison_id')
                before = copy.deepcopy(primary)
                output = root / 'audit'; output.mkdir()
                result = numerical_audit(design, primary, clocks, output, packet=packet)
                self.assertEqual(result['comparison_count'], 2)
                self.assertEqual(result['defined_comparison_count'], 2)
                self.assertEqual(result['attempted_clock_results'], 2)
                self.assertEqual(result['unresolved_request_count'], 1)
                self.assertEqual(result['unresolved_requests'][0]['case_ids'], ['GSM5'] if partial else [])
                self.assertEqual(result['unresolved_requests'][0]['evidence'], request['evidence'])
                self.assertEqual(result['unresolved_requests'][0]['associated_comparison_ids'],
                                 [primary['comparisons'][1]['comparison_id']])
                self.assertEqual(result['accepted_results'], 1)
                self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
                self.assertTrue(all(check['primary_analysis_state']['state'] == 'COMPUTED' for check in result['checks']))
                self.assertIsNone(result['multiplicity']['family_denominator'])
                self.assertEqual(result['multiplicity']['denominator_status'], 'CANNOT VERIFY')
                self.assertTrue(all(row['p_value_by_attempted_family'] is None for row in result['final_outputs']))
                self.assertEqual(json.loads((output / 'unresolved-requests.json').read_text()), result['unresolved_requests'])
                self.assertEqual(primary, before)
                report_dir = root / 'report'; report_dir.mkdir()
                self.assertEqual(report(root, design, primary, result, clocks, report_dir)['unresolved_request_count'], 1)

    def test_unlinked_request_keeps_known_results_but_withholds_full_family_adjustment(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            request = {**design['comparisons'][0], 'target': 'An additional source question',
                       'case_ids': [], 'control_ids': [], 'excluded_ids': [], 'status': 'CANNOT VERIFY',
                       'unresolved': ['Additional groups are undefined'], 'proposed_action': 'obtain_source_evidence'}
            request.pop('comparison_id'); design['comparisons'].append(request)
            output = root / 'audit'; output.mkdir()
            result = numerical_audit(design, primary, clocks, output, packet=packet)
            self.assertEqual(result['comparison_count'], 2)
            self.assertEqual(result['accepted_results'], 2)
            self.assertEqual(result['unresolved_requests'][0]['associated_comparison_ids'], [])
            self.assertEqual(result['study_decision_status'], 'PASS')
            self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
            self.assertTrue(all(row['p_value_by_attempted_family'] is None for row in result['final_outputs']))

    def test_undefined_request_retains_findings_linked_to_a_producer_result(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            request = {**design['comparisons'][0], 'target': 'Additional context with an original-result concern',
                       'case_ids': [], 'control_ids': [], 'excluded_ids': [], 'status': 'CANNOT VERIFY',
                       'proposed_action': 'obtain_source_evidence', 'unresolved': ['Control identity unresolved']}
            request.pop('comparison_id')
            request['findings'] = [{'category': 'SOURCE_DESIGN_CONCERN', 'statement': 'Unresolved source identity',
                'evidence': request['evidence'][:1],
                'affected_result': primary['comparisons'][0]['comparison_id'] + '/Hannum', 'condition': None}]
            design['comparisons'].append(request)
            output = root / 'audit'; output.mkdir()
            result = numerical_audit(design, primary, clocks, output, packet=packet)
            self.assertEqual(result['accepted_results'], 1)
            self.assertEqual(result['checks'][0]['original_scientific_status'], 'CANNOT VERIFY')
            self.assertEqual(result['checks'][0]['model_findings'], request['findings'])
            self.assertEqual(result['checks'][1]['original_scientific_status'], 'PASS')

    def test_only_undefined_requests_have_no_computed_or_missing_result_records(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            primary['rows'] = []; primary['comparisons'] = []
            design['comparisons'] = design['comparisons'][:1]
            request = design['comparisons'][0]
            request.update(case_ids=[], control_ids=[], excluded_ids=[], status='CANNOT VERIFY',
                           proposed_action='obtain_source_evidence', unresolved=['Groups unknown'])
            request.pop('comparison_id')
            output = root / 'audit'; output.mkdir()
            result = numerical_audit(design, primary, clocks, output, packet=packet)
            self.assertEqual(result['comparison_count'], 0)
            self.assertEqual(result['final_outputs'], [])
            self.assertEqual(result['attempted_clock_results'], 0)
            self.assertEqual(result['unresolved_request_count'], 1)
            self.assertIsNone(result['multiplicity']['family_denominator'])
            self.assertEqual(result['study_decision_status'], 'CANNOT VERIFY')

    def test_producer_record_with_undefined_groups_is_preserved(self):
        from audit_system.agent_workflow_2026_09_10 import comparison_identity
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            primary['rows'] = primary['rows'][:1]; primary['comparisons'] = primary['comparisons'][:1]
            cid = comparison_identity(packet['study_id'], [], [], primary['comparisons'][0]['target'])
            for item in primary['rows'] + primary['comparisons']:
                item.update(comparison_id=cid, case_ids=[], control_ids=[])
            before = copy.deepcopy(primary)
            design['comparisons'] = []
            output = root / 'audit'; output.mkdir()
            result = numerical_audit(design, primary, clocks, output, packet=packet)
            self.assertEqual(result['comparison_count'], 1)
            self.assertEqual(result['defined_comparison_count'], 0)
            self.assertEqual(result['attempted_clock_results'], 1)
            self.assertFalse(result['final_outputs'][0].get('absence_record', False))
            self.assertEqual(result['final_outputs'][0]['mean_difference'], before['rows'][0]['mean_difference'])
            self.assertIsNone(result['multiplicity']['family_denominator'])
            self.assertEqual(primary, before)

    def test_different_source_membership_cannot_replace_a_producer_comparison_by_name(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            design['comparisons'] = design['comparisons'][:1]
            child = design['comparisons'][0]
            child['case_ids'], child['control_ids'] = child['control_ids'], child['case_ids']
            child.pop('comparison_id')
            directory = root / 'audit'; directory.mkdir()
            result = numerical_audit(design, primary, clocks, directory, packet=packet)
            self.assertEqual(result['comparison_count'], 3)
            self.assertEqual(result['verified_declared_calculations'], 2)
            self.assertEqual(result['accepted_results'], 0)
            self.assertEqual(sum(row['producer_present'] for row in result['checks']), 2)
            self.assertTrue(any(row['source_design_verification']['status'] == 'FAIL' for row in result['checks']))

    def test_new_child_target_does_not_authorize_a_different_study_question(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            packet['target'] = 'An unrelated phenotype question'
            directory = root / 'audit'; directory.mkdir()
            result = numerical_audit(design, primary, clocks, directory, packet=packet)
            self.assertEqual(result['verified_declared_calculations'], 2)
            self.assertEqual(result['accepted_results'], 0)
            self.assertTrue(all(row['source_design_verification']['study_target_scope']['status'] == 'CANNOT VERIFY'
                                for row in result['checks']))

    def test_declared_source_discovery_scope_does_not_require_generic_target_wording(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            packet.update(target='Source-defined intervention contrasts at the same age',
                          comparison_scope='source_defined_interventions')
            directory = root / 'audit'; directory.mkdir()
            result = numerical_audit(design, primary, clocks, directory, packet=packet)
            self.assertEqual(result['accepted_results'], 2)
            self.assertTrue(all(row['source_design_verification']['study_target_scope']['status'] == 'PASS'
                                for row in result['checks']))
            # Declaring study-wide discovery does not authorize source-inconsistent sample assignments.
            packet['samples'][0]['characteristics_ch1.1'] = 'treatment: unrelated'
            changed = root / 'source-conflict'; changed.mkdir()
            result = numerical_audit(design, primary, clocks, changed, packet=packet)
            self.assertLess(result['accepted_results'], 2)

    def test_unrecognized_scope_cannot_authorize_a_comparison(self):
        from src.pipeline_2026_09_10 import comparison_target_scope
        with self.assertRaisesRegex(ValueError, 'Unknown comparison scope'):
            comparison_target_scope({'target': 'drugX versus vehicle', 'comparison_scope': 'approve_everything'},
                                    'drugX versus vehicle')

    def test_source_discovery_cannot_relabel_an_intervention_by_its_observed_response(self):
        from src.pipeline_2026_09_10 import comparison_target_scope
        with tempfile.TemporaryDirectory() as folder:
            packet, _, _, primary, _ = self.measured_fixture(Path(folder))
            cases, controls = primary['rows'][0]['case_ids'], primary['rows'][0]['control_ids']
            for sample in packet['samples']:
                sample['characteristics_ch1.response'] = 'post-treatment response: ' + (
                    'improved' if sample['sample_id'] in cases else 'unchanged')
            packet.update(comparison_scope='source_defined_interventions', target='improved versus unchanged')
            scope = comparison_target_scope(packet, packet['target'], cases, controls)
            self.assertEqual(scope['status'], 'CANNOT VERIFY')
            self.assertEqual(scope['intervention_binding']['status'], 'PASS')
            self.assertEqual(verify_source_design(packet, cases, controls)['status'], 'CANNOT VERIFY')
            for sample in packet['samples']:
                sample['characteristics_ch1.response'] = 'post-treatment response: retained'
            scope = comparison_target_scope(packet, 'drugX versus vehicle within post-treatment response: retained',
                                            cases, controls)
            self.assertEqual(scope['intervention_binding']['status'], 'PASS')
            self.assertEqual(scope['status'], 'CANNOT VERIFY')
            self.assertEqual(scope['unsupported_source_strata'], ['post treatment response'])

    def test_legacy_single_comparison_cannot_bypass_intervention_scope(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary, _ = self.measured_fixture(root)
            packet['comparison_scope'] = 'source_defined_interventions'
            allowed = root / 'declared-intervention'; allowed.mkdir()
            self.assertEqual(numerical_audit(design, primary, clocks, allowed, packet=packet)['accepted_results'], 1)
            packet['target'] = 'drugX versus vehicle within post-treatment response: retained'
            for sample in packet['samples']:
                sample['characteristics_ch1.response'] = 'post-treatment response: retained'
            blocked = root / 'response-stratum'; blocked.mkdir()
            result = numerical_audit(design, primary, clocks, blocked, packet=packet)
            self.assertEqual(result['accepted_results'], 0)
            self.assertEqual(result['study_decision_status'], 'CANNOT VERIFY')
            self.assertEqual(result['verified_declared_calculations'], 1)
            check = result['checks'][0]
            for key in ('original_reviewed_design_verification', 'candidate_source_design_verification'):
                self.assertEqual(check[key]['study_target_scope']['status'], 'CANNOT VERIFY')
                self.assertEqual(check[key]['study_target_scope']['unsupported_source_strata'], ['post treatment response'])

    def test_legacy_targets_are_bound_separately_for_producer_and_proposal(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary, _ = self.measured_fixture(root)
            unexpected = 'drugX versus vehicle within post-treatment response: retained'
            for sample in packet['samples']:
                sample['characteristics_ch1.response'] = 'post-treatment response: retained'
            primary['rows'][0]['target'] = unexpected
            design['proposed_action'] = 'obtain_source_evidence'
            blocked = root / 'changed-producer-target'; blocked.mkdir()
            result = numerical_audit(design, primary, clocks, blocked, packet=packet)
            self.assertEqual(result['accepted_results'], 0)
            self.assertEqual(result['checks'][0]['original_reviewed_design_verification']['study_target_scope']['comparison_target'], unexpected)
            primary['rows'][0]['target'] = packet['target']
            design.update(target=unexpected, proposed_action='recompute_welch')
            separate = root / 'changed-proposal-target'; separate.mkdir()
            result = numerical_audit(design, primary, clocks, separate, packet=packet)
            self.assertEqual(result['accepted_results'], 1)
            self.assertEqual(result['checks'][0]['selected_output'], 'original')
            self.assertEqual(result['checks'][0]['candidate_source_design_verification']['study_target_scope']['status'], 'CANNOT VERIFY')

    def test_legacy_selected_proposal_preserves_its_checked_target(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary, _ = self.measured_fixture(root)
            packet['comparison_scope'] = 'source_defined_interventions'
            design['target'] = 'drugX versus vehicle within age y: 50'
            directory = root / 'distinct-proposal-target'; directory.mkdir()
            result = numerical_audit(design, primary, clocks, directory, packet=packet)
            self.assertEqual(result['accepted_results'], 1)
            self.assertEqual(result['checks'][0]['selected_output'], 'candidate')
            self.assertEqual(result['candidates'][0]['target'], design['target'])
            self.assertEqual(result['final_outputs'][0]['target'], design['target'])
            self.assertNotEqual(result['final_outputs'][0]['target'], packet['target'])
            self.assertNotIn('target', primary['rows'][0])

    def test_root_ambiguity_blocks_children_but_neutral_observation_does_not(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); packet, clocks, design, primary = self.multiple_comparison_fixture(root)
            design['findings'] = [{'category': 'OBSERVATION', 'statement': 'Shared controls', 'evidence': [], 'affected_result': None}]
            first = root / 'observation'; first.mkdir()
            self.assertEqual(numerical_audit(design, primary, clocks, first, packet=packet)['accepted_results'], 2)
            design['unresolved'] = ['Biological independence requires source adjudication']
            second = root / 'ambiguity'; second.mkdir()
            self.assertEqual(numerical_audit(design, primary, clocks, second, packet=packet)['accepted_results'], 0)

    def test_primary_manifest_requires_all_comparison_clock_rows_and_exact_groups(self):
        for change in ('duplicate', 'omit', 'swap', 'unused', 'repeat_comparison'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); packet, clocks, _, primary = self.multiple_comparison_fixture(root)
                path = Path(primary['path']); manifest_path = Path(primary['manifest_path'])
                with path.open() as handle:
                    reader = csv.DictReader(handle); fields = reader.fieldnames; rows = list(reader)
                if change == 'duplicate': rows.append(rows[0])
                if change == 'omit': rows.pop()
                if change == 'swap': rows[0]['case_ids'] = json.dumps(['GSM5', 'GSM6'])
                manifest = json.loads(manifest_path.read_text())
                if change == 'unused': manifest['unassigned_ids'] = ['GSM1']
                if change == 'repeat_comparison':
                    manifest['comparisons'].append({**manifest['comparisons'][0], 'comparison_name': 'another'})
                manifest_path.write_text(json.dumps(manifest))
                with path.open('w', newline='') as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
                with self.assertRaises(ValueError): read_primary_outputs(path, clocks, packet, manifest_path)

    def test_dependent_adjustment_matches_scipy_with_unresolved_endpoints_in_denominator(self):
        from scipy.stats import false_discovery_control
        values = [.2, .001, .01, .001]
        expected = false_discovery_control(values + [1., 1.], method='by')[:len(values)]
        self.assertEqual(adjust_dependent_family([], 2), [])
        for calculated, independent in zip(adjust_dependent_family(values, 6), expected):
            self.assertAlmostEqual(calculated, independent, places=14)
        self.assertTrue(all(a >= b for a, b in zip(adjust_dependent_family(values, 6), adjust_dependent_family(values, 4))))
        with self.assertRaises(ValueError): adjust_dependent_family(values, 3)

    def test_unbound_target_cannot_be_approved_merely_by_mentioning_the_drug(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);packet,clocks,design,primary,_=self.measured_fixture(root)
            packet['target']='Does drugX reverse radiation resistance?'
            result=numerical_audit(design,primary,clocks,root,packet=packet)
            self.assertEqual(result['checks'][0]['calculation_status'],'PASS')
            self.assertEqual(result['checks'][0]['source_design_verification']['target_binding']['status'],'CANNOT VERIFY')
            self.assertEqual(result['accepted_results'],0)

    def test_analysis_after_abstention_is_not_a_corrected_computed_result(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);packet,clocks,design,primary,_=self.measured_fixture(root)
            primary['rows'][0].update(status='CANNOT VERIFY',p_value='',mean_difference='',effect_size='',statistical_test='')
            result=numerical_audit(design,primary,clocks,root,packet=packet)
            self.assertEqual(result['accepted_results'],1)
            self.assertEqual(result['checks'][0]['primary_analysis_state']['state'],'ABSTAINED')
            self.assertEqual(result['checks'][0]['deterministic_errors'],[])
            self.assertEqual(result['final_outputs'][0]['repair_kind'],'NEW_ANALYSIS_AFTER_ABSTENTION')
            self.assertEqual(result['final_outputs'][0]['historical_release_status'],'CANNOT VERIFY')

    def test_declared_abstention_with_numbers_retains_a_schema_contradiction(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);packet,clocks,design,primary,_=self.measured_fixture(root)
            primary['rows'][0]['status']='CANNOT VERIFY'
            design.update(proposed_action='obtain_source_evidence',unresolved=['fixture unresolved design'])
            result=numerical_audit(design,primary,clocks,root,packet=packet)
            self.assertEqual(result['checks'][0]['primary_analysis_state']['state'],'INCONSISTENT')
            self.assertEqual(result['checks'][0]['calculation_status'],'FAIL')
            self.assertEqual(result['accepted_results'],0)

    def test_verified_arithmetic_repair_preserves_original_and_updates_final_output(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); packet,clocks,design,primary,_=self.measured_fixture(root)
            primary['rows'][0]['mean_difference']+=100000
            before=copy.deepcopy(primary)
            result=numerical_audit(design,primary,clocks,root,packet=packet)
            self.assertEqual(result['checks'][0]['calculation_status'],'FAIL')
            self.assertEqual(result['checks'][0]['repair_status'],'PROPOSED')
            self.assertTrue(result['checks'][0]['producer_error_allegation'])
            self.assertEqual(result['accepted_results'],1)
            self.assertNotEqual(result['final_outputs'][0]['mean_difference'],primary['rows'][0]['mean_difference'])
            self.assertEqual(primary,before)
            self.assertEqual(json.loads((root/'original-analysis.json').read_text()),before)
            self.assertEqual(result['final_outputs'][0]['selected_output'],'candidate')
            self.assertEqual(result['checks'][0]['repair_decision'],'ACCEPTED')
            self.assertEqual(result['checks'][0]['original_scientific_status'], 'FAIL')
            self.assertEqual(result['checks'][0]['candidate_scientific_status'], 'PASS')
            self.assertTrue(result['final_outputs'][0]['repair_applied'])
            report_dir=root/'report'; report_dir.mkdir()
            summary=report(root,design,primary,result,clocks,report_dir)
            self.assertEqual(summary['study_decision_status'],'PASS')
            self.assertFalse(summary['independent_reference_available'])
            self.assertIn('ACCEPTED',(report_dir/'report.html').read_text())

    def test_receipt_feature_count_cannot_overrule_selected_measured_values(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); packet,clocks,_,_,beta=self.measured_fixture(root)
            clocks['checks'][0]['observed_features_min']=0
            packet['samples'][0]['characteristics_ch1.2']='age (y): 50 years'
            result=verify_clock_measurement(packet,clocks,'Hannum',['GSM1','GSM2','GSM3','GSM4'])
            self.assertEqual(result['status'],'PASS')
            beta['GSM4']=float('nan'); beta.loc['not-a-clock-feature']=.3
            source=clocks['checks'][0]['measurement_sources']['raw_beta']
            beta.to_csv(source['path']); source['sha256']=sha(source['path'])
            result=verify_clock_measurement(packet,clocks,'Hannum',['GSM1','GSM2','GSM3'])
            self.assertEqual(result['status'],'PASS')
            result=verify_clock_measurement(packet,clocks,'Hannum',['GSM1','GSM2','GSM3','GSM4'])
            self.assertEqual(result['status'],'FAIL')
            self.assertIn('no measured required CpGs',' '.join(result['failures']))

    def write_series_fixture(self, root, beta, packet, ending=True):
        """Encode the labeled software fixture in GEO's table format."""
        matrix=root/'software-fixture-matrix.txt.gz'
        with gzip.open(matrix,'wt') as handle:
            handle.write('!Series_geo_accession\t"GSE1"\n!Sample_geo_accession\t'+'\t'.join(beta.columns)+'\n')
            handle.write('!series_matrix_table_begin\n')
            beta.rename_axis('ID_REF').to_csv(handle,sep='\t')
            if ending:handle.write('!series_matrix_table_end\n')
        metadata=root/'software-fixture-metadata.csv'
        pd.DataFrame([{**sample,'GSM':sample['sample_id']} for sample in packet['samples']]).to_csv(metadata,index=False)
        return matrix,metadata

    def test_streamed_clock_execution_retains_missing_samples_and_all_source_rows(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); packet,_,_,_,beta=self.measured_fixture(root)
            beta.iloc[0,0]=float('nan')
            beta.loc['unselected-probe']=.2
            matrix,metadata=self.write_series_fixture(root,beta,packet)
            result=check_measurement_profiles(matrix,metadata,root/'validation')
            hannum=result['rows'][0]
            self.assertEqual(hannum['source_samples'],4)
            self.assertEqual(hannum['complete_measurement_samples'],3)
            self.assertEqual(hannum['incomplete_measurement_samples'],1)
            self.assertEqual(hannum['measurement_status'],'CANNOT VERIFY')
            self.assertEqual(result['accepted_intervention_results'],0)
            self.assertIsNone(result['accepted_result_error_rate'])
            saved=pd.read_csv(root/'validation/clock_predictions.csv').set_index('sample_id')
            self.assertTrue(pd.isna(saved.loc['GSM1','Hannum']))
            extraction=json.loads((root/'validation/extraction.json').read_text())
            self.assertEqual(extraction['source_rows'],len(beta))
            self.assertEqual(extraction['selected_sites'],len(beta)-1)

    def test_streamed_matrix_rejects_missing_end_duplicate_probe_and_membership(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); packet,_,_,_,beta=self.measured_fixture(root)
            output=root/'out'; output.mkdir()
            matrix,metadata=self.write_series_fixture(root,beta,packet,ending=False)
            with self.assertRaisesRegex(ValueError,'end marker'):
                extract_clock_sites(matrix,metadata,beta.index,output)
            duplicate=pd.concat([beta,beta.iloc[:1]])
            matrix,metadata=self.write_series_fixture(root,duplicate,packet)
            with self.assertRaisesRegex(ValueError,'duplicate GEO probe'):
                extract_clock_sites(matrix,metadata,beta.index,output)
            packet['samples']=packet['samples'][:-1]
            matrix,metadata=self.write_series_fixture(root,beta,packet)
            with self.assertRaisesRegex(ValueError,'membership differs'):
                extract_clock_sites(matrix,metadata,beta.index,output)

    def test_missing_measured_cpg_cannot_be_hidden_by_reported_full_coverage(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); packet,clocks,_,_,beta=self.measured_fixture(root)
            source=clocks['checks'][0]['measurement_sources']['raw_beta']
            beta.iloc[1:].to_csv(source['path']); source['sha256']=sha(source['path'])
            result=verify_clock_measurement(packet,clocks,'Hannum',['GSM1','GSM2','GSM3','GSM4'])
            self.assertEqual(result['status'],'CANNOT VERIFY')
            self.assertTrue(all(v==len(beta)-1 for v in result['observed_features_by_sample'].values()))

    def test_changed_coefficients_fail_despite_a_matching_declared_hash(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); packet,clocks,_,_,_=self.measured_fixture(root)
            source=clocks['checks'][0]['measurement_sources']['coefficients']
            values=pd.read_csv(source['path'],index_col=0); values.iloc[0,0]+=1
            values.to_csv(source['path']); source['sha256']=sha(source['path'])
            result=verify_clock_measurement(packet,clocks,'Hannum',['GSM1','GSM2','GSM3','GSM4'])
            self.assertEqual(result['status'],'FAIL')
            self.assertIn('pinned model',' '.join(result['failures']))

    def test_outside_measurement_domain_and_missing_source_remain_unresolved(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); packet,clocks,_,_,_=self.measured_fixture(root)
            packet['samples'][0]['characteristics_ch1.5']='tissue: brain'
            result=verify_clock_measurement(packet,clocks,'Hannum',['GSM1','GSM2','GSM3','GSM4'])
            self.assertEqual(result['status'],'CANNOT VERIFY')
            source=clocks['checks'][0]['measurement_sources']['raw_beta']
            Path(source['path']).unlink()
            result=verify_clock_measurement(packet,clocks,'Hannum',['GSM1','GSM2','GSM3','GSM4'])
            self.assertEqual(result['status'],'CANNOT VERIFY')

    def test_age_confounded_fixture_is_rejected_despite_correct_arithmetic(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); packet,clocks,design,primary,_=self.measured_fixture(root)
            for sample in packet['samples'][:2]:sample['characteristics_ch1.2']='age (y): 70'
            result=numerical_audit(design,primary,clocks,root,packet=packet)
            self.assertEqual(result['checks'][0]['calculation_status'],'PASS')
            self.assertEqual(result['checks'][0]['source_design_verification']['status'],'FAIL')
            self.assertEqual(result['accepted_results'],0)
            self.assertEqual(result['study_decision_status'],'FAIL')

    def test_missing_or_repeated_unit_declarations_do_not_assume_independence(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); packet,_,design,_,_=self.measured_fixture(root)
            packet['samples'][2]['characteristics_ch1.0']='donor id: GSM1'
            result=verify_source_design(packet,design['case_ids'],design['control_ids'])
            self.assertEqual(result['status'],'CANNOT VERIFY')
            del packet['samples'][2]['characteristics_ch1.0']
            result=verify_source_design(packet,design['case_ids'],design['control_ids'])
            self.assertEqual(result['status'],'CANNOT VERIFY')

    def test_unit_id_punctuation_is_preserved_and_unknown_values_are_unresolved(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); packet,_,design,_,_=self.measured_fixture(root)
            packet['samples'][0]['characteristics_ch1.0']='donor id: A-1'
            packet['samples'][1]['characteristics_ch1.0']='donor id: A_1'
            self.assertEqual(verify_source_design(packet,design['case_ids'],design['control_ids'])['status'],'PASS')
            packet['samples'][0]['characteristics_ch1.0']='donor id: UNKNOWN'
            self.assertEqual(verify_source_design(packet,design['case_ids'],design['control_ids'])['status'],'CANNOT VERIFY')

    def test_publication_failure_retains_geo_context_and_missing_evidence(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); series=root/'input.soft'; destination=root/'context'
            series.write_text('^SERIES = GSE1\n!Series_title = Software fixture\n'
                              '!Series_overall_design = Declared comparison\n!Series_sample_id = GSM1\n')
            with patch('src.pipeline_inputs_2026_09_10.fetch', side_effect=RuntimeError('download failed')) as call, \
                 patch('src.pipeline_inputs_2026_09_10.time.sleep'):
                path=acquire_study_context('GSE1',series,destination)
            documents,acquisition=read_study_context(path,'GSE1')
            self.assertEqual(call.call_count,1)
            self.assertEqual(len(documents),1)
            self.assertTrue(all(r['status']=='CANNOT VERIFY' for r in acquisition))
            self.assertIn('download failed',json.dumps(acquisition))

    def test_context_offline_mode_skips_publication_network_and_preserves_original(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); series=root/'input.soft'; destination=root/'context'
            payload='^SERIES = GSE1\n!Series_title = Fixture\n!Series_sample_id = GSM1\n'
            series.write_text(payload)
            with patch('src.pipeline_inputs_2026_09_10.fetch') as call:
                path=acquire_study_context('GSE1',series,destination,publication_methods=False)
                call.assert_not_called()
            self.assertEqual(series.read_text(),payload)
            self.assertEqual(read_study_context(path,'GSE1')[1][0]['status'],'NOT_REQUESTED')
            with self.assertRaises(FileExistsError):
                acquire_study_context('GSE1',series,destination,publication_methods=False)

    def test_download_never_overwrites_an_existing_source(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'source.xml'; path.write_text('preserved')
            with patch('src.pipeline_inputs_2026_09_10.subprocess.run') as call, self.assertRaises(FileExistsError):
                fetch('https://example.test/source',path)
            call.assert_not_called()
            self.assertEqual(path.read_text(),'preserved')

    def test_agent_accepts_injected_client_and_executor(self):
        class Client:
            calls = 0
            def invoke(self, messages):
                self.calls += 1
                return AIMessage(content='<execute>print(1)</execute>' if self.calls == 1 else '<solution>done</solution>')
        executed = []
        agent = CodingAgent(llm=Client(), executor=lambda code: executed.append(code) or '1', strict_tags=True, printout=False)
        agent.configure(); _, final = agent.go('software fixture')
        self.assertEqual(executed, ['print(1)']); self.assertIn('done', final)

    def test_incomplete_code_is_not_repaired_or_executed(self):
        class Client:
            def invoke(self, messages): return AIMessage(content='<execute>print(1)')
        agent = CodingAgent(llm=Client(), executor=lambda code: self.fail('must not execute'), strict_tags=True, printout=False)
        agent.configure()
        with self.assertRaisesRegex(ValueError, 'complete'): agent.go('software fixture')

    @unittest.skipUnless(sys.platform.startswith('linux') and shutil.which('bwrap'),
                         'Agent isolation requires Linux and bubblewrap on PATH')
    def test_sandbox_cannot_read_key_home_or_network(self):
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {'OPENAI_API_KEY': 'software-fixture-secret'}):
            text = IsolatedExecutor(folder)("import os,socket\nfrom pathlib import Path\nprint(os.environ.get('OPENAI_API_KEY'))\nprint(Path('/home/kying').exists())\ntry:\n socket.create_connection(('1.1.1.1',443),timeout=1)\n print('NETWORK_OPEN')\nexcept OSError:\n print('NETWORK_BLOCKED')")
            self.assertIn('None\nFalse\nNETWORK_BLOCKED', text)
            self.assertNotIn('software-fixture-secret', text)

    def test_stage_resume_checks_hash_and_does_not_rerun(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); store = StageStore(root, {'version': 1})
            def stage(d): (d / 'data.txt').write_text('source'); return {'ok': True}
            self.assertEqual(store.run('s', stage), store.run('s', lambda d: self.fail('reran')))
            (root / 's/attempt-000/data.txt').write_text('changed')
            with self.assertRaisesRegex(ValueError, 'changed'): store.run('s', stage)

    def test_uncertain_billing_prevents_repeat(self):
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {'OPENAI_API_KEY': 'software-fixture-key'}):
            client = BoundedOpenAI(folder, 'test')
            with patch('src.pipeline_runtime_2026_09_10.subprocess.run', side_effect=TimeoutError):
                with self.assertRaises(TimeoutError): client.invoke([HumanMessage(content='fixture')])
            client = BoundedOpenAI(folder, 'test')
            with patch('src.pipeline_runtime_2026_09_10.subprocess.run') as transport:
                with self.assertRaisesRegex(RuntimeError, 'uncertain billing'): client.invoke([HumanMessage(content='fixture')])
                transport.assert_not_called()

    def test_receipt_mutation_cannot_change_resumed_result(self):
        with tempfile.TemporaryDirectory() as folder:
            store = StageStore(folder, {'version': 1})
            store.run('s', lambda d: {'scientific_status': 'CANNOT VERIFY'})
            path = Path(folder) / 's/complete.json'
            record = json.loads(path.read_text()); record['result']['scientific_status'] = 'PASS'
            path.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, 'receipt changed'):
                store.run('s', lambda d: self.fail('must not rerun'))

    def test_cost_cap_checked_before_transport(self):
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {'OPENAI_API_KEY': 'software-fixture-key'}):
            client = BoundedOpenAI(folder, 'test', max_cost=.000001)
            with patch('src.pipeline_runtime_2026_09_10.subprocess.run') as transport:
                with self.assertRaisesRegex(RuntimeError, 'cost limit'): client.invoke([HumanMessage(content='fixture')])
                transport.assert_not_called()

    def test_generic_provider_key_is_never_forwarded_to_openai(self):
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {'AGENT_API_KEY':'software-other-provider-key'}, clear=True):
            client = BoundedOpenAI(folder, 'test')
            with patch('src.pipeline_runtime_2026_09_10.subprocess.run') as transport:
                with self.assertRaisesRegex(RuntimeError, 'OPENAI_API_KEY'):
                    client.invoke([HumanMessage(content='software fixture')])
                transport.assert_not_called()
            self.assertEqual(list(Path(folder).rglob('request.json')), [])

    def test_explicit_effect_contract_requires_definition_and_preserves_legacy(self):
        from src.pipeline_2026_09_10 import EFFECT_DEFINITION
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'analysis.csv'
            fields = ['clock','case_ids','control_ids','statistical_test','p_value','mean_difference',
                      'effect_size','status','reason','effect_size_definition']
            row = dict(zip(fields, ['fixture','["A","B"]','["C","D"]','Welch',.1,1.,.5,'RESULT','fixture',EFFECT_DEFINITION]))
            def save():
                with path.open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerow(row)
            save()
            result = read_primary_outputs(path, {'clocks':['fixture']}, {'samples':[]}, analysis_contract='explicit_effect_v2')
            self.assertEqual(result['rows'][0]['effect_size_definition'], EFFECT_DEFINITION)
            with self.assertRaises(ValueError): read_primary_outputs(path, {'clocks':['fixture']}, {'samples':[]})
            row['effect_size_definition']='';save()
            with self.assertRaisesRegex(ValueError, 'definition'):
                read_primary_outputs(path, {'clocks':['fixture']}, {'samples':[]}, analysis_contract='explicit_effect_v2')

    def test_duplicate_raw_header_and_empty_matrix_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'source.csv'
            for content in ['gene,GSM1,GSM1\nx,1,2\n', 'gene,GSM1,GSM2\n']:
                path.write_text(content)
                with self.assertRaisesRegex(ValueError, 'AUDIT_FAIL'): read_matrix(path)

    def test_multiline_geo_preamble_never_becomes_a_measurement_header(self):
        # Labeled software regression for the observed physical/logical-row bug.
        content = ('!Series_summary\t"first line\nsecond line\nthird line"\n'
                   '!Sample_geo_accession\t"GSM1"\t"GSM2"\n'
                   '!series_matrix_table_begin\n"ID_REF"\t"GSM1"\t"GSM2"\n'
                   '"001"\t0.2\t0.3\n"cg2"\tNA\t0.4\n'
                   '!series_matrix_table_end\n')
        with tempfile.TemporaryDirectory() as folder:
            for suffix, opener in (('.txt', open), ('.txt.gz', gzip.open)):
                path = Path(folder) / ('fixture' + suffix)
                with opener(path, 'wt', newline='') as handle: handle.write(content)
                frame = read_matrix(path, beta=True)
                self.assertEqual(frame.columns.tolist(), ['GSM1', 'GSM2'])
                self.assertEqual(frame.index.tolist(), ['001', 'cg2'])
                self.assertEqual(frame.iloc[0].tolist(), [.2, .3])
                self.assertTrue(pd.isna(frame.iloc[1, 0]))

    def test_geo_header_is_bound_to_canonical_source_sample_order(self):
        content = ('!Sample_geo_accession\tGSM1\tGSM2\n!series_matrix_table_begin\n'
                   'ID_REF\tGSM1\tGSM2\ncg1\t0.2\t0.3\n!series_matrix_table_end\n')
        variants = [content.replace('ID_REF\tGSM1\tGSM2', header) for header in
                    ['ID_REF\tGSM2\tGSM1', 'ID_REF\tGSM1\tGSM1', 'ID_REF\t0.2\t0.3',
                     'ID_REF\tGSM1\t', 'ID_REF\tGSM1\tGSM２', 'other\tGSM1\tGSM2']]
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'fixture.txt'
            for text in variants:
                with self.subTest(header=text.splitlines()[2]):
                    path.write_text(text)
                    with self.assertRaisesRegex(ValueError, 'AUDIT_FAIL'): read_matrix(path, beta=True)

    def test_ragged_records_and_geo_terminators_are_not_padded_or_dropped(self):
        content = ('!series_matrix_table_begin\nID_REF\tGSM1\tGSM2\n'
                   'cg1\t0.2\t0.3\n!series_matrix_table_end\n')
        variants = [content.replace('cg1\t0.2\t0.3', row) for row in
                    ['cg1\t0.2', 'cg1\t0.2\t0.3\t0.4', '\t0.2\t0.3',
                     'cg1\t0.2\t0.3\ncg1\t0.4\t0.5']]
        variants += [content.replace('!series_matrix_table_end\n', ''),
                     content + '!series_matrix_table_end\n', content + 'cg2\t0.4\t0.5\n',
                     content.replace('!series_matrix_table_end', '!series_matrix_table_end\tNA\tNA')]
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'fixture.txt'
            for text in variants:
                path.write_text(text)
                with self.subTest(text=text), self.assertRaisesRegex(ValueError, 'AUDIT_FAIL'):
                    read_matrix(path, beta=True)

    def test_plain_matrix_aliases_and_literal_feature_ids_survive_width_checks(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'fixture.csv'
            path.write_text('probe,subject A,subject B\n001,0.2,0.3\n002,NA,0.4\n')
            frame = read_matrix(path, beta=True)
            self.assertEqual(frame.index.tolist(), ['001', '002'])
            self.assertEqual(frame.columns.tolist(), ['subject A', 'subject B'])
            path.write_text('probe,subject A,subject B\n001,0.2\n')
            with self.assertRaisesRegex(ValueError, 'fields; expected'): read_matrix(path, beta=True)

    def test_reference_reader_versions_caches_and_preserves_prior_artifacts(self):
        from src.pipeline_inputs_2026_09_10 import reference_means, MATRIX_READER_CONTRACT
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); path = root / 'reference.txt'; cache = root / 'cache'; cache.mkdir()
            samples = ['GSM' + str(i + 1) for i in range(656)]
            path.write_text('!Series_geo_accession\tGSE40279\n!Series_summary\t"fixture\ncontinued"\n'
                            '!Sample_geo_accession\t' + '\t'.join(samples) + '\n!series_matrix_table_begin\n'
                            'ID_REF\t' + '\t'.join(samples) + '\n001\t' + '\t'.join(['0.25'] * 656) +
                            '\n002\t' + '\t'.join(['0.5'] * 656) + '\n!series_matrix_table_end\n')
            old = cache / ('GSE40279-means-' + sha(path)[:16] + '.csv')
            old.write_text('id,mean\nfixture-old-value,0.9\n'); before = old.read_bytes()
            values, record = reference_means(path, cache)
            self.assertEqual(values.index.tolist(), ['001', '002'])
            self.assertEqual(values.tolist(), [.25, .5])
            self.assertEqual(record['reader_contract'], MATRIX_READER_CONTRACT)
            self.assertEqual(old.read_bytes(), before)
            resumed, resumed_record = reference_means(path, cache)
            pd.testing.assert_series_equal(resumed, values)
            self.assertEqual(resumed_record, record)

    def test_source_mutation_between_validation_and_read_is_rejected(self):
        from src.pipeline_inputs_2026_09_10 import _inspect_matrix_stream
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'fixture.csv'; path.write_text('probe,GSM1\ncg1,0.2\n')
            def inspect_then_mutate(*args, **kwargs):
                result = _inspect_matrix_stream(*args, **kwargs)
                path.write_text('probe,GSM1\ncg1,0.8\n')
                return result
            with patch('src.pipeline_inputs_2026_09_10._inspect_matrix_stream', side_effect=inspect_then_mutate):
                with self.assertRaisesRegex(PipelineIntegrityError, 'source changed'):
                    read_matrix(path, beta=True)

    def test_source_parser_preserves_characteristics_and_hides_contacts(self):
        rows = parse_samples('^SAMPLE = GSM1\n!Sample_title = subject\n!Sample_characteristics_ch1 = age: 4\n!Sample_characteristics_ch1 = treatment: A\n!Sample_contact_email = no@example.com', 'GSE1')
        self.assertEqual(rows[0]['characteristics_ch1.1'], 'treatment: A')
        self.assertNotIn('contact_email', rows[0])

    def test_source_descriptions_survive_csv_and_both_audit_packet_paths(self):
        from audit_system.agent_workflow_2026_09_10 import source_packet
        from audit_system.corpus_audit_2026_09_20 import evidence_index
        # Software regression fixture. These labels are not empirical findings.
        payload = ('^SAMPLE = GSM1\n!Sample_title = software specimen\n'
                   '!Sample_characteristics_ch1 = tissue: lung\n'
                   '!Sample_description = software compartment, "A"\n'
                   '!Sample_description = software replicate\t2\n'
                   '!Sample_description = donor_1\n'
                   '!Sample_description.pvalue = HIDDEN_DERIVED_OUTCOME\n'
                   '!Sample_description.verdict = HIDDEN_DERIVED_VERDICT\n'
                   '!Sample_contact_email = HIDDEN_CONTACT\n'
                   '^SAMPLE = GSM2\n!Sample_title = software control\n'
                   '!Sample_characteristics_ch1 = tissue: lung\n')
        rows = parse_samples(payload, 'GSE1')
        expected = {'description': 'software compartment, "A"',
                    'description.1': 'software replicate\t2', 'description.2': 'donor_1'}
        self.assertEqual({k: v for k, v in rows[0].items() if k.startswith('description')}, expected)
        self.assertFalse(any(k.startswith('description') for k in rows[1]))
        self.assertNotIn('HIDDEN_', json.dumps(rows))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'source-metadata.csv'
            with path.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=sorted({k for row in rows for k in row}))
                writer.writeheader(); writer.writerows(rows)
            # Exact-target pipeline and source-defined comprehensive audit share
            # the source-packet reader. Both must retain every source description.
            for comparison_scope in ('exact_target', 'source_defined_interventions'):
                packet = source_packet(path, 'GSE1', 'Software regression fixture',
                                       comparison_scope=comparison_scope)
                self.assertEqual({k: packet['samples'][0][k] for k in expected}, expected)
                index = evidence_index(packet)
                for field, value in expected.items():
                    self.assertTrue(any(record.get('field') == field and record.get('value') == value
                        and record.get('sample_ids') == ['GSM1'] for record in index['evidence'].values()))
                self.assertNotIn('HIDDEN_', json.dumps(packet))

    def test_arithmetic_repair_does_not_imply_scientific_acceptance(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); path = root / 'predictions.csv'
            path.write_text('sample_id,clock\nGSM1,1\nGSM2,3\nGSM3,5\nGSM4,8\n')
            design = {'proposed_action':'recompute_welch','status':'DESIGN_PROPOSED','pairing':'independent','case_ids':['GSM1','GSM2'],
                      'control_ids':['GSM3','GSM4'],'findings':[],'unresolved':[]}
            primary = {'rows':[{'clock':'clock','case_ids':['GSM1','GSM2'],'control_ids':['GSM3','GSM4'],
                       'status':'RESULT','statistical_test':'Welch','p_value':.8,'mean_difference':123,'effect_size':5}]}
            result = numerical_audit(design,primary,{'path':str(path)},root)
            self.assertEqual(result['checks'][0]['detection_status'],'FAIL')
            self.assertEqual(result['checks'][0]['repair_status'],'PROPOSED')
            self.assertEqual(result['checks'][0]['numerical_verification'],'PASS')
            self.assertEqual(result['accepted_results'],0)
            self.assertIsNone(result['accepted_result_error_rate'])

    def test_model_group_disagreement_is_advisory_not_a_confirmed_error(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); path=root/'predictions.csv'
            path.write_text('sample_id,clock\nGSM1,1\nGSM2,3\nGSM3,5\nGSM4,8\n')
            frame=pd.read_csv(path).set_index('sample_id')
            original=recompute_declared_contrast(frame,['GSM1','GSM2'],['GSM3','GSM4'],'clock','Welch')
            original['status']='RESULT'
            design={'proposed_action':'recompute_welch','status':'DESIGN_PROPOSED','pairing':'independent',
                    'case_ids':['GSM3','GSM4'],'control_ids':['GSM1','GSM2'],'findings':[],'unresolved':[]}
            result=numerical_audit(design,{'rows':[original]},{'path':str(path)},root)
            check=result['checks'][0]
            self.assertEqual(check['calculation_status'],'PASS')
            self.assertEqual(check['detection_status'],'CANNOT VERIFY')
            self.assertEqual(check['deterministic_errors'],[])
            self.assertTrue(check['advisory_design_findings'])
            self.assertEqual(result['decisions'][0]['calculation_status'],'PASS')
            self.assertEqual(result['decisions'][0]['scientific_status'],'CANNOT VERIFY')
            self.assertEqual(result['accepted_results'],0)

    def test_unresolved_design_does_not_prevent_detection_of_wrong_arithmetic(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); path=root/'predictions.csv'
            path.write_text('sample_id,clock\nGSM1,1\nGSM2,3\nGSM3,5\nGSM4,8\n')
            design={'proposed_action':'obtain_source_evidence','status':'CANNOT VERIFY','pairing':'unknown',
                    'case_ids':[],'control_ids':[],'findings':[],'unresolved':['Biological units missing']}
            primary={'rows':[{'clock':'clock','case_ids':['GSM1','GSM2'],'control_ids':['GSM3','GSM4'],
                     'status':'RESULT','statistical_test':'Welch','p_value':.8,'mean_difference':123,'effect_size':5}]}
            result=numerical_audit(design,primary,{'path':str(path)},root)
            self.assertEqual(result['checks'][0]['calculation_status'],'FAIL')
            self.assertEqual(result['checks'][0]['detection_status'],'FAIL')
            self.assertEqual(result['decisions'][0]['calculation_status'],'FAIL')
            self.assertEqual(result['scientific_status'],'FAIL')
            self.assertEqual(result['candidates'],[])
            self.assertEqual(result['accepted_results'],0)

    def test_unsupported_primary_test_stays_unresolved(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); path=root/'predictions.csv'
            path.write_text('sample_id,clock\nGSM1,1\nGSM2,3\nGSM3,5\nGSM4,8\n')
            design={'proposed_action':'obtain_source_evidence','status':'CANNOT VERIFY','pairing':'paired',
                    'case_ids':[],'control_ids':[],'findings':[],'unresolved':['Pairing map missing']}
            primary={'rows':[{'clock':'clock','case_ids':['GSM1','GSM2'],'control_ids':['GSM3','GSM4'],
                     'status':'RESULT','statistical_test':'paired t-test','p_value':.8,'mean_difference':123,'effect_size':5}]}
            result=numerical_audit(design,primary,{'path':str(path)},root)
            self.assertEqual(result['checks'][0]['calculation_status'],'CANNOT VERIFY')
            self.assertEqual(result['scientific_status'],'CANNOT VERIFY')
            self.assertFalse(result['checks'][0]['deterministic_errors'])

    def test_welch_spelling_variants_keep_the_unequal_variance_test(self):
        from scipy.stats import t
        # Unequal sample sizes/variances distinguish Welch from the Student path.
        frame = pd.DataFrame({'clock': [1., 2., 8., 5., 6., 7., 8.]},
                             index=['A', 'B', 'C', 'D', 'E', 'F', 'G'])
        case, control = ['A', 'B', 'C'], ['D', 'E', 'F', 'G']
        mean_difference = 11 / 3 - 26 / 4
        va, vb = 43 / 3, 5 / 3
        se2 = va / 3 + vb / 4
        df = se2 ** 2 / ((va / 3) ** 2 / 2 + (vb / 4) ** 2 / 3)
        expected_p = float(2 * t.sf(abs(mean_difference) / se2 ** .5, df))
        for name in ("Welch's two-sample t-test", 'Welch independent-samples t-test'):
            with self.subTest(name=name):
                result = recompute_declared_contrast(frame, case, control, 'clock', name)
                self.assertEqual(result['statistical_test'], 'Welch t-test')
                self.assertAlmostEqual(result['p_value'], expected_p, places=13)
        for name in ('one-sided Welch t-test', 'Welch paired t-test',
                     'paired t-test', 'ANCOVA', 'Welch unknown test'):
            with self.subTest(unsupported=name), self.assertRaises(NotImplementedError):
                recompute_declared_contrast(frame, case, control, 'clock', name)

    def test_nonfinite_predictions_cannot_pass_numeric_verification(self):
        frame=pd.DataFrame({'clock':[1,3,float('nan'),8]},index=['GSM1','GSM2','GSM3','GSM4'])
        with self.assertRaisesRegex(ValueError,'nonfinite'):
            recompute_declared_contrast(frame,['GSM1','GSM2'],['GSM3','GSM4'],'clock','Welch')


@unittest.skipUnless(shutil.which('Rscript'), 'RNA measurement fixtures require Rscript on PATH')
class RnaMeasurementTests(unittest.TestCase):
    """Constructed software fixtures, never empirical evaluation observations."""

    @classmethod
    def setUpClass(cls):
        import joblib
        from src.pipeline_inputs_2026_09_10 import ASSETS, preprocess_rna, predict_rna
        from audit_system.agent_workflow_2026_09_10 import source_packet
        require_development_artifact(ASSETS / 'EN_Hazard_Mouse_All_All_scaleddiff.pkl')
        require_development_artifact(ASSETS / 'EN_Hazard_Rodents_All_All_scaleddiff.pkl')
        require_development_artifact(ASSETS / 'Gene_table_mouse.csv')
        cls.directory = tempfile.TemporaryDirectory()
        cls.root = Path(cls.directory.name)
        cls.ids = ['GSM1', 'GSM2', 'GSM3', 'GSM4']
        models = [joblib.load(ASSETS/name) for name in
                  ('EN_Hazard_Mouse_All_All_scaleddiff.pkl', 'EN_Hazard_Rodents_All_All_scaleddiff.pkl')]
        common = set(map(str, models[0].feature_names)) & set(map(str, models[1].feature_names))
        mapping = pd.read_csv(ASSETS/'Gene_table_mouse.csv', index_col=0).dropna(subset=['Entrez'])
        selected = mapping[mapping['Entrez'].map(lambda value: str(int(value)) in common)].drop_duplicates(subset='Entrez').iloc[:4]
        cls.raw = cls.root/'software-fixture-counts.csv'
        pd.DataFrame([[20,30,40,50],[40,30,20,10],[10,12,20,24],[50,20,60,45]],
                     index=selected.index, columns=cls.ids).to_csv(cls.raw)
        cls.primary = cls.root/'software-fixture-samples.soft'
        text = []
        for i, sid in enumerate(cls.ids):
            text += ['^SAMPLE = '+sid, '!Sample_title = constructed software fixture',
                     '!Sample_organism_ch1 = Mus musculus', '!Sample_molecule_ch1 = total RNA',
                     '!Sample_characteristics_ch1 = tissue: adipose tissue',
                     '!Sample_characteristics_ch1 = treatment: '+('drugX' if i < 2 else 'vehicle'),
                     '!Sample_library_strategy = RNA-Seq', '!Sample_library_source = transcriptomic',
                     '!Sample_library_selection = cDNA',
                     '!Sample_extract_protocol_ch1 = RNA was isolated from homogenized adipose tissue.',
                     '!Sample_extract_protocol_ch1 = Libraries were prepared using an mRNA Library Prep Kit.']
        cls.primary.write_text('\n'.join(text)+'\n')
        rows = parse_samples(cls.primary.read_text(), 'GSE90000000')
        cls.metadata = cls.root/'software-fixture-metadata.csv'
        with cls.metadata.open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        cls.packet = source_packet(cls.metadata, 'GSE90000000', 'drugX versus vehicle', sample_source=cls.primary)
        pre=cls.root/'preprocessing';pre.mkdir();pred=cls.root/'predictions';pred.mkdir()
        cls.processed = preprocess_rna(cls.raw, 'mouse', pre, normalization='published_mouse_rle', source=cls.packet)
        cls.clocks = predict_rna(cls.processed, 'mouse', pred)
        cls.clock = cls.clocks['clocks'][0]

    @classmethod
    def tearDownClass(cls): cls.directory.cleanup()

    def check(self, clocks=None, packet=None, selected=None):
        return verify_clock_measurement(packet or self.packet, clocks or self.clocks, self.clock, selected or self.ids)

    def test_source_assay_fields_survive_packet_and_reference_views(self):
        from audit_system.evaluation_2026_09_10 import reference_source_view
        view = reference_source_view(self.packet)
        for row in view['samples']:
            self.assertEqual(row['library_strategy'], 'RNA-Seq')
            self.assertEqual(row['library_source'], 'transcriptomic')
            self.assertEqual(row['library_selection'], 'cDNA')
            self.assertIn('extract_protocol_ch1.1', row)
        self.assertEqual(view['sample_source']['sha256'], sha(self.primary))

    def test_numerical_profile_pass_has_explicit_imputation_and_no_biological_pass(self):
        result=self.check()
        self.assertEqual(result['status'], 'PASS', result)
        self.assertEqual(result['biological_status'], 'CANNOT VERIFY')
        self.assertEqual(result['observed_features_by_sample'],dict.fromkeys(self.ids,4))
        self.assertTrue(all(value>0 for value in result['imputed_features_by_sample'].values()))
        self.assertIsNone(result['coverage_threshold'])
        changed=copy.deepcopy(self.packet);changed['target']='Another source-bound contrast label'
        self.assertEqual(self.check(packet=changed)['status'],'PASS')

    def test_missing_provenance_and_legacy_mode_abstain_without_error_allegations(self):
        for mutation in ('missing_file','missing_declaration','legacy'):
            clocks=copy.deepcopy(self.clocks);check=clocks['checks'][0]
            if mutation=='missing_file': check['measurement_sources']['sample_source']['path']=str(self.root/'absent.soft')
            if mutation=='missing_declaration': del check['measurement_sources']['sample_source']
            if mutation=='legacy': check['normalization_profile']='clockbase_legacy'
            result=self.check(clocks=clocks)
            self.assertEqual(result['status'],'CANNOT VERIFY',result)
            self.assertFalse(result['failures'])

    def test_primary_source_change_is_detected(self):
        before=self.primary.read_bytes()
        try:
            self.primary.write_bytes(before+b'\n')
            result=self.check()
            self.assertEqual(result['status'],'FAIL')
            self.assertIn('identity differs',result['failures'][0])
        finally:self.primary.write_bytes(before)

    def test_raw_counts_change_is_detected(self):
        before=self.raw.read_bytes()
        try:
            data=pd.read_csv(self.raw,index_col=0);data.iloc[0,0]+=1;data.to_csv(self.raw)
            self.assertEqual(self.check()['status'],'FAIL')
        finally:self.raw.write_bytes(before)

    def test_rehashed_wrong_prediction_is_checked_from_counts(self):
        clocks=copy.deepcopy(self.clocks);path=Path(clocks['path']);before=path.read_bytes()
        try:
            frame=pd.read_csv(path);frame.loc[0,self.clock]+=1;frame.to_csv(path,index=False)
            clocks['prediction_sha256']=sha(path)
            result=self.check(clocks=clocks)
            self.assertEqual(result['status'],'FAIL')
            self.assertIn('numerical values differ',result['failures'][0])
        finally:path.write_bytes(before)

    def test_rehashed_feature_mask_corruption_is_detected(self):
        clocks=copy.deepcopy(self.clocks);declaration=clocks['checks'][0]['measurement_sources']['feature_coverage']
        path=Path(declaration['path']);before=path.read_bytes()
        try:
            coverage=json.loads(before);coverage['imputed_ids_by_sample'][self.ids[0]].pop()
            path.write_text(json.dumps(coverage));declaration['sha256']=sha(path)
            result=self.check(clocks=clocks)
            self.assertEqual(result['status'],'FAIL')
            self.assertIn('coverage',result['failures'][0])
        finally:path.write_bytes(before)

    def test_domain_requires_positive_evidence_for_every_reference_sample(self):
        from audit_system.adapters.measurement_pilot import _rna_primary_samples,_rna_domain
        samples=_rna_primary_samples(self.primary)
        self.assertFalse(_rna_domain(samples)[0])
        samples[self.ids[-1]]['library_strategy']='OTHER'
        self.assertTrue(any(self.ids[-1] in value for value in _rna_domain(samples)[0]))
        samples=_rna_primary_samples(self.primary)
        samples[self.ids[-1]]['extract_protocol_ch1']='Adipose tissue total RNA'
        self.assertTrue(_rna_domain(samples)[0])

    def test_measurement_pass_does_not_override_missing_design(self):
        with tempfile.TemporaryDirectory() as folder:
            frame=pd.read_csv(self.clocks['path']).set_index('sample_id')
            row=recompute_declared_contrast(frame,self.ids[:2],self.ids[2:],self.clock,'Welch');row['status']='RESULT'
            design={'proposed_action':'obtain_source_evidence','status':'CANNOT VERIFY','pairing':'unknown',
                    'case_ids':[],'control_ids':[],'findings':[],'unresolved':['Independent biological units unspecified']}
            result=numerical_audit(design,{'rows':[row]},self.clocks,Path(folder),packet=self.packet)
            self.assertEqual(result['checks'][0]['primary_measurement_verification']['status'],'PASS')
            self.assertEqual(result['accepted_results'],0)
            self.assertEqual(result['scientific_status'],'CANNOT VERIFY')
            entries=[json.loads(line) for line in (Path(folder)/'decision-ledger.jsonl').read_text().splitlines()]
            counts=[entry for entry in entries if entry.get('uri')==str(self.raw.resolve())]
            self.assertEqual(len(counts),1)
            self.assertEqual(counts[0]['classification'],'EMPIRICAL')

    def test_rehashed_absent_feature_zero_padding_is_detected(self):
        clocks=copy.deepcopy(self.clocks);sources=clocks['checks'][0]['measurement_sources']
        transformed=Path(sources['transformed_input']['path']);receipt=Path(sources['preprocessing']['path'])
        originals={path:path.read_bytes() for path in (transformed,receipt)}
        try:
            frame=pd.read_csv(transformed,index_col=0)
            absent=next(column for column in frame if frame[column].isna().all())
            frame.loc[:,absent]=0;frame.to_csv(transformed)
            sources['transformed_input']['sha256']=sha(transformed)
            body=json.loads(receipt.read_text());body['output_sha256']['scaled_diff.csv']=sha(transformed)
            receipt.write_text(json.dumps(body));sources['preprocessing']['sha256']=sha(receipt)
            result=self.check(clocks=clocks)
            self.assertEqual(result['status'],'FAIL')
            self.assertIn('missing-feature mask differs',result['failures'][0])
        finally:
            for path,value in originals.items():path.write_bytes(value)

    def test_implementation_snapshots_survive_an_original_path_edit(self):
        from src.pipeline_inputs_2026_09_10 import snapshot_rna_implementations
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'source.py';path.write_text('# software fixture version one\n')
            records,originals=snapshot_rna_implementations({'software_fixture':path},root)
            path.write_text('# software fixture version two\n')
            archived=records['software_fixture']
            self.assertNotEqual(sha(path),archived['sha256'])
            self.assertEqual(sha(archived['path']),archived['sha256'])
            self.assertNotEqual(archived['path'],archived['executed_path'])

    def test_omitted_source_or_prediction_code_binding_cannot_pass(self):
        for missing in (None, {}, {'path': None, 'sha256': None}):
            clocks=copy.deepcopy(self.clocks);clocks['prediction_implementation_sources']['producer']=missing
            result=self.check(clocks=clocks)
            self.assertEqual(result['status'],'CANNOT VERIFY')
            self.assertFalse(result['failures'])
        clocks=copy.deepcopy(self.clocks)
        clocks['implementation_sha256']={str(self.raw):sha(self.raw)}
        result=self.check(clocks=clocks)
        self.assertEqual(result['status'],'CANNOT VERIFY')
        self.assertFalse(result['failures'])
        clocks=copy.deepcopy(self.clocks);source=clocks['checks'][0]['measurement_sources']['preprocessing']
        path=Path(source['path']);before=path.read_bytes()
        try:
            body=json.loads(before);del body['measurement_sources']['raw_counts']
            path.write_text(json.dumps(body));source['sha256']=sha(path)
            result=self.check(clocks=clocks)
            self.assertEqual(result['status'],'CANNOT VERIFY')
            self.assertFalse(result['failures'])
        finally:path.write_bytes(before)

    def test_missing_executed_prediction_column_is_an_axis_failure(self):
        clocks=copy.deepcopy(self.clocks);path=Path(clocks['path']);before=path.read_bytes()
        try:
            frame=pd.read_csv(path).drop(columns=self.clock);frame.to_csv(path,index=False)
            clocks['prediction_sha256']=sha(path)
            result=self.check(clocks=clocks)
            self.assertEqual(result['status'],'FAIL')
            self.assertIn('clock axis',result['failures'][0])
        finally:path.write_bytes(before)

    def test_changed_bound_implementation_snapshot_fails_identity(self):
        path=Path(self.clocks['checks'][0]['measurement_sources']['preprocessing_python']['path'])
        before=path.read_bytes()
        try:
            path.write_bytes(before+b'\n# changed bound snapshot\n')
            self.assertEqual(self.check()['status'],'FAIL')
        finally:path.write_bytes(before)

    def test_empty_bound_raw_file_is_reported_without_uncaught_exception(self):
        clocks=copy.deepcopy(self.clocks);before=self.raw.read_bytes()
        try:
            self.raw.write_bytes(b'');clocks['checks'][0]['measurement_sources']['raw_counts']['sha256']=sha(self.raw)
            result=self.check(clocks=clocks)
            self.assertEqual(result['status'],'FAIL')
            self.assertIn('absent or duplicate column',result['failures'][0])
        finally:self.raw.write_bytes(before)

    def test_wholly_imputed_score_is_unresolved_without_false_numerical_error(self):
        import joblib
        import numpy as np
        from src.pipeline_inputs_2026_09_10 import ASSETS,preprocess_rna
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            model=joblib.load(ASSETS/'EN_Hazard_Mouse_All_All_scaleddiff.pkl')
            features=list(map(str,model.feature_names));registry=set((ASSETS/'Gene_list_all_4.6.txt').read_text().splitlines())
            mapping=pd.read_csv(ASSETS/'Gene_table_mouse.csv',index_col=0).dropna(subset=['Entrez'])
            unmeasured_registry=registry-set(features)
            other=mapping[mapping['Entrez'].map(lambda value:str(int(value)) in unmeasured_registry)].drop_duplicates(subset='Entrez').iloc[:4]
            raw=root/'software-fixture-unmeasured-counts.csv'
            pd.DataFrame([[20,30,40,50],[40,30,20,10],[10,12,20,24],[50,20,60,45]],
                         index=other.index,columns=self.ids).to_csv(raw)
            pre=root/'preprocessing';pre.mkdir()
            processed=preprocess_rna(raw,'mouse',pre,normalization='published_mouse_rle',source=self.packet)
            imputer,scaler,selector,estimator=[step for _,step in model.steps]
            self.assertEqual(len(estimator.coef_),len(features))
            score=float((imputer.statistics_-scaler.mean_)@estimator.coef_+estimator.intercept_)
            path=root/'software-fixture-wholly-imputed-scores.csv'
            pd.DataFrame({'sample_id':self.ids,self.clock:[score]*len(self.ids)}).to_csv(path,index=False)
            coverage={'feature_ids':features,'training_retained_ids':features,'training_dropped_ids':[],
                'selected_ids':features,'nonzero_coefficient_ids':np.asarray(features)[estimator.coef_!=0].tolist(),
                'imputation_strategy':'median','missing_ids_by_sample':dict.fromkeys(self.ids,features),
                'imputed_ids_by_sample':dict.fromkeys(self.ids,features),'observed_selected_by_sample':dict.fromkeys(self.ids,0),
                'observed_nonzero_by_sample':dict.fromkeys(self.ids,0)}
            coverage_path=root/'software-fixture-feature-coverage.json';coverage_path.write_text(json.dumps(coverage))
            clocks=copy.deepcopy(self.clocks);clocks.update(path=str(path),prediction_sha256=sha(path),clocks=[self.clock])
            check=clocks['checks'][0];clocks['checks']=[check]
            check['observed_features_min']=0
            check['measurement_sources']={**processed['measurement_sources'],'model':check['measurement_sources']['model'],
                'preprocessing':{'path':str(pre/'preprocessing.json'),'sha256':sha(pre/'preprocessing.json')},
                'transformed_input':{'path':processed['scaled_diff'],'sha256':sha(processed['scaled_diff'])},
                'feature_coverage':{'path':str(coverage_path),'sha256':sha(coverage_path)}}
            result=self.check(clocks=clocks)
            self.assertEqual(result['status'],'CANNOT VERIFY',result)
            self.assertFalse(result['failures'])
            self.assertTrue(any('wholly imputed' in reason for reason in result['unresolved']))
            self.assertEqual(result['measured_input_prediction_check'],'CANNOT VERIFY')
            self.assertEqual(result['raw_to_score_verification'],'PASS')


class CodingActionValidationRetryTests(unittest.TestCase):
    """Disclosed software fixtures modeled on three observed action-format shapes."""

    class Client:
        def __init__(self, responses): self.responses, self.calls, self.messages = list(responses), 0, []
        def invoke(self, messages):
            self.calls += 1; self.messages.append(list(messages))
            response = self.responses.pop(0)
            if isinstance(response, Exception): raise response
            return AIMessage(content=response)

    def agent(self, responses, executed):
        client = self.Client(responses)
        validator = CodingActionValidationClient(client, 1)
        agent = CodingAgent(llm=validator, executor=lambda code: executed.append(code) or 'software fixture',
                            strict_tags=True, printout=False)
        agent.configure()
        return client, validator, agent

    def test_mixed_execute_solution_is_not_executed_and_valid_retry_recovers(self):
        rejected = '<execute>print(1)</execute>\n<solution>also done</solution>'
        executed = []
        client, validator, agent = self.agent([rejected, '<solution>done</solution>'], executed)
        _, final = agent.go('software fixture')
        self.assertEqual(final, '<solution>done</solution>')
        self.assertEqual(executed, [])
        self.assertEqual(client.calls, 2)
        self.assertEqual(validator.events[0]['rejected_response'], rejected)
        self.assertEqual(validator.events[0]['retry_status'], 'FORMAT_PASS')
        self.assertIn('FORMAT_REJECTED', client.messages[1][-1].content)

    def test_valid_execute_then_mixed_block_retries_without_reexecuting_rejected_code(self):
        executed = []
        mixed = '<execute>raise AssertionError("must not run")</execute><solution>mixed</solution>'
        _, validator, agent = self.agent([
            '<execute>print("first")</execute>', mixed, '<solution>done</solution>'], executed)
        _, final = agent.go('software fixture')
        self.assertEqual(executed, ['print("first")'])
        self.assertEqual(final, '<solution>done</solution>')
        self.assertEqual(validator.events[0]['rejected_response'], mixed)

    def test_missing_tags_then_repeated_malformed_stops_closed_without_execution(self):
        executed = []
        client, validator, agent = self.agent([
            'progress note without an action tag',
            '<execute>print(1)</execute><solution>mixed</solution>'], executed)
        with self.assertRaisesRegex(ValueError, 'exactly one complete'):
            agent.go('software fixture')
        self.assertEqual(executed, [])
        self.assertEqual(client.calls, 2)
        self.assertEqual(validator.events[0]['retry_status'], 'FORMAT_REJECTED')
        self.assertFalse(validator.events[0]['rejected_response_executed'])
        self.assertFalse(validator.events[0]['retry_response_executed'])

    def test_retry_allowance_is_global_across_the_full_agent_stage(self):
        executed = []
        client, validator, agent = self.agent([
            'first malformed response', '<execute>print("accepted")</execute>',
            '<execute>raise AssertionError("must not run")</execute><solution>mixed</solution>'], executed)
        with self.assertRaisesRegex(ValueError, 'exactly one complete'):
            agent.go('software fixture')
        self.assertEqual(client.calls, 3)
        self.assertEqual(validator.retries_used, 1)
        self.assertEqual(len(validator.events), 2)
        self.assertEqual(validator.events[1]['retry_status'], 'RETRY_ALLOWANCE_EXHAUSTED')
        self.assertIsNone(validator.events[1]['controller_feedback'])
        self.assertIn('must not run', validator.events[1]['rejected_response'])
        self.assertEqual(executed, ['print("accepted")'])

    def test_transport_error_is_not_retried(self):
        client = self.Client([RuntimeError('software fixture HTTP failure')])
        validator = CodingActionValidationClient(client, 1)
        with self.assertRaisesRegex(RuntimeError, 'HTTP failure'):
            validator.invoke([HumanMessage(content='software fixture')])
        self.assertEqual(client.calls, 1)
        self.assertEqual(validator.events, [])

    def test_budget_failure_on_format_retry_is_preserved_and_not_retried(self):
        client = self.Client(['no action tag', RuntimeError('Run API call/cost limit reached')])
        validator = CodingActionValidationClient(client, 1)
        with self.assertRaisesRegex(RuntimeError, 'call/cost limit'):
            validator.invoke([HumanMessage(content='software fixture')])
        self.assertEqual(client.calls, 2)
        self.assertEqual(validator.events[0]['retry_status'], 'TRANSPORT_OR_BUDGET_FAILURE')
        self.assertIn('call/cost limit', validator.events[0]['retry_error'])

    def test_failed_primary_stage_trace_retains_both_invalid_responses_and_feedback(self):
        client = self.Client(['no action tag', '<execute>print(1)</execute><solution>mixed</solution>'])
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); predictions = root / 'predictions.csv'
            predictions.write_text('sample_id,clock\nGSM1,1\n')
            stage = root / 'stage'; stage.mkdir()
            with self.assertRaisesRegex(ValueError, 'exactly one complete'):
                primary_analysis({'study_id': 'GSE1', 'target': 'software fixture'},
                                 {'path': str(predictions), 'clocks': ['clock']}, client, stage,
                                 coding_validation_retries=1)
            receipt = json.loads((stage / 'trace.json').read_text())
            event = receipt['format_validation']['events'][0]
            self.assertEqual(event['rejected_response'], 'no action tag')
            self.assertIn('<solution>mixed</solution>', event['retry_response'])
            self.assertIn('exactly one complete', event['error'])
            self.assertEqual(receipt['request_check']['status'], 'CANNOT VERIFY')
            self.assertFalse(list((stage / 'workspace').glob('execution-*.json')))

    def test_default_is_zero_and_invalid_direct_values_are_rejected(self):
        args = parser().parse_args(['--study', 'GSE1', '--modality', 'human_rna', '--output', 'unused'])
        self.assertEqual(args.coding_validation_retries, 0)
        with self.assertRaisesRegex(ValueError, 'must be 0 or 1'):
            CodingActionValidationClient(object(), 2)


class CodingOutputCompletionTests(unittest.TestCase):
    """Controller software fixtures; no model calls or empirical statistics."""

    abstention_code = '''import csv, json
from pathlib import Path
metadata = json.loads(Path('metadata.json').read_text())
with Path('clock_predictions.csv').open() as handle:
    predictions = list(csv.DictReader(handle))
assert predictions[0]['sample_id'] == metadata['samples'][0]['sample_id']
manifest = {'schema': 'clockbase-primary-comparisons-v1',
            'comparisons': [{'comparison_name': 'unresolved', 'target': metadata['target'],
                             'case_ids': [], 'control_ids': []}],
            'unassigned_ids': [s['sample_id'] for s in metadata['samples']]}
Path('comparisons.json').write_text(json.dumps(manifest))
fields = ['comparison_name', 'clock', 'case_ids', 'control_ids', 'statistical_test',
          'p_value', 'mean_difference', 'effect_size', 'status', 'reason']
with Path('aging_analysis.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerow({'comparison_name': 'unresolved', 'clock': 'clock', 'case_ids': '[]',
                     'control_ids': '[]', 'status': 'CANNOT VERIFY',
                     'reason': 'Software fixture: insufficient source support'})
print(Path('comparisons.json').read_text(), Path('aging_analysis.csv').read_text())
'''

    class Client(CodingActionValidationRetryTests.Client):
        def invoke(self, messages):
            if self.responses and callable(self.responses[0]):
                self.responses[0] = self.responses[0]()
            return super().invoke(messages)

    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.workspace = self.root / 'workspace'; self.workspace.mkdir()
        self.packet = {'study_id': 'GSE1', 'target': 'Software fixture unresolved target',
                       'samples': [{'sample_id': 'GSM1'}]}
        (self.workspace / 'metadata.json').write_text(json.dumps(self.packet))
        self.predictions = self.workspace / 'clock_predictions.csv'
        self.predictions.write_text('sample_id,clock\nGSM1,1\n')
        self.guard = CodingOutputCompletionGuard(IsolatedExecutor(self.workspace))

    def validator(self, responses, retries=1):
        client = self.Client(responses)
        return client, CodingActionValidationClient(client, retries, completion_guard=self.guard)

    @unittest.skipUnless(sys.platform.startswith('linux') and shutil.which('bwrap'),
                         'Agent execution requires Linux and bubblewrap on PATH')
    def test_premature_solution_corrects_through_execution_and_retains_abstention(self):
        stage = self.root / 'stage'; stage.mkdir()
        client = self.Client(['<solution>Files were inaccessible; no outputs.</solution>',
                             '<execute>' + self.abstention_code + '</execute>',
                             '<solution>Saved an unresolved analysis.</solution>'])
        result = primary_analysis(self.packet, {'path': str(self.predictions), 'clocks': ['clock']},
                                  client, stage, coding_validation_retries=1,
                                  coding_completion_policy='outputs_v1')
        self.assertEqual(client.calls, 3)
        self.assertEqual(result['rows'][0]['status'], 'CANNOT VERIFY')
        self.assertEqual(result['rows'][0]['p_value'], '')
        self.assertEqual(result['rows'][0]['case_ids'], [])
        self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
        trace = json.loads((stage / 'trace.json').read_text())
        self.assertEqual(trace['request_check']['status'], 'PASS')
        self.assertEqual(trace['format_validation']['retries_used'], 1)
        self.assertIn('OUTPUT_COMPLETION_REJECTED', client.messages[1][-1].content)
        self.assertIn('CANNOT VERIFY is valid', client.messages[1][-1].content)
        completion = trace['completion_validation']
        self.assertEqual([c['status'] for c in completion['solution_checks']], ['CANNOT VERIFY', 'PASS'])
        self.assertEqual(len(completion['executions']), 1)
        for name, expected in trace['input_sha256'].items():
            self.assertEqual(sha(stage / 'workspace' / name), expected)

    @unittest.skipUnless(sys.platform.startswith('linux') and shutil.which('bwrap'),
                         'Agent execution requires Linux and bubblewrap on PATH')
    def test_completed_solution_passes_without_correction(self):
        self.guard(self.abstention_code)
        client, validator = self.validator(['<solution>CANNOT VERIFY retained.</solution>'])
        response = validator.invoke([HumanMessage(content='Software fixture')])
        self.assertEqual(response.content, '<solution>CANNOT VERIFY retained.</solution>')
        self.assertEqual(client.calls, 1)
        self.assertEqual(validator.retries_used, 0)
        self.assertEqual(validator.events, [])

    def test_output_files_and_a_model_written_receipt_do_not_establish_host_execution(self):
        for name in self.guard.output_names:
            (self.workspace / name).write_text('Software fixture output')
        (self.workspace / 'execution-000.json').write_text('{"software_fixture":true}')
        client, validator = self.validator(['<solution>Done.</solution>'], retries=0)
        with self.assertRaisesRegex(ValueError, 'execution/output completion'):
            validator.invoke([HumanMessage(content='Software fixture')])
        self.assertEqual(client.calls, 1)
        self.assertEqual(self.guard.checks[-1]['host_execution_receipts'], 0)
        self.assertEqual(self.guard.checks[-1]['missing_output_files'], [])

    def test_repeated_premature_solutions_stop_without_automatic_output_creation(self):
        client, validator = self.validator(['<solution>Done.</solution>'] * 2)
        with self.assertRaisesRegex(ValueError, 'shared correction allowance exhausted'):
            validator.invoke([HumanMessage(content='Software fixture')])
        self.assertEqual(client.calls, 2)
        self.assertEqual(validator.retries_used, 1)
        self.assertEqual(validator.events[0]['retry_status'], 'COMPLETION_REJECTED')
        self.assertTrue(all(not (self.workspace / name).exists() for name in self.guard.output_names))
        self.assertFalse(list(self.workspace.glob('execution-*.json')))

    @unittest.skipUnless(sys.platform.startswith('linux') and shutil.which('bwrap'),
                         'Agent execution requires Linux and bubblewrap on PATH')
    def test_missing_each_named_file_blocks_a_solution_after_execution(self):
        self.guard(self.abstention_code)
        for name in self.guard.output_names:
            with self.subTest(name=name):
                path = self.workspace / name; previous = path.read_bytes(); path.unlink()
                client, validator = self.validator(['<solution>Done.</solution>'], retries=0)
                with self.assertRaisesRegex(ValueError, name):
                    validator.invoke([HumanMessage(content='Software fixture')])
                self.assertEqual(client.calls, 1)
                path.write_bytes(previous)

    @unittest.skipUnless(sys.platform.startswith('linux') and shutil.which('bwrap'),
                         'Agent execution requires Linux and bubblewrap on PATH')
    def test_prior_format_retry_consumes_the_completion_allowance(self):
        client, validator = self.validator(['no action tags', '<execute>print("Software fixture")</execute>',
                                           '<solution>Done without outputs.</solution>'])
        agent = CodingAgent(llm=validator, executor=self.guard, strict_tags=True, printout=False)
        agent.configure()
        with self.assertRaisesRegex(ValueError, 'shared correction allowance exhausted'):
            agent.go('Software fixture')
        self.assertEqual(client.calls, 3)
        self.assertEqual(validator.retries_used, 1)
        self.assertEqual(len(self.guard.executions), 1)
        self.assertEqual(validator.events[-1]['retry_status'], 'RETRY_ALLOWANCE_EXHAUSTED')
        self.assertIsNone(validator.events[-1]['controller_feedback'])

    @unittest.skipUnless(sys.platform.startswith('linux') and shutil.which('bwrap'),
                         'Agent execution requires Linux and bubblewrap on PATH')
    def test_prior_completion_retry_consumes_the_format_allowance(self):
        client, validator = self.validator(['<solution>Premature.</solution>',
                                           '<execute>print("Software fixture")</execute>', 'no action tags'])
        agent = CodingAgent(llm=validator, executor=self.guard, strict_tags=True, printout=False)
        agent.configure()
        with self.assertRaisesRegex(ValueError, 'exactly one complete'):
            agent.go('Software fixture')
        self.assertEqual(client.calls, 3)
        self.assertEqual(validator.retries_used, 1)
        self.assertEqual(validator.events[-1]['retry_status'], 'RETRY_ALLOWANCE_EXHAUSTED')

    def test_source_change_is_fatal_before_the_next_call(self):
        self.predictions.write_text('Software fixture mutated source')
        client, validator = self.validator(['<solution>Do not call.</solution>'])
        with self.assertRaisesRegex(PipelineIntegrityError, 'source input changed'):
            validator.invoke([HumanMessage(content='Software fixture')])
        self.assertEqual(client.calls, 0)
        self.assertEqual(validator.retries_used, 0)

    def test_source_change_during_response_is_fatal_without_correction(self):
        def change_source():
            self.predictions.write_text('Software fixture mutated source')
            return '<solution>Premature.</solution>'
        client, validator = self.validator([change_source])
        with self.assertRaisesRegex(PipelineIntegrityError, 'source input changed'):
            validator.invoke([HumanMessage(content='Software fixture')])
        self.assertEqual(client.calls, 1)
        self.assertEqual(validator.retries_used, 0)

    def test_source_drift_on_correction_preserves_response_and_integrity_reason(self):
        def change_source():
            self.predictions.write_text('Software fixture mutated during correction')
            return '<solution>Response remains preserved.</solution>'
        client, validator = self.validator(['<solution>Premature.</solution>', change_source])
        with self.assertRaisesRegex(PipelineIntegrityError, 'source input changed'):
            validator.invoke([HumanMessage(content='Software fixture')])
        self.assertEqual(client.calls, 2)
        self.assertEqual(validator.retries_used, 1)
        self.assertEqual(validator.events[0]['retry_status'], 'INTEGRITY_FAILURE')
        self.assertEqual(validator.events[0]['retry_response'], '<solution>Response remains preserved.</solution>')
        self.assertIn('source input changed', validator.events[0]['retry_error'])

    @unittest.skipUnless(sys.platform.startswith('linux') and shutil.which('bwrap'),
                         'Agent execution requires Linux and bubblewrap on PATH')
    def test_changed_receipt_or_executed_code_is_fatal_without_a_model_call(self):
        self.guard(self.abstention_code)
        for name in self.guard.executions[0]['files']:
            with self.subTest(name=name):
                path = self.workspace / name; previous = path.read_bytes()
                path.write_bytes(previous + b'\nSoftware fixture mutation')
                client, validator = self.validator(['<solution>Do not call.</solution>'])
                with self.assertRaisesRegex(PipelineIntegrityError, 'receipt or code changed'):
                    validator.invoke([HumanMessage(content='Software fixture')])
                self.assertEqual(client.calls, 0)
                path.write_bytes(previous)

    @unittest.skipUnless(sys.platform.startswith('linux') and shutil.which('bwrap'),
                         'Agent execution requires Linux and bubblewrap on PATH')
    def test_script_cannot_replace_its_saved_code_before_host_receipt_binding(self):
        code = "from pathlib import Path\nPath(__file__).write_text('print(\"Software fixture changed code\")')\n"
        with self.assertRaisesRegex(PipelineIntegrityError, 'executed code differs'):
            self.guard(code)
        self.assertEqual(self.guard.executions, [])

    @unittest.skipUnless(sys.platform.startswith('linux') and shutil.which('bwrap'),
                         'Agent execution requires Linux and bubblewrap on PATH')
    def test_next_script_symlink_cannot_overwrite_an_outside_target(self):
        outside = self.root / 'outside-marker.txt'
        original = 'Software fixture: outside target must remain unchanged'
        outside.write_text(original)
        first = "from pathlib import Path\nPath('agent-step-001.py').symlink_to(" + repr(str(outside)) + ")\n"
        self.guard(first)
        self.assertTrue((self.workspace / 'agent-step-001.py').is_symlink())
        with self.assertRaisesRegex(PipelineIntegrityError, 'exists before host execution'):
            self.guard('print("Software fixture second action must not execute")')
        self.assertEqual(outside.read_text(), original)
        self.assertEqual(self.guard.executor.index, 1)
        self.assertEqual(len(self.guard.executions), 1)
        self.assertFalse((self.workspace / 'execution-001.json').exists())
        self.assertFalse((self.workspace / 'execution-001.log').exists())

    def test_broken_script_or_receipt_symlink_is_rejected_before_host_write(self):
        outside = self.root / 'missing-outside-target'
        for name in ('agent-step-000.py', 'execution-000.json'):
            with self.subTest(name=name):
                link = self.workspace / name; link.symlink_to(outside)
                with self.assertRaisesRegex(PipelineIntegrityError, 'exists before host execution'):
                    self.guard('print("Software fixture must not execute")')
                self.assertFalse(outside.exists())
                self.assertTrue(link.is_symlink())
                self.assertEqual(self.guard.executor.index, 0)
                self.assertEqual(self.guard.executions, [])
                link.unlink()

    def test_existing_regular_script_or_receipt_is_preserved_without_execution(self):
        for name in ('agent-step-000.py', 'execution-000.json'):
            with self.subTest(name=name):
                path = self.workspace / name
                original = 'Software fixture: preserve existing entry'
                path.write_text(original)
                with self.assertRaisesRegex(PipelineIntegrityError, 'exists before host execution'):
                    self.guard('print("Software fixture must not execute")')
                self.assertEqual(path.read_text(), original)
                self.assertEqual(self.guard.executor.index, 0)
                self.assertEqual(self.guard.executions, [])
                path.unlink()

    @unittest.skipUnless(sys.platform.startswith('linux') and shutil.which('bwrap'),
                         'Agent execution requires Linux and bubblewrap on PATH')
    def test_output_symlink_is_fatal_and_does_not_trigger_correction(self):
        self.guard(self.abstention_code)
        outside = self.root / 'outside.csv'; outside.write_text('Software fixture outside output')
        output = self.workspace / 'aging_analysis.csv'; output.unlink(); output.symlink_to(outside)
        client, validator = self.validator(['<solution>Done.</solution>'])
        with self.assertRaisesRegex(PipelineIntegrityError, 'output symlink'):
            validator.invoke([HumanMessage(content='Software fixture')])
        self.assertEqual(client.calls, 1)
        self.assertEqual(validator.retries_used, 0)
        self.assertEqual(outside.read_text(), 'Software fixture outside output')

    def test_transport_policy_budget_and_integrity_failures_never_trigger_another_call(self):
        errors = [TimeoutError('Software fixture transport failure'),
                  PermissionError('Software fixture provider policy failure'),
                  RuntimeError('Run API call/cost limit reached'),
                  PipelineIntegrityError('Software fixture model mismatch')]
        for failure in errors:
            for prior_rejection in (False, True):
                with self.subTest(failure=type(failure).__name__, prior_rejection=prior_rejection):
                    responses = (['<solution>Premature.</solution>'] if prior_rejection else []) + [failure]
                    client, validator = self.validator(responses)
                    with self.assertRaises(type(failure)):
                        validator.invoke([HumanMessage(content='Software fixture')])
                    self.assertEqual(client.calls, 2 if prior_rejection else 1)
                    self.assertEqual(validator.retries_used, int(prior_rejection))

    @unittest.skipUnless(sys.platform.startswith('linux') and shutil.which('bwrap'),
                         'Agent execution requires Linux and bubblewrap on PATH')
    def test_output_content_remains_subject_to_existing_validator(self):
        stage = self.root / 'invalid-stage'; stage.mkdir()
        code = self.abstention_code + "\nPath('aging_analysis.csv').write_text('invalid,software,fixture\\n')\n"
        client = self.Client(['<execute>' + code + '</execute>', '<solution>Saved.</solution>'])
        with self.assertRaisesRegex(ValueError, 'output schema differs'):
            primary_analysis(self.packet, {'path': str(self.predictions), 'clocks': ['clock']},
                             client, stage, coding_validation_retries=1, coding_completion_policy='outputs_v1')
        self.assertEqual(client.calls, 2)
        trace = json.loads((stage / 'trace.json').read_text())
        self.assertEqual(trace['completion_validation']['solution_checks'][-1]['status'], 'PASS')
        self.assertEqual(trace['format_validation']['retries_used'], 0)

    def test_legacy_solution_and_receipt_remain_unchanged(self):
        client = self.Client(['<solution>Legacy final response.</solution>'])
        validator = CodingActionValidationClient(client, 1)
        self.assertEqual(validator.invoke([]).content, '<solution>Legacy final response.</solution>')
        self.assertEqual(validator.receipt(), {'schema': 'clockbase-coding-action-validation-v1',
            'configured_retries': 1, 'retries_used': 0, 'events': [],
            'scope': 'Action-tag syntax only; no scientific content or output acceptance validation'})
        args = parser().parse_args(['--study', 'GSE1', '--modality', 'human_rna', '--output', 'unused'])
        self.assertEqual(args.coding_completion_policy, 'legacy')
        self.assertEqual(args.coding_validation_retries, 0)
        with self.assertRaisesRegex(ValueError, 'Unknown coding completion policy'):
            primary_analysis({}, {}, client, self.root / 'not-created', coding_completion_policy='unknown')



class IDATSourceAcquisitionTests(unittest.TestCase):
    """Software fixtures exercise ownership/resume, without scientific measurements."""
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.series=self.root/'series.soft'; self.samples=self.root/'samples.soft'
        self.series.write_text('^SERIES = GSE1234\n!Series_sample_id = GSM1234\n')
        self.sample_text=('^SAMPLE = GSM1234\n!Sample_platform_id = GPL23976\n'+''.join(
            '!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1nnn/GSM1234/suppl/'
            'GSM1234_software-fixture_'+ch+'.idat.gz\n' for ch in ('Grn','Red')))
        self.samples.write_text(self.sample_text)
        from src import pipeline_inputs_2026_09_10 as adapter
        self.adapter=adapter
        self.plan=self.root/'plan.json'
        self.adapter.prepare_dnam_idat_plan('GSE1234',self.series,self.samples,self.plan)
        self.calls=[]

    def download_fixture(self,url,destination,limit):
        self.calls.append(url)
        # Binary magic only: no empirical array and no claim that an IDAT decoder accepts it.
        destination.write_bytes(gzip.compress(b'IDAT software regression fixture'))
        write_new(Path(str(destination)+'.source.json'),{'url':url,'sha256':sha(destination),
            'size_bytes':destination.stat().st_size})
        return destination

    def test_literal_source_pairing_never_infers_biological_units(self):
        plan=json.loads(self.plan.read_text())
        self.assertFalse(plan['biological_units_inferred'])
        self.assertEqual(plan['scientific_status'],'CANNOT VERIFY')
        self.assertEqual(plan['pairs'][0]['sample_id'],'GSM1234')
        self.assertEqual(set(plan['pairs'][0]['channels']),{'Grn','Red'})
        self.assertTrue(plan['pairs'][0]['channels']['Grn']['https_url'].startswith('https://ftp.ncbi.nlm.nih.gov/'))

    def test_source_ownership_platform_and_pair_failures(self):
        broken=[self.sample_text.replace('/GSM1234/suppl/','/GSM9999/suppl/'),
                self.sample_text.replace('GPL23976','GPL13534'),
                self.sample_text.replace('_Red.idat.gz','_Grn.idat.gz'),
                self.sample_text.replace('GSM1234_software-fixture_Red','GSM1234_different_Red'),
                self.sample_text.replace('ftp.ncbi.nlm.nih.gov','other.example')]
        for text in broken:
            with self.subTest(text=text):
                self.samples.write_text(text)
                with self.assertRaises(ValueError):
                    self.adapter.dnam_idat_source_pairs('GSE1234',self.series,self.samples)

    def test_changed_source_or_profile_stops_before_network(self):
        self.samples.write_text(self.sample_text+'!Sample_title = Changed fixture\n')
        with patch.object(self.adapter,'fetch') as fetch_call, self.assertRaisesRegex(ValueError,'source changed'):
            self.adapter.acquire_dnam_idat(self.plan,self.root/'cache')
        fetch_call.assert_not_called()
        self.samples.write_text(self.sample_text)
        plan=json.loads(self.plan.read_text());plan['profile']['beta']['offset']=0
        self.plan.write_text(json.dumps(plan))
        with patch.object(self.adapter,'fetch') as fetch_call, self.assertRaisesRegex(ValueError,'profile'):
            self.adapter.acquire_dnam_idat(self.plan,self.root/'cache')
        fetch_call.assert_not_called()

    def test_completed_acquisition_resumes_without_download(self):
        with patch.object(self.adapter,'fetch',side_effect=self.download_fixture):
            first=self.adapter.acquire_dnam_idat(self.plan,self.root/'cache')
            second=self.adapter.acquire_dnam_idat(self.plan,self.root/'cache')
        self.assertEqual(first,second);self.assertEqual(len(self.calls),2)
        self.assertEqual(len(first['files']),8)
        self.assertEqual(first['new_clock_runs'],0)

    def test_changed_raw_bytes_and_missing_roles_reject_resume(self):
        with patch.object(self.adapter,'fetch',side_effect=self.download_fixture):
            value=self.adapter.acquire_dnam_idat(self.plan,self.root/'cache')
        path=Path(value['files']['raw:GSM1234:Grn']['path'])
        before=path.read_bytes();path.write_bytes(before+b'changed')
        with self.assertRaisesRegex(ValueError,'changed'):
            self.adapter.acquire_dnam_idat(self.plan,self.root/'cache')
        path.write_bytes(before)
        receipt=self.root/'cache/acquisition.json'
        value['files'].pop('raw:GSM1234:Red');receipt.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError,'required files'):
            self.adapter.acquire_dnam_idat(self.plan,self.root/'cache')

    def test_failed_download_is_retained_and_next_attempt_is_separate(self):
        def interrupted(url,destination,limit):
            Path(str(destination)+'.partial').write_text('software fixture interrupted download')
            raise RuntimeError('fixture transport stop')
        with patch.object(self.adapter,'fetch',side_effect=interrupted),self.assertRaisesRegex(RuntimeError,'transport'):
            self.adapter.acquire_dnam_idat(self.plan,self.root/'cache')
        partials=list((self.root/'cache').glob('**/*.partial'))
        self.assertEqual(len(partials),1)
        original=partials[0].read_bytes()
        with patch.object(self.adapter,'fetch',side_effect=self.download_fixture):
            value=self.adapter.acquire_dnam_idat(self.plan,self.root/'cache')
        self.assertEqual(partials[0].read_bytes(),original)
        self.assertIn('attempt-002',value['files']['compressed:GSM1234:Grn']['path'])

    def test_decompression_bound_preserves_failed_bytes(self):
        plan=json.loads(self.plan.read_text());plan['limits']['decompressed_bytes_per_file']=8
        self.plan.write_text(json.dumps(plan))
        with patch.object(self.adapter,'fetch',side_effect=self.download_fixture),self.assertRaisesRegex(ValueError,'decompression byte'):
            self.adapter.acquire_dnam_idat(self.plan,self.root/'cache')
        self.assertTrue(list((self.root/'cache').glob('uncompressed/*.partial')))
        self.assertFalse((self.root/'cache/acquisition.json').exists())

    def test_old_and_new_command_line_dispatch_remain_available(self):
        import subprocess,sys
        for command,expected in [([], '--series'),(['prepare-idat'],'--samples'),
                                 (['acquire-idat'],'--plan'),(['process-idat'],'--acquisition')]:
            with self.subTest(command=command):
                result=subprocess.run([sys.executable,'-B','-m','src.pipeline_inputs_2026_09_10',
                    *command,'--help'],capture_output=True,text=True,timeout=30)
                self.assertEqual(result.returncode,0,result.stderr)
                self.assertIn(expected,result.stdout)

    def test_existing_plan_is_not_overwritten(self):
        before=self.plan.read_bytes()
        with self.assertRaises(FileExistsError):
            self.adapter.prepare_dnam_idat_plan('GSE1234',self.series,self.samples,self.plan)
        self.assertEqual(before,self.plan.read_bytes())

    def test_supported_profile_is_explicit_and_strict_v1_remains_selectable(self):
        self.assertEqual(json.loads(self.plan.read_text())['profile'],self.adapter.DNAM_IDAT_SUPPORTED_PROFILE)
        strict=self.root/'strict-plan.json'
        self.adapter.prepare_dnam_idat_plan('GSE1234',self.series,self.samples,strict,
            profile_name='minfi_ssnoob_beta100_detection01_v1')
        self.assertEqual(json.loads(strict.read_text())['profile'],self.adapter.DNAM_IDAT_PROFILE)
        with self.assertRaisesRegex(ValueError,'Unknown IDAT'):
            self.adapter.prepare_dnam_idat_plan('GSE1234',self.series,self.samples,self.root/'unknown.json',profile_name='invented')
        self.assertFalse((self.root/'unknown.json').exists())

    def strict_acquisition_fixture(self):
        strict=self.root/'strict-plan.json'
        self.adapter.prepare_dnam_idat_plan('GSE1234',self.series,self.samples,strict,
            profile_name='minfi_ssnoob_beta100_detection01_v1')
        with patch.object(self.adapter,'fetch',side_effect=self.download_fixture):
            self.adapter.acquire_dnam_idat(strict,self.root/'strict-cache')
        return self.root/'strict-cache/acquisition.json'

    def test_new_profile_reuses_verified_raw_files_without_download_or_overwrite(self):
        prior=self.strict_acquisition_fixture();before=prior.read_bytes()
        with patch.object(self.adapter,'fetch') as fetch_call:
            value=self.adapter.acquire_dnam_idat(self.plan,self.root/'reuse-cache',reuse_acquisition=prior)
            resumed=self.adapter.acquire_dnam_idat(self.plan,self.root/'reuse-cache',reuse_acquisition=prior)
        fetch_call.assert_not_called()
        self.assertEqual(value,resumed);self.assertEqual(prior.read_bytes(),before)
        self.assertEqual(value['new_downloads'],0)
        self.assertEqual(value['reused_acquisition_sha256'],sha(prior))
        self.assertEqual(value['profile'],self.adapter.DNAM_IDAT_SUPPORTED_PROFILE)
        self.assertEqual(value['files'],json.loads(before)['files'])
        self.assertEqual(len(list((self.root/'reuse-cache').iterdir())),1)

    def test_reuse_rejects_changed_source_bytes_before_network(self):
        prior=self.strict_acquisition_fixture()
        old=json.loads(prior.read_text())
        Path(old['files']['raw:GSM1234:Grn']['path']).write_bytes(b'IDAT changed software fixture')
        with patch.object(self.adapter,'fetch') as fetch_call,self.assertRaisesRegex(ValueError,'changed'):
            self.adapter.acquire_dnam_idat(self.plan,self.root/'reuse-cache',reuse_acquisition=prior)
        fetch_call.assert_not_called()
        self.assertFalse((self.root/'reuse-cache/acquisition.json').exists())

    def test_reuse_enforces_new_resource_limits(self):
        prior=self.strict_acquisition_fixture()
        plan=json.loads(self.plan.read_text());plan['limits']['compressed_bytes_per_file']=1
        self.plan.write_text(json.dumps(plan))
        with patch.object(self.adapter,'fetch') as fetch_call,self.assertRaisesRegex(ValueError,'new plan bound'):
            self.adapter.acquire_dnam_idat(self.plan,self.root/'reuse-cache',reuse_acquisition=prior)
        fetch_call.assert_not_called()

    def test_reuse_rejects_a_different_source_snapshot(self):
        prior=self.strict_acquisition_fixture()
        other=self.root/'samples-new.soft';other.write_text(self.sample_text+'!Sample_title = New fixture snapshot\n')
        plan=self.root/'other-plan.json'
        self.adapter.prepare_dnam_idat_plan('GSE1234',self.series,other,plan)
        with patch.object(self.adapter,'fetch') as fetch_call,self.assertRaisesRegex(ValueError,'sources differ'):
            self.adapter.acquire_dnam_idat(plan,self.root/'reuse-cache',reuse_acquisition=prior)
        fetch_call.assert_not_called()

    def test_partial_cached_channels_obey_each_new_per_file_bound(self):
        for kind, limit_key in [('compressed', 'compressed_bytes_per_file'), ('raw', 'decompressed_bytes_per_file')]:
            with self.subTest(kind=kind):
                output = self.root / ('partial-' + kind)
                value = self.adapter.acquire_dnam_idat(self.plan, output, fetch_callback=self.download_fixture)
                (output/'acquisition.json').rename(self.root / ('preserved-' + kind + '.json'))
                before = {str(p):sha(p) for p in output.rglob('*') if p.is_file()}
                limits = dict(self.adapter.DNAM_IDAT_LIMITS)
                limits[limit_key] = Path(value['files'][kind + ':GSM1234:Grn']['path']).stat().st_size - 1
                smaller = self.root / ('smaller-' + kind + '.json')
                self.adapter.prepare_dnam_idat_plan('GSE1234', self.series, self.samples, smaller, limits=limits)
                with self.assertRaisesRegex(ValueError, 'Cached ' + kind + ' IDAT.*byte bound'):
                    self.adapter.acquire_dnam_idat(smaller, output,
                        fetch_callback=lambda *args, **kwargs: self.fail('Download on retained-channel rejection'))
                self.assertEqual(before, {str(p):sha(p) for p in output.rglob('*') if p.is_file()})
                self.assertFalse((output/'acquisition.json').exists())

    def test_completed_receipt_cannot_authorize_oversized_cached_files(self):
        for kind, key in [('compressed','compressed_bytes_per_file'), ('raw','decompressed_bytes_per_file')]:
            with self.subTest(kind=kind):
                output = self.root / ('completed-' + kind)
                value = self.adapter.acquire_dnam_idat(self.plan, output, fetch_callback=self.download_fixture)
                limits = dict(self.adapter.DNAM_IDAT_LIMITS)
                limits[key] = Path(value['files'][kind + ':GSM1234:Grn']['path']).stat().st_size - 1
                smaller = self.root / ('claimed-' + kind + '.json')
                self.adapter.prepare_dnam_idat_plan('GSE1234', self.series, self.samples, smaller, limits=limits)
                # Adversarial SOFTWARE FIXTURE: a receipt claims the smaller
                # frozen plan, while its hashed physical files exceed it.
                value.update(plan_path=str(smaller), plan_sha256=sha(smaller))
                (output/'acquisition.json').write_text(json.dumps(value))
                with self.assertRaisesRegex(ValueError, 'Completed ' + kind + ' IDAT.*byte bound'):
                    self.adapter.acquire_dnam_idat(smaller, output,
                        fetch_callback=lambda *args, **kwargs: self.fail('Download on completed-cache rejection'))

    def test_reuse_receipt_respects_output_cap_and_preserves_failed_bytes(self):
        prior = self.strict_acquisition_fixture(); before = prior.read_bytes()
        output = self.root / 'bounded-reuse'
        with patch.object(self.adapter, 'fetch', side_effect=AssertionError('Network during reuse')) as download:
            with self.assertRaisesRegex(ValueError, 'output.*byte cap'):
                self.adapter.acquire_dnam_idat(self.plan, output, reuse_acquisition=prior, max_output_bytes=32)
            retained = (output/'acquisition.json').read_bytes()
            self.assertGreater(len(retained), 32)
            with self.assertRaisesRegex(ValueError, 'output.*byte cap'):
                self.adapter.acquire_dnam_idat(self.plan, output, reuse_acquisition=prior, max_output_bytes=32)
        self.assertEqual((output/'acquisition.json').read_bytes(), retained)
        self.assertEqual(prior.read_bytes(), before); download.assert_not_called()

    def test_refresh_callback_owns_downloads_and_completed_acquisition_has_no_calls(self):
        with patch.object(self.adapter, 'fetch', side_effect=AssertionError('Unaccounted network')) as fallback:
            value = self.adapter.acquire_dnam_idat(self.plan, self.root/'callback-cache', fetch_callback=self.download_fixture)
            resumed = self.adapter.acquire_dnam_idat(self.plan, self.root/'callback-cache',
                fetch_callback=lambda *args, **kwargs: self.fail('Download on completed acquisition'))
        self.assertEqual(value, resumed); self.assertEqual(len(self.calls), 2); fallback.assert_not_called()

    def test_smaller_limits_keep_exact_frozen_keys_and_reject_expansion(self):
        limits = {key: min(value, 10000) for key, value in self.adapter.DNAM_IDAT_LIMITS.items()}
        limits['samples'] = 1
        path = self.root/'smaller-plan.json'
        self.adapter.prepare_dnam_idat_plan('GSE1234', self.series, self.samples, path, limits=limits)
        self.assertEqual(self.adapter._checked_idat_plan(path)['limits'], limits)
        for changed in ({**limits, 'extra': 1}, {**limits, 'samples': True}, {**limits, 'samples': 65}):
            with self.subTest(limits=changed), self.assertRaisesRegex(ValueError, 'resource limits'):
                self.adapter.prepare_dnam_idat_plan('GSE1234', self.series, self.samples,
                                                   self.root/'rejected-plan.json', limits=changed)
        self.assertFalse((self.root/'rejected-plan.json').exists())

    def test_retained_acquisition_byte_cap_stops_decompression_and_preserves_partial(self):
        def expanded_fixture(url, destination, limit):
            destination.write_bytes(gzip.compress(b'IDAT software fixture '+b'x'*20000))
            write_new(Path(str(destination)+'.source.json'), {'url':url, 'sha256':sha(destination),
                'size_bytes':destination.stat().st_size})
            return destination
        with patch.object(self.adapter, 'fetch', side_effect=AssertionError('Unaccounted network')):
            with self.assertRaisesRegex(ValueError, 'output.*byte cap'):
                self.adapter.acquire_dnam_idat(self.plan, self.root/'bounded-cache',
                    fetch_callback=expanded_fixture, max_output_bytes=1000)
        self.assertTrue(list((self.root/'bounded-cache').glob('uncompressed/*.partial')))
        self.assertFalse((self.root/'bounded-cache/acquisition.json').exists())

    def test_processing_completed_resume_never_calls_acquisition_or_worker(self):
        acquisition = self.root/'acquisition.json'; acquisition.write_text('{"software_fixture":true}')
        output = self.root/'completed-processing'; output.mkdir()
        matrix = output/'software-matrix.csv'; matrix.write_text('SOFTWARE FIXTURE; no methylation values\n')
        write_new(output/'provenance.json', {'files': {
            'protocol': {'sha256':sha(self.plan)}, 'acquisition': {'sha256':sha(acquisition)},
            'matrix': {'path':str(matrix), 'sha256':sha(matrix)}}})
        with patch.object(self.adapter, 'acquire_dnam_idat', side_effect=AssertionError('Acquisition on completed processing')):
            with patch.object(self.adapter, '_run_idat_worker', side_effect=AssertionError('Worker on completed processing')):
                with patch('audit_system.report_2026_09_10.verify_dnam_idat_provenance', return_value={'software_fixture':True}):
                    result = self.adapter.process_dnam_idat(self.plan, acquisition, output)
        self.assertTrue(result['resumed'])

    def test_processing_missing_acquisition_never_downloads(self):
        with patch.object(self.adapter, 'fetch', side_effect=AssertionError('Processing attempted network')) as download:
            with self.assertRaisesRegex(FileNotFoundError, 'completed IDAT acquisition'):
                self.adapter.process_dnam_idat(self.plan, self.root/'absent/acquisition.json', self.root/'processing')
        download.assert_not_called(); self.assertFalse((self.root/'processing').exists())

    def test_processing_output_cap_preserves_failure_before_any_worker(self):
        self.adapter.acquire_dnam_idat(self.plan, self.root/'cache', fetch_callback=self.download_fixture)
        with patch.object(self.adapter.subprocess, 'Popen', side_effect=AssertionError('Unexpected worker')) as launch:
            with self.assertRaisesRegex(ValueError, 'output.*byte cap'):
                self.adapter.process_dnam_idat(self.plan, self.root/'cache/acquisition.json',
                                              self.root/'processing', max_output_bytes=1)
        launch.assert_not_called()
        self.assertEqual(json.loads((self.root/'processing/failure.json').read_text())['status'], 'FAIL')

    def test_worker_output_and_time_guards_stop_software_child(self):
        import sys
        for mode in ('bytes', 'seconds'):
            with self.subTest(mode=mode):
                output=self.root/mode; output.mkdir()
                code='import time; time.sleep(10)'
                if mode=='bytes':
                    code='from pathlib import Path; Path('+repr(str(output/'software-only.bin'))+').write_bytes(b"x"*8192); '+code
                with (output/'stdout.txt').open('x') as stdout, (output/'stderr.txt').open('x') as stderr:
                    with self.assertRaises(ValueError if mode=='bytes' else TimeoutError):
                        self.adapter._run_idat_worker([sys.executable,'-c',code],stdout,stderr,
                            2 if mode=='bytes' else .05,output,max_output_bytes=1000)
                if mode=='bytes': self.assertEqual((output/'software-only.bin').stat().st_size,8192)

    def test_nested_refresh_deadline_cannot_extend_enclosing_preparation_budget(self):
        import signal, time
        from src.update_pipeline_2026_09_20 import source_deadline
        with self.assertRaises(TimeoutError):
            with source_deadline(.02, 'Software fixture preparation'):
                with source_deadline(1, 'Software fixture source'):
                    time.sleep(.04)
        self.assertEqual(signal.getitimer(signal.ITIMER_REAL)[0], 0)

    @unittest.skipUnless(os.path.isdir('/proc/self'), 'Process-group software fixture requires Linux procfs')
    def test_worker_guards_stop_descendants_after_the_leader_exits(self):
        import ctypes, signal, subprocess, sys, time
        libc = ctypes.CDLL(None, use_errno=True)
        prior = ctypes.c_int()
        self.assertEqual(libc.prctl(37, ctypes.byref(prior), 0, 0, 0), 0)
        self.assertEqual(libc.prctl(36, 1, 0, 0, 0), 0)
        try:
            for mode in ('seconds', 'bytes', 'leader_exit'):
                with self.subTest(mode=mode):
                    output = self.root / ('descendant-' + mode); output.mkdir()
                    pid_file, pulse, leader_file = output/'child.pid', output/'pulse.txt', output/'leader.pid'
                    child_code = ('import os,signal,time; from pathlib import Path; '
                        'signal.signal(signal.SIGTERM,signal.SIG_IGN); '
                        'Path('+repr(str(pid_file))+').write_text(str(os.getpid())); '
                        'p=Path('+repr(str(pulse))+'); end=time.monotonic()+10; '
                        + ('Path('+repr(str(output/'oversized.bin'))+').write_bytes(b"x"*8192); ' if mode=='bytes' else '')
                        + '\nwhile time.monotonic()<end:\n p.write_text(str(time.monotonic())); time.sleep(.02)\n')
                    parent_code = ('import os,subprocess,sys,time; from pathlib import Path; '
                        'Path('+repr(str(leader_file))+').write_text(str(os.getpid())); '
                        'subprocess.Popen([sys.executable,"-c",'+repr(child_code)+']); '
                        + ('\nwhile not Path('+repr(str(pulse))+').exists(): time.sleep(.01)\n' if mode=='leader_exit' else 'time.sleep(10)'))
                    try:
                        with (output/'stdout.txt').open('x') as stdout, (output/'stderr.txt').open('x') as stderr:
                            if mode == 'leader_exit':
                                result = self.adapter._run_idat_worker([sys.executable,'-c',parent_code], stdout,stderr,2,output,max_output_bytes=1000)
                                self.assertEqual(result.returncode, 0)
                            else:
                                with self.assertRaises(ValueError if mode=='bytes' else TimeoutError):
                                    self.adapter._run_idat_worker([sys.executable,'-c',parent_code],stdout,stderr,
                                        2 if mode=='bytes' else .5,output,max_output_bytes=1000)
                        self.assertTrue(pid_file.is_file()); self.assertTrue(pulse.is_file())
                        frozen = pulse.read_bytes(); time.sleep(.15)
                        self.assertEqual(pulse.read_bytes(), frozen, 'Owned descendant continued writing after guard return')
                        child = int(pid_file.read_text()); deadline = time.monotonic()+2
                        while time.monotonic() < deadline:
                            waited, status = os.waitpid(child, os.WNOHANG)
                            if waited:
                                self.assertTrue(os.WIFSIGNALED(status)); self.assertEqual(os.WTERMSIG(status), signal.SIGKILL)
                                break
                            time.sleep(.02)
                        else:
                            self.fail('Owned descendant was not killed/reaped')
                    finally:
                        if leader_file.exists():
                            try: os.killpg(int(leader_file.read_text()), signal.SIGKILL)
                            except ProcessLookupError: pass
                        if pid_file.exists():
                            try: os.waitpid(int(pid_file.read_text()), 0)
                            except ChildProcessError: pass
        finally:
            self.assertEqual(libc.prctl(36, prior.value, 0, 0, 0), 0)

    @unittest.skipUnless(os.path.isfile('/proc/self/stat'), 'Thread-state software fixture requires Linux procfs')
    def test_worker_cleanup_retains_zombie_leader_with_live_thread(self):
        """Software fixture: pthread_exit of the leader does not end its process."""
        import signal, subprocess, sys, time
        output = self.root / 'thread-leader-cleanup'; output.mkdir()
        pulse = output / 'software-pulse.bin'
        code = ('import ctypes, signal, threading, time; from pathlib import Path; '
            'signal.signal(signal.SIGTERM, signal.SIG_IGN); p=Path('+repr(str(pulse))+'); '
            'ready=threading.Event(); '
            '\ndef worker():\n'
            ' ready.set(); end=time.monotonic()+5\n'
            ' with p.open("ab", buffering=0) as handle:\n'
            '  while time.monotonic()<end:\n'
            '   handle.write(b"x"); time.sleep(.002)\n'
            'threading.Thread(target=worker).start(); ready.wait()\n'
            'libc=ctypes.CDLL(None); libc.pthread_exit.argtypes=[ctypes.c_void_p]; '
            'libc.pthread_exit.restype=None; libc.pthread_exit(None)\n')
        with (output/'stdout.txt').open('x') as stdout, (output/'stderr.txt').open('x') as stderr:
            process = subprocess.Popen([sys.executable, '-c', code], start_new_session=True,
                stdout=stdout, stderr=stderr, env={'PATH':'/usr/bin:/bin', 'PYTHONDONTWRITEBYTECODE':'1'})
            try:
                deadline = time.monotonic()+3
                while time.monotonic()<deadline:
                    fields = (Path('/proc')/str(process.pid)/'stat').read_text().rsplit(')',1)[1].split()
                    if fields[0]=='Z' and int(fields[17])>1 and pulse.exists(): break
                    self.assertIsNone(process.poll(), 'Threaded fixture exited before the checked state')
                    time.sleep(.005)
                else: self.fail('Threaded fixture did not establish a zombie leader with live threads')
                self.assertEqual(int(fields[2]), process.pid)
                before = pulse.stat().st_size; time.sleep(.03)
                self.assertGreater(pulse.stat().st_size, before, 'Software thread must still be live')
                self.assertIn(process.pid, self.adapter._idat_running_group_members(process.pid))
            finally:
                try: os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError: pass
                process.wait(timeout=3)
            self.assertEqual(self.adapter._idat_running_group_members(process.pid), [])



class SourceFieldReferenceTests(unittest.TestCase):
    """Labeled software mutations of archived responses, never new model evidence."""
    def setUp(self):
        from src import pipeline_2026_09_10 as pipeline
        self.pipeline = pipeline
        self.policy = pipeline.SOURCE_FIELD_REFERENCE_POLICY
        base = pipeline.ROOT / ('paper/revision_2026-09-07/evidence/'
            'continuous-audit-pipeline-2026-09-20/update-controller/'
            'dnam-xlsx-ingestion-followup-2026-09-23/raw-idat-adapter-2026-09-24/'
            'development-pipeline-001/stopped-run-001/archive/studies/human_dnam/GSE228867')
        require_development_artifact(base)
        versions = list(base.iterdir())
        self.assertEqual(len(versions), 1)
        source = versions[0] / 'pipeline/03_blinded_design/attempt-000'
        self.packet = json.loads((source / 'source-packet.json').read_text())
        self.raw = [(source / f'response-validation/{i:03d}-response.txt').read_text() for i in (0, 1)]
        self.small = {'samples': [{'sample_id': 'GSM1', 'title': '  HumanMethylatio850\nα\t',
                                   'characteristics_ch1.1': 'replicate: A'}],
                      'study_documents': [{'source_id': 'GEO:GSE1',
                                           'fields': {'overall_design': 'Same literal\nß'}}]}

    @staticmethod
    def refs(value):
        """Explicit fixture mutation selecting the field of each original citation."""
        if isinstance(value, dict):
            if set(value) in ({'sample_id', 'field', 'quote'}, {'source_id', 'field', 'quote'}):
                return {key: val for key, val in value.items() if key != 'quote'} | {'quote_from_field': True}
            return {key: SourceFieldReferenceTests.refs(val) for key, val in value.items()}
        if isinstance(value, list):
            return [SourceFieldReferenceTests.refs(val) for val in value]
        return value

    def normalize(self, value, packet=None):
        raw = value if isinstance(value, str) else json.dumps(value)
        return self.pipeline.source_contract_candidate(raw, self.small if packet is None else packet, self.policy)

    def test_full_field_literals_preserve_unicode_whitespace_and_numbered_fields(self):
        from audit_system.agent_workflow_2026_09_10 import _check_evidence
        entries = [{'sample_id': 'GSM1', 'field': name, 'quote_from_field': True}
                   for name in ('title', 'characteristics_ch1.1')]
        entries += [{'source_id': 'GEO:GSE1', 'field': 'overall_design', 'quote_from_field': True}]
        value = {'evidence': entries}; original = copy.deepcopy(value)
        normalized, changes = self.normalize(value)
        expected = [self.small['samples'][0]['title'], self.small['samples'][0]['characteristics_ch1.1'],
                    self.small['study_documents'][0]['fields']['overall_design']]
        self.assertEqual([item['quote'] for item in normalized['evidence']], expected)
        self.assertEqual(_check_evidence(normalized['evidence'], self.small)[1], {'GSM1'})
        self.assertEqual(value, original)
        import hashlib
        self.assertEqual([item['inserted_sha256'] for item in changes],
                         [hashlib.sha256(text.encode()).hexdigest() for text in expected])

    def test_malformed_ambiguous_or_nonexistent_references_are_rejected(self):
        good = {'sample_id': 'GSM1', 'field': 'title', 'quote_from_field': True}
        cases = [good | {'quote_from_field': v} for v in (False, 1, 'true', None)]
        cases += [good | {'quote': 'extra'}, good | {'source_id': 'GEO:GSE1'}, good | {'extra': 0},
                  good | {'sample_id': 'GSM2'}, good | {'sample_id': ['GSM1']},
                  good | {'field': 'Title'}, good | {'field': 'sample_id'}, good | {'field': ''},
                  good | {'field': ['title']}, {'source_id': 'GSM1', 'field': 'title', 'quote_from_field': True},
                  {'sample_selector': {'where': []}, 'field': 'title', 'quote_from_field': True}]
        for entry in cases:
            with self.subTest(entry=entry), self.assertRaises(ValueError):
                self.normalize({'evidence': [entry]})
        with self.assertRaises(ValueError): self.normalize(good)
        with self.assertRaises(ValueError): self.normalize({'evidence': [good | {'quote_from_field': True}], 'other': good})
        duplicate = '{"evidence":[{"sample_id":"GSM1","sample_id":"GSM1","field":"title","quote_from_field":true}]}'
        with self.assertRaisesRegex(ValueError, 'Duplicate'): self.normalize(duplicate)

    def test_source_duplicates_empty_and_nonstring_fields_are_rejected(self):
        for namespace in ('samples', 'study_documents'):
            packet = copy.deepcopy(self.small); packet[namespace] += copy.deepcopy(packet[namespace])
            with self.subTest(namespace=namespace), self.assertRaisesRegex(ValueError, 'unique'):
                self.normalize({'evidence': []}, packet)
        for literal in ('', ' \n', None, 1, ['text']):
            packet = copy.deepcopy(self.small); packet['samples'][0]['title'] = literal
            with self.subTest(literal=literal), self.assertRaisesRegex(ValueError, 'nonempty string'):
                self.normalize({'evidence': [{'sample_id': 'GSM1', 'field': 'title', 'quote_from_field': True}]}, packet)

    def test_visible_plain_compact_and_self_review_sources_match_exactly(self):
        from audit_system.source_packet import compact_source_packet
        source = copy.deepcopy(self.small)
        for transport in (source, compact_source_packet(source), {'source_packet': compact_source_packet(source),
                                                                  'preceding_analysis': {'rows': []}}):
            self.assertEqual(self.pipeline.source_contract_visible_packet(source, transport), source)
        for changed in ('hidden', 'altered'):
            transport = copy.deepcopy(source)
            if changed == 'hidden': del transport['samples'][0]['title']
            else: transport['study_documents'][0]['fields']['overall_design'] += 'x'
            with self.assertRaisesRegex(ValueError, 'differs'):
                self.pipeline.source_contract_visible_packet(source, transport)
        damaged = compact_source_packet(source); damaged['sample_table']['packet_sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'digest differs'):
            self.pipeline.source_contract_visible_packet(source, damaged)

    def test_expansion_and_saved_json_byte_limits_stop_before_admission(self):
        value = {'evidence': [{'sample_id': 'GSM1', 'field': 'title', 'quote_from_field': True}] * 2}
        for key in self.pipeline.SOURCE_FIELD_REFERENCE_LIMITS:
            with self.subTest(limit=key), patch.dict(self.pipeline.SOURCE_FIELD_REFERENCE_LIMITS, {key: 1}), \
                    self.assertRaisesRegex(ValueError, 'byte|bytes'):
                self.normalize(value)

    def test_original_bad_quote_and_missing_sample_evidence_remain_rejected(self):
        for index, message in ((0, 'condition sample lacks'), (1, 'Unsupported source quotation')):
            candidate, _ = self.normalize(self.raw[index], self.packet)
            with self.subTest(index=index), self.assertRaisesRegex(ValueError, message):
                parse_design(json.dumps(candidate), self.packet, review_profile='source_domains_v4')
        candidate, _ = self.normalize(self.refs(json.loads(self.raw[0])), self.packet)
        with self.assertRaisesRegex(ValueError, 'condition sample lacks'):
            parse_design(json.dumps(candidate), self.packet, review_profile='source_domains_v4')

    def test_default_and_older_policy_do_not_admit_field_references(self):
        raw = json.dumps(self.refs(json.loads(self.raw[1])))
        for policy in ('strict', 'conservative_status_v1'):
            candidate, changes = self.pipeline.source_contract_candidate(raw, self.packet, policy)
            self.assertTrue(any('quote_from_field' in e for e in candidate['comparisons'][0]['evidence']))
            with self.subTest(policy=policy), self.assertRaisesRegex(ValueError, 'Evidence schema mismatch'):
                parse_design(json.dumps(candidate), self.packet, review_profile='source_domains_v4')
        args = parser().parse_args(['--study', 'GSE1', '--modality', 'human_dnam', '--output', 'unused'])
        self.assertEqual(args.source_contract_policy, 'strict')

    def test_nested_references_leave_scientific_fields_and_uncertainty_unchanged(self):
        value = self.refs(json.loads(self.raw[1]))
        candidate, changes = self.normalize(value, self.packet)
        checked = parse_design(json.dumps(candidate), self.packet, review_profile='source_domains_v4')
        self.assertEqual(checked['scientific_status'], 'CANNOT VERIFY')
        self.assertTrue(changes)
        expected, _ = self.pipeline.conservative_source_status(json.dumps(value))
        self.assertEqual(self.refs(candidate), expected)
        self.assertTrue(any('.covariate_plan[' in x['path'] for x in changes))
        self.assertTrue(any('.audit_domains.' in x['path'] for x in changes))
        self.assertTrue(any('.findings[' in x['path'] for x in changes))
        direct, direct_changes = self.normalize(candidate, self.packet)
        self.assertEqual(direct, candidate); self.assertEqual(direct_changes, [])

    def test_assessment_retains_original_failed_and_normalized_responses_without_extra_calls(self):
        from audit_system.agent_workflow_2026_09_10 import write_new
        replies = [json.dumps(self.refs(json.loads(raw)), ensure_ascii=False) for raw in self.raw]
        calls = []
        def invoke(messages):
            calls.append(copy.deepcopy(messages)); return AIMessage(content=replies[len(calls) - 1])
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            result = self.pipeline.assess_design(self.packet, SimpleNamespace(invoke=invoke), root,
                review_profile='source_domains_v4', source_binding_profile='literal_source_v1',
                source_contract_policy=self.policy, validation_retries=1)
            self.assertEqual(len(calls), 2)
            feedback = json.loads(calls[1][-1].content)
            self.assertEqual(feedback['coverage_issue_count'], 1)
            self.assertEqual(feedback['citation_issue_count'], 0)
            self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
            self.assertEqual(json.loads((root/'model-design.json').read_text()), json.loads(replies[1]))
            for index in (0, 1):
                receipt = json.loads((root/f'contract-adjustments/{index:03d}-receipt.json').read_text())
                self.assertEqual((root/f'contract-adjustments/{index:03d}-original-response.txt').read_text(), replies[index])
                self.assertEqual(receipt['source_packet_sha256'], sha(root/'source-packet.json'))
                self.assertEqual(receipt['source_transport_sha256'], sha(root/'source-transport.json'))
                self.assertEqual(receipt['implementation_sha256'], sha(self.pipeline.__file__))
                self.assertFalse(receipt['accepted_result'])
                self.assertEqual(receipt['contract_status'], 'REJECTED' if index == 0 else 'PARTIALLY_ADMISSIBLE')
            self.assertEqual(json.loads((root/'source-contract-admissibility.json').read_text())['scientific_status'], 'CANNOT VERIFY')

    def test_rejected_reference_retains_original_and_failure_receipt(self):
        value = self.refs(json.loads(self.raw[1])); value['comparisons'][0]['evidence'][0]['quote'] = 'ambiguous'
        raw = json.dumps(value)
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            with self.assertRaisesRegex(ValueError, 'exclusive exact'):
                self.pipeline.assess_design(self.packet, SimpleNamespace(invoke=lambda _: AIMessage(content=raw)), root,
                    review_profile='source_domains_v4', source_binding_profile='literal_source_v1',
                    source_contract_policy=self.policy)
            receipt = json.loads((root/'contract-adjustments/000-receipt.json').read_text())
            self.assertEqual(receipt['contract_status'], 'REJECTED')
            self.assertEqual((root/'contract-adjustments/000-original-response.txt').read_text(), raw)
            self.assertFalse((root/'contract-adjustments/000-candidate.json').exists())
            self.assertFalse((root/'design.json').exists())

    def test_update_controller_requires_new_version_and_forwards_policy(self):
        from audit_system.tests.test_update_pipeline_2026_09_20 import UpdatePipelineTests
        from src import update_pipeline_2026_09_20 as update
        fixture = UpdatePipelineTests(methodName='runTest'); fixture.setUp(); self.addCleanup(fixture.doCleanups)
        plan = fixture.plan([fixture.study()])
        previous = fixture.execute(plan)['studies'][0]['version_id']
        plan['pipeline']['source_contract_policy'] = self.policy
        blocked = fixture.execute(plan)
        self.assertIn('supersedes_version', blocked['studies'][0]['error'])
        plan['studies'][0]['supersedes_version'] = previous
        current = fixture.execute(plan)['studies'][0]['version_id']
        self.assertNotEqual(current, previous)
        request = copy.deepcopy(fixture.calls[-1]); request['result_path'] = str(fixture.root/'field-policy-child.json')
        with patch('src.pipeline_2026_09_10.main', return_value={'status': 'SOFTWARE_FIXTURE'}) as main:
            update.child(request)
        argv = main.call_args.args[0]
        self.assertEqual(argv[argv.index('--source-contract-policy') + 1], self.policy)
        for invalid in (None, True, [], 'automatic_quote_repair'):
            plan['pipeline']['source_contract_policy'] = invalid
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, 'Unknown source contract policy'):
                update.validate_plan(plan)

    def make_replay_fixture(self):
        """Relocated and intentionally rewritten SOFTWARE FIXTURE; no API call."""
        from audit_system import agent_workflow_2026_09_10 as workflow, source_packet
        from audit_system.tests.test_update_pipeline_2026_09_20 import SourceDesignReplayTests
        from src import update_pipeline_2026_09_20 as update
        fixture = SourceDesignReplayTests(methodName='runTest'); fixture.setUp(); self.addCleanup(fixture.doCleanups)
        source = fixture.origin/'pipeline/03_blinded_design'
        api = fixture.origin/'pipeline/api/design-000'
        original = json.loads((api/'response.json').read_text())
        rewritten = json.dumps(self.refs(json.loads(original['choices'][0]['message']['content'])))
        destination = source/'field-reference-software-fixture'; destination.mkdir()
        captured = []
        def invoke(messages):
            captured.extend(messages); return AIMessage(content=rewritten)
        result = self.pipeline.assess_design(fixture.packet, SimpleNamespace(invoke=invoke), destination,
            review_profile='source_domains_v4', source_binding_profile='literal_source_v1', source_contract_policy=self.policy)
        config_path = fixture.origin/'pipeline/config.json'; config = json.loads(config_path.read_text())
        config['source_contract_policy'] = self.policy
        producing = {'audit_system/source_packet.py': sha(source_packet.__file__),
                     'audit_system/agent_workflow_2026_09_10.py': sha(workflow.__file__),
                     'src/pipeline_2026_09_10.py': sha(self.pipeline.__file__)}
        config['implementation'].update(producing)
        config_path.write_text(json.dumps(config))
        # This fixture regenerates diagnostics and policy with current code;
        # retain coherent provenance instead of claiming the archived producer.
        version_path = fixture.origin/'version.json'; version = json.loads(version_path.read_text())
        for declaration in (version, version['identity']):
            declaration['implementation'].update(producing)
            declaration['pipeline_options']['source_contract_policy'] = self.policy
        version['version_id'] = update.identity(version['identity'])
        version_path.write_text(json.dumps(version))
        (fixture.origin/'version.json.sha256').write_text(sha(version_path)+'\n')
        self.assertEqual(json.loads((destination/'domain-verification.json').read_text())['source_diagnostic_profile'],
                         source_packet.SOURCE_DIAGNOSTIC_CURRENT)
        request_path = api/'request.json'; request = json.loads(request_path.read_text())
        request['messages'] = [{'role': 'system' if i == 0 else 'user', 'content': msg.content}
                               for i, msg in enumerate(captured)]
        request_path.write_text(json.dumps(request))
        original['choices'][0]['message']['content'] = rewritten; (api/'response.json').write_text(json.dumps(original))
        receipt_path = api/'receipt.json'; receipt = json.loads(receipt_path.read_text())
        receipt['response_sha256'] = sha(api/'response.json'); receipt_path.write_text(json.dumps(receipt))
        completed = {'result': result, 'outputs': {str(p.relative_to(source)): sha(p)
                      for p in destination.rglob('*') if p.is_file()}, 'elapsed_seconds': 0}
        (source/'complete.json').write_text(json.dumps(completed))
        (source/'complete.sha256').write_text(sha(source/'complete.json')+'\n')
        fixture.parameters['source_contract_policy'] = self.policy
        fixture.freeze()
        return fixture, destination

    def test_replay_validates_field_reference_provenance_and_rejects_bound_mutations(self):
        fixture, source = self.make_replay_fixture()
        result, replay = fixture.replay()
        self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(json.loads((replay/'source-design-replay.json').read_text())['new_model_calls'], 0)
        path = source/'contract-adjustments/000-receipt.json'; original = path.read_text()
        for key in ('source_packet_sha256', 'source_transport_sha256', 'implementation_sha256',
                    'original_response_sha256', 'candidate_sha256', 'normalization_limits', 'schema',
                    'scope', 'contract_status', 'accepted_result'):
            with self.subTest(field=key):
                receipt = json.loads(original); receipt[key] = '0'*64
                path.write_text(json.dumps(receipt)); fixture.rebind_source_fixture()
                with self.assertRaisesRegex(ValueError, 'field reference .*differs'):
                    fixture.replay()
                path.write_text(original); fixture.rebind_source_fixture()


    def test_replay_requires_candidate_hash_and_exact_final_admissibility_alias(self):
        fixture, source = self.make_replay_fixture()
        path = source/'contract-adjustments/000-receipt.json'; original = path.read_text()
        receipt = json.loads(original); del receipt['candidate_sha256']
        path.write_text(json.dumps(receipt)); fixture.rebind_source_fixture()
        with self.assertRaisesRegex(ValueError, 'field reference .*differs'):
            fixture.replay()
        path.write_text(original); fixture.rebind_source_fixture()
        alias = source/'source-contract-admissibility.json'; alias_original = alias.read_text()
        for field, value in (('scientific_status', 'PASS'), ('accepted_result', True), ('contract_status', 'REJECTED')):
            receipt = json.loads(alias_original); receipt[field] = value
            alias.write_text(json.dumps(receipt)); fixture.rebind_source_fixture()
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, 'final admissibility differs'):
                fixture.replay()
            alias.write_text(alias_original); fixture.rebind_source_fixture()


class ServiceTierTransportTests(unittest.TestCase):
    """Intercepted software transport fixtures; no credentials or API calls."""
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.messages = [HumanMessage(content='Software fixture')]

    def client(self, root, **options):
        return BoundedOpenAI(root, 'test', **{
            'model':'gpt-6-luna', 'max_calls':2, 'max_cost':1, 'output_tokens':20,
            'prompt_bytes':10000, 'input_per_million':.2, 'output_per_million':.75,
            'require_model_match':True, **options})

    @staticmethod
    def response(**options):
        body = {'model':'gpt-6-luna', 'usage':{'prompt_tokens':10,'completion_tokens':3},
                'choices':[{'finish_reason':'stop','message':{'content':'<solution>fixture</solution>'}}],
                **options}
        return SimpleNamespace(returncode=0, stdout=json.dumps(body).encode()+b'\n200')

    def test_invalid_explicit_tiers_stop_before_directory_or_transport(self):
        from src.pipeline_runtime_2026_09_10 import openai_request
        with patch('src.pipeline_runtime_2026_09_10.subprocess.run') as transport:
            for i, value in enumerate((None, False, 1, [], {}, '', 'auto', 'priority', 'flex')):
                with self.subTest(value=value), self.assertRaisesRegex(ValueError,'service_tier must be default'):
                    self.client(self.root/str(i), service_tier=value)
                self.assertFalse((self.root/str(i)).exists())
                with self.assertRaisesRegex(ValueError,'service_tier must be default'):
                    openai_request('gpt-6-luna', [], 20, service_tier=value)
            transport.assert_not_called()

    def test_omitted_wire_bytes_and_receipt_schema_remain_unchanged(self):
        for json_mode in (False, True):
            root = self.root/str(json_mode)
            with patch.dict(os.environ, {'OPENAI_API_KEY':'software-fixture-key'}, clear=True), \
                    patch('src.pipeline_runtime_2026_09_10.subprocess.run', return_value=self.response()) as transport:
                self.assertIn('fixture', self.client(root,json_response=json_mode).invoke(self.messages).content)
                self.client(root,json_response=json_mode).invoke(self.messages)
                self.assertEqual(transport.call_count, 1)
            # Literal pre-option envelope and ordering, including its saved newline.
            expected = {'model':'gpt-6-luna','messages':[{'role':'user','content':'Software fixture'}],
                        'max_completion_tokens':20,'reasoning_effort':'low','store':False}
            if json_mode:
                expected['response_format'] = {'type':'json_object'}
            self.assertEqual((root/'test-000/request.json').read_bytes(),
                             (json.dumps(expected,indent=2,allow_nan=False)+'\n').encode())
            self.assertEqual(set(json.loads((root/'test-000/receipt.json').read_text())),
                {'usage','estimated_cost_usd','billed_cost_usd','response_sha256',
                 'requested_model','returned_model','elapsed_seconds'})
            self.assertFalse((root/'test-000/service-tier-error.json').exists())

    def test_cli_omission_and_source_and_producer_forwarding(self):
        argv = ['--study','GSE1','--modality','human_rna','--output',str(self.root),
                '--model','gpt-6-luna','--audit-model','gpt-6-luna',
                '--input-usd-per-million','.2','--output-usd-per-million','.75',
                '--audit-input-usd-per-million','.2','--audit-output-usd-per-million','.75',
                '--audit-profile','source_domains_v4','--source-contract-policy','conservative_status_v1']
        self.assertNotIn('service_tier', vars(parser().parse_args(argv)))
        args = parser().parse_args(argv+['--service-tier','default'])
        for actor in ('design', 'primary'):
            root = self.root/actor
            with patch.dict(os.environ, {'OPENAI_API_KEY':'software-fixture-key'}, clear=True), \
                    patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                          return_value=self.response(service_tier='default')) as transport:
                build_client(args,root,actor).invoke(self.messages)
                self.assertEqual(transport.call_count,1)
            request = json.loads((root/'api'/f'{actor}-000/request.json').read_text())
            self.assertEqual(request['service_tier'],'default')
            self.assertEqual(request['model'],'gpt-6-luna')
            self.assertEqual('response_format' in request,actor=='design')

    def test_live_tier_rejection_is_integrity_failure_before_content_or_retry(self):
        for i, options in enumerate(({}, {'service_tier':None}, {'service_tier':'priority'},
                                     {'service_tier':False}, {'service_tier':[]}, {'service_tier':{}})):
            root = self.root/str(i)
            client = self.client(root,service_tier='default')
            correction = CodingActionValidationClient(client,1)
            with patch.dict(os.environ, {'OPENAI_API_KEY':'software-fixture-key'}, clear=True), \
                    patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                          return_value=self.response(**options)) as transport, \
                    patch('langchain_core.messages.AIMessage',side_effect=AssertionError('No rejected content use')):
                with self.assertRaisesRegex(PipelineIntegrityError,'Returned service tier'):
                    correction.invoke(self.messages)
                with self.assertRaisesRegex(PipelineIntegrityError,'Prior service-tier rejection'):
                    self.client(root,service_tier='default').invoke(self.messages)
                self.assertEqual(transport.call_count,1)
            slot = root/'test-000'
            self.assertTrue((slot/'response.json').is_file())
            self.assertTrue((slot/'reservation.json').is_file())
            self.assertFalse((slot/'receipt.json').exists())
            diagnostic = json.loads((slot/'service-tier-error.json').read_text())
            self.assertEqual(diagnostic['category'],'OPENAI_SERVICE_TIER_MISMATCH')
            self.assertEqual(diagnostic['billing_status'],'uncertain')
            self.assertFalse(diagnostic['automatic_retry'])
            self.assertTrue(diagnostic['raw_response_retained'])
            self.assertEqual(diagnostic['request_sha256'],sha(slot/'request.json'))
            self.assertEqual(diagnostic['response_sha256'],sha(slot/'response.json'))
            self.assertEqual(correction.retries_used,0)
            self.assertEqual(correction.events,[])

    def test_cached_tier_rejection_prevents_content_and_keeps_child_reservation(self):
        for i, tier in enumerate((None, 'priority')):
            root = self.root/str(i)
            with patch.dict(os.environ, {'OPENAI_API_KEY':'software-fixture-key'}, clear=True), \
                    patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                          return_value=self.response(service_tier='default')):
                self.client(root,service_tier='default').invoke(self.messages)
            slot = root/'test-000'
            response = json.loads((slot/'response.json').read_text())
            if tier is None:
                del response['service_tier']
            else:
                response['service_tier'] = tier
            # Deliberate software tampering also rehashes the receipt; semantic checks must reject it.
            (slot/'response.json').write_text(json.dumps(response))
            receipt = json.loads((slot/'receipt.json').read_text())
            receipt['response_sha256'] = sha(slot/'response.json')
            (slot/'receipt.json').write_text(json.dumps(receipt))
            before = (slot/'receipt.json').read_bytes()
            with patch('src.pipeline_runtime_2026_09_10.subprocess.run') as transport, \
                    patch('langchain_core.messages.AIMessage',side_effect=AssertionError('No cached content use')):
                with self.assertRaisesRegex(PipelineIntegrityError,'Returned service tier'):
                    self.client(root,service_tier='default').invoke(self.messages)
                transport.assert_not_called()
            self.assertEqual((slot/'receipt.json').read_bytes(),before)
            self.assertTrue((slot/'service-tier-error.json').is_file())
            reserve = json.loads((slot/'reservation.json').read_text())['reserved_usd']
            another = self.client(root,service_tier='default',max_cost=reserve*1.01)
            another.actor = 'another'
            with patch.dict(os.environ, {'OPENAI_API_KEY':'software-fixture-key'}, clear=True), \
                    patch('src.pipeline_runtime_2026_09_10.subprocess.run') as transport:
                with self.assertRaisesRegex(RuntimeError,'call/cost limit'):
                    another.invoke(self.messages)
                transport.assert_not_called()

    def test_exact_model_remains_required_with_explicit_default(self):
        with patch.dict(os.environ, {'OPENAI_API_KEY':'software-fixture-key'}, clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                      return_value=self.response(model='gpt-5.6-luna',service_tier='default')) as transport, \
                patch('langchain_core.messages.AIMessage',side_effect=AssertionError('No mismatched model content')):
            for _ in range(2):
                with self.assertRaisesRegex(PipelineIntegrityError,'Returned model differs'):
                    self.client(self.root,service_tier='default').invoke(self.messages)
            self.assertEqual(transport.call_count,1)

    def test_interrupted_cache_without_receipt_validates_tier_without_retry(self):
        from src.pipeline_runtime_2026_09_10 import openai_request
        for i,tier in enumerate((None,'priority','default')):
            root=self.root/str(i);slot=root/'test-000';slot.mkdir(parents=True)
            request=openai_request('gpt-6-luna',[{'role':'user','content':'Software fixture'}],20,
                                   service_tier='default')
            write_new(slot/'request.json',request)
            write_new(slot/'reservation.json',{'reserved_usd':.01,'started_unix':1,
                      'model':'gpt-6-luna','automatic_retries':0,'concurrency':1})
            response={'model':'gpt-6-luna','choices':[{'message':{'content':'must not be consumed'}}]}
            if tier is not None:response['service_tier']=tier
            write_new(slot/'response.json',response)
            error=RuntimeError if tier=='default' else PipelineIntegrityError
            with self.subTest(tier=tier), patch('src.pipeline_runtime_2026_09_10.subprocess.run') as transport, \
                    patch('langchain_core.messages.AIMessage',side_effect=AssertionError('No unsettled content')):
                with self.assertRaises(error):
                    self.client(root,service_tier='default').invoke(self.messages)
                transport.assert_not_called()
            self.assertFalse((slot/'receipt.json').exists())
            self.assertEqual((slot/'service-tier-error.json').exists(),tier!='default')

    def test_tier_rejection_on_format_correction_stops_without_third_call(self):
        bad = self.response(service_tier='default',choices=[{'finish_reason':'stop','message':{'content':'invalid action'}}])
        correction = CodingActionValidationClient(self.client(self.root,service_tier='default'),1)
        with patch.dict(os.environ, {'OPENAI_API_KEY':'software-fixture-key'}, clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=[bad,self.response()]) as transport:
            with self.assertRaisesRegex(PipelineIntegrityError,'Returned service tier'):
                correction.invoke(self.messages)
            self.assertEqual(transport.call_count,2)
        self.assertEqual(correction.retries_used,1)
        self.assertTrue((self.root/'test-000/receipt.json').is_file())
        self.assertFalse((self.root/'test-001/receipt.json').exists())
        self.assertNotIn('retry_response',correction.events[0])

    def test_explicit_tier_does_not_expand_correction_call_budget(self):
        bad = self.response(service_tier='default',choices=[{'finish_reason':'stop','message':{'content':'invalid action'}}])
        correction = CodingActionValidationClient(self.client(self.root,service_tier='default',max_calls=1),1)
        with patch.dict(os.environ, {'OPENAI_API_KEY':'software-fixture-key'}, clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',return_value=bad) as transport:
            with self.assertRaisesRegex(RuntimeError,'call/cost limit'):
                correction.invoke(self.messages)
            self.assertEqual(transport.call_count,1)
        self.assertEqual(len(list(self.root.glob('*/reservation.json'))),1)


if __name__ == '__main__': unittest.main()
