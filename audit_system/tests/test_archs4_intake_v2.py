import csv
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from audit_system import archs4_update_2026_09_10 as intake
from audit_system.tests import require_development_artifact


class HardenedIntakeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.release = self.root / 'release.h5'
        with h5py.File(self.release, 'w') as h5:
            h5.create_dataset('data/expression', data=np.ones((2, 2), dtype=np.int32))
            h5.create_dataset('meta/genes/ensembl_gene', data=np.array([b'ENSG1', b'ENSG2']))
            h5.create_dataset('meta/samples/geo_accession', data=np.array([b'GSM1', b'GSM2']))
            h5.create_dataset('meta/samples/series_id', data=np.array([b'GSE1', b'GSE1']))
            h5.create_dataset('meta/samples/organism_ch1', data=np.array([b'Homo sapiens', b'Homo sapiens']))
        self.metadata = self.root / 'metadata.csv'
        with self.metadata.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['GSM', 'gse_id', 'organism_ch1', 'molecule_ch1', 'title', 'contact_email'])
            writer.writeheader()
            writer.writerows([
                {'GSM': 'GSM1', 'gse_id': 'GSE1', 'organism_ch1': 'Homo sapiens', 'molecule_ch1': 'total RNA', 'title': 'a', 'contact_email': 'private@example.org'},
                {'GSM': 'GSM2', 'gse_id': 'GSE1', 'organism_ch1': 'Homo sapiens', 'molecule_ch1': 'total RNA', 'title': 'b', 'contact_email': 'private@example.org'},
            ])
        self.context = self.root / 'context.json'
        self.context.write_text('{}')
        self.asset = self.root / 'dependency.txt'; self.asset.write_text('dependency')
        self.model = self.root / 'model.pkl'; self.model.write_bytes(b'model')
        self.module = Path(intake.__file__).resolve()

    @staticmethod
    def digest(path):
        h = hashlib.sha256()
        with Path(path).open('rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                h.update(block)
        return h.hexdigest()

    def declaration(self, path, **extra):
        value = {'path': str(Path(path).resolve()), 'sha256': self.digest(path), 'size_bytes': Path(path).stat().st_size}
        value.update(extra)
        return value

    def plan(self, output=None, max_bytes=100000):
        frozen = [
            {**self.declaration(self.module), 'kind': 'implementation'},
            {**self.declaration(self.asset), 'kind': 'dependency_asset'},
            {**self.declaration(self.model), 'kind': 'model_asset'},
        ]
        release = self.declaration(self.release, expression_shape=[2, 2], sample_count=2)
        return {
            'schema': intake.HARDENED_INTAKE_SCHEMA,
            'output': str(Path(output or (self.root / 'v2-run')).resolve()),
            'runtime': {'python': platform.python_version(),
                        'packages': {'numpy': importlib.metadata.version('numpy')}},
            'frozen_files': frozen,
            'releases': {'human': release, 'mouse': release},
            'limits': {'max_studies': 1, 'max_samples_per_study': 2,
                       'max_extracted_bytes': max_bytes, 'concurrency': 1,
                       'max_numeric_chars': 64, 'max_identifier_chars': 128},
            'studies': [{'study_id': 'GSE1', 'species': 'human',
                         'target': 'source design audit', 'normalization': 'clockbase_legacy',
                         'source_files': {'metadata': self.declaration(self.metadata),
                                          'context': self.declaration(self.context)}}],
        }

    def packet_stub(self, metadata, study_id, target, *, context=None):
        self.assertEqual(Path(metadata).resolve(), self.metadata.resolve())
        self.assertEqual(Path(context).resolve(), self.context.resolve())
        return {'study_id': study_id, 'target': target,
                'samples': [{'sample_id': 'GSM1', 'organism_ch1': 'Homo sapiens', 'molecule_ch1': 'total RNA', 'title': 'a'},
                            {'sample_id': 'GSM2', 'organism_ch1': 'Homo sapiens', 'molecule_ch1': 'total RNA', 'title': 'b'}],
                'withheld': ['contact information'], 'source_sha256': self.digest(metadata)}

    def extract_stub(self, hdf5, study, species, source_samples, output):
        path = Path(output) / 'counts.csv'
        with path.open('w', newline='') as handle:
            writer = csv.writer(handle); writer.writerow(['ensembl_gene', *source_samples]); writer.writerow(['ENSG1', 1, 2])
        return {'path': str(path), 'sha256': self.digest(path), 'sample_count': len(source_samples)}

    def test_tampered_release_rejected_by_full_sha(self):
        plan = self.plan(max_bytes=100000)
        self.release.write_bytes(self.release.read_bytes() + b'tamper')
        with self.assertRaisesRegex(ValueError, 'changed before hardened intake'):
            intake.run_hardened_intake(plan)

    def test_conservative_cap_rejected_before_extraction(self):
        plan = self.plan(max_bytes=1)
        with patch('src.pipeline_inputs_2026_09_10.extract_rna', side_effect=AssertionError('extraction must not run')):
            with self.assertRaisesRegex(ValueError, 'bound exceeds cap'):
                intake.run_hardened_intake(plan)

    def test_existing_v1_style_directory_is_rejected(self):
        v1 = self.root / 'run'; v1.mkdir(); (v1 / 'old-receipt.json').write_text('{}')
        plan = self.plan(output=v1, max_bytes=100000)
        with self.assertRaisesRegex(ValueError, 'new empty output directory'):
            intake.run_hardened_intake(plan)

    def test_failed_stage_is_terminal_and_resume_does_not_retry(self):
        plan = self.plan(max_bytes=100000)
        preprocess_calls = []

        def fail_preprocess(*args, **kwargs):
            preprocess_calls.append(1)
            raise ValueError('declared fixture failure')

        with patch('audit_system.agent_workflow_2026_09_10.source_packet', side_effect=self.packet_stub), \
             patch('src.pipeline_inputs_2026_09_10.extract_rna', side_effect=self.extract_stub), \
             patch('src.pipeline_inputs_2026_09_10.preprocess_rna', side_effect=fail_preprocess):
            first = intake.run_hardened_intake(plan)
            second = intake.run_hardened_intake(plan)
        self.assertEqual(first, second)
        self.assertEqual(len(preprocess_calls), 1)
        self.assertEqual(first['results'][0]['preprocess']['status'], 'CANNOT VERIFY')
        self.assertFalse(first['results'][0]['preprocess']['automatic_retry'])

    def test_source_packet_is_compact_model_view_without_contact_fields(self):
        plan = self.plan(max_bytes=100000)
        with patch('audit_system.agent_workflow_2026_09_10.source_packet', side_effect=self.packet_stub), \
             patch('src.pipeline_inputs_2026_09_10.extract_rna', side_effect=self.extract_stub), \
             patch('src.pipeline_inputs_2026_09_10.preprocess_rna', side_effect=ValueError('stop after packet')):
            result = intake.run_hardened_intake(plan)
        packet_path = Path(result['results'][0]['source_packet']['value']['path'])
        packet = json.loads(packet_path.read_text())
        self.assertEqual(packet['sample_table']['sample_ids'], ['GSM1', 'GSM2'])
        self.assertNotIn('contact_email', json.dumps(packet))
        self.assertEqual(packet['study_id'], 'GSE1')

    def test_reuse_matrix_is_labeled_and_skips_fresh_extraction(self):
        reused = self.root / 'verified-counts.csv'
        reused.write_text('ensembl_gene,GSM1,GSM2\nENSG1,1,2\n')
        plan = self.plan(max_bytes=100000)
        plan['studies'][0]['reuse_extraction'] = self.declaration(reused)
        with patch('audit_system.agent_workflow_2026_09_10.source_packet', side_effect=self.packet_stub), \
             patch('src.pipeline_inputs_2026_09_10.extract_rna', side_effect=AssertionError('fresh extraction must not run')), \
             patch('src.pipeline_inputs_2026_09_10.preprocess_rna', side_effect=ValueError('stop after reuse')):
            result = intake.run_hardened_intake(plan)
        value = result['results'][0]['extract']['value']
        self.assertEqual(value['input_mode'], 'REUSED_VERIFIED_INGESTED_MATRIX')
        self.assertEqual(value['sample_count'], 2)

    def test_reuse_matrix_byte_cap_is_checked_before_stage(self):
        reused = self.root / 'verified-counts.csv'
        reused.write_text('ensembl_gene,GSM1,GSM2\nENSG1,1,2\n')
        plan = self.plan(max_bytes=reused.stat().st_size - 1)
        plan['studies'][0]['reuse_extraction'] = self.declaration(reused)
        with self.assertRaisesRegex(ValueError, 'Reused extraction exceeds byte cap before stage'):
            intake.run_hardened_intake(plan)


class HardenedV21ProfileTest(HardenedIntakeTest):
    def strict_plan(self):
        plan = self.plan(max_bytes=100000)
        plan['schema'] = intake.HARDENED_INTAKE_SCHEMA_V21
        plan['limits']['max_numeric_chars'] = intake.HARDENED_V21_NUMERIC_CHARS
        plan['limits']['max_identifier_chars'] = intake.HARDENED_V21_IDENTIFIER_CHARS
        required = intake._intake_v21_expected_assets('human')
        frozen = []
        for name, kind in sorted(required):
            path = self.root / name
            path.write_text(name)
            frozen.append({**self.declaration(path), 'kind': kind})
        plan['frozen_files'] = frozen
        return plan

    def test_v21_rejects_duplicate_study_ids(self):
        plan = self.strict_plan()
        plan['limits']['max_studies'] = 2
        plan['studies'] = [plan['studies'][0], dict(plan['studies'][0])]
        with self.assertRaisesRegex(ValueError, 'study IDs must be unique'):
            intake.validate_hardened_intake_plan_v21(plan)

    def test_v21_rejects_weakened_conservative_character_bounds(self):
        plan = self.strict_plan()
        plan['limits']['max_numeric_chars'] = 1
        with self.assertRaisesRegex(ValueError, 'numeric bound cannot be weakened'):
            intake.validate_hardened_intake_plan_v21(plan)

    def test_v21_rejects_complex_expression_dtype(self):
        with h5py.File(self.release, 'w') as h5:
            h5.create_dataset('data/expression', data=np.ones((2, 2), dtype=np.complex64))
            h5.create_dataset('meta/genes/ensembl_gene', data=np.array([b'ENSG1', b'ENSG2']))
            h5.create_dataset('meta/samples/geo_accession', data=np.array([b'GSM1', b'GSM2']))
            h5.create_dataset('meta/samples/series_id', data=np.array([b'GSE1', b'GSE1']))
            h5.create_dataset('meta/samples/organism_ch1', data=np.array([b'Homo sapiens', b'Homo sapiens']))
        plan = self.strict_plan()
        with self.assertRaisesRegex(ValueError, 'rejects complex'):
            intake.validate_hardened_intake_plan_v21(plan)

    def test_v21_checks_observed_identifier_length(self):
        with self.metadata.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['GSM', 'gse_id', 'organism_ch1', 'molecule_ch1'])
            writer.writeheader()
            for index in range(2):
                writer.writerow({'GSM': 'G' * 128 + str(index), 'gse_id': 'GSE1',
                                 'organism_ch1': 'Homo sapiens', 'molecule_ch1': 'total RNA'})
        plan = self.strict_plan()
        with self.assertRaisesRegex(ValueError, 'metadata identifier exceeds'):
            intake.validate_hardened_intake_plan_v21(plan)

    def test_v21_requires_every_executed_asset(self):
        plan = self.strict_plan()
        plan['frozen_files'] = [item for item in plan['frozen_files']
                                if Path(item['path']).name != 'Gene_table_human.csv']
        with self.assertRaisesRegex(ValueError, 'execution assets are incomplete'):
            intake.validate_hardened_intake_plan_v21(plan)

    def test_v21_final_hash_recheck_catches_mid_preflight_tamper(self):
        plan = self.strict_plan()
        original = intake._intake_v21_release_preflight
        def preflight_and_tamper(*args, **kwargs):
            value = original(*args, **kwargs)
            path = self.root / 'FUN.Map_to_orthologs.R'
            payload = path.read_bytes()
            path.write_bytes(bytes([payload[0] ^ 1]) + payload[1:])
            return value
        with patch.object(intake, '_intake_v21_release_preflight', side_effect=preflight_and_tamper):
            with self.assertRaisesRegex(ValueError, 'File SHA256 changed before hardened intake'):
                intake.validate_hardened_intake_plan_v21(plan)

    def test_v21_run_lock_is_exclusive(self):
        destination = self.root / 'v21-run'
        with intake._intake_v21_run_lock(destination):
            with self.assertRaisesRegex(RuntimeError, 'already running'):
                with intake._intake_v21_run_lock(destination):
                    pass


class DiscoveryAvailabilityTest(unittest.TestCase):
    """Small software fixtures; these matrices are not empirical evidence."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        verifier_path = Path(__file__).resolve().parents[2] / 'src/verify_ingestion_2026_09_20.py'  # released copy
        require_development_artifact(verifier_path)
        spec = importlib.util.spec_from_file_location('discovered_input_verifier', verifier_path)
        self.verifier = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.verifier)
        self.release = self.root / 'software-fixture.h5'
        with h5py.File(self.release, 'w') as h5:
            h5.create_dataset('data/expression', data=np.zeros((2, 4), dtype=np.int32))
            h5.create_dataset('meta/samples/geo_accession', data=np.array([b'GSM1', b'GSM2', b'GSM3', b'GSM4']))
            h5.create_dataset('meta/samples/series_id', data=np.array([b'GSE99, GSE1', b'GSE1', b'GSE2', b'GSE3']))
            h5.create_dataset('meta/samples/organism_ch1', data=np.array([b'Homo sapiens']*4))
        self.discovery = self.root / 'discovery'
        self.discovery.mkdir()
        self.input_hashes = {}
        self.rows = []
        for uid, study, samples, status in [
            ('1', 'GSE1', ['GSM1', 'GSM2'], 'UNRESOLVED_NO_ASSAY_INPUT'),
            ('2', 'GSE2', ['GSM3', 'GSM4'], 'UNRESOLVED_NO_ASSAY_INPUT'),
            ('3', 'GSE3', ['GSM4'], 'SHARED_NEW_CANDIDATE_SAMPLES'),
            ('4', 'GSE4', ['GSM6'], 'UNRESOLVED_NO_ASSAY_INPUT'),
        ]:
            source = self.discovery / (study + '.soft')
            source.write_text('^SERIES = ' + study + '\n!Series_title = software fixture\n'
                '!Series_status = Public on Feb 20 2025\n' + ''.join('!Series_sample_id = ' + s + '\n' for s in samples))
            digest = self.verifier.sha(source)
            self.rows.append({'geo_uid': uid, 'accession': study, 'status': status,
                'first_public_date': '2025-02-20', 'sample_ids': samples,
                'series_source_sha256': digest,
                'source_receipt': {'status': 'PASS', 'response_path': str(source), 'sha256': digest}})
        self.seal(self.discovery / 'config.json', {'plan': {'modality': 'human_rna',
            'date_from': '2025-01-01', 'date_to': '2025-12-31'}})
        self.refresh_summary()

    def seal(self, path, value):
        path.write_text(json.dumps(value))
        path.with_suffix('.json.sha256').write_text(self.verifier.sha(path))

    def refresh_summary(self):
        (self.discovery / 'candidate-ledger.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in self.rows))
        outputs = {str(p): self.verifier.sha(p) for p in self.discovery.iterdir()
                   if p.name not in {'summary.json', 'summary.json.sha256'}}
        self.seal(self.discovery / 'summary.json', {
            'schema': 'clockbase-update-discovery-result-v1', 'status': 'COMPLETED',
            'enumerated_uids': len(self.rows), 'output_sha256': outputs, 'input_sha256': self.input_hashes})

    def run_readback(self):
        return self.verifier.discovery_availability(self.discovery/'summary.json', self.release,
                                                   self.verifier.sha(self.release))

    def test_secondary_membership_exact_only_and_overlap_exclusion(self):
        result = self.run_readback()
        self.assertEqual(set(result['molecular_sources']), {'GSE1'})
        self.assertEqual(result['counts']['by_availability'], {
            'EXACT_GEO_MEMBERSHIP': 2, 'SAMPLE_MEMBERSHIP_MISMATCH': 1, 'STUDY_ABSENT_FROM_RELEASE': 1})
        self.assertEqual(result['rows'][1]['missing_sample_ids'], ['GSM4'])
        self.assertFalse(result['rows'][2]['declared_as_available_input'])
        self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(result['original_paper_2025_eligibility'], 'QUARANTINED')

    def test_extra_archs4_sample_does_not_admit_partial_source(self):
        row = self.rows[0]
        row['sample_ids'] = ['GSM1']
        source = Path(row['source_receipt']['response_path'])
        source.write_text(source.read_text().replace('!Series_sample_id = GSM2\n', ''))
        row['series_source_sha256'] = row['source_receipt']['sha256'] = self.verifier.sha(source)
        self.refresh_summary()
        result = self.run_readback()
        self.assertEqual(result['rows'][0]['extra_sample_ids'], ['GSM2'])
        self.assertEqual(result['molecular_sources'], {})

    def test_source_date_disagreement_rejected(self):
        self.rows[0]['first_public_date'] = '2025-02-21'
        self.refresh_summary()
        with self.assertRaisesRegex(ValueError, 'first-public date differs'):
            self.run_readback()

    def test_source_output_tamper_rejected(self):
        (self.discovery/'GSE1.soft').write_text('changed software fixture')
        with self.assertRaisesRegex(ValueError, 'Discovery output changed'):
            self.run_readback()

    def test_release_hash_and_duplicate_sample_rejected(self):
        with self.assertRaisesRegex(ValueError, 'SHA256 differs'):
            self.verifier.discovery_availability(self.discovery/'summary.json', self.release, '0'*64)
        with h5py.File(self.release, 'r+') as h5:
            h5['meta/samples/geo_accession'][1] = b'GSM1'
        with self.assertRaisesRegex(ValueError, 'sample identities'):
            self.run_readback()

    def test_source_organism_mismatch_does_not_admit(self):
        with h5py.File(self.release, 'r+') as h5:
            h5['meta/samples/organism_ch1'][0] = b'Mus musculus'
        result = self.run_readback()
        self.assertEqual(result['rows'][0]['availability'], 'ORGANISM_MISMATCH')
        self.assertEqual(result['molecular_sources'], {})

    def test_duplicate_accession_rejected_instead_of_overwritten(self):
        self.rows[1]['accession'] = self.rows[0]['accession']
        self.refresh_summary()
        with self.assertRaisesRegex(ValueError, 'accessions.*duplicated'):
            self.run_readback()

    def test_stale_candidate_status_cannot_hide_shared_new_samples(self):
        self.rows[2]['status'] = 'UNRESOLVED_NO_ASSAY_INPUT'
        self.refresh_summary()
        result = self.run_readback()
        self.assertEqual(result['rows'][2]['operational_exclusions'], ['SHARED_CANDIDATE_SAMPLES'])
        self.assertNotIn('GSE3', result['molecular_sources'])

    def test_stale_status_cannot_hide_seen_accession_or_sample_overlap(self):
        for seen_id, expected_exclusion in [('GSE1', 'SEEN_ACCESSION'), ('GSE99', 'SHARED_SEEN_SAMPLES')]:
            with self.subTest(seen_id=seen_id):
                worklist, samples = self.root/'seen.jsonl', self.root/'seen-samples.json'
                worklist.write_text(json.dumps({'study_id': seen_id})+'\n')
                samples.write_text(json.dumps([{'accession': seen_id, 'sample_ids': ['GSM1']}]))
                config_path = self.discovery/'config.json'
                config = json.loads(config_path.read_text())
                for role, path in [('seen_worklist', worklist), ('seen_sample_sources', samples)]:
                    digest = self.verifier.sha(path)
                    config['plan'][role] = {'path': str(path), 'sha256': digest}
                    self.input_hashes[str(path)] = digest
                self.seal(config_path, config)
                self.refresh_summary()
                result = self.run_readback()
                self.assertIn(expected_exclusion, result['rows'][0]['operational_exclusions'])
                self.assertEqual(result['molecular_sources'], {})

    def test_unbound_configuration_and_incomplete_summary_rejected(self):
        path = self.discovery/'summary.json'
        summary = json.loads(path.read_text())
        summary['output_sha256'].pop(str(self.discovery/'config.json'))
        self.seal(path, summary)
        with self.assertRaisesRegex(ValueError, 'configuration is not bound'):
            self.run_readback()
        self.refresh_summary()
        summary = json.loads(path.read_text())
        summary['status'] = 'RUNNING'
        self.seal(path, summary)
        with self.assertRaisesRegex(ValueError, 'completed discovery summary'):
            self.run_readback()

    def test_duplicate_ledger_sample_is_not_silently_deduplicated(self):
        self.rows[0]['sample_ids'].append('GSM1')
        self.refresh_summary()
        with self.assertRaisesRegex(ValueError, 'GEO source membership differs'):
            self.run_readback()

    def test_summary_changed_during_release_scan_cannot_be_bound_as_original(self):
        original = self.verifier.sha
        release_hash = original(self.release)
        path = self.discovery/'summary.json'
        mutated = False
        def mutate_after_release_hash(file):
            nonlocal mutated
            result = original(file)
            if Path(file) == self.release and not mutated:
                mutated = True
                changed = json.loads(path.read_text())
                changed['unexpected_update'] = 'disclosed software fixture mutation'
                path.write_text(json.dumps(changed))
                path.with_suffix('.json.sha256').write_text(original(path))
            return result
        with patch.object(self.verifier, 'sha', side_effect=mutate_after_release_hash):
            with self.assertRaisesRegex(ValueError, 'summary, configuration or verifier changed'):
                self.verifier.discovery_availability(path, self.release, release_hash)


class DiscoverySourceProfilesTest(unittest.TestCase):
    """SOFT software fixtures, not observed cohort evidence."""

    def setUp(self):
        fixture = DiscoveryAvailabilityTest()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.fixture = fixture
        self.verifier = fixture.verifier
        self.discovery = fixture.discovery
        fixture.rows = [fixture.rows[0], fixture.rows[3]]
        self.candidates = []
        for row in fixture.rows:
            row['status'] = 'READY_FOR_INGESTION_ATTEMPT'
            sample_path = self.discovery / (row['accession'] + '-samples.soft')
            sample_path.write_text(''.join(
                '^SAMPLE = ' + sample + '\n!Sample_title = software fixture\n'
                '!Sample_organism_ch1 = Homo sapiens\n!Sample_molecule_ch1 = total RNA\n'
                for sample in row['sample_ids']))
            self.candidates.append({'study_id': row['accession'], 'modality': 'human_rna',
                'sources': {'series': {'path': row['source_receipt']['response_path'],
                                      'sha256': row['series_source_sha256']},
                            'samples': {'path': str(sample_path), 'sha256': self.verifier.sha(sample_path)},
                            'archs4': {'path': str(fixture.release), 'sha256': self.verifier.sha(fixture.release)}}})
        self.refresh()

    def refresh(self):
        for row in self.candidates:
            for record in row['sources'].values():
                record['sha256'] = self.verifier.sha(record['path'])
        (self.discovery / 'candidates.json').write_text(json.dumps(
            {'schema': 'clockbase-update-candidates-v1', 'studies': self.candidates}))
        self.fixture.refresh_summary()
        path = self.discovery / 'summary.json'
        summary = json.loads(path.read_text())
        summary['ready_for_ingestion_attempt'] = len(self.candidates)
        self.fixture.seal(path, summary)

    def profile(self):
        return self.verifier.discovery_source_profiles(self.discovery / 'summary.json')

    def mutate_samples(self, transform):
        path = Path(self.candidates[0]['sources']['samples']['path'])
        path.write_text(transform(path.read_text()))
        self.refresh()

    def test_literal_supported_profile_is_not_scientific_acceptance(self):
        result = self.profile()
        self.assertEqual(result['supported_study_ids'], ['GSE1', 'GSE4'])
        self.assertEqual(result['counts']['enrolled_samples'], 3)
        self.assertEqual(result['expression_cells_read'], 0)
        self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
        self.assertTrue(all(r['scientific_status'] == 'CANNOT VERIFY' for r in result['studies']))

    def test_mixed_non_rna_is_retained_without_dropping_other_study(self):
        self.mutate_samples(lambda s: s.replace('!Sample_molecule_ch1 = total RNA',
                                               '!Sample_molecule_ch1 = protein', 1))
        result = self.profile()
        self.assertEqual(result['counts']['enrolled_studies'], 2)
        self.assertEqual(result['supported_study_ids'], ['GSE4'])
        self.assertEqual(result['stopped_study_ids'], ['GSE1'])
        self.assertEqual(result['studies'][0]['molecule_label_counts'], {'protein': 1, 'total RNA': 1})
        self.assertEqual(result['studies'][0]['scientific_status'], 'CANNOT VERIFY')

    def test_polya_is_outside_profile_not_called_erroneous(self):
        self.mutate_samples(lambda s: s.replace('total RNA', 'polyA RNA'))
        result = self.profile()
        self.assertEqual(result['studies'][0]['profile_status'], 'UNSUPPORTED_SOURCE_PROFILE')
        self.assertIn('unsupported does not mean erroneous', result['scope'])
        self.assertEqual(result['counts']['molecule_label_occurrences'], {'polyA RNA': 2, 'total RNA': 1})

    def test_duplicate_and_missing_annotations_remain_ambiguous(self):
        self.mutate_samples(lambda s: s.replace('!Sample_molecule_ch1 = total RNA',
            '!Sample_molecule_ch1 = total RNA\n!Sample_molecule_ch1 = total RNA', 1))
        result = self.profile()
        self.assertEqual(result['studies'][0]['profile_status'], 'UNSUPPORTED_SOURCE_PROFILE')
        self.assertTrue(any('ambiguous molecule_ch1' in value for value in result['studies'][0]['limitations']))

    def test_platform_section_cannot_supply_missing_sample_annotation(self):
        self.mutate_samples(lambda s: s.replace('!Sample_molecule_ch1 = total RNA',
            '^PLATFORM = GPL1\n!Sample_molecule_ch1 = total RNA', 1))
        result = self.profile()
        self.assertEqual(result['studies'][0]['profile_status'], 'UNSUPPORTED_SOURCE_PROFILE')
        self.assertEqual(result['studies'][0]['samples'][0]['molecule_ch1'], [])

    def test_sample_membership_mismatch_is_not_reconciled_by_position(self):
        self.mutate_samples(lambda s: s.replace('^SAMPLE = GSM2', '^SAMPLE = GSM8'))
        with self.assertRaisesRegex(ValueError, 'sample SOFT membership differs'):
            self.profile()

    def test_foreign_organism_retained_as_unresolved(self):
        self.mutate_samples(lambda s: s.replace('Homo sapiens', 'Mus musculus', 1))
        result = self.profile()
        self.assertEqual(result['studies'][0]['profile_status'], 'UNSUPPORTED_SOURCE_PROFILE')
        self.assertTrue(any('organism' in issue for issue in result['studies'][0]['limitations']))

    def test_candidate_omission_from_ready_ledger_is_rejected(self):
        self.candidates.pop()
        self.refresh()
        with self.assertRaisesRegex(ValueError, 'differ from the ready ledger'):
            self.profile()

    def test_changed_sample_payload_is_rejected(self):
        path = Path(self.candidates[0]['sources']['samples']['path'])
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'Discovery output changed'):
            self.profile()


class IndependentBatchSourceReadbackTest(unittest.TestCase):
    """Tiny deterministic software fixtures; none are empirical measurements."""

    def setUp(self):
        from audit_system.tests.test_archs4_batch_extraction_2026_09_21 import Archs4BatchExtractionTest
        self.fixture = Archs4BatchExtractionTest()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.summary = self.fixture.run_basic()
        self.path = self.fixture.output / 'batch-summary.json'
        path = Path(__file__).resolve().parents[2] / 'src/verify_ingestion_2026_09_20.py'  # released copy
        require_development_artifact(path)
        spec = importlib.util.spec_from_file_location('independent_batch_readback', path)
        self.verifier = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.verifier)

    def seal_summary(self):
        self.summary['artifacts'] = {
            str(path.relative_to(self.fixture.output)): self.fixture.sha(path)
            for path in self.fixture.output.rglob('*') if path.is_file()
            and path.name not in {'batch-summary.json', 'batch-summary.json.sha256', '.intake.lock'}}
        self.path.write_text(json.dumps(self.summary))
        self.path.with_suffix('.json.sha256').write_text(self.fixture.sha(self.path))

    def rewrite_matrix(self, transform):
        """Reseal producer records so only an independent source read detects error."""
        row = self.summary['results'][0]
        matrix = self.fixture.output / row['matrix']['path']
        matrix.write_text(transform(matrix.read_text()))
        row['matrix'].update(sha256=self.fixture.sha(matrix), size_bytes=matrix.stat().st_size)
        extraction = self.fixture.output / row['extraction_receipt']['path']
        payload = json.loads(extraction.read_text())
        payload['matrix'] = row['matrix']
        extraction.write_text(json.dumps(payload))
        row['extraction_receipt']['sha256'] = self.fixture.sha(extraction)
        row['extraction_receipt']['size_bytes'] = extraction.stat().st_size
        stage_path = self.fixture.output / row['stage_receipt']['path']
        stage = json.loads(stage_path.read_text())
        stage['result']['value'].update(matrix=row['matrix'],
                                        receipt_sha256=row['extraction_receipt']['sha256'])
        terminal_path = stage_path.parent / next(name for name in stage['outputs']
                                                 if Path(name).name == 'terminal-attempt.json')
        terminal_path.write_text(json.dumps(stage['result']))
        for path in (matrix, extraction, terminal_path):
            stage['outputs'][str(path.relative_to(stage_path.parent))] = self.fixture.sha(path)
        stage_path.write_text(json.dumps(stage))
        stage_path.with_suffix('.sha256').write_text(self.fixture.sha(stage_path))
        row['stage_receipt']['sha256'] = self.fixture.sha(stage_path)
        self.summary['counts']['extracted_bytes'] = sum(
            row['matrix']['size_bytes'] for row in self.summary['results'] if row['matrix'])
        self.seal_summary()

    def readback(self):
        return self.verifier.verify_batch_extraction(self.path)

    def test_every_cell_matches_and_source_stops_remain_unresolved_read_only(self):
        before = {str(p): (self.fixture.sha(p), p.stat().st_mtime_ns)
                  for p in self.fixture.root.rglob('*') if p.is_file()}
        result = self.readback()
        after = {str(p): (self.fixture.sha(p), p.stat().st_mtime_ns)
                 for p in self.fixture.root.rglob('*') if p.is_file()}
        self.assertEqual(before, after)
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['study_status_counts'], {'CANNOT VERIFY': 1, 'PASS': 2})
        self.assertEqual(result['enrolled_studies'], 3)
        self.assertEqual(result['cells_compared'], 6)
        self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(result['original_paper_2025_eligibility'], 'QUARANTINED')
        self.assertEqual(result['results'][2]['sample_ids'], ['GSM4'])
        self.assertEqual(result['model_calls'], 0)

    def test_source_cell_change_fails_even_with_resealed_producer_hashes(self):
        self.rewrite_matrix(lambda s: s.replace('ENSG1,1,2', 'ENSG1,9,2'))
        result = self.readback()
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(result['results'][0]['status'], 'FAIL')
        self.assertIn('count cell differs', result['results'][0]['reason'])
        self.assertEqual(result['results'][1]['status'], 'PASS')
        self.assertEqual(result['results'][2]['status'], 'CANNOT VERIFY')

    def test_sample_order_change_cannot_pass_on_matching_dimensions(self):
        self.rewrite_matrix(lambda s: s.replace('GSM1,GSM2', 'GSM2,GSM1'))
        result = self.readback()
        self.assertEqual(result['results'][0]['status'], 'FAIL')
        self.assertIn('sample order', result['results'][0]['reason'])

    def test_feature_order_change_cannot_pass_on_matching_values(self):
        self.rewrite_matrix(lambda s: s.replace('ENSG1,1,2\nENSG2,5,6', 'ENSG2,1,2\nENSG1,5,6'))
        result = self.readback()
        self.assertEqual(result['results'][0]['status'], 'FAIL')
        self.assertIn('feature identity', result['results'][0]['reason'])

    def test_negative_count_cannot_wrap_to_unsigned(self):
        self.rewrite_matrix(lambda s: s.replace('ENSG1,1,2', 'ENSG1,-1,2'))
        result = self.readback()
        self.assertEqual(result['results'][0]['status'], 'FAIL')
        self.assertIn('unsigned integer literal', result['results'][0]['reason'])

    def test_decimal_count_cannot_be_rounded_to_an_integer(self):
        self.rewrite_matrix(lambda s: s.replace('ENSG1,1,2', 'ENSG1,1.0,2'))
        self.assertEqual(self.readback()['results'][0]['status'], 'FAIL')

    def test_omitted_unresolved_case_rejected_before_cell_read(self):
        self.summary['results'].pop()
        self.seal_summary()
        with self.assertRaisesRegex(ValueError, 'complete source cohort'):
            self.readback()

    def test_denominator_drift_rejected(self):
        self.summary['counts']['extracted_samples'] += 1
        self.seal_summary()
        with self.assertRaisesRegex(ValueError, 'counts or terminal summary'):
            self.readback()

    def test_ambiguous_path_scope_rejected(self):
        self.summary.pop('path_base')
        self.seal_summary()
        with self.assertRaisesRegex(ValueError, 'completed scoped extraction batch'):
            self.readback()

    def test_unsealed_payload_change_rejected(self):
        path = self.fixture.output / self.summary['results'][0]['matrix']['path']
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'artifact hash differs'):
            self.readback()

    def test_origin_release_change_rejected(self):
        with h5py.File(self.fixture.release, 'r+') as handle:
            handle['data/expression'][0, 0] = 9
        with self.assertRaisesRegex(ValueError, 'Independent release hash differs'):
            self.readback()

    def rewrite_extraction(self, transform):
        row = self.summary['results'][0]
        path = self.fixture.output / row['extraction_receipt']['path']
        value = json.loads(path.read_text())
        transform(value)
        path.write_text(json.dumps(value))
        # This also refreshes the extraction receipt and terminal result bindings.
        self.rewrite_matrix(lambda text: text)

    def test_secondary_study_membership_must_match_hdf5(self):
        self.rewrite_extraction(lambda r: r['all_study_memberships']['GSM1'].append('GSE999'))
        result = self.readback()
        self.assertEqual(result['results'][0]['status'], 'FAIL')
        self.assertIn('membership', result['results'][0]['reason'])

    def test_metadata_hash_must_resolve_to_completed_metadata_stage(self):
        self.rewrite_extraction(lambda r: r['source'].update(metadata_sha256='0' * 64))
        self.assertEqual(self.readback()['results'][0]['status'], 'FAIL')

    def test_context_hash_must_resolve_to_completed_metadata_stage(self):
        self.rewrite_extraction(lambda r: r['source'].update(context_sha256='0' * 64))
        self.assertEqual(self.readback()['results'][0]['status'], 'FAIL')

    def test_stage_status_cannot_disagree_with_terminal_attempt(self):
        row = self.summary['results'][0]
        path = self.fixture.output / row['stage_receipt']['path']
        value = json.loads(path.read_text())
        value['result'] = {'status': 'CANNOT VERIFY', 'value': None}
        path.write_text(json.dumps(value))
        row['stage_receipt']['sha256'] = self.fixture.sha(path)
        path.with_suffix('.sha256').write_text(self.fixture.sha(path))
        self.seal_summary()
        result = self.readback()
        self.assertEqual(result['results'][0]['status'], 'FAIL')
        self.assertIn('result and retained attempt disagree', result['results'][0]['reason'])

    def test_unresolved_summary_cannot_hide_successful_extraction(self):
        row = self.summary['results'][0]
        row.update(status='EXTRACTION_UNRESOLVED', matrix=None, extraction_receipt=None,
                   error='software fixture changed terminal label')
        counts = self.summary['counts']
        counts.update(extracted_studies=1, extraction_unresolved_studies=1,
                      extracted_samples=self.summary['results'][1]['matrix']['sample_count'],
                      extracted_bytes=self.summary['results'][1]['matrix']['size_bytes'])
        self.seal_summary()
        result = self.readback()
        self.assertEqual(result['results'][0]['status'], 'FAIL')
        self.assertIn('disagrees with its terminal result', result['results'][0]['reason'])

    def test_supported_profile_with_failed_extraction_retains_unresolved_stage(self):
        output = self.fixture.root / 'bounded-failure'
        plan_path = self.fixture.root / 'bounded-failure-plan.json'
        self.fixture.write_plan(plan_path, output, lambda p: p['limits'].update(
            max_extracted_bytes_per_study=1))
        self.fixture.run_basic(output=output, plan=plan_path)
        result = self.verifier.verify_batch_extraction(output / 'batch-summary.json')
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['study_status_counts'], {'CANNOT VERIFY': 3})
        self.assertEqual(result['cells_compared'], 0)
        self.assertIn('per-study cap', result['results'][0]['reason'])


if __name__ == '__main__':
    unittest.main()
