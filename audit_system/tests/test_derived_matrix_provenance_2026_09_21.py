"""Software fixtures for provenance integrity; no empirical results or API calls."""
import copy
import importlib.util
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from src import update_pipeline_2026_09_20 as update
from audit_system.pipeline_measurement_evidence_2026_09_20 import build_measurement_evidence
from audit_system.tests.test_pipeline_measurement_evidence_2026_09_20 import _fixture
from audit_system.tests import require_development_artifact


class DerivedMatrixFixture:
    def __init__(self, root, matrix=None):
        self.root = root
        self.batch = root / 'batch'; self.batch.mkdir()
        self.stage = self.batch / 'studies/GSE1/extract'; self.stage.mkdir(parents=True)
        attempt = self.stage / 'attempt-000'; attempt.mkdir()
        self.matrix = attempt / 'counts.csv'
        if matrix is None:
            self.matrix.write_text(',GSM1,GSM2\nf1,1,2\nf2,3,4\n')
        else:
            shutil.copyfile(matrix, self.matrix)
        self.series = root / 'series.soft'
        self.series.write_text('^SERIES = GSE1\n!Series_status = Public on Jan 03 2025\n!Series_sample_id = GSM1\n!Series_sample_id = GSM2\n')
        self.samples = root / 'samples.soft'
        self.samples.write_text(''.join(f'^SAMPLE = {sample}\n!Sample_title = Software fixture\n!Sample_organism_ch1 = Homo sapiens\n!Sample_molecule_ch1 = total RNA\n' for sample in ('GSM1', 'GSM2')))
        self.profile_script = root / 'profile-script.py'; self.profile_script.write_text('# Software fixture\n')
        self.verifier = root / 'independent-verifier.py'; self.verifier.write_text('# Separate software verifier fixture\n')
        self.implementation = self.batch / 'implementation/batch.py'; self.implementation.parent.mkdir()
        self.implementation.write_text('# Frozen batch software fixture\n')
        self.profile_row = {'study_id': 'GSE1', 'modality': 'human_rna', 'sample_ids': ['GSM1', 'GSM2'],
                            'profile_status': 'SUPPORTED_SOURCE_PROFILE', 'sources': {'series': self.decl(self.series), 'samples': self.decl(self.samples)}}
        self.profile = root / 'source-profiles.json'
        self.write(self.profile, {'schema': 'clockbase-rna-source-profile-census-v1', 'studies': [self.profile_row],
            'supported_study_ids': ['GSE1'], 'stopped_study_ids': [], 'scientific_status': 'CANNOT VERIFY',
            'verification_script_sha256': update.sha(self.profile_script),
            'software_blinding_marker': 'DO_NOT_INCLUDE_DERIVATION_RECEIPTS_IN_SOURCE_DESIGN'})
        self.extraction = attempt / 'extraction.json'
        self.write(self.extraction, {'schema': 'clockbase-archs4-batch-study-extraction-v1', 'path_base': 'batch_root',
            'study_id': 'GSE1', 'modality': 'human_rna', 'species': 'human', 'status': 'PASS',
            'input_mode': 'FRESH_ARCHS4_BATCH_EXTRACTION', 'scientific_status': 'CANNOT VERIFY',
            'matrix': {**self.relative(self.matrix), 'size_bytes': self.matrix.stat().st_size, 'feature_count': 2, 'sample_count': 2},
            'sample_ids_in_matrix_order': ['GSM1', 'GSM2'], 'all_study_memberships': {'GSM1': ['GSE1'], 'GSM2': ['GSE1']},
            'source': {'release_sha256': 'a' * 64, 'release_size_bytes': 62000000000,
                       'series_sha256': update.sha(self.series), 'samples_sha256': update.sha(self.samples),
                       'source_profile_descriptor_sha256': update.sha(self.profile),
                       'source_profile_record_sha256': update.identity(self.profile_row)}})
        self.stage_receipt = self.stage / 'complete.json'
        self.write(self.stage_receipt, {'outputs': {str(p.relative_to(self.stage)): update.sha(p) for p in (self.matrix, self.extraction)},
                                       'result': {'status': 'PASS'}})
        self.stage_receipt.with_suffix('.sha256').write_text(update.sha(self.stage_receipt) + '\n')
        self.plan = root / 'batch-plan.json'
        execution_inputs = {'selected_study_ids': ['GSE1'], 'source_profile_descriptor': self.decl(self.profile),
            'release': {'path': '/archive/not-opened-by-child.h5', 'sha256': 'a' * 64, 'size_bytes': 62000000000},
            'runtime': {'python': 'software-fixture'}, 'limits': {'max_studies': 1},
            'implementation': {'batch': self.decl(self.implementation)}}
        self.write(self.plan, {'schema': 'clockbase-archs4-batch-extraction-plan-v1', 'output': str(self.batch), **execution_inputs})
        self.config = self.batch / 'config.json'
        self.seal(self.config, {'schema': 'clockbase-archs4-batch-extraction-config-v1',
            'plan_path': str(self.plan), 'plan_sha256': update.sha(self.plan), **execution_inputs})
        self.releases = {}
        for phase in ('before', 'after'):
            path = self.batch / ('release-' + phase + '.json')
            self.seal(path, {'schema': 'clockbase-archs4-batch-release-check-v1', 'phase': phase,
                'path': '/archive/not-opened-by-child.h5', 'sha256': 'a' * 64, 'size_bytes': 62000000000,
                'expected_sha256': 'a' * 64, 'status': 'PASS', 'hash_scope': 'FULL_HDF5',
                'stat': {'device': 1, 'inode': 2, 'size_bytes': 62000000000, 'mtime_ns': 3, 'ctime_ns': 4}})
            self.releases[phase] = self.relative(path)
        self.execution = self.batch / 'execution-attempt.json'
        self.seal(self.execution, {'schema': 'clockbase-archs4-batch-execution-attempt-v1', 'status': 'COMPLETED',
            'plan_sha256': update.sha(self.plan), 'release_before': self.releases['before'], 'release_after': self.releases['after'],
            'result_count': 1, 'model_api_calls': 0, 'preprocessing_runs': 0, 'clock_runs': 0})
        self.summary = self.batch / 'batch-summary.json'
        self.summary_data = {'schema': 'clockbase-archs4-batch-extraction-summary-v1', 'path_base': 'batch_root',
            'status': 'COMPLETED', 'execution_status': 'COMPLETED', 'plan_sha256': update.sha(self.plan),
            'config_sha256': update.sha(self.config), 'model_api_calls': 0, 'preprocessing_runs': 0, 'clock_runs': 0,
            'source_profile_descriptor': self.decl(self.profile),
            'source_hdf5': {'sha256': 'a' * 64, 'size_bytes': 62000000000},
            'release_before': self.releases['before'], 'release_after': self.releases['after'],
            'implementation': {'batch': self.relative(self.implementation)},
            'execution_attempt': self.relative(self.execution), 'scientific_status': 'CANNOT VERIFY',
            'results': [{'study_id': 'GSE1', 'modality': 'human_rna', 'status': 'EXTRACTED',
                'source_profile_status': 'SUPPORTED_SOURCE_PROFILE', 'source_profile_record_sha256': update.identity(self.profile_row),
                'series': self.decl(self.series), 'samples': self.decl(self.samples), 'matrix': self.relative(self.matrix),
                'extraction_receipt': self.relative(self.extraction), 'stage_receipt': self.relative(self.stage_receipt)}],
            'artifacts': {str(p.relative_to(self.batch)): update.sha(p) for p in self.batch.rglob('*') if p.is_file()}}
        self.seal(self.summary, self.summary_data)
        self.verification = root / 'independent-every-cell-readback.json'
        self.verification_data = {'schema': 'clockbase-archs4-independent-every-cell-v1', 'status': 'PASS',
            'batch_summary_sha256': update.sha(self.summary), 'release_sha256': 'a' * 64,
            'release_sha256_before': 'a' * 64, 'release_sha256_after': 'a' * 64,
            'source_profile_sha256': update.sha(self.profile), 'verification_script': self.decl(self.verifier),
            'verification_script_sha256': update.sha(self.verifier), 'scientific_status': 'CANNOT VERIFY',
            'model_calls': 0, 'preprocessing_runs': 0, 'clock_runs': 0,
            'results': [{'study_id': 'GSE1', 'modality': 'human_rna', 'status': 'PASS',
                'matrix_sha256': update.sha(self.matrix), 'extraction_receipt_sha256': update.sha(self.extraction),
                'sample_ids': ['GSM1', 'GSM2'], 'feature_count': 2, 'sample_count': 2, 'cells_compared': 4,
                'feature_ids_equal': True, 'sample_ids_equal': True, 'all_cells_equal': True}]}
        self.write(self.verification, self.verification_data)
        self.wrapper = root / 'matrix-provenance.json'
        self.wrapper_data = {'schema': update.MATRIX_PROVENANCE_SCHEMA, 'study_id': 'GSE1', 'modality': 'human_rna',
            'matrix': self.decl(self.matrix), 'series': self.decl(self.series), 'samples': self.decl(self.samples),
            'batch_summary': self.decl(self.summary), 'extraction_receipt': self.decl(self.extraction),
            'verification_receipt': self.decl(self.verification), 'source_profile_descriptor': self.decl(self.profile),
            'source_profile_script': self.decl(self.profile_script), 'source_profile_record_sha256': update.identity(self.profile_row),
            'source_hdf5': {'sha256': 'a' * 64, 'size_bytes': 62000000000, 'hash_scope': 'historical_batch_pre_and_post'},
            'implementation': {'batch': self.decl(self.implementation)},
            'scientific_status': 'CANNOT VERIFY', 'original_paper_2025_eligibility': 'QUARANTINED'}
        self.write(self.wrapper, self.wrapper_data)

    @staticmethod
    def write(path, value):
        path.write_text(json.dumps(value, indent=2) + '\n')

    @staticmethod
    def decl(path):
        return {'path': str(path.resolve()), 'sha256': update.sha(path)}

    def relative(self, path):
        return {'path': str(path.relative_to(self.batch)), 'sha256': update.sha(path)}

    def seal(self, path, value):
        self.write(path, value)
        path.with_suffix('.json.sha256').write_text(update.sha(path) + '\n')

    def entry(self):
        return {'study_id': 'GSE1', 'modality': 'human_rna', 'sources': {
            'matrix': self.wrapper_data['matrix'], 'series': self.wrapper_data['series'], 'samples': self.wrapper_data['samples'],
            'matrix_provenance': self.decl(self.wrapper)}}

    def update_wrapper(self):
        self.write(self.wrapper, self.wrapper_data)

    def update_verification(self):
        self.write(self.verification, self.verification_data)
        self.wrapper_data['verification_receipt'] = self.decl(self.verification)
        self.update_wrapper()

    def update_summary(self):
        self.seal(self.summary, self.summary_data)
        self.verification_data['batch_summary_sha256'] = update.sha(self.summary)
        self.wrapper_data['batch_summary'] = self.decl(self.summary)
        self.update_verification()


class DerivedMatrixProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fixture = DerivedMatrixFixture(self.root)
        self.limits = {'max_input_bytes': 2000000, 'allow_network': False}

    def collect(self, entry=None):
        return update.collect_sources(entry or self.fixture.entry(), self.root / 'acquisition', self.limits, {})

    def test_nested_sources_and_historical_origin_without_hdf5_access(self):
        original_sha = update.sha
        def protected_sha(path):
            self.assertNotEqual(Path(path).suffix, '.h5', 'Child must not rehash historical HDF5')
            return original_sha(path)
        with patch.object(update, 'sha', side_effect=protected_sha):
            sources = self.collect()
            binding = update.read_matrix_provenance(self.fixture.wrapper, original_sha(self.fixture.wrapper), sources, 'GSE1', 'human_rna')
        self.assertIn('matrix_provenance:verification_script', sources)
        self.assertIn('matrix_provenance:release_after:seal', sources)
        self.assertEqual(binding['measurement']['cells_compared'], 4)
        self.assertFalse(binding['measurement']['hdf5_rehashed_for_this_child'])
        self.assertEqual(binding['measurement']['scientific_status'], 'CANNOT VERIFY')
        self.assertNotIn(str(self.root), json.dumps(binding['measurement']))

    def test_arbitrary_matrix_without_derivation_keeps_existing_behavior(self):
        entry = self.fixture.entry(); del entry['sources']['matrix_provenance']
        self.assertEqual(set(self.collect(entry)), {'series', 'samples', 'matrix'})

    def test_provenance_requires_local_rna_matrix(self):
        for alteration in ('url', 'dnam', 'archs4'):
            with self.subTest(alteration=alteration):
                entry = self.fixture.entry()
                if alteration == 'url': entry['sources']['matrix_provenance'] = {'url': 'https://example.com/provenance.json'}
                elif alteration == 'dnam': entry['modality'] = 'human_dnam'
                else: entry['sources']['archs4'] = entry['sources'].pop('matrix')
                with self.assertRaisesRegex(ValueError, 'Derived matrix provenance'):
                    self.collect(entry)

    def test_missing_and_changed_nested_receipts_stop(self):
        self.fixture.verification.write_text('{}')
        with self.assertRaisesRegex(ValueError, 'nested SHA256 differs: verification_receipt'):
            self.collect()
        self.fixture.verification.unlink()
        with self.assertRaisesRegex(ValueError, 'missing nested file: verification_receipt'):
            self.collect()

    def test_foreign_same_content_matrix_path_is_rejected(self):
        copied = self.root / 'foreign-counts.csv'; shutil.copyfile(self.fixture.matrix, copied)
        self.fixture.wrapper_data['matrix'] = self.fixture.decl(copied); self.fixture.update_wrapper()
        with self.assertRaisesRegex(ValueError, 'foreign or changed matrix lineage'):
            self.collect()

    def test_foreign_source_metadata_is_rejected_after_updated_outer_hash(self):
        foreign = self.root / 'foreign-series.soft'; shutil.copyfile(self.fixture.series, foreign)
        self.fixture.wrapper_data['series'] = self.fixture.decl(foreign); self.fixture.update_wrapper()
        with self.assertRaisesRegex(ValueError, 'foreign series source lineage'):
            self.collect()

    def test_unsupported_schema_foreign_study_and_scientific_overclaim_rejected(self):
        original = copy.deepcopy(self.fixture.wrapper_data)
        for field, value, error in [('schema', 'unsupported-v9', 'unsupported wrapper schema'),
                                    ('study_id', 'GSE2', 'study or RNA modality differs'),
                                    ('scientific_status', 'PASS', 'scientific or temporal acceptance')]:
            with self.subTest(field=field):
                self.fixture.wrapper_data = dict(original, **{field: value}); self.fixture.update_wrapper()
                with self.assertRaisesRegex(ValueError, error): self.collect()

    def test_independent_readback_cannot_omit_cells_or_duplicate_study(self):
        original = copy.deepcopy(self.fixture.verification_data)
        for change, error in [
            (lambda x: x['results'][0].update(cells_compared=3), 'every-cell verification dimensions'),
            (lambda x: x['results'][0].update(all_cells_equal=False), 'every-cell verification dimensions'),
            (lambda x: x['results'][0].update(sample_ids=['GSM2', 'GSM1']), 'matrix sample identity/order differs'),
            (lambda x: x['results'].append(copy.deepcopy(x['results'][0])), 'ambiguous verification study inventory')]:
            with self.subTest(error=error):
                self.fixture.verification_data = copy.deepcopy(original); change(self.fixture.verification_data)
                self.fixture.update_verification()
                with self.assertRaisesRegex(ValueError, error): self.collect()

    def test_changed_implementation_is_bound_on_resume(self):
        sources = self.collect()
        self.fixture.implementation.write_text('# changed software fixture\n')
        with self.assertRaisesRegex(ValueError, 'Frozen source changed: matrix_provenance:implementation:batch'):
            update.verify_sources(sources)
        with self.assertRaisesRegex(ValueError, 'nested SHA256 differs: implementation:batch'):
            self.collect()

    def test_batch_relative_path_base_and_completion_seals_are_required(self):
        original = copy.deepcopy(self.fixture.summary_data)
        self.fixture.summary_data['path_base'] = 'receipt_parent'
        self.fixture.seal(self.fixture.summary, self.fixture.summary_data)
        self.fixture.wrapper_data['batch_summary'] = self.fixture.decl(self.fixture.summary); self.fixture.update_wrapper()
        with self.assertRaisesRegex(ValueError, 'unsupported batch summary schema or path base'):
            self.collect()
        self.fixture.seal(self.fixture.summary, original)
        self.fixture.wrapper_data['batch_summary'] = self.fixture.decl(self.fixture.summary); self.fixture.update_wrapper()
        self.fixture.summary.with_suffix('.json.sha256').unlink()
        with self.assertRaisesRegex(ValueError, 'missing or changed completion seal: batch_summary'):
            self.collect()

    def test_same_bytes_unbound_nested_source_does_not_enter_measurement(self):
        sources = self.collect()
        duplicate = self.root / 'duplicate-verifier.py'; shutil.copyfile(self.fixture.verifier, duplicate)
        sources['matrix_provenance:verification_script']['path'] = str(duplicate)
        # Plain hash verification is insufficient to substitute a source role.
        update.verify_sources(sources)
        binding = update.read_matrix_provenance(self.fixture.wrapper, update.sha(self.fixture.wrapper), sources, 'GSE1', 'human_rna')
        self.assertNotEqual(binding['files']['verification_script']['path'], str(duplicate))

    def test_batch_matrix_membership_and_zero_execution_scope_cannot_be_relaxed(self):
        original = copy.deepcopy(self.fixture.summary_data)
        for mutation, error in [
            (lambda x: x['artifacts'].pop(str(self.fixture.matrix.relative_to(self.fixture.batch))), 'matrix is absent from batch artifact inventory'),
            (lambda x: x.update(clock_runs=1), 'batch scope is not extraction-only')]:
            with self.subTest(error=error):
                self.fixture.summary_data = copy.deepcopy(original); mutation(self.fixture.summary_data)
                self.fixture.update_summary()
                with self.assertRaisesRegex(ValueError, error): self.collect()

    def test_resealed_release_stat_change_or_omission_is_rejected(self):
        path = self.fixture.batch / 'release-after.json'
        original = json.loads(path.read_text())
        for mutation, error in [
            (lambda x: x.pop('stat'), 'historical release stat receipt is incomplete'),
            (lambda x: x['stat'].update(inode=99), 'historical HDF5 identity/stat changed')]:
            with self.subTest(error=error):
                value = copy.deepcopy(original); mutation(value); self.fixture.seal(path, value)
                declared = self.fixture.relative(path); self.fixture.summary_data['release_after'] = declared
                execution = json.loads(self.fixture.execution.read_text()); execution['release_after'] = declared
                self.fixture.seal(self.fixture.execution, execution)
                self.fixture.summary_data['execution_attempt'] = self.fixture.relative(self.fixture.execution)
                for changed in (path, path.with_suffix('.json.sha256'), self.fixture.execution, self.fixture.execution.with_suffix('.json.sha256')):
                    self.fixture.summary_data['artifacts'][str(changed.relative_to(self.fixture.batch))] = update.sha(changed)
                self.fixture.update_summary()
                with self.assertRaisesRegex(ValueError, error): self.collect()

    def test_resealed_plan_config_disagreement_is_rejected(self):
        config = json.loads(self.fixture.config.read_text()); config['limits']['max_studies'] = 99
        self.fixture.seal(self.fixture.config, config)
        self.fixture.summary_data['config_sha256'] = update.sha(self.fixture.config)
        for changed in (self.fixture.config, self.fixture.config.with_suffix('.json.sha256')):
            self.fixture.summary_data['artifacts'][str(changed.relative_to(self.fixture.batch))] = update.sha(changed)
        self.fixture.update_summary()
        with self.assertRaisesRegex(ValueError, 'plan/config release, runtime or limits differ'):
            self.collect()

    def test_resealed_plan_cannot_claim_another_batch_output(self):
        plan = json.loads(self.fixture.plan.read_text()); plan['output'] = str(self.root / 'another-batch')
        self.fixture.write(self.fixture.plan, plan)
        config = json.loads(self.fixture.config.read_text()); config['plan_sha256'] = update.sha(self.fixture.plan)
        self.fixture.seal(self.fixture.config, config)
        execution = json.loads(self.fixture.execution.read_text()); execution['plan_sha256'] = update.sha(self.fixture.plan)
        self.fixture.seal(self.fixture.execution, execution)
        self.fixture.summary_data.update(plan_sha256=update.sha(self.fixture.plan), config_sha256=update.sha(self.fixture.config),
                                         execution_attempt=self.fixture.relative(self.fixture.execution))
        for path in (self.fixture.config, self.fixture.config.with_suffix('.json.sha256'),
                     self.fixture.execution, self.fixture.execution.with_suffix('.json.sha256')):
            self.fixture.summary_data['artifacts'][str(path.relative_to(self.fixture.batch))] = update.sha(path)
        self.fixture.update_summary()
        with self.assertRaisesRegex(ValueError, 'batch plan output does not identify this archived batch'):
            self.collect()

    def test_source_only_packet_does_not_receive_derivation(self):
        from audit_system.agent_workflow_2026_09_10 import source_packet
        sources = self.collect(); directory = self.root / 'metadata'; directory.mkdir()
        acquired = update.frozen_metadata('GSE1', directory, sources)
        packet = source_packet(Path(acquired['metadata']), 'GSE1', 'Source-defined intervention comparisons',
            context=acquired['study_context'], comparison_scope='source_defined_interventions',
            sample_source=acquired['source_sample_records'])
        content = json.dumps(packet)
        for marker in ('DO_NOT_INCLUDE_DERIVATION_RECEIPTS_IN_SOURCE_DESIGN', 'matrix_provenance', 'cells_compared', 'historical_batch_pre_and_post'):
            self.assertNotIn(marker, content)

    def test_version_identity_and_zero_call_resume_include_nested_files(self):
        requests = []
        def executor(request_path, log, timeout):
            request = json.loads(request_path.read_text()); requests.append(request)
            directory = Path(request['pipeline_root']); directory.mkdir()
            (directory / 'input.txt').write_text('SOFTWARE FIXTURE ONLY')
            log.write_text('No model, clock, or preprocessing execution')
            update.write_new(Path(request['result_path']), {'execution_status': 'COMPLETED'})
            return 0
        plan = {'schema': update.SCHEMA, 'studies': [self.fixture.entry()], 'limits': {
            'max_studies': 1, 'max_input_bytes': 2000000, 'max_samples_per_study': 10, 'max_child_seconds': 10,
            'concurrency': 1, 'allow_network': False, 'max_total_model_calls': 10, 'max_total_cost_usd': 2},
            'comprehensive_audit': {'model': 'gpt-5.6-luna', 'max_calls': 2, 'max_cost_usd': 0.5,
                'input_usd_per_million': 0.5, 'output_usd_per_million': 1.8}}
        plan_path = self.root / 'plan.json'; self.fixture.write(plan_path, plan)
        with patch.object(update, 'implementation_hashes', return_value={'src/update_pipeline_2026_09_20.py': update.sha(update.__file__)}):
            first = update.run(plan_path, self.root / 'registry', 'input', executor=executor)
            second = update.run(plan_path, self.root / 'registry', 'input', executor=executor)
        self.assertEqual(first['studies'][0]['version_id'], second['studies'][0]['version_id'])
        self.assertEqual(len(requests), 1)
        version = update.read_sealed(Path(first['studies'][0]['version_root']) / 'version.json')
        self.assertEqual(version['identity']['source_sha256']['matrix_provenance:verification_script'], update.sha(self.fixture.verifier))
        self.assertIn('matrix_provenance:stage_receipt:seal', version['sources'])


class DerivedMatrixMeasurementTests(unittest.TestCase):
    def test_late_path_free_evidence_and_exact_nested_binding(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            pipeline, index, original = _fixture(root / 'pipeline-fixture')
            fixture = DerivedMatrixFixture(root, Path(original['matrix']['path']))
            sources = update.collect_sources(fixture.entry(), root / 'acquisition', {'max_input_bytes': 2000000, 'allow_network': False}, {})
            sources['fixture_asset'] = original['fixture_asset']
            index['study_id'] = 'GSE1'
            bundle = build_measurement_evidence(pipeline, 'GSE1', 'human_rna', index, sources)
            record = bundle['evidence']['M00004']
            self.assertEqual(record['kind'], 'historical_archs4_matrix_derivation')
            self.assertEqual(record['value']['cells_compared'], 4)
            self.assertEqual(bundle['scientific_status'], 'CANNOT VERIFY')
            self.assertNotIn(str(root), json.dumps(bundle['evidence']))
            omitted = dict(sources); omitted.pop('matrix_provenance:verification_script')
            with self.assertRaisesRegex(ValueError, 'nested sources differ from sealed version inventory'):
                build_measurement_evidence(pipeline, 'GSE1', 'human_rna', index, omitted)
            duplicate = root / 'duplicate-verifier.py'; shutil.copyfile(fixture.verifier, duplicate)
            foreign = copy.deepcopy(sources)
            foreign['matrix_provenance:verification_script']['path'] = str(duplicate)
            with self.assertRaisesRegex(ValueError, 'file lacks exact sealed version binding'):
                build_measurement_evidence(pipeline, 'GSE1', 'human_rna', index, foreign)

    def test_producer_independent_reader_and_wrapper_handoff(self):
        from audit_system.tests.test_archs4_batch_extraction_2026_09_21 import Archs4BatchExtractionTest
        fixture = Archs4BatchExtractionTest(); fixture.setUp(); self.addCleanup(fixture.doCleanups)
        profile_script = fixture.root / 'frozen-profile-script.py'; profile_script.write_text('# Software source-profile fixture\n')
        descriptor = json.loads(fixture.descriptor.read_text())
        descriptor['verification_script_sha256'] = update.sha(profile_script)
        fixture.descriptor.write_text(json.dumps(descriptor, sort_keys=True) + '\n')
        fixture.write_plan(fixture.plan_path, fixture.output)
        batch = fixture.run_basic()
        repository = Path(__file__).resolve().parents[2]
        source_verifier = repository / 'src/verify_ingestion_2026_09_20.py'  # released copy
        require_development_artifact(source_verifier)
        verifier_path = fixture.root / 'frozen-independent-verifier.py'; shutil.copyfile(source_verifier, verifier_path)
        spec = importlib.util.spec_from_file_location('independent_count_reader_fixture', verifier_path)
        verifier = importlib.util.module_from_spec(spec); spec.loader.exec_module(verifier)
        verified = verifier.verify_batch_extraction(fixture.output / 'batch-summary.json', fixture.release)
        self.assertEqual(verified['status'], 'PASS')
        self.assertEqual([row['status'] for row in verified['results']], ['PASS', 'PASS', 'CANNOT VERIFY'])
        verification_path = fixture.root / 'independent-readback.json'
        verification_path.write_text(json.dumps(verified, indent=2) + '\n')
        row = next(row for row in batch['results'] if row['study_id'] == 'GSE1')
        absolute = lambda declaration: {'path': str((fixture.output / declaration['path']).resolve()), 'sha256': declaration['sha256']}
        sources = {name: {'path': row[name]['path'], 'sha256': row[name]['sha256']} for name in ('series', 'samples')}
        sources['matrix'] = absolute(row['matrix'])
        wrapper = {'schema': update.MATRIX_PROVENANCE_SCHEMA, 'study_id': 'GSE1', 'modality': 'human_rna', **sources,
            'batch_summary': DerivedMatrixFixture.decl(fixture.output / 'batch-summary.json'),
            'extraction_receipt': absolute(row['extraction_receipt']),
            'verification_receipt': DerivedMatrixFixture.decl(verification_path),
            'source_profile_descriptor': DerivedMatrixFixture.decl(fixture.descriptor),
            'source_profile_script': DerivedMatrixFixture.decl(profile_script),
            'source_profile_record_sha256': row['source_profile_record_sha256'],
            'source_hdf5': {'sha256': batch['source_hdf5']['sha256'], 'size_bytes': batch['source_hdf5']['size_bytes'],
                            'hash_scope': 'historical_batch_pre_and_post'},
            'implementation': {role: absolute(value) for role, value in batch['implementation'].items()},
            'scientific_status': 'CANNOT VERIFY', 'original_paper_2025_eligibility': 'QUARANTINED'}
        wrapper_path = fixture.root / 'matrix-provenance.json'; wrapper_path.write_text(json.dumps(wrapper, indent=2) + '\n')
        sources['matrix_provenance'] = DerivedMatrixFixture.decl(wrapper_path)
        entry = {'study_id': 'GSE1', 'modality': 'human_rna', 'sources': sources}
        collected = update.collect_sources(entry, fixture.root / 'acquisition', {'max_input_bytes': 4000000, 'allow_network': False}, {})
        self.assertEqual(collected['matrix']['sha256'], row['matrix']['sha256'])
        self.assertEqual(collected['matrix_provenance:verification_script']['sha256'], verified['verification_script_sha256'])


class BatchUpdateHandoffTests(unittest.TestCase):
    """Tiny software batch exercises the producer and separate cell reader."""
    def setUp(self):
        from audit_system.tests.test_archs4_batch_extraction_2026_09_21 import Archs4BatchExtractionTest
        fixture = self.fixture = Archs4BatchExtractionTest()
        fixture.setUp(); self.addCleanup(fixture.doCleanups)
        self.root = fixture.root
        self.profile_script = self.root / 'frozen-profile-script.py'
        self.profile_script.write_text('# Software source-profile fixture\n')
        descriptor = json.loads(fixture.descriptor.read_text())
        descriptor['verification_script_sha256'] = update.sha(self.profile_script)
        fixture.descriptor.write_text(json.dumps(descriptor, sort_keys=True) + '\n')
        fixture.write_plan(fixture.plan_path, fixture.output)
        self.batch = fixture.run_basic()
        repository = Path(__file__).resolve().parents[2]
        source_verifier = repository / 'src/verify_ingestion_2026_09_20.py'  # released copy
        require_development_artifact(source_verifier)
        self.verifier_path = self.root / 'frozen-independent-verifier.py'
        shutil.copyfile(source_verifier, self.verifier_path)
        spec = importlib.util.spec_from_file_location('independent_batch_handoff_fixture', self.verifier_path)
        verifier = importlib.util.module_from_spec(spec); spec.loader.exec_module(verifier)
        self.checked = verifier.verify_batch_extraction(fixture.output / 'batch-summary.json', fixture.release)
        self.assertEqual(self.checked['status'], 'PASS')
        self.verification = self.root / 'independent-readback.json'
        self.verification.write_text(json.dumps(self.checked, indent=2) + '\n')
        self.plan = {'schema': update.BATCH_HANDOFF_SCHEMA, 'sources': {
            'batch_summary': DerivedMatrixFixture.decl(fixture.output / 'batch-summary.json'),
            'verification_receipt': DerivedMatrixFixture.decl(self.verification),
            'source_profile_script': DerivedMatrixFixture.decl(self.profile_script)},
            'limits': {'max_studies': 3, 'max_total_matrix_bytes': 1000000}}
        self.plan_path = self.root / 'handoff-plan.json'
        self.write_plan()
        self.output = self.root / 'handoff'

    def write_plan(self):
        self.plan_path.write_text(json.dumps(self.plan, indent=2) + '\n')

    def write_verification(self):
        self.verification.write_text(json.dumps(self.checked, indent=2) + '\n')
        self.plan['sources']['verification_receipt'] = DerivedMatrixFixture.decl(self.verification)
        self.write_plan()

    def run_handoff(self):
        # Empirical HDF5 release is unavailable to this operation by contract.
        original_sha = update.sha
        def protected_sha(path):
            self.assertNotEqual(Path(path).suffix.lower(), '.h5', 'Handoff must not reopen the historical HDF5')
            return original_sha(path)
        with patch.object(update, 'sha', side_effect=protected_sha), \
                patch.object(update, 'child', side_effect=AssertionError('No child execution')), \
                patch.object(update.urllib.request, 'urlopen', side_effect=AssertionError('No network')):
            return update.build_batch_handoff(self.plan_path, self.output)

    def test_end_to_end_handoff_retains_cohort_and_validates_candidates(self):
        result = self.run_handoff()
        self.assertEqual(result['counts']['enrolled_studies'], 3)
        self.assertEqual(result['counts']['candidate_studies'], 2)
        self.assertEqual(result['counts']['candidate_samples'], 3)
        self.assertEqual(result['counts']['candidate_cells_compared'], 6)
        self.assertEqual(result['counts']['dispositions'], {'READY_FOR_UPDATE': 2, 'SOURCE_PROFILE_UNRESOLVED': 1})
        self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(result['original_paper_2025_eligibility'], 'QUARANTINED')
        self.assertEqual(result['hdf5_files_opened'], 0)
        self.assertEqual(result['model_calls'], 0)
        ledger = update.read_sealed(self.output / 'candidate-ledger.json')['studies']
        self.assertEqual([row['study_id'] for row in ledger], ['GSE1', 'GSE2', 'GSE3'])
        self.assertEqual([row['audit_status'] for row in ledger], ['NOT_YET_AUDITED'] * 3)
        self.assertIsNone(ledger[2]['candidate'])
        self.assertTrue(ledger[2]['reason'])
        candidates = update.read_sealed(self.output / 'candidates.json')
        self.assertEqual(candidates['schema'], 'clockbase-update-candidates-v1')
        self.assertEqual([row['study_id'] for row in candidates['studies']], ['GSE1', 'GSE2'])
        for row in candidates['studies']:
            sources = update.collect_sources(row, self.root / 'unused-acquisition',
                                             {'max_input_bytes': 4000000, 'allow_network': False}, {})
            self.assertIn('matrix_provenance:verification_script', sources)
            self.assertNotIn('archs4', sources)
            self.assertEqual(set(row['sources']), {'series', 'samples', 'matrix', 'matrix_provenance'})
        self.assertFalse((self.root / 'unused-acquisition').exists())
        from audit_system.agent_workflow_2026_09_10 import source_packet
        directory = self.root / 'blinded-metadata'; directory.mkdir()
        acquired = update.frozen_metadata('GSE2', directory, sources)
        packet = source_packet(Path(acquired['metadata']), 'GSE2', 'Source-defined intervention comparisons',
            context=acquired['study_context'], comparison_scope='source_defined_interventions',
            sample_source=acquired['source_sample_records'])
        for marker in ('matrix_provenance', 'cells_compared', 'historical_batch_pre_and_post', 'READY_FOR_UPDATE'):
            self.assertNotIn(marker, json.dumps(packet))

    def test_resume_verifies_unchanged_complete_inventory_without_rewriting(self):
        first = self.run_handoff()
        inventory = {str(p.relative_to(self.output)): (update.sha(p), p.stat().st_mtime_ns)
                     for p in self.output.rglob('*') if p.is_file()}
        second = self.run_handoff()
        self.assertEqual(first, second)
        self.assertEqual(inventory, {str(p.relative_to(self.output)): (update.sha(p), p.stat().st_mtime_ns)
                                   for p in self.output.rglob('*') if p.is_file()})

    def test_resume_rejects_matrix_nested_source_and_implementation_tampering(self):
        self.run_handoff()
        paths = [self.fixture.output / self.batch['results'][0]['matrix']['path'],
                 self.fixture.output / self.batch['results'][2]['stop_receipt']['path'], self.verifier_path]
        for path in paths:
            with self.subTest(path=path.name):
                original = path.read_bytes()
                path.write_bytes(original + b'\n')
                with self.assertRaisesRegex(ValueError, 'byte count|SHA256 differs'): self.run_handoff()
                path.write_bytes(original)

    def test_resume_rejects_changed_candidate_wrapper_and_extra_output(self):
        self.run_handoff()
        for name in ('candidates.json', 'wrappers/human_rna-GSE1.json'):
            with self.subTest(name=name):
                path = self.output / name; original = path.read_bytes()
                path.write_bytes(original + b'\n')
                with self.assertRaisesRegex(ValueError, 'completed handoff output changed'): self.run_handoff()
                path.write_bytes(original)
        (self.output / 'extra.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'completed output inventory changed'): self.run_handoff()

    def test_resealed_resume_cannot_add_foreign_source_or_scientific_acceptance(self):
        self.run_handoff()
        path = self.output / 'handoff.json'; original = json.loads(path.read_text())
        foreign = self.root / 'not-consumed.json'; foreign.write_text('{}')
        for mutation, error in [('source', 'source inventory changed'), ('science', 'scope changed')]:
            with self.subTest(mutation=mutation):
                changed = copy.deepcopy(original)
                if mutation == 'source': changed['input_sha256'][str(foreign)] = update.sha(foreign)
                else: changed['scientific_status'] = 'PASS'
                path.write_text(json.dumps(changed))
                path.with_suffix('.json.sha256').write_text(update.sha(path) + '\n')
                with self.assertRaisesRegex(ValueError, error): self.run_handoff()

    def test_whole_readback_pass_is_required_and_no_partial_admission(self):
        self.checked['status'] = 'FAIL'; self.write_verification()
        with self.assertRaisesRegex(ValueError, 'whole independent every-cell readback must PASS'): self.run_handoff()
        self.assertFalse(self.output.exists())

    def test_missing_duplicate_reordered_and_foreign_verification_rows_fail(self):
        original = copy.deepcopy(self.checked)
        for mutation in ('missing', 'duplicate', 'order', 'foreign'):
            with self.subTest(mutation=mutation):
                self.checked = copy.deepcopy(original)
                if mutation == 'missing': self.checked['results'].pop()
                elif mutation == 'duplicate': self.checked['results'].append(copy.deepcopy(self.checked['results'][0]))
                elif mutation == 'order': self.checked['results'].reverse()
                else: self.checked['batch_summary']['path'] = str(self.root / 'foreign-summary.json')
                self.write_verification()
                with self.assertRaisesRegex(ValueError, 'inventory|different batch'): self.run_handoff()
                self.assertFalse(self.output.exists())

    def test_unresolved_readback_cannot_change_sample_identity(self):
        self.checked['results'][2]['sample_ids'] = ['GSM999']
        self.write_verification()
        with self.assertRaisesRegex(ValueError, 'unresolved readback sample inventory differs'): self.run_handoff()
        self.assertFalse(self.output.exists())

    def test_nonlocal_source_and_wrong_hash_fail_before_output(self):
        original = copy.deepcopy(self.plan)
        for change in ('url', 'hash', 'hdf5'):
            with self.subTest(change=change):
                self.plan = copy.deepcopy(original)
                source = self.plan['sources']['verification_receipt']
                if change == 'url': source['path'] = 'https://example.com/source.json'
                elif change == 'hash': source['sha256'] = '0' * 64
                else: source['path'] = str(self.fixture.release)
                self.write_plan()
                with self.assertRaisesRegex(ValueError, 'local path|SHA256 differs|must not open'): self.run_handoff()
                self.assertFalse(self.output.exists())

    def test_matrix_and_study_limits_fail_without_partial_candidate_output(self):
        for name, value in [('max_studies', 2), ('max_total_matrix_bytes', 1)]:
            with self.subTest(name=name):
                previous = self.plan['limits'][name]
                self.plan['limits'][name] = value; self.write_plan()
                with self.assertRaisesRegex(ValueError, 'exceeds limit|exceeds aggregate limit'): self.run_handoff()
                self.assertFalse(self.output.exists())
                self.plan['limits'][name] = previous

    def test_unfinished_output_is_preserved(self):
        self.output.mkdir(); retained = self.output / 'operator-note.txt'; retained.write_text('retain me\n')
        with self.assertRaisesRegex(ValueError, 'unfinished handoff retained'): self.run_handoff()
        self.assertEqual(retained.read_text(), 'retain me\n')
        self.assertFalse((self.output / 'candidates.json').exists())

    def test_provenance_rejection_is_retained_without_candidate(self):
        original = update.read_matrix_provenance
        def reject_one(path, expected, sources, study, modality):
            if study == 'GSE1': raise update.UpdateIntegrityError('Deliberate software contract rejection')
            return original(path, expected, sources, study, modality)
        with patch.object(update, 'read_matrix_provenance', side_effect=reject_one):
            result = self.run_handoff()
        self.assertEqual(result['counts']['dispositions'], {
            'PROVENANCE_UNRESOLVED': 1, 'READY_FOR_UPDATE': 1, 'SOURCE_PROFILE_UNRESOLVED': 1})
        candidates = update.read_sealed(self.output / 'candidates.json')['studies']
        self.assertEqual([row['study_id'] for row in candidates], ['GSE2'])
        ledger = update.read_sealed(self.output / 'candidate-ledger.json')['studies']
        self.assertIn('Deliberate software contract rejection', ledger[0]['reason'])
        self.assertEqual(ledger[0]['scientific_status'], 'CANNOT VERIFY')

    def test_cli_handoff_route_and_mode_exclusion(self):
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(update.main(['--batch-handoff', str(self.plan_path), '--output', str(self.output)]), 0)
            self.assertEqual(update.main(['--batch-handoff', str(self.plan_path), '--plan', str(self.plan_path),
                                          '--output', str(self.output)]), 1)
        self.assertTrue((self.output / 'candidates.json').is_file())


if __name__ == '__main__':
    unittest.main()
