"""Bounded bootstrap software fixtures; no empirical inputs or performance claims."""
import copy
import csv
import fcntl
import gzip
import hashlib
import io
import json
from pathlib import Path
import socket
import shutil
import threading
import unittest
from unittest.mock import patch

import h5py

from audit_system.tests import test_future_refresh_2026_09_21 as refresh_tests
from src import pipeline_inputs_2026_09_10 as sources
from src import update_pipeline_2026_09_20 as update


class RefreshBootstrapTests(unittest.TestCase):
    def setUp(self):
        self.fixture = refresh_tests.FutureRefreshTests(); self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root; self.plan = self.fixture.plan
        self.original_release = self.fixture.fixture.release
        self.payload = self.original_release.read_bytes()
        with h5py.File(self.original_release, 'r') as handle:
            samples = int(handle['data/expression'].shape[1])
        self.catalog = [{'id': 1, 'checksum': 'software-opaque-checksum-v1',
                         'species': 'human', 'data_level': 'gene', 'version_major': 2,
                         'version_minor': 'latest', 'timestamp': 'Mon, 01 Jan 2024 00:00:00 GMT',
                         'file_size': len(self.payload), 'samples': samples}]
        self.plan.pop('archs4'); self.calls = []
        self.plan['network_budget'] = {'max_requests': 100, 'max_bytes': 10000000}
        self.plan['source_acquisition'] = {'policy': 'official_archs4_v2_latest',
            'cache_dir': str(self.root/'immutable-cache'),
            'limits': {'max_requests': 4, 'max_bytes': 1000000, 'max_seconds': 30,
                       'max_source_bytes': 1000000, 'max_metadata_bytes': 100000, 'min_free_bytes': 1},
            'budget': {'max_seconds': 90, 'max_cache_bytes': 10000000}}
        self.fixture.write_plan()

    def transport(self, method, url, headers, limit):
        self.calls.append((method, url, headers, limit))
        if url == 'https://archs4.org/api/versionfile':
            self.assertEqual(method, 'GET')
            return {'Content-Type': 'application/json'}, json.dumps(self.catalog).encode()
        self.assertEqual(url, 'https://s3.k8s.maayanlab.cloud/archs4/files/human_gene_v2.latest.h5')
        result = {'ETag': '"software-fixture-etag"', 'Content-Length': str(len(self.payload))}
        if method == 'HEAD': return result, b''
        self.assertEqual(headers, {'If-Match': '"software-fixture-etag"'})
        return result, self.payload

    def execute(self, until='input', acquisition_transport=None, executor=None):
        self.assertTrue((update.ROOT / update.REFRESH_VERIFIER).is_file())
        return update.refresh(self.fixture.path, self.fixture.output, until,
            transport=self.fixture.transport, executor=executor or self.fixture.software_child,
            acquisition_transport=acquisition_transport or self.transport)

    def terminal(self):
        return update.read_sealed(self.fixture.output/'source-acquisitions'/update.identity(self.plan)/'result.json')

    @staticmethod
    def dnam_fixture(*, extra_relation=False):
        """Disclosed two-sample software fixture; no empirical measurements."""
        ids = ['GSM1', 'GSM2']
        def row(name, values): return name+'\t'+'\t'.join('"'+v+'"' for v in values)+'\n'
        metadata = [row('!Series_geo_accession', ['GSE40279']),
                    row('!Series_relation', ['BioProject: software-fixture'])]
        if extra_relation: metadata.append(row('!Series_relation', ['Added: software-fixture']))
        metadata += ['\n', row('!Sample_geo_accession', ids),
                     row('!Sample_platform_id', ['GPL13534']*2),
                     row('!Sample_organism_ch1', ['Homo sapiens']*2),
                     row('!Sample_molecule_ch1', ['genomic DNA']*2)]
        table = row('ID_REF', ids)+row('cg1', ['0.1', '0.2'])+row('cg2', ['0.3', '0.4'])
        text = ''.join(metadata)+'!series_matrix_table_begin\n'+table+'!series_matrix_table_end\n'
        parsed = list(csv.reader(io.StringIO(''.join(metadata), newline=''), delimiter='\t', strict=True))
        canonical = [record for record in parsed if not record or record[0] != '!Series_relation']
        header = ['ID_REF']+ids
        profile = {**sources.DNAM_REFERENCE_PROFILE,
            'ordered_canonical_metadata_record_count': len(canonical),
            'ordered_canonical_metadata_sha256': hashlib.sha256(json.dumps(
                canonical, ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode()).hexdigest(),
            'table_bytes_including_header_excluding_end_marker': len(table.encode()),
            'table_sha256': hashlib.sha256(table.encode()).hexdigest(),
            'feature_rows': 2,
            'ordered_header_sha256': hashlib.sha256(json.dumps(
                header, ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode()).hexdigest(),
            'sample_count': 2}
        return gzip.compress(text.encode(), mtime=0), profile, text

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_orchestrated_input_pipeline_and_report_resume_use_no_new_transport(self):
        def real_input(request_path, log, timeout):
            request = json.loads(request_path.read_text()); self.assertEqual(request['until'], 'input')
            update.child(request); return 0
        result = self.execute(executor=real_input)
        self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 2, 'SOURCE_PROFILE_UNRESOLVED': 1})
        self.assertEqual(len(self.calls), 4)
        self.assertEqual(result['reserved_network_requests'], 24)
        self.assertEqual(result['reserved_network_bytes'], 2000000)
        saved = self.terminal(); self.assertEqual(saved['status'], 'COMPLETED')
        self.assertEqual(saved['requests_used'], 4)
        self.assertEqual(saved['response_bytes'], len(self.payload)+len(json.dumps(self.catalog).encode()))
        self.assertEqual(saved['resolved_source']['sha256'], update.sha(self.original_release))
        first = {str(p): update.sha(p) for p in self.fixture.output.rglob('*') if p.is_file()}
        with patch.object(socket.socket, 'connect', side_effect=AssertionError('No network on resume')):
            resumed = update.refresh(self.fixture.path, self.fixture.output, 'input',
                transport=lambda *x: self.fail('Discovery on resume'), executor=lambda *x: self.fail('Child on resume'),
                acquisition_transport=lambda *x: self.fail('Source acquisition on resume'))
        self.assertEqual(result, resumed)
        self.assertEqual(first, {str(p): update.sha(p) for p in self.fixture.output.rglob('*') if p.is_file()})
        report = self.execute('report')
        self.assertEqual(len(self.calls), 4)
        self.assertEqual(report['dispositions']['AUDIT_EXECUTION_COMPLETED'], 2)

    def test_new_revision_reuses_verified_object_and_preserves_origin_receipt(self):
        self.execute(); old = self.terminal(); old_count = len(self.calls)
        self.plan['revision'] = 'explicit-second-release-check'; self.fixture.write_plan(); self.execute()
        new = self.terminal()
        self.assertEqual(new['requests_used'], 2)
        self.assertEqual(len(self.calls), old_count+2)
        self.assertEqual(new['resolved_source'], old['resolved_source'])
        self.assertTrue(any('/source-acquisitions/' in path and '/003/receipt.json' in path
                            for path in new['external_sha256']))

    def test_catalog_id_and_opaque_checksum_separately_change_selected_cache_object(self):
        self.execute(); first = self.terminal()
        paths = {first['resolved_source']['path']}
        original = copy.deepcopy(self.catalog[0])
        for field, value in [('id', 2), ('checksum', 'software-opaque-checksum-v2')]:
            with self.subTest(changed_catalog_field=field):
                self.catalog[0] = {**original, field: value}
                self.plan['revision'] = 'catalog-change-' + field
                self.fixture.write_plan(); count = len(self.calls); self.execute()
                saved = self.terminal(); source = Path(saved['resolved_source']['path'])
                self.assertEqual(len(self.calls) - count, 4)
                self.assertEqual(saved['requests_used'], 4)
                self.assertNotIn(str(source), paths); paths.add(str(source))
                self.assertEqual(saved['resolved_source']['sha256'], first['resolved_source']['sha256'])
                selection = json.loads((self.fixture.output/'source-acquisitions'/update.identity(self.plan)
                                       /'selection/human-selection.json').read_text())
                receipt = json.loads((source.parent/'human_download.json').read_text())
                self.assertEqual(selection, receipt['declaration'])
                self.assertEqual(selection['id'], self.catalog[0]['id'])
                self.assertEqual(selection['checksum'], self.catalog[0]['checksum'])
                self.assertEqual(selection['catalog_payload_checksum_verification']['status'], 'CANNOT VERIFY')

    def test_missing_catalog_identity_stops_before_head_without_implicit_retry(self):
        self.catalog[0].pop('checksum')
        with self.assertRaisesRegex(RuntimeError, 'lacks usable id/checksum'): self.execute()
        self.assertEqual([method for method, *_ in self.calls], ['GET'])
        before = len(self.calls)
        with self.assertRaisesRegex(RuntimeError, 'no implicit retry'): self.execute()
        self.assertEqual(len(self.calls), before)
        self.assertNotIn('resolved_source', self.terminal())

    def test_receipts_exclude_cookie_and_session_headers_but_keep_verification_metadata(self):
        excluded = {'Set-Cookie': 'software-private-cookie-sentinel',
                    'Authorization': 'software-private-auth-sentinel',
                    'X-Session-Token': 'software-private-session-sentinel',
                    'Server': 'software-server-sentinel'}
        def transport(*args):
            headers, payload = self.transport(*args)
            return {**headers, **excluded, 'lAsT-mOdIfIeD': 'Mon, 01 Jan 2024 00:00:00 GMT',
                    'cOnTeNt-TyPe': 'application/octet-stream'}, payload
        self.execute(acquisition_transport=transport)
        request_root = self.fixture.output/'source-acquisitions'/update.identity(self.plan)/'requests'
        receipts = sorted(request_root.glob('*/receipt.json'))
        self.assertEqual(len(receipts), 4)
        for path in receipts:
            value = update.read_sealed(path); headers = value['response_headers']
            self.assertEqual(value['response_header_policy'], 'source-verification-allowlist-v1')
            self.assertTrue(set(headers) <= {'etag', 'content-length', 'last-modified', 'content-type'})
            self.assertEqual(headers['last-modified'], 'Mon, 01 Jan 2024 00:00:00 GMT')
            self.assertEqual(headers['content-type'], 'application/octet-stream')
            if value['url'] != 'https://archs4.org/api/versionfile':
                self.assertEqual(headers['etag'], '"software-fixture-etag"')
                self.assertEqual(headers['content-length'], str(len(self.payload)))
            for sentinel in excluded.values(): self.assertNotIn(sentinel, path.read_text())
        with patch.object(socket.socket, 'connect', side_effect=AssertionError('No network on resume')):
            self.execute(acquisition_transport=lambda *args: self.fail('No source transport on resume'))

    def test_changed_cache_and_origin_receipt_are_rejected_before_transport(self):
        self.execute(); saved = self.terminal(); source = Path(saved['resolved_source']['path'])
        source.write_bytes(source.read_bytes()+b'changed')
        with self.assertRaisesRegex(update.UpdateIntegrityError, 'Cached source'):
            self.execute(acquisition_transport=lambda *x: self.fail('Network after cache change'))

    def test_origin_receipt_tamper_is_not_hidden_by_cache_hit(self):
        self.execute(); saved = self.terminal()
        origin = next(Path(path) for path in saved['external_sha256'] if '/003/receipt.json' in path and not path.endswith('.sha256'))
        origin.write_text('{}')
        with self.assertRaisesRegex(update.UpdateIntegrityError, 'Acquisition output changed|Cached source'):
            self.execute(acquisition_transport=lambda *x: self.fail('Network after origin change'))

    def test_new_revision_rejects_resealed_cache_receipt_against_original_acquisition(self):
        self.execute(); source = Path(self.terminal()['resolved_source']['path'])
        path = next(source.parent.glob('*_download.json')); value = json.loads(path.read_text())
        value['source_etag_verified'] = False; path.write_text(json.dumps(value))
        Path(str(path)+'.sha256').write_text(update.sha(path)+'\n')
        self.plan['revision'] = 'after-cache-receipt-drift'; self.fixture.write_plan()
        with self.assertRaisesRegex(RuntimeError, 'original completed acquisition'): self.execute()

    def test_ambiguous_catalog_retains_failure_and_never_retries(self):
        self.catalog.append(copy.deepcopy(self.catalog[0]))
        with self.assertRaisesRegex(RuntimeError, 'ambiguous'): self.execute()
        self.assertEqual(self.terminal()['requests_used'], 1)
        self.assertEqual(self.terminal()['status'], 'UNRESOLVED')
        before = len(self.calls)
        with self.assertRaisesRegex(RuntimeError, 'no implicit retry'): self.execute()
        self.assertEqual(len(self.calls), before)
        reservation = next((self.fixture.output/'network-reservations').glob('*-source.json'))
        self.assertEqual(update.read_sealed(reservation)['max_bytes'], 1000000)

    def test_changed_post_download_etag_stops_without_source_admission(self):
        def transport(*args):
            headers, payload = self.transport(*args)
            if len(self.calls) == 4: headers['ETag'] = '"changed"'
            return headers, payload
        with self.assertRaisesRegex(RuntimeError, 'source changed'): self.execute(acquisition_transport=transport)
        self.assertNotIn('resolved_source', self.terminal())
        self.assertTrue(list((self.root/'immutable-cache').rglob('*.part')))

    def test_wrong_hdf_sample_dimension_is_rejected(self):
        self.catalog[0]['samples'] += 1
        with self.assertRaisesRegex(RuntimeError, 'dimensions differ'): self.execute()
        self.assertEqual(self.terminal()['status'], 'UNRESOLVED')

    def test_metadata_overflow_and_unknown_transport_failure_keep_reservations(self):
        self.plan['source_acquisition']['limits']['max_metadata_bytes'] = 4; self.fixture.write_plan()
        with self.assertRaisesRegex(RuntimeError, 'response cap'): self.execute()
        saved = self.terminal()
        self.assertEqual(saved['response_bytes'], 5)
        self.assertEqual(saved['requests_used'], 1)
        self.assertEqual(saved['status'], 'UNRESOLVED')

    def test_uncertain_failure_accounts_full_request_cap(self):
        def stopped(*args): raise TimeoutError('Software transport failure')
        with self.assertRaisesRegex(RuntimeError, 'Software transport failure'): self.execute(acquisition_transport=stopped)
        self.assertEqual(self.terminal()['response_bytes'], 100001)
        self.assertEqual(self.terminal()['requests_used'], 1)

    def test_declared_release_cap_stops_before_head_or_download(self):
        self.plan['source_acquisition']['limits']['max_source_bytes'] = 1; self.fixture.write_plan()
        with self.assertRaisesRegex(RuntimeError, 'source-byte limit'): self.execute()
        self.assertEqual(len(self.calls), 1)

    def test_storage_stop_precedes_payload_download(self):
        self.plan['source_acquisition']['budget']['max_cache_bytes'] = 10000; self.fixture.write_plan()
        with self.assertRaisesRegex(RuntimeError, 'storage cap'): self.execute()
        self.assertEqual([m for m, *_ in self.calls], ['GET', 'HEAD'])

    def test_free_space_and_output_lock_stop_before_network(self):
        usage = update.shutil.disk_usage(self.root)
        with patch.object(update.shutil, 'disk_usage', return_value=usage._replace(free=0)):
            with self.assertRaisesRegex(RuntimeError, 'storage/free-space'): self.execute()
        self.assertEqual(self.calls, [])
        with (self.fixture.output/'.refresh.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(RuntimeError, 'Another process'): self.execute()
        self.assertEqual(self.calls, [])

    def test_global_network_and_time_caps_survive_new_revisions(self):
        self.plan['source_acquisition']['budget']['max_seconds'] = 30; self.fixture.write_plan(); self.execute()
        self.plan['revision'] = 'second'; self.fixture.write_plan(); before = len(self.calls)
        with self.assertRaisesRegex(RuntimeError, 'Lifetime source acquisition time'): self.execute()
        self.assertEqual(len(self.calls), before)

    def test_source_network_reservation_precedes_transport(self):
        self.plan['network_budget']['max_bytes'] = 999999; self.fixture.write_plan()
        with self.assertRaisesRegex(RuntimeError, 'network reservation'): self.execute()
        self.assertEqual(self.calls, [])

    def test_cumulative_source_and_discovery_byte_reservations_are_not_reset(self):
        self.plan['network_budget']['max_bytes'] = 4000000; self.fixture.write_plan()
        self.execute(); self.plan['revision'] = 'second'; self.fixture.write_plan(); self.execute()
        self.plan['revision'] = 'third'; self.fixture.write_plan(); before = len(self.calls)
        with self.assertRaisesRegex(RuntimeError, 'network reservation'): self.execute()
        self.assertEqual(len(self.calls), before)

    def test_budget_or_code_drift_on_completed_acquisition_rejects_without_network(self):
        self.execute(); before = len(self.calls)
        implementations = update.implementation_hashes()
        implementations['src/pipeline_inputs_2026_09_10.py'] = '0'*64
        with patch.object(update, 'implementation_hashes', return_value=implementations):
            with self.assertRaisesRegex(update.UpdateIntegrityError, 'input/code/runtime'): self.execute()
        self.assertEqual(len(self.calls), before)
        self.plan['source_acquisition']['budget']['max_seconds'] += 1; self.fixture.write_plan()
        with self.assertRaisesRegex(update.UpdateIntegrityError, 'lifetime limits'): self.execute()
        self.assertEqual(len(self.calls), before)

    def test_incomplete_acquisition_is_terminal_and_keeps_reservation(self):
        self.catalog.append(copy.deepcopy(self.catalog[0]))
        with self.assertRaisesRegex(RuntimeError, 'ambiguous'): self.execute()
        terminal = self.fixture.output/'source-acquisitions'/update.identity(self.plan)/'result.json'
        terminal.unlink(); Path(str(terminal)+'.sha256').unlink()
        before = len(self.calls)
        with self.assertRaisesRegex(RuntimeError, 'Interrupted source acquisition'): self.execute()
        self.assertEqual(len(self.calls), before)
        self.assertEqual(self.terminal()['status'], 'UNRESOLVED')

    def test_missing_release_cap_fails_before_any_acquisition(self):
        self.plan.pop('max_release_bytes'); self.fixture.write_plan()
        with self.assertRaises(ValueError): self.execute()
        self.assertEqual(self.calls, [])

    def test_live_deadline_stops_transport_without_retry(self):
        self.plan['source_acquisition']['limits']['max_seconds'] = 1; self.fixture.write_plan()
        def slow(*args): threading.Event().wait(3); self.fail('Acquisition deadline was not enforced')
        with self.assertRaisesRegex(RuntimeError, 'wall-time'): self.execute(acquisition_transport=slow)
        self.assertEqual(self.terminal()['requests_used'], 1)

    def test_local_and_acquisition_sources_cannot_be_combined(self):
        self.plan['archs4'] = {'path': str(self.original_release), 'sha256': update.sha(self.original_release)}
        self.fixture.write_plan()
        with self.assertRaisesRegex(ValueError, 'exclusive'): self.execute()
        self.assertEqual(self.calls, [])

    def test_dnam_bootstrap_verifies_profile_identity_without_means(self):
        payload, profile, _ = self.dnam_fixture(extra_relation=True)
        checksum = hashlib.sha256(payload).hexdigest()
        self.plan['modality'] = 'human_dnam'; self.plan.pop('max_release_bytes')
        self.plan['source_acquisition']['policy'] = 'geo_gse40279'; self.fixture.ids = []; self.fixture.query_total = 0
        self.plan['source_acquisition']['limits']['max_source_bytes'] = 100000
        self.fixture.write_plan()
        def transport(method, url, headers, limit):
            self.assertEqual(url, sources.DNAM_REFERENCE_URL); self.assertEqual(method, 'GET')
            self.calls.append((method, url, headers, limit)); return {}, payload
        with patch.object(sources, 'DNAM_REFERENCE_PROFILE', profile), \
                patch.object(sources, 'reference_means', side_effect=AssertionError('No means in bootstrap')):
            self.execute(acquisition_transport=transport)
            self.assertEqual(len(self.calls), 1); self.assertFalse(self.terminal()['means_computed'])
            self.assertEqual(self.terminal()['role'], 'dnam_reference')
            reference = Path(self.terminal()['resolved_source']['path'])
            receipt = json.loads((reference.parent/'reference.json').read_text())
            self.assertEqual(receipt['schema'], 'clockbase-bounded-dnam-reference-v2')
            self.assertEqual(receipt['source_sha256'], checksum)
            self.assertEqual(len(receipt['relation_records']), 2)
            self.assertEqual(receipt['profile_sha256'], sources._dnam_reference_profile_sha256())
            update.refresh(self.fixture.path, self.fixture.output, 'input', transport=lambda *x: self.fail('No discovery on resume'),
                acquisition_transport=lambda *x: self.fail('No reference on resume'))

    def test_dnam_malformed_gzip_is_rejected_before_means(self):
        self.plan['modality'] = 'human_dnam'; self.plan.pop('max_release_bytes')
        self.plan['source_acquisition']['policy'] = 'geo_gse40279'
        self.plan['source_acquisition']['limits']['max_source_bytes'] = 100000; self.fixture.write_plan()
        with patch.object(sources, 'reference_means', side_effect=AssertionError('No means in bootstrap')):
            with self.assertRaisesRegex(RuntimeError, 'gzip/table is malformed'):
                self.execute(acquisition_transport=lambda *x: ({}, b'WRNG'))
        self.assertEqual(self.terminal()['requests_used'], 1)
        self.assertNotIn('resolved_source', self.terminal())

    def test_dnam_header_identity_is_checked_independently_of_profile_hash(self):
        _, profile, text = self.dnam_fixture()
        payload = gzip.compress(text.replace('"GSE40279"', '"GSE999"').encode(), mtime=0)
        self.plan['modality'] = 'human_dnam'; self.plan.pop('max_release_bytes')
        self.plan['source_acquisition']['policy'] = 'geo_gse40279'
        self.plan['source_acquisition']['limits']['max_source_bytes'] = 100000; self.fixture.write_plan()
        with patch.object(sources, 'DNAM_REFERENCE_PROFILE', profile):
            with self.assertRaisesRegex(RuntimeError, 'identity/sample/platform'):
                self.execute(acquisition_transport=lambda *x: ({}, payload))

    def test_dnam_relation_only_drift_preserves_complete_records(self):
        original, profile, _ = self.dnam_fixture()
        changed, _, _ = self.dnam_fixture(extra_relation=True)
        with patch.object(sources, 'DNAM_REFERENCE_PROFILE', profile):
            observations = []
            for index, payload in enumerate((original, changed)):
                path = self.root/('reference-'+str(index)+'.gz'); path.write_bytes(payload)
                observations.append(sources.verify_dnam_reference_source(path))
        self.assertNotEqual(observations[0]['source_sha256'], observations[1]['source_sha256'])
        self.assertEqual(observations[0]['table_sha256'], observations[1]['table_sha256'])
        self.assertEqual(len(observations[0]['relation_records']), 1)
        self.assertEqual(len(observations[1]['relation_records']), 2)

    def test_dnam_profile_rejects_changed_values_ids_metadata_and_table_bounds(self):
        _, profile, text = self.dnam_fixture()
        variants = {
            'value': text.replace('"0.4"', '"0.5"'),
            'feature': text.replace('cg2', 'cg9'),
            'sample': text.replace('"GSM2"', '"GSM3"', 1),
            'header_order': text.replace('ID_REF\t"GSM1"\t"GSM2"', 'ID_REF\t"GSM2"\t"GSM1"'),
            'row_order': text.replace('cg1\t"0.1"\t"0.2"\ncg2\t"0.3"\t"0.4"',
                                      'cg2\t"0.3"\t"0.4"\ncg1\t"0.1"\t"0.2"'),
            'scientific_metadata': text.replace('"Homo sapiens"', '"Mus musculus"', 1),
            'blank_metadata_removed': text.replace('\n\n!Sample_geo_accession', '\n!Sample_geo_accession'),
            'blank_metadata_changed': text.replace('\n\n!Sample_geo_accession', '\n!Changed_metadata\n!Sample_geo_accession'),
            'extra_row': text.replace('!series_matrix_table_end', 'cg3\t"0.1"\t"0.2"\n!series_matrix_table_end'),
            'missing_end': text.replace('!series_matrix_table_end\n', ''),
            'changed_end_marker': text.replace('!series_matrix_table_end\n', '!series_matrix_table_end \n'),
            'trailing': text+'extra\n',
        }
        with patch.object(sources, 'DNAM_REFERENCE_PROFILE', profile):
            for name, body in variants.items():
                with self.subTest(name=name):
                    path = self.root/(name+'.gz'); path.write_bytes(gzip.compress(body.encode(), mtime=0))
                    with self.assertRaises(ValueError): sources.verify_dnam_reference_source(path)

    def test_dnam_profile_enforces_decompression_metadata_and_line_bounds(self):
        payload, profile, _ = self.dnam_fixture()
        path = self.root/'bounded.gz'; path.write_bytes(payload)
        with patch.object(sources, 'DNAM_REFERENCE_PROFILE', profile):
            for field, value in [('DNAM_REFERENCE_MAX_DECOMPRESSED_BYTES', 20),
                                 ('DNAM_REFERENCE_MAX_METADATA_BYTES', 20),
                                 ('DNAM_REFERENCE_MAX_LINE_BYTES', 20)]:
                with self.subTest(field=field), patch.object(sources, field, value):
                    with self.assertRaisesRegex(ValueError, 'bound|limit'): sources.verify_dnam_reference_source(path)

    def test_dnam_profile_reports_missing_source_and_invalid_cap_separately(self):
        missing = self.root/'missing-GSE40279.gz'
        with self.assertRaisesRegex(FileNotFoundError, 'Cannot proceed - required data file not found: '+str(missing)):
            sources.verify_dnam_reference_source(missing)
        payload, profile, _ = self.dnam_fixture()
        path = self.root/'present-GSE40279.gz'; path.write_bytes(payload)
        with patch.object(sources, 'DNAM_REFERENCE_PROFILE', profile):
            with self.assertRaisesRegex(ValueError, 'positive source-byte limit'):
                sources.verify_dnam_reference_source(path, max_source_bytes=0)

    def test_dnam_v2_cache_hit_and_legacy_local_path_use_same_profile(self):
        payload, profile, _ = self.dnam_fixture(extra_relation=True)
        cache = self.root/'dnam-cache'; calls = []
        def request(method, url, destination, maximum, headers):
            calls.append((method, url, maximum)); destination.write_bytes(payload)
            return {'receipt': {'path': str(self.root/'fixture-request.json'), 'sha256': '0'*64}}
        with patch.object(sources, 'DNAM_REFERENCE_PROFILE', profile):
            reference = sources.acquire_dnam_reference(cache, request=request, max_source_bytes=1000000)
            self.assertEqual(len(calls), 1)
            self.assertIn('gse40279-profile-v2-', str(reference))
            self.assertEqual(reference, sources.acquire_dnam_reference(cache,
                request=lambda *args: self.fail('No request on cache hit'), max_source_bytes=1000000))
            receipt = reference.parent/'reference.json'; original = json.loads(receipt.read_text())
            for field, value in [('sha256', '0'*64), ('reference_study', 'GSE999')]:
                with self.subTest(resealed_alias=field):
                    changed = {**original, field: value}; receipt.write_text(json.dumps(changed)+'\n')
                    Path(str(receipt)+'.sha256').write_text(sources.sha(receipt)+'\n')
                    with self.assertRaises(ValueError): sources.acquire_dnam_reference(cache,
                        request=lambda *args: self.fail('No request after resealed alias'), max_source_bytes=1000000)
            receipt.write_text(json.dumps(original)+'\n')
            Path(str(receipt)+'.sha256').write_text(sources.sha(receipt)+'\n')
            local = self.root/'legacy'/reference.name; local.parent.mkdir(); local.write_bytes(payload)
            self.assertEqual(local, sources.acquire_dnam_reference(local.parent, max_source_bytes=1000000))
            local.write_bytes(gzip.compress(b'changed', mtime=0))
            with self.assertRaises(ValueError): sources.acquire_dnam_reference(local.parent, max_source_bytes=1000000)
            reference.write_bytes(gzip.compress(b'changed', mtime=0))
            with self.assertRaises(ValueError): sources.acquire_dnam_reference(cache,
                request=lambda *args: self.fail('No request after tamper'), max_source_bytes=1000000)


if __name__ == '__main__': unittest.main()
