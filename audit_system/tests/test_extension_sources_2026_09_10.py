"""Software fixtures test release/identity handling; no extension discoveries."""
from pathlib import Path
import json
import gzip
import tempfile
import unittest
from unittest.mock import patch
import urllib.error
import urllib.parse
import io
import signal
import time

import h5py
from audit_system.archs4_update_2026_09_10 import sample_inventory,remaining_storage
from audit_system.extension_eligibility_2026_09_10 import (series_record, parse_family_source, collect_family_sources,
    build_family_manifest, prepare_family_acquisition, FamilyAcquirer, acquire_family_sources, ACQUISITION_VERSION)
from audit_system.evaluation_2026_09_10 import build_families
from audit_system.report_2026_09_10 import verify_study_families
from audit_system.agent_workflow_2026_09_10 import sha


def acquisition_fixture(root, *, max_requests=30, max_series=5, old_relative=False, ambiguous=False):
    """Software-only identifiers/CSV values; never empirical evaluation inputs."""
    records = [{'accession': 'GSE1', 'source_id': 'fixture:1', 'source_kind': 'sample_metadata',
                'sample_ids': ['GSM1'], 'unresolved': []}]
    if old_relative:
        records.append({'accession': 'GSE9', 'source_id': 'fixture:old-primary', 'source_kind': 'geo_series_primary',
                        'sample_ids': ['GSM9'], 'pubmed_ids': ['10'], 'unresolved': []})
    registry = build_families(records, ['GSE9'])
    def pin(name, data):
        path = root / name
        path.write_text(json.dumps(data))
        return {'path': str(path), 'sha256': sha(path)}
    atlas = root / 'software-fixture.csv'
    atlas.write_text('gse_id,effect_size,verdict\nGSE1,99999,SOFTWARE_FIXTURE\nGSE9,-99999,SOFTWARE_FIXTURE\n')
    if ambiguous:
        with atlas.open('a') as handle: handle.write('GSE1/GSE9,0,SOFTWARE_FIXTURE\n')
    plan = {'schema': ACQUISITION_VERSION, 'registry': pin('registry.json', registry),
            'source_records': pin('records.json', records), 'exposure_addenda': [],
            'atlas_sources': [{'path': str(atlas), 'sha256': sha(atlas), 'modality': 'software_fixture',
                              'clock_root': str(root / 'absent-clocks')}],
            'family_limit': 1, 'max_series': max_series, 'max_requests': max_requests, 'max_response_bytes': 2000000,
            'max_reverse_ids': 256, 'gsm_query_batch': 64, 'timeout_seconds': 30, 'request_interval_seconds': .35}
    path = root / 'plan.json'
    path.write_text(json.dumps(plan))
    result = prepare_family_acquisition(path, root / 'selection')
    return Path(result['selection']), plan


class FixtureResponse(io.BytesIO):
    status = 200
    url = 'https://www.ncbi.nlm.nih.gov/software-fixture'
    headers = {}


class FamilyAcquisitionTests(unittest.TestCase):
    def test_selection_retains_missing_clock_input_and_excludes_exposed_families(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            selected, plan = acquisition_fixture(root)
            value = json.loads(selected.read_text())
            self.assertEqual(value['seed_studies'], ['GSE1'])
            self.assertFalse(value['clock_availability_used_for_selection'])
            self.assertFalse(value['candidate_studies'][0]['clock_files'][0]['exists'])
            self.assertEqual(value['candidate_studies'][1]['screen'], 'EXPOSED')
            self.assertNotIn('99999', selected.read_text())
            (root / 'software-fixture.csv').write_text('gse_id,effect_size\nGSE1,0\n')
            with self.assertRaisesRegex(ValueError, 'Atlas source changed'):
                prepare_family_acquisition(root / 'plan.json', root / 'changed')

    def test_search_truncation_warnings_and_summary_identity_are_not_admission(self):
        obj = object.__new__(FamilyAcquirer)
        obj.plan = {'max_reverse_ids': 256}
        for result in ({'count': '2', 'idlist': ['1']},
                       {'count': '1', 'idlist': ['1'], 'warninglist': {'phrasesnotfound': ['fixture']}},
                       {'count': '2', 'idlist': ['1', '1']}):
            with patch.object(obj, 'json_request', return_value=({'esearchresult': result}, {})):
                with self.assertRaises(ValueError): obj.search('gds', 'software fixture')
        obj.plan['max_reverse_ids'] = 1
        with patch.object(obj, 'json_request', return_value=({'esearchresult': {'count': '2', 'idlist': ['1', '2']}}, {})):
            with self.assertRaises(ValueError): obj.search('gds', 'software fixture')
        good = {'result': {'uids': ['123'], '123': {'uid': '123', 'accession': 'GDS7', 'gse': '88'}}}
        with patch.object(obj, 'json_request', return_value=(good, {})):
            self.assertEqual(obj.summaries(['123'])[0], ['GSE88'])
        good['result']['123']['gse'] = '0'
        with patch.object(obj, 'json_request', return_value=(good, {})):
            with self.assertRaisesRegex(ValueError, 'explicit Series identity'): obj.summaries(['123'])
        good['result']['123'] = {'uid': '124', 'accession': 'GSE88'}
        with patch.object(obj, 'json_request', return_value=(good, {})):
            with self.assertRaisesRegex(ValueError, 'identity/error'): obj.summaries(['123'])

    def test_failed_requests_resume_without_retry_and_enforce_transport_budget(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            selected, _ = acquisition_fixture(root, max_requests=1)
            obj = FamilyAcquirer(root / 'run', selected)
            failure = urllib.error.HTTPError('https://fixture', 503, 'fixture', {}, io.BytesIO(b'software fixture failure'))
            with patch('urllib.request.urlopen', side_effect=failure) as transport:
                first = obj.fetch('series', 'https://fixture/one')
                second = obj.fetch('series', 'https://fixture/one')
                exhausted = obj.fetch('series', 'https://fixture/two')
            self.assertEqual(transport.call_count, 1)
            self.assertEqual(first, second)
            self.assertEqual(first['status'], 'CANNOT VERIFY')
            self.assertIn('allocation exhausted', exhausted['error'])
            Path(first['response_path']).write_bytes(b'changed fixture')
            with self.assertRaisesRegex(ValueError, 'Completed stage output changed'):
                obj.fetch('series', 'https://fixture/one')

    def test_interrupted_request_is_retained_without_automatic_repeat(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            selected, _ = acquisition_fixture(root)
            obj = FamilyAcquirer(root / 'run', selected)
            with patch('urllib.request.urlopen', side_effect=KeyboardInterrupt):
                with self.assertRaises(KeyboardInterrupt): obj.fetch('series', 'https://fixture/interrupted')
            with patch('urllib.request.urlopen', side_effect=AssertionError('network repeat')):
                result = obj.fetch('series', 'https://fixture/interrupted')
            self.assertEqual(result['status'], 'CANNOT VERIFY')
            self.assertIn('Interrupted prior request', result['error'])

    def test_truncated_http_body_cannot_pass_as_a_complete_series(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            selected, _ = acquisition_fixture(root)
            obj = FamilyAcquirer(root / 'run', selected)
            response = FixtureResponse(b'^SERIES = GSE1\n!Series_status = Public on Jan 01 2020\n!Series_sample_id = GSM1\n')
            response.headers = {'Content-Length': str(len(response.getvalue()) + 10)}
            with patch('urllib.request.urlopen', return_value=response):
                result = obj.fetch('series', 'https://fixture/truncated')
            self.assertEqual(result['status'], 'CANNOT VERIFY')
            self.assertIn('Content-Length', result['error'])

    def test_pending_relative_with_old_primary_cannot_complete_current_closure(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            selected, _ = acquisition_fixture(root, max_series=1, old_relative=True)
            payload = b'^SERIES = GSE1\n!Series_status = Public on Jan 01 2020\n!Series_sample_id = GSM1\n!Series_pubmed_id = 10\n'
            def reverse(kind, ids):
                return (['GSE1', 'GSE9'] if kind == 'PubMed' else ['GSE1']), []
            with patch('urllib.request.urlopen', return_value=FixtureResponse(payload)), patch.object(FamilyAcquirer, 'reverse', side_effect=reverse):
                result = acquire_family_sources(selected, root / 'run')
            self.assertEqual(result['pending_series'], ['GSE9'])
            self.assertEqual(result['selected_family_linkage'], {'CANNOT VERIFY': 1})

    def test_changed_selection_or_family_dependency_rejects_before_network(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            selected, _ = acquisition_fixture(root)
            with patch('audit_system.extension_eligibility_2026_09_10.acquisition_dependencies', return_value={'changed': 'fixture'}):
                with self.assertRaisesRegex(ValueError, 'implementation changed'):
                    FamilyAcquirer(root / 'run', selected)
            value = json.loads(selected.read_text())
            value['plan']['max_requests'] = 999999
            selected.write_text(json.dumps(value))
            with self.assertRaisesRegex(ValueError, 'selection artifact changed'):
                FamilyAcquirer(root / 'run', selected)

    def test_registry_and_pinned_source_graph_must_agree_before_selection(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            _, plan = acquisition_fixture(root)
            registry = root / 'registry.json'
            data = json.loads(registry.read_text())
            data['studies'].append(dict(data['studies'][0], study_id='GSE999'))
            registry.write_text(json.dumps(data))
            plan['registry']['sha256'] = sha(registry)
            (root / 'changed-plan.json').write_text(json.dumps(plan))
            with self.assertRaisesRegex(ValueError, 'differs from its pinned'):
                prepare_family_acquisition(root / 'changed-plan.json', root / 'other-selection')

    def test_missing_and_compound_atlas_ids_remain_in_screening_without_invented_links(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            _, plan = acquisition_fixture(root)
            atlas = root / 'software-fixture.csv'
            with atlas.open('a') as handle:
                handle.write(',0,SOFTWARE_FIXTURE\nGSE1/GSE9,0,SOFTWARE_FIXTURE\n')
            plan['atlas_sources'][0]['sha256'] = sha(atlas)
            (root / 'identity-plan.json').write_text(json.dumps(plan))
            result = prepare_family_acquisition(root / 'identity-plan.json', root / 'identity-selection')
            value = json.loads(Path(result['selection']).read_text())
            self.assertEqual(value['seed_studies'], ['GSE1'])
            self.assertEqual(result['unresolved_atlas_identity_rows'], 2)
            self.assertEqual(value['unresolved_atlas_identity_rows'][1]['mentioned_accessions'], ['GSE1', 'GSE9'])
            self.assertEqual(len(value['candidate_studies'][0]['unresolved_identity_rows']), 1)

    def test_total_transport_deadline_restores_signal_state(self):
        from audit_system.extension_eligibility_2026_09_10 import acquisition_deadline
        previous = signal.getsignal(signal.SIGALRM)
        with self.assertRaisesRegex(TimeoutError, 'total wall-clock deadline'):
            with acquisition_deadline(.02):
                time.sleep(.2)
        self.assertEqual(signal.getsignal(signal.SIGALRM), previous)
        self.assertEqual(signal.getitimer(signal.ITIMER_REAL), (0., 0.))

    def test_compound_identity_blocks_admission_without_creating_a_family_edge(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            selected, _ = acquisition_fixture(root, ambiguous=True)
            payload = b'^SERIES = GSE1\n!Series_status = Public on Jan 01 2020\n!Series_sample_id = GSM1\n'
            with patch('urllib.request.urlopen', return_value=FixtureResponse(payload)), patch.object(FamilyAcquirer, 'reverse', return_value=(['GSE1'], [])):
                result = acquire_family_sources(selected, root / 'run')
            self.assertEqual(result['selected_family_linkage'], {'CANNOT VERIFY': 1})
            self.assertEqual(result['selected_family_splits'], {'UNRESOLVED': 1})
            selected_family = json.loads(next((root / 'run').glob('zz_manifest/attempt-*/selected-families.json')).read_text())
            self.assertEqual(selected_family[0]['members'], ['GSE1'])

    def test_primary_and_reverse_index_closure_propagates_exposure_and_resumes(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            selected, _ = acquisition_fixture(root)
            def transport(request, **_):
                url = urllib.parse.urlsplit(request.full_url)
                params = urllib.parse.parse_qs(url.query)
                if url.path.endswith('acc.cgi'):
                    study = params['acc'][0]
                    sample = 'GSM1' if study == 'GSE1' else 'GSM9'
                    body = f'^SERIES = {study}\n!Series_status = Public on Jan 01 2020\n!Series_sample_id = {sample}\n!Series_pubmed_id = 10\n'.encode()
                elif url.path.endswith('elink.fcgi'):
                    body = json.dumps({'linksets': [{'dbfrom': 'pubmed', 'ids': ['10'], 'linksetdbs': [
                        {'dbto': 'gds', 'linkname': 'pubmed_gds', 'links': ['101', '109']}]}]}).encode()
                elif url.path.endswith('esearch.fcgi'):
                    uid = '101' if 'GSM1[ACCN]' in params['term'][0] else '109'
                    body = json.dumps({'esearchresult': {'count': '1', 'idlist': [uid]}}).encode()
                else:
                    ids = params['id'][0].split(',')
                    data = {'uids': ids}
                    for uid in ids: data[uid] = {'uid': uid, 'accession': 'GSE1' if uid == '101' else 'GSE9'}
                    body = json.dumps({'result': data}).encode()
                return FixtureResponse(body)
            with patch('urllib.request.urlopen', side_effect=transport), patch('time.sleep'):
                result = acquire_family_sources(selected, root / 'run')
            self.assertEqual(result['parsed_primary_series'], 2)
            self.assertEqual(result['selected_family_splits'], {'development': 1})
            self.assertEqual(result['selected_family_linkage'], {'PASS': 1})
            self.assertEqual(result['evaluation_release'], 'CLOSED')
            before = {str(p): sha(p) for p in (root / 'run').rglob('*') if p.is_file()}
            with patch('urllib.request.urlopen', side_effect=AssertionError('network repeat')):
                after = acquire_family_sources(selected, root / 'run')
            self.assertEqual(result, after)
            self.assertEqual(before, {str(p): sha(p) for p in (root / 'run').rglob('*') if p.is_file()})


class ExtensionSourceTests(unittest.TestCase):
    def test_metadata_links_use_exact_ids_and_preserve_failed_sources(self):
        payload=b'GSM,geo_accession,title\nGSM1,GSM1,fixture\nGSM2,GSM2,fixture\n'
        record=parse_family_source(Path('GSE1.csv'),'sample_metadata',payload)
        self.assertEqual(record['sample_ids'],['GSM1','GSM2'])
        self.assertNotIn('title',record)
        for bad in (payload.replace(b'GSM2,GSM2',b'GSM2,GSM3'), payload+b'GSM3\n', b'GSM\n'):
            with self.assertRaises(ValueError):parse_family_source(Path('GSE1.csv'),'sample_metadata',bad)
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); metadata=root/'metadata';metadata.mkdir()
            (metadata/'GSE1.csv').write_bytes(payload)
            (metadata/'GSE2.csv').write_text('title\nmissing sample IDs\n')
            records,sources,execution=collect_family_sources([metadata],[],[],root/'cache')
            self.assertEqual(execution['parsed'],2)
            self.assertEqual(sum(s['status']=='CANNOT VERIFY' for s in sources),1)
            self.assertEqual({r['accession'] for r in records},{'GSE1','GSE2'})

    def test_family_source_cache_resume_checks_input_and_cached_content(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'GSE1.csv';source.write_text('GSM\nGSM1\n')
            before=collect_family_sources([source],[],[],root/'cache')
            after=collect_family_sources([source],[],[],root/'cache')
            self.assertEqual(before[:2],after[:2]);self.assertEqual(after[2]['resumed'],1)
            source.write_text('GSM\nGSM2\n')
            with self.assertRaisesRegex(ValueError,'source changed during resume'):
                collect_family_sources([source],[],[],root/'cache')
            source.write_text('GSM\nGSM1\n')
            cached=next((root/'cache').glob('*/*.json'))
            changed=json.loads(cached.read_text());changed['payload']['record']['sample_ids']=['GSM2']
            cached.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError,'Cached family source changed'):
                collect_family_sources([source],[],[],root/'cache')

    def test_family_manifest_blocks_partial_sources_and_keeps_exposure_provenance(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);metadata=root/'metadata';metadata.mkdir()
            for accession in ('GSE1','GSE2'):(metadata/(accession+'.csv')).write_text('GSM\nGSM1\n')
            exposed=root/'software-fixture-exposure.json';exposed.write_text('{"study":"GSE1"}')
            result=build_family_manifest([metadata],[],[],[exposed],root/'cache',root/'out')
            self.assertEqual(result['directly_exposed_studies'],1)
            self.assertEqual(result['family_exposed_studies'],2)
            self.assertEqual(result['split_counts'],{'development':2})
            self.assertFalse(result['source_coverage_complete'])
            self.assertEqual(result['evaluation_release'],'CLOSED')
            verified=verify_study_families(root/'out',root/'verified')
            self.assertEqual(verified['independent_graph_components'],1)
            self.assertEqual(verified['additional_family_exclusions'],1)
            family_path=root/'out/families.json.gz'
            with gzip.open(family_path,'rt') as handle:family=json.load(handle)
            family['edges'][0]['identifier']='GSM999'
            with gzip.open(family_path,'wt') as handle:json.dump(family,handle)
            result['output_sha256']['families.json.gz']=sha(family_path)
            (root/'out/manifest.json').write_text(json.dumps(result))
            with self.assertRaisesRegex(ValueError,'witness is not present'):
                verify_study_families(root/'out',root/'wrong-witness')

    def test_primary_and_archived_series_links_have_different_source_assurance(self):
        payload=b'^SERIES = GSE1\n!Series_status = Public on Jan 01 2025\n!Series_sample_id = GSM1\n!Series_pubmed_id = 12345\n!Series_relation = BioProject: PRJNA10\n!Series_relation = SubSeries of: GSE2\n'
        primary=parse_family_source(Path('GSE1.soft'),'geo_series_primary',payload)
        self.assertEqual(primary['linked_accessions'],['GSE2'])
        self.assertEqual(primary['bioproject_ids'],['PRJNA10'])
        self.assertEqual(primary['pubmed_ids'],['12345'])
        self.assertEqual(primary['unresolved'],[])
        parsed=series_record('GSE1',payload)
        receipt=parse_family_source(Path('GSE1.json'),'geo_series_receipt',json.dumps(parsed).encode())
        self.assertTrue(receipt['unresolved'])
    def test_geo_series_release_is_distinct_from_modification(self):
        payload=b'^SERIES = GSE123\n!Series_status = Public on Jan 15 2025\n!Series_last_update_date = Sep 09 2026\n!Series_sample_id = GSM1\n!Series_sample_id = GSM2\n!Series_relation = BioProject: PRJNA123\n'
        record=series_record('GSE123',payload)
        self.assertEqual(record['series_first_public_date'],'2025-01-15')
        self.assertEqual(record['sample_ids'],['GSM1','GSM2'])
        self.assertEqual(record['relations'],['BioProject: PRJNA123'])
        for bad in [payload.replace(b'Public on Jan 15 2025',b'Private'),
                    payload+b'!Series_sample_id = GSM1\n',payload.replace(b'GSE123',b'GSE124')]:
            with self.assertRaises(ValueError):series_record('GSE123',bad)

    def test_archs4_preserves_all_series_memberships_and_rejects_duplicate_samples(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'software_fixture.h5'
            with h5py.File(path,'w') as h5:
                h5.create_dataset('meta/samples/geo_accession',data=[b'GSM1',b'GSM2'])
                h5.create_dataset('meta/samples/series_id',data=[b'GSE1,GSE2',b'GSE2'])
            self.assertEqual(sample_inventory(path)['GSM1'],['GSE1','GSE2'])
            with h5py.File(path,'r+') as h5:h5['meta/samples/geo_accession'][1]=b'GSM1'
            with self.assertRaises(ValueError):sample_inventory(path)

    def test_sparse_partial_size_does_not_erase_download_storage_reserve(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'human_gene_v2.latest.h5.part'
            size=4*1024*1024
            with path.open('wb') as handle:handle.seek(size-1);handle.write(b'0')
            remaining=remaining_storage({'human':{'file_size':size}},Path(directory))
            self.assertEqual(remaining,size-min(size,path.stat().st_blocks*512))
            if path.stat().st_blocks*512<size:self.assertGreater(remaining,0)


if __name__=='__main__':unittest.main()
