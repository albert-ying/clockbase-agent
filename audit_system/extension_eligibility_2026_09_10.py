"""Verify GEO release evidence and overlap without evaluating biological results.

Retains one source receipt per accession. Acquisition dates are provenance,
never temporal eligibility. Unknown family overlap prevents final admission.
"""
import argparse
import fcntl
from collections import Counter
import csv
from contextlib import contextmanager
import gzip
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import re
import signal
import time
import urllib.request
import urllib.parse
import urllib.error
import sys
import tempfile

from .agent_workflow_2026_09_10 import sha, write_new

FAMILY_SOURCE_VERSION = 'family-source-2026-09-16-v1'
ACQUISITION_VERSION = 'family-acquisition-2026-09-17-v1'


def acquisition_dependencies():
    root = Path(__file__).resolve().parents[1]
    return {name: sha(root / name) for name in ('audit_system/extension_eligibility_2026_09_10.py',
        'audit_system/evaluation_2026_09_10.py', 'audit_system/agent_workflow_2026_09_10.py',
        'src/pipeline_runtime_2026_09_10.py')}


@contextmanager
def acquisition_deadline(seconds):
    """Bound the complete HTTP exchange on the supported Linux main thread."""
    if signal.getitimer(signal.ITIMER_REAL) != (0., 0.):
        raise RuntimeError('Source acquisition cannot replace an existing process deadline')
    def expired(_signum, _frame):
        raise TimeoutError('Source HTTP exchange exceeded its total wall-clock deadline')
    previous = signal.signal(signal.SIGALRM, expired)
    try:
        signal.setitimer(signal.ITIMER_REAL, seconds)
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.)
        signal.signal(signal.SIGALRM, previous)


def pinned_json(item):
    path = Path(item['path']).resolve()
    if not path.is_file():
        raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
    if sha(path) != item['sha256']:
        raise ValueError('Pinned source changed: ' + str(path))
    if path.suffix == '.gz':
        with gzip.open(path, 'rt') as handle:
            return json.load(handle)
    return json.loads(path.read_text())


def prepare_family_acquisition(plan_path, output):
    """Select by atlas IDs and existing family hashes, never by numerical results."""
    from src.pipeline_runtime_2026_09_10 import StageStore
    from .evaluation_2026_09_10 import build_families
    plan_path = Path(plan_path).resolve()
    plan = json.loads(plan_path.read_text())
    if plan.get('schema') != ACQUISITION_VERSION:
        raise ValueError('Unrecognized family acquisition plan')
    for key in ('family_limit', 'max_series', 'max_requests', 'max_response_bytes', 'max_reverse_ids',
                'gsm_query_batch', 'timeout_seconds'):
        if type(plan.get(key)) is not int or plan[key] <= 0:
            raise ValueError('Positive integer acquisition limit required: ' + key)
    if not .35 <= plan.get('request_interval_seconds', 0) <= 60:
        raise ValueError('Sequential NCBI requests require at least 0.35 seconds between starts')
    if plan['max_reverse_ids'] > 10000 or plan['gsm_query_batch'] > 100:
        raise ValueError('Reverse-query request exceeds supported retrieval bounds')
    registry = pinned_json(plan['registry'])
    records = pinned_json(plan['source_records'])
    if registry.get('schema') != 'source-study-families-v1' or registry['evaluation_release'] != 'CLOSED':
        raise ValueError('Source acquisition requires a closed family registry')
    if build_families(records, registry['exposed_studies'], registry['seed']) != registry:
        raise ValueError('Family registry differs from its pinned membership/linkage records')
    exposed = set(registry['exposed_studies'])
    for item in plan.get('exposure_addenda', []):
        values = pinned_json(item)[item['field']]
        if not isinstance(values, list) or any(not re.fullmatch(r'GSE\d+', value) for value in values):
            raise ValueError('Exposure addendum must name exact study identities')
        exposed.update(values)
    families = {family['family_id']: family for family in registry['families']}
    owners = {row['study_id']: row['family_id'] for row in registry['studies']}
    if len(owners) != len(registry['studies']):
        raise ValueError('Duplicate registry study identity')
    # Validate the family identity instead of treating a stored split as a decision.
    for key, family in families.items():
        expected = hashlib.sha256((str(registry['seed']) + '|' + ';'.join(sorted(family['members']))).encode()).hexdigest()
        if key != 'family:' + expected or any(owners.get(s) != key for s in family['members']):
            raise ValueError('Registry family identity or membership is inconsistent')
    candidates = {}
    atlas_counts = []
    identity_failures = []
    for item in plan['atlas_sources']:
        path = Path(item['path']).resolve()
        if sha(path) != item['sha256']:
            raise ValueError('Atlas source changed: ' + str(path))
        ids = set()
        # Only gse_id is retained. Clock, score, effect, P value and verdict
        # fields are neither interpreted nor used for selection.
        with path.open(newline='') as handle:
            reader = csv.reader(handle, strict=True)
            fields = next(reader)
            if fields.count('gse_id') != 1:
                raise ValueError('Atlas needs one exact gse_id column')
            column = fields.index('gse_id')
            for number, row in enumerate(reader, 2):
                if len(row) != len(fields):
                    raise ValueError('Malformed atlas row width; identity projection cannot safely continue')
                if not re.fullmatch(r'GSE\d+', row[column]):
                    identity_failures.append({'source_path': str(path), 'source_sha256': item['sha256'],
                        'csv_record_number': number, 'gse_id': row[column], 'status': 'CANNOT VERIFY',
                        'reason': 'Missing or non-single study identity; no biological-family relation inferred',
                        'mentioned_accessions': sorted(set(re.findall(r'(?<![A-Za-z0-9])GSE\d+(?![A-Za-z0-9])', row[column])))})
                    continue
                ids.add(row[column])
        if sha(path) != item['sha256']:
            raise ValueError('Atlas source changed during identity projection')
        available = 0
        for accession in sorted(ids):
            clock = Path(item['clock_root']).resolve() / (accession + '.csv')
            exists = clock.is_file()
            available += exists
            candidates.setdefault(accession, {'study_id': accession, 'modalities': [], 'clock_files': []})
            candidates[accession]['modalities'].append(item['modality'])
            candidates[accession]['clock_files'].append({'path': str(clock), 'exists': exists})
        atlas_counts.append({'modality': item['modality'], 'study_ids': len(ids), 'clock_files_present': available})
    eligible = set()
    for row in candidates.values():
        key = owners.get(row['study_id'])
        row['prior_family_id'] = key
        row['family_exposed'] = bool(key and set(families[key]['members']) & exposed)
        row['screen'] = 'CANNOT VERIFY' if key is None else 'EXPOSED' if row['family_exposed'] else 'UNEXPOSED'
        row['unresolved_identity_rows'] = [failure for failure in identity_failures if row['study_id'] in failure['mentioned_accessions']]
        if row['screen'] == 'UNEXPOSED':
            eligible.add(key)
    chosen = sorted(eligible)[:plan['family_limit']]
    seeds = sorted({s for key in chosen for s in families[key]['members']})
    if not chosen or len(seeds) > plan['max_series']:
        raise ValueError('Selected source families are empty or exceed the series allocation')
    config = {'schema': ACQUISITION_VERSION, 'plan_sha256': sha(plan_path), 'plan': plan,
              'implementation_sha256': sha(__file__), 'dependency_sha256': acquisition_dependencies()}
    store = StageStore(output, config)

    def save(directory):
        result = {'schema': ACQUISITION_VERSION, 'plan': plan, 'plan_sha256': sha(plan_path),
                  'implementation_sha256': sha(__file__), 'exposed_studies': sorted(exposed),
                  'dependency_sha256': config['dependency_sha256'],
                  'selected_prior_family_ids': chosen, 'seed_studies': seeds,
                  'selection_rule': 'First unexposed observed families in ascending SHA-256 family order; all members retained',
                  'clock_availability_used_for_selection': False,
                  'candidate_studies': sorted(candidates.values(), key=lambda row: row['study_id']),
                  'unresolved_atlas_identity_rows': identity_failures,
                  'atlas_counts': atlas_counts, 'eligible_observed_families': len(eligible),
                  'evaluation_release': 'CLOSED', 'numerical_results_inspected': False}
        write_new(directory / 'selection.json', result)
        return {'selection': str(directory / 'selection.json'), 'selection_sha256': sha(directory / 'selection.json'),
                'selected_prior_families': len(chosen), 'seed_studies': len(seeds),
                'candidate_studies': len(candidates), 'eligible_observed_families': len(eligible),
                'unresolved_atlas_identity_rows': len(identity_failures),
                'atlas_counts': atlas_counts, 'evaluation_release': 'CLOSED'}
    return store.run('00_selection', save)


class FamilyAcquirer:
    """Sequential primary-source requests with immutable receipts and bounded closure."""
    def __init__(self, output, selection_path):
        from src.pipeline_runtime_2026_09_10 import StageStore
        self.selection_path = Path(selection_path).resolve()
        stage = self.selection_path.parent.parent
        complete = stage / 'complete.json'
        if sha(complete) != (stage / 'complete.sha256').read_text().strip():
            raise ValueError('Selection preparation receipt changed')
        prepared = json.loads(complete.read_text())
        if (prepared['outputs'].get(str(self.selection_path.relative_to(stage))) != sha(self.selection_path) or
                prepared['result']['selection_sha256'] != sha(self.selection_path)):
            raise ValueError('Prepared selection artifact changed')
        self.selection = json.loads(self.selection_path.read_text())
        self.plan = self.selection['plan']
        if (self.selection.get('schema') != ACQUISITION_VERSION or self.selection['implementation_sha256'] != sha(__file__) or
                self.selection.get('dependency_sha256') != acquisition_dependencies()):
            raise ValueError('Selection version or implementation changed; prepare a new acquisition')
        self.root = Path(output).resolve()
        self.store = StageStore(self.root, {'schema': ACQUISITION_VERSION,
            'selection_path': str(self.selection_path), 'selection_sha256': sha(self.selection_path),
            'implementation_sha256': sha(__file__), 'dependency_sha256': acquisition_dependencies(),
            'paid_model_calls': 0, 'evaluation_release': 'CLOSED'})
        self.last_request = 0.
        self.request_log = []

    def fetch(self, kind, url):
        key = hashlib.sha256((kind + '|' + url).encode()).hexdigest()
        name = 'request_' + key
        folder = self.root / name

        def acquire(directory):
            record = {'kind': kind, 'url': url, 'status': 'CANNOT VERIFY', 'automatic_retries': 0}
            prior = [p for p in folder.glob('attempt-*/request.json') if p.parent != directory]
            if prior:
                record['error'] = 'Interrupted prior request retained; no automatic repeat'
                record['prior_requests'] = [str(p) for p in prior]
                write_new(directory / 'receipt.json', record)
                return record
            count = sum(1 for _ in self.root.glob('request_*/attempt-*/request.json'))
            if count >= self.plan['max_requests']:
                record['error'] = 'Frozen request allocation exhausted before transport'
                write_new(directory / 'receipt.json', record)
                return record
            time.sleep(max(0., self.plan['request_interval_seconds'] - (time.monotonic() - self.last_request)))
            record['retrieved_at_utc'] = datetime.now(timezone.utc).isoformat()
            write_new(directory / 'request.json', {'kind': kind, 'url': url, 'started_at_utc': record['retrieved_at_utc']})
            self.last_request = time.monotonic()
            started = self.last_request
            body = None
            try:
                request = urllib.request.Request(url, headers={'User-Agent': 'ClockBaseAudit-source-acquisition/1'})
                with acquisition_deadline(self.plan['timeout_seconds']):
                    try:
                        response = urllib.request.urlopen(request, timeout=self.plan['timeout_seconds'])
                    except urllib.error.HTTPError as error:
                        response = error
                    with response:
                        record['http_status'] = response.status
                        record['final_url'] = response.url
                        declared = response.headers.get('Content-Length')
                        record['declared_response_bytes'] = int(declared) if declared is not None else None
                        body = response.read(self.plan['max_response_bytes'] + 1)
                if len(body) > self.plan['max_response_bytes']:
                    raise ValueError('Response exceeds frozen byte ceiling; retained prefix is incomplete')
                if declared is not None and len(body) != record['declared_response_bytes']:
                    raise ValueError('Response ended before its declared Content-Length; retained prefix is incomplete')
                if record['http_status'] != 200:
                    raise ValueError('HTTP ' + str(record['http_status']))
                record['status'] = 'PASS'
            except Exception as error:
                record['error'] = f'{type(error).__name__}: {error}'
            if body is not None:
                path = directory / ('response.soft' if kind == 'series' else 'response.json')
                with path.open('xb') as handle:
                    handle.write(body)
                record.update(response_path=str(path), response_sha256=sha(path), response_bytes=len(body))
            record['elapsed_seconds'] = time.monotonic() - started
            write_new(directory / 'receipt.json', record)
            return record
        record = self.store.run(name, acquire)
        self.request_log.append(record)
        return record

    def json_request(self, endpoint, parameters):
        parameters = dict(parameters, retmode='json', tool='ClockBaseAudit')
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' + endpoint + '?' + urllib.parse.urlencode(parameters)
        receipt = self.fetch(endpoint, url)
        if receipt['status'] != 'PASS':
            raise ValueError(receipt['error'])
        data = json.loads(Path(receipt['response_path']).read_text())
        if data.get('error') or data.get('ERROR'):
            raise ValueError('NCBI query error: ' + str(data.get('error') or data.get('ERROR')))
        return data, receipt

    def search(self, db, term):
        data, receipt = self.json_request('esearch.fcgi', {'db': db, 'term': term, 'retmax': self.plan['max_reverse_ids']})
        result = data['esearchresult']
        ids = result.get('idlist', [])
        if result.get('errorlist') or result.get('warninglist') or result.get('ERROR'):
            raise ValueError('NCBI search has an error or warning; no completeness claim')
        if (int(result['count']) != len(ids) or len(ids) > self.plan['max_reverse_ids'] or len(ids) != len(set(ids)) or
                any(not re.fullmatch(r'\d+', str(uid)) for uid in ids)):
            raise ValueError('NCBI search is truncated or has invalid/duplicate identifiers')
        return ids, [receipt]

    def links(self, db, uid):
        data, receipt = self.json_request('elink.fcgi', {'dbfrom': db, 'db': 'gds', 'id': uid, 'linkname': db + '_gds'})
        sets = data['linksets']
        if len(sets) != 1 or sets[0].get('dbfrom') != db or list(map(str, sets[0].get('ids', []))) != [str(uid)]:
            raise ValueError('NCBI reverse-link source identity is ambiguous')
        if sets[0].get('error') or sets[0].get('ERROR'):
            raise ValueError('NCBI reverse-link error')
        groups = sets[0].get('linksetdbs', [])
        if any(g.get('dbto') != 'gds' or g.get('linkname') != db + '_gds' for g in groups) or len(groups) > 1:
            raise ValueError('Unexpected reverse-link namespace')
        ids = list(map(str, groups[0].get('links', []))) if groups else []
        if len(ids) > self.plan['max_reverse_ids'] or len(set(ids)) != len(ids) or any(not re.fullmatch(r'\d+', uid) for uid in ids):
            raise ValueError('Reverse-link set is invalid or exceeds the frozen retrieval allocation')
        return ids, [receipt]

    def summaries(self, ids):
        if not ids:
            return [], []
        data, receipt = self.json_request('esummary.fcgi', {'db': 'gds', 'id': ','.join(sorted(ids))})
        result = data['result']
        if len(result['uids']) != len(set(ids)) or set(result['uids']) != set(ids):
            raise ValueError('GEO summary response does not match requested identifiers')
        accessions = set()
        for uid in ids:
            row = result[uid]
            if str(row.get('uid')) != uid or row.get('error'):
                raise ValueError('GEO summary record identity/error')
            accession = row.get('accession', '')
            if re.fullmatch(r'GSE\d+', accession):
                accessions.add(accession)
            elif re.fullmatch(r'GDS\d+', accession) and re.fullmatch(r'[1-9]\d*', str(row.get('gse', ''))):
                # This is the explicit Series field in the source record,
                # never an accession inferred from UID order or dates.
                accessions.add('GSE' + str(row['gse']))
            else:
                raise ValueError('GEO summary has no supported explicit Series identity')
        return sorted(accessions), [receipt]

    def reverse(self, kind, identifiers):
        receipts = []
        if kind == 'GSM':
            term = '(' + ' OR '.join(sid + '[ACCN]' for sid in identifiers) + ') AND gse[ETYP]'
            ids, used = self.search('gds', term)
            receipts.extend(used)
        else:
            identifier = identifiers[0]
            db, uid = 'pubmed', identifier
            if kind == 'BioProject':
                db = 'bioproject'
                ids, used = self.search(db, identifier + '[PRJA]')
                receipts.extend(used)
                if len(ids) != 1:
                    raise ValueError('BioProject accession does not resolve to exactly one indexed record')
                uid = ids[0]
            ids, used = self.links(db, uid)
            receipts.extend(used)
        accessions, used = self.summaries(ids)
        receipts.extend(used)
        return accessions, receipts


def acquire_family_sources(selection_path, output):
    """Expand a frozen selection without inspecting clock values or study results."""
    from .evaluation_2026_09_10 import build_families
    acquirer = FamilyAcquirer(output, selection_path)
    with (acquirer.root / '.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another acquisition is using this output directory') from None
        records = pinned_json(acquirer.plan['source_records'])
        registry = pinned_json(acquirer.plan['registry'])
        owners = {sid: family['members'] for family in registry['families'] for sid in family['members']}
        seen, queue = set(), set(acquirer.selection['seed_studies'])
        new_records, series, queries = [], [], []
        for failure in acquirer.selection['unresolved_atlas_identity_rows']:
            source_id = 'atlas-identity:' + hashlib.sha256(json.dumps(failure, sort_keys=True).encode()).hexdigest()
            for accession in failure['mentioned_accessions']:
                new_records.append({'accession': accession, 'source_kind': 'atlas_identity_ambiguity',
                    'source_id': source_id, 'sample_ids': [], 'source_identity_failure': failure,
                    'unresolved': ['Atlas row has an unresolved compound study identity; related-output eligibility is unverified']})
        reverse_cache = {}
        while queue and len(seen) < acquirer.plan['max_series']:
            accession = min(queue, key=lambda sid: (sid not in acquirer.selection['seed_studies'], sid))
            queue.remove(accession)
            if accession in seen:
                continue
            seen.add(accession)
            queue.update(set(owners.get(accession, [])) - seen)
            url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&form=text&view=full'
            receipt = acquirer.fetch('series', url)
            entry = {'study_id': accession, 'receipt': receipt, 'status': 'CANNOT VERIFY'}
            source_id = 'acquired-series:' + accession
            record = {'accession': accession, 'source_id': source_id, 'sample_ids': [],
                      'source_kind': 'geo_series_primary', 'unresolved': []}
            try:
                if receipt['status'] != 'PASS':
                    raise ValueError(receipt['error'])
                payload = Path(receipt['response_path']).read_bytes()
                identity = series_record(accession, payload)
                record = parse_family_source(Path(receipt['response_path']), 'geo_series_primary', payload)
                record['source_id'] = source_id
                entry.update(status='PASS', first_public_date=identity['series_first_public_date'],
                             source_status_quote=identity['source_status_quote'], sample_count=len(record['sample_ids']))
                queue.update(set(record['linked_accessions']) - seen)
                groups = [('PubMed', [pid]) for pid in record['pubmed_ids']]
                groups += [('BioProject', [pid]) for pid in record['bioproject_ids']]
                batch = acquirer.plan['gsm_query_batch']
                groups += [('GSM', record['sample_ids'][i:i+batch]) for i in range(0, len(record['sample_ids']), batch)]
                for kind, ids in groups:
                    key = (kind, tuple(ids))
                    if key not in reverse_cache:
                        result = {'kind': kind, 'identifiers': ids, 'status': 'CANNOT VERIFY'}
                        request_start = len(acquirer.request_log)
                        try:
                            links, receipts = acquirer.reverse(kind, ids)
                            result.update(status='PASS', linked_studies=links, receipts=receipts)
                        except Exception as error:
                            result['error'] = f'{type(error).__name__}: {error}'
                        result['receipts'] = acquirer.request_log[request_start:]
                        reverse_cache[key] = result
                    query = dict(reverse_cache[key], source_study=accession)
                    if query['status'] == 'PASS' and accession not in query['linked_studies']:
                        query.update(status='CANNOT VERIFY', error='Reverse index omits the source study; completeness unverified')
                    queries.append(query)
                    if query['status'] != 'PASS':
                        record['unresolved'].append(kind + ' reverse linkage: ' + query['error'])
                    # Partial/anchor-missing index hits still require inspection;
                    # their absence from a validated family is never assumed.
                    linked = query.get('linked_studies', [])
                    queue.update(set(linked) - seen)
                    if kind in {'PubMed', 'BioProject'}:
                        for related in linked:
                            new_records.append({'accession': related, 'source_kind': 'ncbi_index_link',
                                'source_id': 'index:' + hashlib.sha256(json.dumps(query, sort_keys=True).encode()).hexdigest(),
                                'sample_ids': [], 'pubmed_ids': ids if kind == 'PubMed' else [],
                                'bioproject_ids': ids if kind == 'BioProject' else [], 'unresolved': []})
            except Exception as error:
                record['unresolved'].append(f'{type(error).__name__}: {error}')
                entry['error'] = record['unresolved'][-1]
            new_records.append(record)
            series.append(entry)
            print(f'Primary families: {len(seen)} series attempted; {len(queue)} pending', flush=True)
        for pending in sorted(queue - seen):
            new_records.append({'accession': pending, 'source_kind': 'bounded_acquisition_gap',
                'source_id': 'pending-series:' + pending, 'sample_ids': [],
                'unresolved': ['Discovered relative was not acquired within the frozen series allocation; current reverse linkage remains unverified']})
        # A GSM index hit needs primary membership before it can become an edge.
        # Unread/contradictory hits keep the originating record unresolved.
        primary = {r['accession']: r for r in new_records if r['source_kind'] == 'geo_series_primary'}
        for query in queries:
            if query['kind'] != 'GSM':
                continue
            for linked in query.get('linked_studies', []):
                other = primary.get(linked)
                if other is None or not (set(other['sample_ids']) & set(query['identifiers'])):
                    primary[query['source_study']]['unresolved'].append('GSM index hit lacks primary membership confirmation: ' + linked)
        result = build_families(records + new_records, acquirer.selection['exposed_studies'], registry['seed'])
        selected = set(acquirer.selection['seed_studies'])
        selected_families = [f for f in result['families'] if selected & set(f['members'])]
        summary = {'schema': ACQUISITION_VERSION, 'seed_studies': len(selected),
                   'unresolved_atlas_identity_rows': len(acquirer.selection['unresolved_atlas_identity_rows']),
                   'selected_prior_families': len(acquirer.selection['selected_prior_family_ids']),
                   'attempted_primary_series': len(series), 'parsed_primary_series': sum(s['status'] == 'PASS' for s in series),
                   'pending_series': sorted(queue - seen), 'reverse_query_status': dict(Counter(q['status'] for q in queries)),
                   'distinct_reverse_queries': len(reverse_cache),
                   'selected_current_families': len(selected_families),
                   'selected_family_linkage': dict(Counter(f['source_linkage_status'] for f in selected_families)),
                   'selected_family_splits': dict(Counter(f['split'] for f in selected_families)),
                   'network_requests_started': sum(1 for _ in acquirer.root.glob('request_*/attempt-*/request.json')),
                   'evaluation_release': 'CLOSED', 'held_out_studies_evaluated': 0, 'paid_model_calls': 0,
                   'numerical_results_inspected': False,
                   'scope': 'Primary GEO plus indexed GSM, publication and BioProject closure within the frozen request/series limits; unindexed links remain unknown'}
        artifacts = {'families.json.gz': result, 'acquired-source-records.json.gz': new_records,
                     'primary-series.json': series, 'reverse-queries.json': queries,
                     'selected-families.json': selected_families}
        summary['logical_artifact_sha256'] = {name: hashlib.sha256(json.dumps(content, sort_keys=True,
            separators=(',', ':')).encode()).hexdigest() for name, content in artifacts.items()}
        if acquisition_dependencies() != acquirer.selection['dependency_sha256']:
            raise ValueError('Acquisition implementation changed during execution')

        def save(directory):
            for name, content in (('families.json.gz', result), ('acquired-source-records.json.gz', new_records)):
                with (directory / name).open('xb') as raw:
                    with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as handle:
                        handle.write(json.dumps(content, sort_keys=True, separators=(',', ':')).encode())
            write_new(directory / 'primary-series.json', series)
            write_new(directory / 'reverse-queries.json', queries)
            write_new(directory / 'selected-families.json', selected_families)
            write_new(directory / 'summary.json', summary)
            return summary
        # Compare freshly reconstructed values even on a completed resume.
        complete = acquirer.root / 'zz_manifest' / 'complete.json'
        if complete.exists() and json.loads(complete.read_text())['result'] != summary:
            raise ValueError('Acquisition summary changed on resume')
        return acquirer.store.run('zz_manifest', save)


def parse_family_source(path, kind, payload):
    """Read membership/links only; do not extract numerical study outcomes."""
    record = {'accession': None, 'sample_ids': [], 'pubmed_ids': [], 'bioproject_ids': [],
              'linked_accessions': [], 'source_kind': kind, 'unresolved': []}
    if kind == 'sample_metadata':
        if not re.fullmatch(r'GSE\d+', path.stem):
            raise ValueError('Metadata filename has no unambiguous study accession')
        record['accession'] = path.stem
        csv.field_size_limit(64 * 1024 * 1024)
        reader = csv.reader(io.StringIO(payload.decode('utf-8-sig')), strict=True)
        fields = next(reader)
        keys = [i for i, field in enumerate(fields) if field in {'GSM', 'gsm_id', 'geo_accession'}]
        if not keys:
            raise ValueError('No source GSM identity field')
        sample_ids = set(); duplicate_rows = 0
        for row in reader:
            if len(row) != len(fields):
                raise ValueError('Source metadata row width differs from header')
            identities = {row[i].strip() for i in keys}
            if len(identities) != 1 or not re.fullmatch(r'GSM\d+', next(iter(identities))):
                raise ValueError('Conflicting or invalid source GSM identity')
            if 'gse_id' in fields and row[fields.index('gse_id')] != path.stem:
                raise ValueError('Metadata row and filename study identity differ')
            sid = next(iter(identities))
            duplicate_rows += sid in sample_ids
            sample_ids.add(sid)
        if not sample_ids:
            raise ValueError('Source metadata contains no sample memberships')
        record.update(sample_ids=sorted(sample_ids), repeated_membership_rows=duplicate_rows,
                      study_identity_basis='Local accession filename and explicit GSM columns')
        return record
    if kind == 'geo_series_primary':
        identities = re.findall(r'^\^SERIES = (GSE\d+)\s*$', payload.decode('utf-8'), flags=re.M)
        if len(identities) != 1:
            raise ValueError('Primary series source identity is missing/ambiguous')
        data = series_record(identities[0], payload)
    elif kind == 'geo_series_receipt':
        data = json.loads(payload)
        if not re.fullmatch(r'GSE\d+', data.get('accession', '')):
            raise ValueError('Archived series receipt has no valid study identity')
        record['unresolved'].append('Archived parsed receipt; raw primary series payload is not in this source record')
        if data.get('error'):
            record['unresolved'].append('Archived acquisition failed: ' + str(data['error']))
    else:
        raise ValueError('Unrecognized family source adapter')
    record.update(accession=data['accession'], sample_ids=sorted(set(data.get('sample_ids', []))),
                  pubmed_ids=sorted(set(data.get('pubmed_ids', []))))
    if not record['sample_ids']:
        record['unresolved'].append('Source series sample membership is unavailable')
    relations = data.get('relations', [])
    record['linked_accessions'] = sorted({accession for text in relations
        for accession in re.findall(r'(?<![A-Za-z0-9])GSE\d+(?![A-Za-z0-9])', text)})
    record['bioproject_ids'] = sorted({accession for text in relations
        for accession in re.findall(r'(?<![A-Za-z0-9])PRJ(?:NA|EB|DB)\d+(?![A-Za-z0-9])', text)})
    record['relation_quotes'] = relations
    return record


def collect_family_sources(metadata_roots, series_roots, receipt_roots, cache):
    """Hash each source and resume unchanged parsed records from a local cache."""
    cache = Path(cache).resolve(); cache.mkdir(parents=True, exist_ok=True)
    implementation = {str(path): sha(path) for path in (Path(__file__).resolve(),
        Path(__file__).with_name('evaluation_2026_09_10.py').resolve())}
    config = {'schema': FAMILY_SOURCE_VERSION, 'implementation_sha256': implementation,
              'metadata_roots': [str(Path(p).resolve()) for p in metadata_roots],
              'series_roots': [str(Path(p).resolve()) for p in series_roots],
              'receipt_roots': [str(Path(p).resolve()) for p in receipt_roots]}
    config_path = cache / 'config.json'
    if config_path.exists():
        if json.loads(config_path.read_text()) != config:
            raise ValueError('Family cache implementation/configuration changed; use a new cache')
    else:
        write_new(config_path, config)
    inputs = {}
    for roots, kind, pattern in ((metadata_roots, 'sample_metadata', 'GSE*.csv'),
                                 (series_roots, 'geo_series_primary', '*.soft'),
                                 (receipt_roots, 'geo_series_receipt', 'GSE*.json')):
        for root in roots:
            root = Path(root).resolve()
            if not root.exists():
                raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(root))
            paths = sorted(root.rglob(pattern)) if root.is_dir() else [root]
            if not paths:
                raise ValueError('No family sources matched: ' + str(root))
            for path in paths:
                path = path.resolve()
                if path in inputs and inputs[path] != kind:
                    raise ValueError('Source has conflicting adapter declarations')
                inputs[path] = kind
    if not inputs:
        raise ValueError('Family map requires source records')
    records, sources = [], []
    parsed = resumed = 0
    for index, (path, kind) in enumerate(sorted(inputs.items()), 1):
        key = hashlib.sha256((kind + '|' + str(path)).encode()).hexdigest()
        destination = cache / key[:2] / (key + '.json')
        digest = sha(path)
        if destination.is_file():
            envelope = json.loads(destination.read_text())
            encoded = json.dumps(envelope['payload'], sort_keys=True, separators=(',', ':')).encode()
            if hashlib.sha256(encoded).hexdigest() != envelope['payload_sha256']:
                raise ValueError('Cached family source changed: ' + str(destination))
            item = envelope['payload']
            if item['sha256'] != digest or item['path'] != str(path):
                raise ValueError('Family source changed during resume: ' + str(path))
            resumed += 1
        else:
            before = path.stat()
            item = {'source_id': 'source:' + key, 'path': str(path), 'sha256': digest,
                    'bytes': before.st_size, 'mtime_ns': before.st_mtime_ns, 'kind': kind,
                    'status': 'PASS', 'record': None}
            try:
                payload = path.read_bytes()
                if hashlib.sha256(payload).hexdigest() != digest:
                    raise ValueError('Family source changed while reading')
                item['record'] = parse_family_source(path, kind, payload)
                item['record']['source_id'] = item['source_id']
            except (OSError, ValueError, TypeError, KeyError, UnicodeError, csv.Error, StopIteration) as error:
                item.update(status='CANNOT VERIFY', error=f'{type(error).__name__}: {error}')
                accession = re.findall(r'GSE\d+', path.stem)
                if len(set(accession)) == 1:
                    item['record'] = {'accession': accession[0], 'source_id': item['source_id'],
                                      'sample_ids': [], 'source_kind': kind,
                                      'unresolved': ['Source parsing failed: ' + item['error']]}
            if (path.stat().st_size, path.stat().st_mtime_ns) != (before.st_size, before.st_mtime_ns):
                raise ValueError('Family source changed during extraction: ' + str(path))
            destination.parent.mkdir(parents=True, exist_ok=True)
            encoded = json.dumps(item, sort_keys=True, separators=(',', ':')).encode()
            envelope = {'payload': item, 'payload_sha256': hashlib.sha256(encoded).hexdigest()}
            # Commit the payload and its checksum together. An interrupted
            # temporary file never appears as a completed cache entry.
            with tempfile.NamedTemporaryFile(mode='w', dir=destination.parent, suffix='.partial', delete=False) as handle:
                json.dump(envelope, handle, sort_keys=True, separators=(',', ':'))
                pending = Path(handle.name)
            os.link(pending, destination)  # Exclusive: never replace a competing completed entry.
            pending.unlink()
            parsed += 1
        sources.append({k: value for k, value in item.items() if k != 'record'})
        if item['record'] is not None:
            records.append(item['record'])
        if index % 1000 == 0:
            print(f'Family sources: {index}/{len(inputs)}; parsed {parsed}, resumed {resumed}', file=sys.stderr, flush=True)
    if any(sha(path) != digest for path, digest in implementation.items()):
        raise ValueError('Family implementation changed during source collection')
    return records, sources, {'parsed': parsed, 'resumed': resumed, 'implementation_sha256': implementation}


def build_family_manifest(metadata_roots, series_roots, receipt_roots, exposure_paths, cache, output):
    from .evaluation_2026_09_10 import build_families
    output = Path(output).resolve(); output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    records, sources, execution = collect_family_sources(metadata_roots, series_roots, receipt_roots, cache)
    exposure_sources, exposure = [], {}
    for path in sorted({Path(path).resolve() for path in exposure_paths}):
        payload = path.read_bytes()
        # These files are explicitly nominated development records. Extracting
        # accession mentions is conservative and never confers test eligibility.
        ids = sorted(set(re.findall(r'(?<![A-Za-z0-9])GSE\d+(?![A-Za-z0-9])', payload.decode('utf-8'))))
        if not ids:
            raise ValueError('Nominated exposure source has no accession: ' + str(path))
        identity = {'path': str(path), 'sha256': hashlib.sha256(payload).hexdigest(), 'accessions': ids}
        exposure_sources.append(identity)
        for accession in ids:
            exposure.setdefault(accession, []).append(str(path))
    if not exposure:
        raise ValueError('Explicit development-exposure records are required')
    result = build_families(records, exposure)
    summary = {'created_at_utc': datetime.now(timezone.utc).isoformat(),
               'source_files': len(sources), 'source_parse_failures': sum(s['status'] != 'PASS' for s in sources),
               'studies': len(result['studies']), 'families': len(result['families']), 'edge_witnesses': len(result['edges']),
               'multi_study_families': sum(len(f['members']) > 1 for f in result['families']),
               'directly_exposed_studies': len(exposure),
               'family_exposed_studies': sum(s['family_exposed'] for s in result['studies']),
               'split_counts': dict(Counter(s['split'] for s in result['studies'])),
               'source_linkage_status': dict(Counter(s['source_linkage_status'] for s in result['studies'])),
               'evaluation_release': 'CLOSED', 'held_out_studies_evaluated': 0,
               'scope': result['scope'], 'api_calls': 0, 'clock_execution_performed': False,
               'elapsed_seconds': time.monotonic()-started, 'cache_execution': execution,
               'source_coverage_complete': False}
    for name, value in [('families.json.gz', result), ('source-records.json.gz', records),
                        ('source-catalog.json.gz', sources), ('development-exposure.json.gz', exposure_sources)]:
        # Fixed gzip header time makes equivalent payload bytes reproducible.
        with (output/name).open('xb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as compressed:
            compressed.write((json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode())
    with (output/'study-enrollment.csv').open('x', newline='') as handle:
        columns = ['study_id','family_id','split','directly_exposed','family_exposed','known_sample_count','source_linkage_status']
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator='\n'); writer.writeheader()
        writer.writerows({key: row[key] for key in columns} for row in result['studies'])
    summary['output_sha256'] = {path.name: sha(path) for path in output.iterdir() if path.is_file()}
    write_new(output/'manifest.json', summary)
    return summary


def series_record(accession, payload):
    fields={}
    for line in payload.decode('utf-8').splitlines():
        key,sep,value=line.partition(' = ')
        if sep and key in ('^SERIES','!Series_status','!Series_sample_id','!Series_relation',
                           '!Series_pubmed_id','!Series_supplementary_file','!Series_platform_id'):
            fields.setdefault(key,[]).append(value.strip())
    if fields.get('^SERIES')!=[accession] or len(fields.get('!Series_status',[]))!=1:
        raise ValueError('GEO series identity or public-status evidence is missing/ambiguous')
    status=fields['!Series_status'][0]
    if not status.startswith('Public on '):
        raise ValueError('No public-release status')
    released=datetime.strptime(status.removeprefix('Public on '),'%b %d %Y').date().isoformat()
    samples=fields.get('!Series_sample_id',[])
    if not samples or len(set(samples))!=len(samples) or any(not re.fullmatch(r'GSM\d+',s) for s in samples):
        raise ValueError('Series sample list is absent or invalid')
    return {'accession':accession,'series_first_public_date':released,'source_status_quote':status,
            'sample_ids':samples,'relations':fields.get('!Series_relation',[]),
            'pubmed_ids':fields.get('!Series_pubmed_id',[]),
            'supplementary_files':fields.get('!Series_supplementary_file',[]),
            'platform_ids':fields.get('!Series_platform_id',[]),
            'response_sha256':hashlib.sha256(payload).hexdigest()}


def read_previous(overviews, metadata):
    studies=set(); projects=set(); samples=set(); failures=[]
    sources={}
    for path in sorted(Path(overviews).glob('*.csv')):
        sources[str(path.resolve())]=sha(path)
        with path.open() as handle:
            for row in csv.DictReader(handle):
                studies.update(re.findall(r'GSE\d+',row.get('accession','')))
                projects.update(re.findall(r'PRJ[A-Z]+\d+',row.get('project','')))
    # This screen is scoped to available old DNAm metadata; other modalities
    # and absent old metadata still need the family map before admission.
    for path in sorted(Path(metadata).glob('GSE*.csv')):
        studies.add(path.stem)
        try:
            sources[str(path.resolve())]=sha(path)
            with path.open() as handle:
                reader=csv.DictReader(handle)
                keys=[k for k in ('GSM','gsm_id','geo_accession') if k in (reader.fieldnames or [])]
                if not keys:
                    raise ValueError('No sample identity field')
                for row in reader:
                    values=[row[k] for k in keys]
                    if len(set(values))!=1 or not re.fullmatch(r'GSM\d+',values[0] or ''):
                        raise ValueError('Conflicting or invalid sample identity')
                    samples.add(values[0])
        except Exception as error:
            failures.append({'path':str(path),'error':f'{type(error).__name__}: {error}'})
    if not studies or not samples:
        raise ValueError('Old accession/sample screen has no usable inputs')
    return studies,projects,samples,failures,sources


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates',type=Path)
    parser.add_argument('--old-overviews',type=Path)
    parser.add_argument('--old-metadata',type=Path)
    parser.add_argument('--build-family-map',action='store_true')
    parser.add_argument('--metadata-root',type=Path,action='append',default=[])
    parser.add_argument('--series-root',type=Path,action='append',default=[])
    parser.add_argument('--series-receipts',type=Path,action='append',default=[])
    parser.add_argument('--exposure',type=Path,action='append',default=[])
    parser.add_argument('--family-cache',type=Path)
    parser.add_argument('--prepare-family-acquisition',type=Path)
    parser.add_argument('--acquire-family-sources',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if sum(bool(value) for value in (args.prepare_family_acquisition,args.acquire_family_sources,args.build_family_map,args.candidates)) != 1:
        raise ValueError('Choose exactly one family preparation/acquisition, family-map or release-screen operation')
    if args.prepare_family_acquisition:
        print(json.dumps(prepare_family_acquisition(args.prepare_family_acquisition,args.output))); return
    if args.acquire_family_sources:
        print(json.dumps(acquire_family_sources(args.acquire_family_sources,args.output))); return
    if args.build_family_map:
        if args.family_cache is None:raise ValueError('--family-cache is required')
        result=build_family_manifest(args.metadata_root,args.series_root,args.series_receipts,args.exposure,
                                     args.family_cache,args.output)
        print(json.dumps(result)); return
    if any(value is None for value in (args.candidates,args.old_overviews,args.old_metadata)):
        raise ValueError('Release screen requires --candidates, --old-overviews and --old-metadata')
    config={'implementation_sha256':sha(__file__),'candidates_sha256':sha(args.candidates),
            'old_overviews':str(args.old_overviews.resolve()),'old_metadata':str(args.old_metadata.resolve())}
    args.output.mkdir(parents=True,exist_ok=True)
    lock=args.output/'config.json'
    if lock.exists() and json.loads(lock.read_text())!=config:
        raise ValueError('Screen configuration changed; use a new output directory')
    if not lock.exists():write_new(lock,config)
    previous,projects,samples,failures,sources=read_previous(args.old_overviews,args.old_metadata)
    provenance=args.output/'old-inputs.json'
    if provenance.exists() and json.loads(provenance.read_text())!=sources:
        raise ValueError('Old corpus inputs changed during resume')
    if not provenance.exists():write_new(provenance,sources)
    candidates=json.loads(args.candidates.read_text())['studies']
    records=[]
    for candidate in candidates:
        accession=candidate['accession']
        if not re.fullmatch(r'GSE\d+',accession):raise ValueError('Invalid accession')
        receipt=args.output/(accession+'.json')
        if receipt.exists():
            record=json.loads(receipt.read_text())
        else:
            url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&form=text&view=full'
            record={'accession':accession,'url':url,'retrieved_at_utc':datetime.now(timezone.utc).isoformat(),
                    'eligibility':'CANNOT VERIFY','historical_corpus_boundary_verified':False}
            try:
                time.sleep(.35) # Sequential requests stay below NCBI's unkeyed rate ceiling.
                with urllib.request.urlopen(url,timeout=30) as response:
                    body=response.read(2000001)
                if len(body)>2000000:raise ValueError('Series response exceeds bounded size')
                record.update(series_record(accession,body))
                record['esummary_public_date']=candidate['pdat'].replace('/','-')
                record['release_date_agreement']=record['series_first_public_date']==record['esummary_public_date']
                record['shared_old_sample_ids']=sorted(set(record['sample_ids'])&samples)
                record['accession_seen_before']=accession in previous
                linked_gse={g for line in record['relations'] for g in re.findall(r'GSE\d+',line)}
                linked_projects={p for line in record['relations'] for p in re.findall(r'PRJ[A-Z]+\d+',line)}
                record['linked_old_accessions']=sorted(linked_gse&previous)
                record['linked_old_projects']=sorted(linked_projects&projects)
                record['overlap_detected']=bool(record['shared_old_sample_ids'] or record['accession_seen_before'] or record['linked_old_accessions'] or record['linked_old_projects'])
                record['missing_for_admission']=['Complete old-corpus sample/family/publication map',
                    'Molecular input and per-clock measurement applicability','Independent audit and frozen cohort manifest']
                if record['overlap_detected']:record['eligibility']='FAIL'
            except Exception as error:
                record['error']=f'{type(error).__name__}: {error}'
            write_new(receipt,record)
        records.append(record)
    summary={'candidate_studies':len(records),'date_agreement_studies':sum(r.get('release_date_agreement',False) for r in records),
             'by_verified_release_year':dict(Counter(r['series_first_public_date'][:4] for r in records if r.get('release_date_agreement'))),
             'overlap_detected_studies':sum(r.get('overlap_detected',False) for r in records),
             'retrieval_failures':sum('error' in r for r in records),
             'old_sample_ids_inspected':len(samples),'old_metadata_failures':failures,
             'eligible_studies':None,'audited_completed_studies':0,'biological_results_computed':False}
    if (args.output/'summary.json').exists():
        if json.loads((args.output/'summary.json').read_text())!=summary:raise ValueError('Completed summary changed')
    else:write_new(args.output/'summary.json',summary)
    print(json.dumps(summary))


if __name__=='__main__':
    try:main()
    except Exception as error:raise SystemExit(f'Extension eligibility stopped: {type(error).__name__}: {error}')
