"""Source acquisition and adapters for the existing ClockBase clock producers."""
import csv
from collections import Counter
from contextlib import contextmanager
import gzip
import hashlib
import importlib.util
import io
import json
import os
import posixpath
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

import h5py
import numpy as np
import pandas as pd

from audit_system.agent_workflow_2026_09_10 import sha, write_new
from audit_system.extension_eligibility_2026_09_10 import series_record
from audit_system.source_packet import CONTEXT_VERSION, context_documents, article_methods, series_design
from audit_system.source_packet import read_source_csv, SAMPLE_METHOD_FIELDS
from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError

ASSETS = Path(__file__).resolve().parents[1] / 'pipeline_assets/rna_4_6_2026_09_10'
GEO_METADATA_VIEWS = frozenset({'full', 'brief'})
DNAM_REFERENCE_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz'
# Source-verified original producer reference; changing it requires a new implementation.
DNAM_REFERENCE_SHA256 = 'aa10085d9e53a17a736398078de7af9287e10946915f1ad1b5637e172093f0d6'
DNAM_REFERENCE_BYTES = 1237679690
# Frozen from the historical source, with only complete !Series_relation
# metadata records excluded. The source-derived preparation is recorded under
# update-controller/source-bootstrap-empirical-001/reference-profile-preparation-002.json.
DNAM_REFERENCE_PROFILE = {
    'schema': 'clockbase-gse40279-source-table-profile-v2',
    'ordered_canonical_metadata_record_count': 67,
    'ordered_canonical_metadata_sha256': 'c39676efa3010c56baeb6bdb22b2c03b6dbf3e284a0492bcd4498491b910e58b',
    'table_bytes_including_header_excluding_end_marker': 3134661444,
    'table_sha256': 'f5accdff861b0ef666f95925540a29543ca0530b83b38175e0d0569437447c20',
    'feature_rows': 473034,
    'ordered_header_sha256': '030f4560fb9af912101b5c6707af398a26aa7f87d35bf26029911b7ce28d555d',
    'sample_count': 656,
}
DNAM_REFERENCE_MAX_METADATA_BYTES = 1_000_000
DNAM_REFERENCE_MAX_LINE_BYTES = 1_000_000
DNAM_REFERENCE_MAX_DECOMPRESSED_BYTES = 3_200_000_000


def checked_assets():
    declaration = json.loads((ASSETS / 'manifest.json').read_text())
    for name, record in declaration['files'].items():
        path = ASSETS / name
        if not path.is_file():
            raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
        if path.read_bytes()[:80].startswith(b'version https://git-lfs.github.com/spec/'):
            raise ValueError('Missing model payload; run git lfs pull --include="pipeline_assets/**"')
        if sha(path) != record['sha256']:
            raise PipelineIntegrityError('Pinned pipeline asset changed: ' + name)
    return declaration


def fetch(url, destination, limit=20_000_000):
    if not url.startswith('https://'):
        raise ValueError('Public source acquisition requires HTTPS')
    part = destination.with_name(destination.name + '.partial')
    if destination.exists() or part.exists() or destination.with_name(destination.name + '.source.json').exists():
        raise FileExistsError('Source destination already exists: ' + str(destination))
    result = subprocess.run(['curl', '--fail', '--silent', '--show-error', '--location', '--proto', '=https',
                             '--proto-redir', '=https', '--max-redirs', '3', '--max-time', '180',
                             '--max-filesize', str(limit), '--output', str(part), url], capture_output=True, timeout=190)
    if result.returncode:
        raise RuntimeError('Public source download failed (curl %d): %s' % (result.returncode, url))
    if part.stat().st_size == 0 or part.stat().st_size > limit:
        raise ValueError('Source response is empty or exceeds download limit')
    part.rename(destination)
    write_new(destination.with_name(destination.name + '.source.json'),
              {'url': url, 'retrieved_unix': time.time(), 'sha256': sha(destination),
               'size_bytes': destination.stat().st_size})
    return destination


def parse_samples(payload, study):
    samples = []
    current = None
    allowed = {'title', 'source_name_ch1', 'organism_ch1', 'molecule_ch1', 'platform_id',
               'characteristics_ch1', 'description', 'relation'} | SAMPLE_METHOD_FIELDS
    for line in payload.splitlines():
        key, sep, value = line.partition(' = ')
        if key == '^SAMPLE' and sep:
            if current is not None: samples.append(current)
            current = {'GSM': value.strip(), 'gse_id': study}
        elif current is not None and sep and key.startswith('!Sample_'):
            field = key.removeprefix('!Sample_')
            if field in allowed:
                if field in current:
                    field += '.' + str(sum(k.startswith(field) for k in current))
                current[field] = value.strip()
    if current is not None: samples.append(current)
    if not samples or len({s['GSM'] for s in samples}) != len(samples):
        raise ValueError('GEO source sample identities are absent or duplicated')
    return samples


def acquire_study_context(study, series, output, *, publication_methods=True):
    """Bounded public-source lookup; unavailable Methods remain explicit gaps.

    Up to two Europe PMC searches and three full-text downloads. Papers are
    selected by repository links/accession membership, never their outcomes.
    The response snapshots are retained; only extracted Methods enter packets.
    """
    if not re.fullmatch(r'GSE\d+', study):
        raise ValueError('Study context requires a GEO accession')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    local_series = output / 'series.soft'
    if Path(series).resolve() != local_series.resolve():
        if local_series.exists():
            raise FileExistsError('Source destination already exists: ' + str(local_series))
        shutil.copyfile(series, local_series)
    fields = series_design(local_series.read_bytes(), study)
    pubmed_ids = fields.get('pubmed_id', '').splitlines()
    if any(not re.fullmatch(r'\d+', pid) for pid in pubmed_ids):
        raise ValueError('GEO PubMed identifier is malformed')
    sources = [{'kind': 'geo_series', 'id': study, 'path': local_series.name,
                'sha256': sha(local_series),
                'url': f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={study}&targ=self&form=text&view=full'}]
    acquisition = []
    publications = {}

    def search(query, name):
        url = ('https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&pageSize=100&query='
               + urllib.parse.quote(query, safe=''))
        try:
            time.sleep(.35)
            path = fetch(url, output / (name + '.json'), limit=2_000_000)
            response = json.loads(path.read_text())
            rows = response.get('resultList', {}).get('result', [])
            count = response['hitCount']
            if type(count) is not int or not isinstance(rows, list):
                raise ValueError('Publication search response schema mismatch')
            acquisition.append({'step': name, 'status': 'RETRIEVED', 'query': query,
                                'hit_count': count, 'returned_count': len(rows),
                                'complete_search_response': count == len(rows)})
            for row in rows:
                pmcid = row.get('pmcid', '')
                if re.fullmatch(r'PMC\d+', pmcid or ''):
                    publications[pmcid] = row.get('id')
            if rows and not publications:
                acquisition.append({'step': name, 'status': 'CANNOT VERIFY',
                                    'reason': 'Search records have no PMC full-text identifier'})
        except (OSError, ValueError, RuntimeError, KeyError, subprocess.SubprocessError) as error:
            acquisition.append({'step': name, 'status': 'CANNOT VERIFY',
                                'reason': f'{type(error).__name__}: {error}'})

    if publication_methods:
        if pubmed_ids:
            search('(' + ' OR '.join('EXT_ID:' + pid for pid in sorted(set(pubmed_ids))) + ') AND SRC:MED',
                   'geo-citation-search')
        # A missing GEO citation is not evidence that no associated paper exists.
        if not publications:
            search('BODY:' + study, 'accession-body-search')
        selected = sorted(publications, key=lambda pid: int(pid[3:]))[:3]
        if len(publications) > len(selected):
            acquisition.append({'step': 'publication-limit', 'status': 'CANNOT VERIFY',
                                'reason': 'More linked full texts than the declared three-paper limit',
                                'omitted_ids': sorted(set(publications) - set(selected))})
        for pmcid in selected:
            url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML'
            try:
                time.sleep(.35)
                path = fetch(url, output / (pmcid + '.xml'), limit=10_000_000)
                documents = article_methods(path.read_bytes(), study, pmcid, pubmed_ids)
                sources.append({'kind': 'publication_jats', 'id': pmcid, 'path': path.name,
                                'url': url, 'sha256': sha(path)})
                acquisition.append({'step': 'publication-methods', 'source_id': pmcid,
                                    'status': 'EXTRACTED' if documents else 'CANNOT VERIFY',
                                    'method_sections': len(documents),
                                    'reason': None if documents else 'No separable Methods section in JATS'})
            except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
                acquisition.append({'step': 'publication-methods', 'source_id': pmcid,
                                    'status': 'CANNOT VERIFY', 'reason': f'{type(error).__name__}: {error}'})
        if not selected:
            acquisition.append({'step': 'publication-methods', 'status': 'CANNOT VERIFY',
                                'reason': 'No accessible linked PMC full text identified'})
    else:
        acquisition.append({'step': 'publication-methods', 'status': 'NOT_REQUESTED'})
    context = {'version': CONTEXT_VERSION, 'study_id': study, 'sources': sources,
               'documents': context_documents(study, sources, output), 'acquisition': acquisition,
               'limits': {'searches': 2, 'publication_full_texts': 3,
                          'search_response_bytes': 2_000_000, 'publication_bytes': 10_000_000},
               'source_acquired_unix': time.time(),
               'interpretation': 'Source association and quotation checks; cohort correspondence and design validity need assessment'}
    path = output / 'study-context.json'
    write_new(path, context)
    return path


def acquire_metadata(study, output, *, publication_methods=True, metadata_view='full'):
    if metadata_view not in GEO_METADATA_VIEWS:
        raise ValueError('GEO metadata view must be one of: brief, full')
    series = fetch(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={study}&targ=self&form=text&view=full', output / 'series.soft')
    record = series_record(study, series.read_bytes())
    time.sleep(.35)
    source = fetch(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={study}&targ=gsm&form=text&view={metadata_view}', output / 'samples.soft')
    samples = parse_samples(source.read_text(), study)
    if {s['GSM'] for s in samples} != set(record['sample_ids']):
        raise ValueError('GEO series and sample metadata membership differ')
    path = output / 'metadata.csv'
    fields = sorted({k for sample in samples for k in sample})
    with path.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(samples)
    context = acquire_study_context(study, series, output, publication_methods=publication_methods)
    record.update(metadata=str(path), source_sample_records=str(source), series_source=str(series),
                  study_context=str(context), sample_metadata_view=metadata_view)
    write_new(output / 'source-membership.json', record)
    return record


def select_archs4_release(catalog, species):
    """Select an official release; catalog timestamps are not study-public dates."""
    from email.utils import parsedate_to_datetime
    if species not in {'human', 'mouse'}: raise ValueError('Unsupported ARCHS4 species')
    if isinstance(catalog, dict): catalog = catalog.get('versionfiles', catalog.get('data', catalog.get('files', [])))
    if not isinstance(catalog, list) or any(not isinstance(row, dict) for row in catalog):
        raise ValueError('Malformed official ARCHS4 catalog')
    choices = [c for c in catalog if c.get('species') == species and c.get('data_level') == 'gene'
               and str(c.get('version_major')) == '2' and str(c.get('version_minor')) == 'latest']
    if not choices: raise ValueError('Official ARCHS4 release is missing')
    choices.sort(key=lambda c: parsedate_to_datetime(c['timestamp']), reverse=True)
    if len(choices) > 1 and parsedate_to_datetime(choices[0]['timestamp']) == parsedate_to_datetime(choices[1]['timestamp']):
        raise ValueError('Official ARCHS4 release selection is ambiguous')
    declaration = dict(choices[0])
    declaration['url'] = f'https://s3.k8s.maayanlab.cloud/archs4/files/{species}_gene_v2.latest.h5'
    return declaration


def acquire_archs4(species, cache, *, request=None, selection_dir=None, max_source_bytes=None, metadata_bytes=2000000):
    """Select and pin an official release before downloading; cache locally."""
    from audit_system.archs4_update_2026_09_10 import download, remaining_storage
    if request is not None:
        # Selection is per explicit refresh revision; immutable objects may be shared.
        import hashlib
        selection_dir = Path(selection_dir); selection_dir.mkdir(parents=True, exist_ok=True)
        catalog_path = selection_dir / 'catalog.json'
        request('GET', 'https://archs4.org/api/versionfile', catalog_path, metadata_bytes, {})
        declaration = select_archs4_release(json.loads(catalog_path.read_text()), species)
        catalog_id, catalog_checksum = declaration.get('id'), declaration.get('checksum')
        if (not ((type(catalog_id) is int and catalog_id > 0)
                 or (isinstance(catalog_id, str) and catalog_id.strip()))
                or not isinstance(catalog_checksum, str) or not catalog_checksum.strip()):
            raise ValueError('Official ARCHS4 selected record lacks usable id/checksum')
        declaration = {key: declaration[key] for key in ('id', 'checksum', 'species', 'data_level',
                       'version_major', 'version_minor', 'timestamp', 'file_size', 'samples', 'url')}
        declaration['catalog_payload_checksum_verification'] = {
            'status': 'CANNOT VERIFY',
            'reason': 'The official catalog checksum algorithm and hashed representation are not established'}
        if (type(declaration.get('file_size')) is not int or declaration['file_size'] <= 0
                or declaration['file_size'] > max_source_bytes):
            raise ValueError('Selected ARCHS4 release exceeds the explicit source-byte limit')
        head = request('HEAD', declaration['url'], None, 0, {})
        etag = head['headers'].get('etag')
        if (not etag or etag.startswith('W/') or head['headers'].get('content-length') != str(declaration['file_size'])):
            raise ValueError('Official ARCHS4 HEAD differs from catalog or lacks a strong ETag')
        declaration['http_headers'] = {'Etag': etag}
        write_new(selection_dir / (species+'-selection.json'), declaration)
        object_id = hashlib.sha256(json.dumps(declaration, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        destination = cache / 'objects' / object_id
        result = download(species, declaration, destination, request=request)
        return Path(result['path'])
    cache.mkdir(parents=True, exist_ok=True)
    catalog_path = cache / 'catalog.json'
    if not catalog_path.exists(): fetch('https://archs4.org/api/versionfile', catalog_path)
    selection = cache / (species + '-selection.json')
    if selection.exists():
        declaration = json.loads(selection.read_text())
    else:
        declaration = select_archs4_release(json.loads(catalog_path.read_text()), species)
        with urllib.request.urlopen(urllib.request.Request(declaration['url'], method='HEAD'), timeout=30) as response:
            declaration['http_headers'] = {'Etag': response.headers['Etag']}
        write_new(selection, declaration)
    if shutil.disk_usage(cache).free < remaining_storage({species: declaration}, cache) + 5 * 1024**3:
        raise RuntimeError('Insufficient disk for ARCHS4 species file plus 5 GiB reserve')
    return Path(download(species, declaration, cache)['path'])


def extract_rna(hdf5, study, species, source_samples, output):
    decode = lambda values: [v.decode() if isinstance(v, bytes) else str(v) for v in values]
    before = hdf5.stat()
    with h5py.File(hdf5, 'r') as h5:
        ids = decode(h5['meta/samples/geo_accession'][:])
        memberships = decode(h5['meta/samples/series_id'][:])
        if len(ids) != len(set(ids)) or len(ids) != len(memberships):
            raise ValueError('ARCHS4 sample identities are not unique/aligned')
        indices = [i for i, line in enumerate(memberships) if study in {g.strip() for g in line.split(',')}]
        if not indices:
            raise ValueError('Study has no ARCHS4 samples: ' + study)
        selected = [ids[i] for i in indices]
        if set(selected) != set(source_samples):
            raise ValueError('ARCHS4/GEO sample membership differs; no silent partial-study extraction')
        organism = decode(h5['meta/samples/organism_ch1'][indices])
        expected = {'mouse': 'Mus musculus', 'human': 'Homo sapiens'}[species]
        if set(organism) != {expected}: raise ValueError('ARCHS4 organism differs from requested modality')
        genes = decode(h5['meta/genes/ensembl_gene'][:])
        frame = pd.DataFrame(h5['data/expression'][:, indices], index=genes, columns=selected)
        member_records = {ids[i]: memberships[i].split(',') for i in indices}
    if (before.st_size, before.st_mtime_ns) != (hdf5.stat().st_size, hdf5.stat().st_mtime_ns):
        raise PipelineIntegrityError('ARCHS4 source changed during extraction')
    require_matrix(frame, beta=False)
    path = output / 'counts.csv'
    frame.rename_axis('ensembl_gene').to_csv(path, mode='x')
    receipt = {'path': str(path), 'sha256': sha(path), 'feature_count': len(frame), 'sample_count': len(frame.columns),
               'source_hdf5': str(hdf5), 'source_size': before.st_size, 'source_mtime_ns': before.st_mtime_ns,
               'all_study_memberships': member_records,
               'source_file_sha256': None, 'hash_scope': 'Full extracted matrix hashed; HDF5 identity also checked by run configuration'}
    write_new(output / 'extraction.json', receipt)
    return receipt


def require_matrix(frame, *, beta):
    if frame.empty or frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError('AUDIT_FAIL: empty molecular matrix or duplicate identity')
    if any(pd.isna(v) or not str(v).strip() for v in [*frame.index, *frame.columns]):
        raise ValueError('AUDIT_FAIL: blank molecular identity')
    values = frame.to_numpy(dtype=float)
    finite = np.isfinite(values)
    if np.isinf(values).any() or (finite.sum(axis=0) == 0).any():
        raise ValueError('AUDIT_FAIL: sample has no finite observed measurements')
    if np.any(values[finite] < 0) or (beta and np.any(values[finite] > 1)):
        raise ValueError('AUDIT_FAIL: measurement range conflicts with input contract')
    if not beta and (not finite.all() or np.any(values.sum(axis=0) <= 0)):
        raise ValueError('AUDIT_FAIL: missing or empty RNA counts')


MATRIX_READER_CONTRACT = 'geo_stream_round_trip_v2'


def _inspect_matrix_stream(handle, *, geo, delimiter):
    """Validate logical CSV records and retain a seek cookie after the header.

    Physical line counts cannot be passed as pandas skiprows when the GEO
    preamble contains quoted newlines. The same open stream is used for the
    subsequent value read. No sample or feature identifier is inferred.
    """
    reader = csv.reader(iter(handle.readline, ''), delimiter=delimiter, strict=True)
    accessions, sample_rows = [], []
    if geo:
        for row in reader:
            if row == ['!series_matrix_table_begin']:
                break
            if row and row[0] == '!Series_geo_accession': accessions.extend(row[1:])
            if row and row[0] == '!Sample_geo_accession': sample_rows.append(row[1:])
        else:
            raise ValueError('AUDIT_FAIL: GEO matrix begin marker missing')
    header = next(reader, [])
    if (len(header) < 2 or len(set(header)) != len(header)
            or any(not value.strip() for value in header[1:])):
        raise ValueError('AUDIT_FAIL: missing or duplicate raw header sample IDs')
    if geo:
        if header[0] != 'ID_REF' or any(not re.fullmatch(r'GSM[0-9]+', value) for value in header[1:]):
            raise ValueError('AUDIT_FAIL: GEO header requires ID_REF and canonical GSM identifiers')
        if sample_rows and sample_rows != [header[1:]]:
            raise ValueError('AUDIT_FAIL: GEO matrix and metadata sample order or membership differs')
    data_start = handle.tell()
    features, seen, ended = [], set(), False
    for record, row in enumerate(reader, 1):
        if geo and row == ['!series_matrix_table_end']:
            if ended: raise ValueError('AUDIT_FAIL: repeated GEO matrix end marker')
            ended = True
            continue
        if ended:
            raise ValueError('AUDIT_FAIL: content after GEO matrix end marker')
        if len(row) != len(header):
            raise ValueError(f'AUDIT_FAIL: matrix CSV record {record} has {len(row)} fields; expected {len(header)}')
        if not row[0].strip() or row[0] in seen:
            raise ValueError('AUDIT_FAIL: blank or duplicate molecular feature identity')
        if geo and row[0] in {'!series_matrix_table_begin', '!series_matrix_table_end'}:
            raise ValueError('AUDIT_FAIL: malformed GEO matrix marker')
        features.append(row[0]); seen.add(row[0])
    if geo and not ended:
        raise ValueError('AUDIT_FAIL: GEO matrix end marker missing')
    return {'header': header, 'features': features, 'data_start': data_start,
            'geo': geo, 'delimiter': delimiter, 'accessions': accessions,
            'reader_contract': MATRIX_READER_CONTRACT}


def read_matrix(path, *, beta=False):
    opener = gzip.open if str(path).endswith('.gz') else open
    checksum = sha(path)
    with opener(path, 'rt', newline='') as handle:
        first = handle.readline()
        geo = first.startswith(('!Series_', '^SERIES', '!Sample_', '!series_matrix_table_begin'))
        delimiter = '\t' if geo or '\t' in first else ','
        handle.seek(0)
        source = _inspect_matrix_stream(handle, geo=geo, delimiter=delimiter)
        handle.seek(source['data_start'])
        frame = pd.read_csv(handle, sep=delimiter, index_col=0, names=source['header'],
                            header=None, nrows=len(source['features']), converters={0: str},
                            float_precision='round_trip')
    if sha(path) != checksum:
        raise PipelineIntegrityError('Molecular source changed during validation/read')
    if list(frame.columns) != source['header'][1:] or list(frame.index) != source['features']:
        raise ValueError('AUDIT_FAIL: parsed molecular identities differ from literal source records')
    frame.index = frame.index.astype(str)
    require_matrix(frame, beta=beta)
    return frame


@contextmanager
def _xlsx_beta_stream(path):
    """Read stored XLSX cell tokens without evaluating formulas or guessing samples.

    One worksheet is required because this entry point has no sheet selector.
    Coordinates, not worksheet dimensions, define the observed cell layout.
    """
    ns = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
    rel_ns = '{http://schemas.openxmlformats.org/package/2006/relationships}'
    office_ns = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
    with zipfile.ZipFile(path) as archive:
        entries = archive.infolist()
        names = [entry.filename for entry in entries]
        if len(names) != len(set(names)):
            raise ValueError('XLSX contains duplicate archive members')
        if sum(entry.file_size for entry in entries) > 3_200_000_000:
            raise ValueError('XLSX exceeds the 3.2 GB decompressed input limit')
        workbook = ET.fromstring(archive.read('xl/workbook.xml'))
        sheets = workbook.findall(ns+'sheets/'+ns+'sheet')
        if len(sheets) != 1:
            raise ValueError('XLSX beta recovery requires exactly one worksheet; sheet selection is unresolved')
        relationships = ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        selected = [item for item in relationships.findall(rel_ns+'Relationship')
                    if item.get('Id') == sheets[0].get(office_ns+'id')]
        if (len(selected) != 1 or selected[0].get('TargetMode') == 'External'
                or selected[0].get('Type') != office_ns[1:-1]+'/worksheet'):
            raise ValueError('XLSX worksheet relationship is missing or ambiguous')
        target = selected[0].get('Target', '')
        member = posixpath.normpath(target.lstrip('/') if target.startswith('/') else posixpath.join('xl', target))
        if not member.startswith('xl/worksheets/') or member not in names:
            raise ValueError('XLSX worksheet target is outside the declared archive')
        strings = []
        if 'xl/sharedStrings.xml' in names:
            with archive.open('xl/sharedStrings.xml') as source:
                events = ET.iterparse(source, events=('start', 'end'))
                _, root = next(events)
                for event, element in events:
                    if event == 'end' and element.tag == ns+'si':
                        if element.find(ns+'rPh') is not None:
                            raise ValueError('XLSX phonetic string annotations require explicit interpretation')
                        strings.append(''.join(item.text or '' for item in element.iter(ns+'t')))
                        root.clear()

        def rows():
            width = None
            observed_row = 0
            sheet_data = None
            with archive.open(member) as source:
                for event, element in ET.iterparse(source, events=('start', 'end')):
                    if event == 'start' and element.tag == ns+'mergeCell':
                        raise ValueError('XLSX merged cells require an explicit interpretation and cannot define sample measurements')
                    if event == 'start' and element.tag == ns+'sheetData':
                        sheet_data = element
                    if event != 'end' or element.tag != ns+'row':
                        continue
                    observed_row += 1
                    if element.get('r') != str(observed_row):
                        raise ValueError('XLSX beta table contains missing, repeated or unordered rows')
                    cells = {}
                    for cell in element.findall(ns+'c'):
                        match = re.fullmatch(r'([A-Z]+)([1-9][0-9]*)', cell.get('r', ''))
                        if match is None or int(match[2]) != observed_row:
                            raise ValueError('XLSX cell coordinate does not match its row')
                        column = 0
                        for letter in match[1]: column = column * 26 + ord(letter) - 64
                        if column > 16384 or column in cells:
                            raise ValueError('XLSX has an invalid or repeated column coordinate')
                        if cell.find(ns+'f') is not None:
                            raise ValueError('XLSX formula cells are not measured beta inputs')
                        kind = cell.get('t', 'n')
                        value = cell.findtext(ns+'v', default='')
                        if kind == 's':
                            if not re.fullmatch(r'[0-9]+', value) or int(value) >= len(strings):
                                raise ValueError('XLSX shared-string reference is invalid')
                            value = strings[int(value)]
                        elif kind == 'inlineStr':
                            if cell.find('.//'+ns+'rPh') is not None:
                                raise ValueError('XLSX phonetic string annotations require explicit interpretation')
                            value = ''.join(item.text or '' for item in cell.iter(ns+'t'))
                        elif kind != 'n':
                            raise ValueError('XLSX boolean, error or date cell is not a measured beta value')
                        if (observed_row == 1 or column == 1) and kind not in {'s', 'inlineStr'}:
                            raise ValueError('XLSX sample and feature identifiers must be explicit text cells')
                        cells[column] = value
                    if width is None:
                        width = max(cells, default=0)
                        if width < 2 or set(cells) != set(range(1, width+1)):
                            raise ValueError('XLSX beta header has missing sample identifiers')
                    if not cells or max(cells) > width:
                        raise ValueError('XLSX beta row has no cells or extends beyond the source header')
                    yield [cells.get(column, '') for column in range(1, width+1)]
                    element.clear()
                    if sheet_data is not None: sheet_data.clear()
        stream = rows()
        try:
            yield sheets[0].get('name'), stream
        finally:
            stream.close()


@contextmanager
def _beta_supplement_rows(path, layout):
    if layout.get('format') == 'xlsx':
        with _xlsx_beta_stream(path) as (worksheet, rows):
            if worksheet != layout['worksheet']:
                raise ValueError('XLSX beta worksheet changed')
            yield rows
    else:
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt', newline='') as handle:
            delimiter = layout['delimiter']
            yield ((shlex.split(line, comments=False, posix=True) for line in handle)
                   if delimiter == 'space' else csv.reader(handle, delimiter=delimiter))


def beta_table_header(path):
    """Inspect a plain or GEO beta table without loading its molecular values."""
    if str(path).lower().endswith('.xlsx'):
        with _xlsx_beta_stream(path) as (worksheet, rows):
            header = next(rows, [])
            sample_ids = header[1:]
            if not sample_ids or any(not value for value in sample_ids) or len(set(sample_ids)) != len(sample_ids):
                raise ValueError('XLSX beta measurement columns lack distinct sample identifiers; annotation columns are not samples')
            empty = next(rows, None) is None
        return {'header': header, 'delimiter': None, 'skiprows': 0, 'geo': False,
                'feature_column': True, 'empty': empty, 'format': 'xlsx', 'worksheet': worksheet}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        first = handle.readline()
        geo = first.startswith(('!Series_', '^SERIES'))
        skipped = 0
        if geo:
            skipped = 1
            for line in handle:
                skipped += 1
                if line.strip() == '!series_matrix_table_begin':
                    break
            else:
                raise ValueError('GEO matrix begin marker is missing')
            first = handle.readline()
        if '\t' in first:
            delimiter = '\t'
            header = next(csv.reader([first], delimiter=delimiter))
            feature_column = True
        elif ',' in first:
            delimiter = ','
            header = next(csv.reader([first], delimiter=delimiter))
            feature_column = True
        else:
            # Some GEO-linked beta tables omit the feature-column label and
            # quote each sample token in a whitespace-delimited header.  Use
            # shell-style tokenization only for this explicitly detected
            # layout; it preserves coordinate labels verbatim.
            delimiter = 'space'
            header = shlex.split(first, comments=False, posix=True)
            feature_column = False
        sample_ids = header[1:] if feature_column else header
        if len(sample_ids) < 1 or any(not s for s in sample_ids) or len(set(sample_ids)) != len(sample_ids):
            raise ValueError('Missing or duplicate beta-table sample identifiers')
        first_data_raw = handle.readline()
        first_data = first_data_raw.strip()
        empty = first_data == '!series_matrix_table_end' if geo else not first_data
        if empty and any(line.strip() for line in handle):
            raise ValueError('Content follows an apparently empty matrix')
    return {'header': header, 'delimiter': delimiter, 'skiprows': skipped, 'geo': geo,
            'feature_column': feature_column, 'empty': empty}


def recovery_primary_samples(series, study, sample_source=None):
    """Read assay identity evidence from a GEO Matrix header or family SOFT.

    Stop before Matrix values and skip SOFT sample tables. Repeated source
    characteristics remain separate strings, without inferred column aliases.
    """
    def header_lines():
        paths = [series] + ([sample_source] if sample_source is not None else [])
        for path in paths:
            opener = gzip.open if str(path).endswith('.gz') else open
            with opener(path, 'rt', newline='') as handle:
                for line in handle:
                    if line.strip() == '!series_matrix_table_begin': break
                    yield line

    studies, members, matrix, soft = [], [], {}, []
    current, in_table = None, False
    for line in header_lines():
        stripped = line.strip()
        if stripped == '!sample_table_begin': in_table = True; continue
        if stripped == '!sample_table_end': in_table = False; continue
        if in_table: continue
        parts = next(csv.reader([line], delimiter='\t'))
        if len(parts) > 1 and re.fullmatch(r'!Sample_[A-Za-z0-9_]+', parts[0]):
            matrix.setdefault(parts[0][8:], []).append(parts[1:])
            continue
        key, separator, value = line.partition(' = ')
        if not separator:
            if len(parts) != 2: continue
            key, value = parts
        value = value.strip().strip('"')
        if key.startswith('^') and current is not None:
            soft.append(current); current = None
        if key in {'^SERIES', '!Series_geo_accession'}: studies.append(value)
        if key == '!Series_sample_id': members.extend(value.split())
        if key == '^SAMPLE': current = {'geo_accession':[value]}
        elif key.startswith('!Sample_') and current is not None:
            current.setdefault(key[8:], []).append(value)
    if current is not None: soft.append(current)
    if set(studies) != {study}: raise ValueError('Recovery primary study identity differs')
    if matrix and soft: raise ValueError('Recovery primary source mixes sample layouts')
    if matrix:
        identity_rows = matrix.get('geo_accession', [])
        if len(identity_rows) != 1: raise ValueError('Primary Series GSM row is missing/duplicated')
        ids = identity_rows[0]
        if any(len(row) != len(ids) for rows in matrix.values() for row in rows):
            raise ValueError('Primary Series sample field widths differ')
        soft = [{key:[row[i] for row in rows] for key,rows in matrix.items()} for i in range(len(ids))]
    samples = {}
    for sample in soft:
        ids = sample.get('geo_accession', [])
        if (not ids or len(set(ids)) != 1 or not re.fullmatch(r'GSM\d+', ids[0])
                or ids[0] in samples):
            raise ValueError('Primary Series sample identities are missing/conflicting/duplicated')
        samples[ids[0]] = sample
    if not samples: raise ValueError('Primary Series sample identity evidence is missing')
    if sample_source is not None and not members:
        raise ValueError('Separate sample records require declared primary Series membership')
    if members and (len(set(members)) != len(members) or set(members) != set(samples)):
        raise ValueError('Primary Series declared membership differs from sample records')
    return samples


def reconcile_recovery_metadata(metadata, series, study, sample_source=None):
    """Require local identity labels to agree with primary GEO sample fields."""
    rows, _ = read_source_csv(metadata)
    primary = recovery_primary_samples(series, study, sample_source)
    by_gsm = {}
    for row in rows:
        ids = [row[k] for k in ('GSM','gsm_id','geo_accession') if k in row]
        if not ids or len(set(ids)) != 1 or not re.fullmatch(r'GSM\d+', ids[0]) or ids[0] in by_gsm:
            raise ValueError('Metadata has missing, conflicting or duplicate GSM identities')
        if row.get('gse_id', study) != study: raise ValueError('Metadata mixes studies')
        by_gsm[ids[0]] = row
    if set(by_gsm) != set(primary): raise ValueError('Metadata membership differs from primary Series')
    for gsm, row in by_gsm.items():
        if 'title' in row and primary[gsm].get('title') != [row['title']]:
            raise ValueError('Metadata title differs from primary Series for '+gsm)
        for field, value in row.items():
            if field.startswith('characteristics_ch1'):
                label, separator, _ = value.partition(':')
                if separator and label.strip().casefold() in {'sample id','sample identifier','sample code'}:
                    if value not in primary[gsm].get('characteristics_ch1', []):
                        raise ValueError('Metadata sample label differs from primary Series for '+gsm)
    return by_gsm, primary


def recovery_alias_owners(primary):
    """Retain every possible supported alias, including labels absent locally."""
    owners = {}
    for gsm, sample in primary.items():
        tokens = {gsm}
        for title in sample.get('title', []):
            match = re.search(r'\[([^\[\]]+)\]\s*$', title)
            if match: tokens.add(match[1])
        for value in sample.get('characteristics_ch1', []):
            label, separator, token = value.partition(':')
            if separator and label.strip().casefold() in {'sample id','sample identifier','sample code'}:
                tokens.add(token.strip())
        for value in sample.get('supplementary_file', []):
            name = Path(urllib.parse.urlparse(value).path).name
            match = re.fullmatch(r'(GSM\d+)_([^_]+(?:_[^_]+)*)_(?:Grn|Red)\.idat\.gz', name)
            if match and match[1] == gsm:
                tokens.add(match[2])
        for token in tokens: owners.setdefault(token, set()).add(gsm)
    return owners


def propose_dnam_recovery(study, failed_matrix, supplement, metadata, series, *, sample_source=None):
    """Propose a lossless replacement only for a verified empty source table.

    Source title tokens identify assay columns, never biological independence.
    All samples must map uniquely; there is no positional or intersection join.
    """
    if not re.fullmatch(r'GSE\d+', study): raise ValueError('Invalid recovery study identity')
    paths = {key: Path(value).resolve() for key, value in {
        'failed_matrix': failed_matrix, 'supplement': supplement, 'metadata': metadata, 'series': series}.items()}
    if sample_source is not None: paths['sample_records'] = Path(sample_source).resolve()
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
    if paths['supplement'] in {paths['failed_matrix'], paths['metadata'], paths['series']}:
        raise ValueError('Recovery requires a distinct supplementary source')
    original, candidate = (beta_table_header(paths[name]) for name in ('failed_matrix', 'supplement'))
    if not original['empty'] or candidate['empty'] or candidate['geo']:
        raise ValueError('Recovery requires an empty original and a nonempty plain supplementary beta table')
    opener = gzip.open if str(series).endswith('.gz') else open
    with opener(series, 'rt') as handle:
        declarations = []
        identities = []
        for line in handle:
            if line.strip() == '!series_matrix_table_begin': break
            key, separator, value = line.partition(' = ')
            if not separator:
                parts = next(csv.reader([line], delimiter='\t'))
                if len(parts) == 2: key, value = parts
                else: continue
            value = value.strip().strip('"')
            if key in {'^SERIES', '!Series_geo_accession'}: identities.append(value)
            if key == '!Series_supplementary_file': declarations.append(value)
    if set(identities) != {study}:
        raise ValueError('GEO source does not establish the nominated study identity')
    expected_prefix = '/geo/series/' + study[:-3] + 'nnn/' + study + '/suppl/'
    urls = []
    for value in declarations:
        parsed = urllib.parse.urlparse(value)
        if (parsed.scheme in {'https','ftp'} and parsed.hostname == 'ftp.ncbi.nlm.nih.gov'
                and parsed.path == expected_prefix + paths['supplement'].name
                and not parsed.query and not parsed.fragment):
            urls.append('https://ftp.ncbi.nlm.nih.gov' + parsed.path)
    if len(set(urls)) != 1:
        raise ValueError('Supplementary file is not uniquely linked by this GEO source')
    by_gsm, primary = reconcile_recovery_metadata(paths['metadata'], paths['series'], study, paths.get('sample_records'))
    primary_owners = recovery_alias_owners(primary)
    aliases = {}
    for gsm, row in by_gsm.items():
        options = [(gsm, 'GSM', gsm, 'explicit_GSM')]
        title = row.get('title', '')
        match = re.search(r'\[([^\[\]]+)\]\s*$', title)
        if match: options.append((match[1], 'title', title, 'exact_bracketed_source_title_token'))
        for field, value in row.items():
            if field.startswith('characteristics_ch1'):
                label, sep, token = value.partition(':')
                if sep and label.strip().casefold() in {'sample id','sample identifier','sample code'}:
                    options.append((token.strip(), field, value, 'explicit_source_sample_label'))
        for value in primary[gsm].get('supplementary_file', []):
            name = Path(urllib.parse.urlparse(value).path).name
            match = re.fullmatch(r'(GSM\d+)_([^_]+(?:_[^_]+)*)_(?:Grn|Red)\.idat\.gz', name)
            if match and match[1] == gsm:
                options.append((match[2], 'supplementary_file', value,
                                'exact_sample_supplementary_idat_coordinate'))
        for token, field, quote, rule in options:
            aliases.setdefault(token, []).append({'matrix_sample_id': token, 'gsm_id': gsm,
                                                  'field': field, 'quote': quote, 'mapping_rule': rule})
    original_ids = original['header'][1:] if original['feature_column'] else original['header']
    candidate_ids = candidate['header'][1:] if candidate['feature_column'] else candidate['header']
    if set(original_ids) != set(by_gsm):
        raise ValueError('Empty input membership differs from source metadata')
    crosswalk = []
    for token in candidate_ids:
        matches = aliases.get(token, [])
        if len({m['gsm_id'] for m in matches}) != 1:
            raise ValueError('Supplement sample has no unique source binding: ' + token)
        if primary_owners.get(token) != {matches[0]['gsm_id']}:
            raise ValueError('Supplement sample alias is ambiguous in primary Series: ' + token)
        crosswalk.append(matches[0])
    if len(crosswalk) != len(by_gsm) or {m['gsm_id'] for m in crosswalk} != set(by_gsm):
        raise ValueError('Supplement does not cover each GEO sample exactly once')
    return {'schema': 'source-linked-dnam-recovery-v4', 'study_id': study,
            'detection_status': 'FAIL', 'detected_error': 'Original molecular table has zero feature rows',
            'proposed_action': 'reconstruct_beta_input_from_linked_supplement',
            'original_feature_count': 0, 'source_url': urls[0], 'crosswalk': crosswalk,
            'source_identity_verification': {'status':'PASS', 'samples':len(by_gsm),
                'scope':'GSM membership, metadata titles and explicit sample labels reconciled to primary Series'},
            'source_files': {key: {'path': str(path), 'sha256': sha(path), 'size_bytes': path.stat().st_size,
                                   'mtime_ns': path.stat().st_mtime_ns} for key, path in paths.items()},
            'source_layout': candidate, 'scientific_status': 'CANNOT VERIFY',
            'scope': 'Input recovery only; feature identifiers retained without asserting all are CpGs; no imputation, biological-unit inference, clock validity or result release'}


def reconstruct_dnam_input(proposal, output):
    """Stream and validate every source value while preserving its text token."""
    output = Path(output)
    for record in proposal['source_files'].values():
        if sha(record['path']) != record['sha256']: raise ValueError('Recovery source changed')
    source = Path(proposal['source_files']['supplement']['path'])
    target = output / 'recovered-beta.tsv.gz'
    header = proposal['source_layout']['header']
    ids = [row['gsm_id'] for row in proposal['crosswalk']]
    seen, finite_counts, missing, count = set(), np.zeros(len(ids), dtype=np.int64), 0, 0
    categories = Counter()
    with _beta_supplement_rows(source, proposal['source_layout']) as reader, target.open('xb') as raw_out:
        if next(reader) != header: raise ValueError('Supplement header changed')
        with gzip.GzipFile(fileobj=raw_out, filename='', mode='wb', compresslevel=1, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding='utf-8', newline='') as stream:
                writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
                writer.writerow(['ID_REF', *ids])
                for row in reader:
                    if len(row) != len(header) + (0 if proposal['source_layout']['feature_column'] else 1):
                        expected = len(header) + (0 if proposal['source_layout']['feature_column'] else 1)
                        raise ValueError(f'Supplement row {count+1} has {len(row)} fields; expected {expected}')
                    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.:+-]*', row[0]) or row[0] in seen:
                        raise ValueError(f'Supplement row {count+1} has invalid/duplicate feature identity: {row[0]!r}')
                    # Source assay tables can contain non-cg features. Retain all
                    # identities; the clock adapter selects its exact required CpGs.
                    category = 'cg_digits' if re.fullmatch(r'cg\d+', row[0]) else 'ch_prefix' if row[0].startswith('ch.') else 'other'
                    categories[category] += 1
                    values = np.array([np.nan if value.strip().lower() in {'','na','nan','null'} else float(value)
                                       for value in row[1:]], dtype=float)
                    finite = np.isfinite(values)
                    if np.isinf(values).any() or np.any((values[finite] < 0) | (values[finite] > 1)):
                        raise ValueError('Supplement has nonfinite or out-of-range beta measurements')
                    finite_counts += finite; missing += int((~finite).sum())
                    seen.add(row[0]); count += 1; writer.writerow(row)
                    if count % 100000 == 0: print(f'Recovery validated {count} feature rows', flush=True)
    if not count or (finite_counts == 0).any():
        raise ValueError('Recovered input has no measured features for one or more samples')
    for record in proposal['source_files'].values():
        if sha(record['path']) != record['sha256']: raise ValueError('Recovery source changed during reconstruction')
    result = {'path': str(target), 'sha256': sha(target), 'feature_count': count, 'sample_count': len(ids),
              'missing_values': missing, 'finite_values_per_sample': dict(zip(ids, finite_counts.tolist())),
              'feature_id_categories': dict(categories), 'feature_annotation_status': 'CANNOT VERIFY',
              'value_changes': 0, 'imputation': 'none', 'normalization': 'none',
              'sample_selection': 'All declared samples, source column order preserved',
              'repair_status': 'PROPOSED', 'scientific_status': 'CANNOT VERIFY'}
    write_new(output / 'reconstruction.json', result)
    return result


def acquire_dnam(study, source_samples, matrix, output, *, metadata=None, series=None, supplement=None, sample_source=None, idat_provenance=None):
    idat_binding = None
    if idat_provenance is not None:
        if matrix is None or supplement is not None:
            raise ValueError('IDAT provenance requires its derived matrix and cannot be combined with supplement recovery')
        from audit_system.report_2026_09_10 import verify_dnam_idat_provenance
        verified = verify_dnam_idat_provenance(idat_provenance, matrix_path=matrix)
        if verified['study_id'] != study or set(verified['sample_ids']) != set(source_samples):
            raise ValueError('IDAT provenance study/sample membership differs')
        for role, path in [('series', series), ('samples', sample_source)]:
            if path is None or sha(path) != verified['files'][role]['sha256']:
                raise ValueError('IDAT provenance physical source differs from the current '+role)
        idat_binding = {'path': str(Path(idat_provenance).resolve()), 'sha256': sha(idat_provenance)}
    if matrix is None:
        prefix = study[:-3] + 'nnn'
        base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{study}/matrix/'
        listing = fetch(base, output / 'matrix-index.html')
        matches = sorted(set(re.findall(r'href="([^"/]+\.txt\.gz)"', listing.read_text())))
        if len(matches) != 1:
            raise ValueError('GEO Series Matrix is absent or ambiguous; supply a source-linked beta matrix with --matrix')
        matrix = fetch(base + matches[0], output / matches[0], limit=2_000_000_000)
    if supplement is not None:
        if metadata is None or series is None:
            raise ValueError('Supplement recovery requires primary metadata and GEO series sources')
        from audit_system.report_2026_09_10 import verify_dnam_recovery
        proposal = propose_dnam_recovery(study, matrix, supplement, metadata, series, sample_source=sample_source)
        if {row['gsm_id'] for row in proposal['crosswalk']} != set(source_samples):
            raise ValueError('Recovered source membership differs from requested GEO samples')
        write_new(output / 'recovery-proposal.json', proposal)
        result = reconstruct_dnam_input(proposal, output)
        verification = verify_dnam_recovery(proposal, result)
        write_new(output / 'recovery-verification.json', verification)
        result.update(source=str(supplement), source_sha256=proposal['source_files']['supplement']['sha256'],
                      repair_status=verification['status'], repaired_original=str(matrix))
        write_new(output / 'molecular-input.json', result)
        return result
    frame = read_matrix(matrix, beta=True)
    if set(frame.columns) != set(source_samples): raise ValueError('DNAm/GEO sample membership differs')
    if any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.:+-]*', str(p)) for p in frame.index):
        raise ValueError('DNAm adapter requires unambiguous source feature identifiers')
    target = output / 'beta.csv'; frame.to_csv(target, mode='x')
    result = {'path': str(target), 'sha256': sha(target), 'source': str(matrix),
              'source_sha256': sha(matrix), 'feature_count': len(frame), 'sample_count': len(frame.columns)}
    if idat_binding is not None:
        result['idat_provenance'] = idat_binding
    write_new(output / 'molecular-input.json', result)
    return result


def _dnam_reference_profile_sha256():
    payload = json.dumps(DNAM_REFERENCE_PROFILE, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    return hashlib.sha256(payload).hexdigest()


def verify_dnam_reference_source(path, *, max_source_bytes=2000000000):
    """Check bounded GSE40279 source identity, metadata and every table byte.

    The sole permitted source drift is in complete !Series_relation metadata
    records. Their exact values are retained in the returned evidence.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError('Cannot proceed - required data file not found: '+str(path))
    if max_source_bytes <= 0:
        raise ValueError('GSE40279 reference requires a positive source-byte limit')
    original_size = path.stat().st_size
    if original_size == 0 or original_size > max_source_bytes:
        raise ValueError('GSE40279 reference exceeds the explicit source-byte limit')
    original_sha = sha(path)
    metadata_bytes = bytearray(); table_sha = hashlib.sha256()
    table_bytes = decompressed_bytes = feature_rows = 0
    header = None; stage = 'metadata'
    try:
        with gzip.open(path, 'rb') as handle:
            while True:
                line = handle.readline(DNAM_REFERENCE_MAX_LINE_BYTES+1)
                if not line: break
                if len(line) > DNAM_REFERENCE_MAX_LINE_BYTES:
                    raise ValueError('GSE40279 source record exceeds the bounded line size')
                decompressed_bytes += len(line)
                if decompressed_bytes > DNAM_REFERENCE_MAX_DECOMPRESSED_BYTES:
                    raise ValueError('GSE40279 source exceeds the bounded decompressed size')
                if stage == 'metadata':
                    metadata_bytes.extend(line)
                    if len(metadata_bytes) > DNAM_REFERENCE_MAX_METADATA_BYTES:
                        raise ValueError('GSE40279 source metadata exceeds its byte limit')
                    if line.strip() == b'!series_matrix_table_begin': stage = 'header'
                elif stage == 'header':
                    header = next(csv.reader([line.decode('utf-8')], delimiter='\t', strict=True))
                    table_sha.update(line); table_bytes += len(line); stage = 'table'
                elif stage == 'table':
                    if line == b'!series_matrix_table_end\n': stage = 'done'
                    else:
                        table_sha.update(line); table_bytes += len(line); feature_rows += 1
                else:
                    raise ValueError('GSE40279 source has bytes after its table end marker')
    except (EOFError, OSError, UnicodeError, csv.Error) as error:
        raise ValueError('GSE40279 source gzip/table is malformed: '+str(error)) from error
    if stage != 'done' or header is None:
        raise ValueError('GSE40279 source lacks a complete GEO matrix table')
    try:
        records = list(csv.reader(io.StringIO(metadata_bytes.decode('utf-8'), newline=''),
                                  delimiter='\t', strict=True))
    except (UnicodeError, csv.Error) as error:
        raise ValueError('GSE40279 source metadata is malformed: '+str(error)) from error
    if not records or records[-1] != ['!series_matrix_table_begin']:
        raise ValueError('GSE40279 source metadata/table boundary differs')
    records = records[:-1]
    fields = {}
    for row in records:
        if row: fields.setdefault(row[0], []).append(row[1:])
    ids = fields.get('!Sample_geo_accession', [[]])
    if (fields.get('!Series_geo_accession') != [['GSE40279']] or len(ids) != 1
            or len(ids[0]) != DNAM_REFERENCE_PROFILE['sample_count']
            or len(set(ids[0])) != len(ids[0]) or any(not re.fullmatch(r'GSM[0-9]+', v) for v in ids[0])
            or header != ['ID_REF']+ids[0]
            or fields.get('!Sample_platform_id') != [['GPL13534']*len(ids[0])]
            or fields.get('!Sample_organism_ch1') != [['Homo sapiens']*len(ids[0])]
            or fields.get('!Sample_molecule_ch1') != [['genomic DNA']*len(ids[0])]):
        raise ValueError('GSE40279 source identity/sample/platform header differs')
    canonical = [row for row in records if not row or row[0] != '!Series_relation']
    relations = [row for row in records if row and row[0] == '!Series_relation']
    actual = {
        'ordered_canonical_metadata_record_count': len(canonical),
        'ordered_canonical_metadata_sha256': hashlib.sha256(json.dumps(
            canonical, ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode()).hexdigest(),
        'table_bytes_including_header_excluding_end_marker': table_bytes,
        'table_sha256': table_sha.hexdigest(),
        'feature_rows': feature_rows,
        'ordered_header_sha256': hashlib.sha256(json.dumps(
            header, ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode()).hexdigest(),
        'sample_count': len(ids[0]),
    }
    for key, value in actual.items():
        if value != DNAM_REFERENCE_PROFILE[key]:
            raise ValueError('GSE40279 source profile differs: '+key)
    if path.stat().st_size != original_size or sha(path) != original_sha:
        raise PipelineIntegrityError('GSE40279 source changed during profile verification')
    return {'profile_schema': DNAM_REFERENCE_PROFILE['schema'],
            'profile_sha256': _dnam_reference_profile_sha256(),
            'source_sha256': original_sha, 'size_bytes': original_size,
            'decompressed_bytes': decompressed_bytes, **actual,
            'source_sample_ids': ids[0], 'source_platform': 'GPL13534',
            'relation_records': relations}


def acquire_dnam_reference(cache, *, request=None, max_source_bytes=2000000000):
    """Obtain profile-verified GSE40279 without calculating means or clocks."""
    cache = Path(cache); cache.mkdir(parents=True, exist_ok=True)
    if request is None:
        reference = cache / 'GSE40279_series_matrix.txt.gz'
        if not reference.exists(): fetch(DNAM_REFERENCE_URL, reference, limit=max_source_bytes)
        verify_dnam_reference_source(reference, max_source_bytes=max_source_bytes)
        return reference
    profile_sha = _dnam_reference_profile_sha256()
    directory = cache / 'objects' / ('gse40279-profile-v2-'+profile_sha[:24]); directory.mkdir(parents=True, exist_ok=True)
    reference = directory / 'GSE40279_series_matrix.txt.gz'; receipt = directory / 'reference.json'
    if receipt.exists():
        sidecar = Path(str(receipt)+'.sha256')
        if not sidecar.is_file() or sha(receipt) != sidecar.read_text().strip():
            raise ValueError('Cached GSE40279 reference receipt changed')
        saved = json.loads(receipt.read_text())
        if (saved.get('schema') != 'clockbase-bounded-dnam-reference-v2'
                or saved.get('reference_study') != 'GSE40279'
                or saved.get('url') != DNAM_REFERENCE_URL or saved.get('profile_sha256') != profile_sha
                or saved.get('path') != str(reference.resolve()) or not reference.is_file()):
            raise ValueError('Cached GSE40279 reference or receipt changed')
        observed = verify_dnam_reference_source(reference, max_source_bytes=max_source_bytes)
        if (saved.get('sha256') != observed['source_sha256']
                or any(saved.get(key) != observed[key] for key in observed)):
            raise ValueError('Cached GSE40279 reference profile or receipt changed')
        return reference
    part = reference.with_name(reference.name+'.partial')
    if reference.exists() or part.exists(): raise FileExistsError('Unverified reference retained; no implicit retry')
    response = request('GET', DNAM_REFERENCE_URL, part, max_source_bytes, {})
    observed = verify_dnam_reference_source(part, max_source_bytes=max_source_bytes)
    part.rename(reference)
    value = {'schema': 'clockbase-bounded-dnam-reference-v2', 'reference_study': 'GSE40279',
             'path': str(reference.resolve()), 'url': DNAM_REFERENCE_URL,
             'sha256': observed['source_sha256'], 'original_producer_sha256': DNAM_REFERENCE_SHA256,
             'original_producer_size_bytes': DNAM_REFERENCE_BYTES,
             'response_receipt': response['receipt'], 'means_computed': False,
             'scientific_status': 'CANNOT VERIFY', **observed}
    write_new(receipt, value)
    with Path(str(receipt)+'.sha256').open('x') as sidecar:
        sidecar.write(sha(receipt)+'\n')
    return reference


def reference_means(reference, cache):
    """Reproduce the inspected producer's GSE40279 population-mean reference.

    Read all reference measurements in chunks; no stand-in imputation vector.
    """
    import gzip
    cache.mkdir(parents=True, exist_ok=True)
    if reference is None:
        reference = acquire_dnam_reference(cache)
    checksum = sha(reference)
    target = cache / ('GSE40279-means-' + checksum[:16] + '-' + MATRIX_READER_CONTRACT + '.csv')
    receipt = target.with_suffix('.json')
    if receipt.exists():
        record = json.loads(receipt.read_text())
        if (record.get('reader_contract') != MATRIX_READER_CONTRACT
                or record['source_sha256'] != checksum or sha(target) != record['means_sha256']):
            raise ValueError('DNAm reference cache changed')
        return pd.read_csv(target, index_col=0, converters={0: str}, float_precision='round_trip')['mean'], record
    opener = gzip.open if str(reference).endswith('.gz') else open
    all_means = []; seen = set(); observed_features = []
    with opener(reference, 'rt', newline='') as handle:
        source = _inspect_matrix_stream(handle, geo=True, delimiter='\t')
        header = source['header']
        if source['accessions'] != ['GSE40279'] or len(header[1:]) != 656:
            raise ValueError('Reference accession or 656-sample membership dimension differs from Biolearn GSE40279 definition')
        handle.seek(source['data_start'])
        for chunk in pd.read_csv(handle, sep='\t', names=header, header=None, index_col=0,
                                 converters={0: str}, nrows=len(source['features']), chunksize=10000,
                                 float_precision='round_trip'):
            if list(chunk.columns) != header[1:]:
                raise ValueError('AUDIT_FAIL: parsed reference sample identities differ from source')
            require_matrix(chunk, beta=True)
            if seen.intersection(chunk.index): raise ValueError('Duplicate CpGs in reference')
            seen.update(chunk.index); observed_features.extend(chunk.index)
            all_means.append(chunk.mean(axis=1))
    if not all_means or observed_features != source['features']:
        raise ValueError('Reference is incomplete or parsed feature order differs from source')
    if sha(reference) != checksum:
        raise PipelineIntegrityError('Reference source changed during validation/read')
    means = pd.concat(all_means).rename('mean')
    means.to_csv(target, mode='x')
    record = {'reference_study': 'GSE40279', 'source': str(reference), 'source_sha256': checksum,
              'source_mtime_ns': reference.stat().st_mtime_ns, 'means_sha256': sha(target),
              'reference_samples': len(header)-1, 'reference_features': len(means), 'created_unix': time.time(),
              'reader_contract': MATRIX_READER_CONTRACT}
    write_new(receipt, record)
    return means, record


def preprocess_dnam(raw, reference, cache, output, *, idat_provenance=None):
    from biolearn.imputation import impute_from_standard
    frame = read_matrix(raw, beta=True)
    if idat_provenance is not None:
        from audit_system.report_2026_09_10 import verify_dnam_idat_provenance
        if not isinstance(idat_provenance, dict) or set(idat_provenance) != {'path','sha256'}:
            raise ValueError('IDAT preprocessing requires its acquisition path/hash binding')
        if sha(idat_provenance['path']) != idat_provenance['sha256']:
            raise ValueError('IDAT provenance changed after acquisition')
        verified_idat = verify_dnam_idat_provenance(idat_provenance['path'])
        original = read_matrix(verified_idat['files']['matrix']['path'], beta=True)
        if not frame.equals(original):
            raise ValueError('IDAT stage-02 beta values differ from the validated derived matrix')
    means, reference_record = reference_means(reference, cache)
    joined = frame.reindex(means.index.union(frame.index))
    imputed = impute_from_standard(joined, means)
    verified = joined.where(joined.notna(), means, axis=0)
    if not imputed.index.equals(verified.index) or not imputed.columns.equals(verified.columns):
        raise ValueError('Population-imputation sample/feature identities changed')
    if not np.allclose(imputed, verified, equal_nan=True, rtol=0, atol=0):
        raise ValueError('Independent population-imputation calculation differs')
    path = output / 'imputed_beta.csv'; imputed.to_csv(path, mode='x')
    record = {'path': str(path), 'raw': str(raw), 'reference': reference_record,
              'status': 'PASS', 'raw_methylation_preprocessing_recomputed': True,
              'independent_imputation_verification': 'PASS', 'idat_processing': False}
    if idat_provenance is not None:
        record.update(idat_processing=True, idat_provenance=idat_provenance,
                      idat_measurement=verified_idat['measurement'])
    write_new(output / 'preprocessing.json', record)
    return record


def predict_dnam(preprocessed, output):
    from biolearn.data_library import GeoData
    from biolearn.model_gallery import ModelGallery
    frame = read_matrix(Path(preprocessed['path']), beta=True)
    raw = read_matrix(Path(preprocessed['raw']), beta=True)
    predictions = pd.DataFrame(index=frame.columns); checks = []
    for name in ['Horvathv1', 'Hannum', 'PhenoAge']:
        model = ModelGallery().get(name, imputation_method='none')
        coefficients = model.coefficients.iloc[:, 0]
        sites = [s for s in coefficients.index if s != 'intercept']
        observed = raw.reindex(sites).notna().sum(axis=0)
        if (observed == 0).any(): raise ValueError('Sample has no observed sites for clock: ' + name)
        values = frame.reindex(sites)
        if not np.isfinite(values.to_numpy(float)).all():
            raise ValueError('CpG values unavailable after source-reference imputation: ' + name)
        computed = model.predict(GeoData(metadata=pd.DataFrame(index=frame.columns), dnam=frame.copy()))['Predicted'].reindex(frame.columns).to_numpy(float)
        independent = coefficients.loc[sites].to_numpy(float) @ values.to_numpy(float) + float(coefficients.get('intercept', 0))
        if name == 'Horvathv1':
            transformed = independent + .696
            independent = np.where(transformed < 0, 21*np.exp(transformed)-1, 21*transformed+20)
        if not np.isfinite(computed).all() or not np.allclose(computed, independent, rtol=1e-10, atol=1e-10):
            raise ValueError('Independent methylation prediction failed: ' + name)
        predictions[name] = computed
        coefficient_path = output / (name + '-coefficients.csv')
        coefficients.to_csv(coefficient_path, mode='x')
        checks.append({'clock': name, 'numerical_status': 'PASS', 'required_features': len(sites),
                       'observed_features_min': int(observed.min()), 'measurement_applicability': 'CANNOT VERIFY',
                       'measurement_sources': {
                           'raw_beta': {'path': str(Path(preprocessed['raw']).resolve()), 'sha256': sha(preprocessed['raw'])},
                           'coefficients': {'path': str(coefficient_path.resolve()), 'sha256': sha(coefficient_path)}},
                       'max_abs_independent_difference': float(np.max(np.abs(computed-independent)))})
        if preprocessed.get('idat_provenance') is not None:
            absent = int((~pd.Index(sites).isin(raw.index)).sum())
            missing = raw.reindex(sites).isna().sum(axis=0)
            checks[-1]['idat_feature_coverage'] = {
                str(sample): {'required': len(sites), 'observed_after_qc': int(observed[sample]),
                    'absent_from_array': absent, 'failed_qc_or_nonfinite': int(missing[sample])-absent,
                    'filled_from_reference': int(missing[sample]), 'remaining_nonfinite': 0}
                for sample in frame.columns}
    path = output / 'clock_predictions.csv'; predictions.rename_axis('sample_id').to_csv(path, mode='x')
    result = {'path': str(path), 'clocks': list(predictions.columns), 'sample_count': len(predictions),
              'checks': checks, 'clock_predictions_recomputed': True, 'scientific_status': 'CANNOT VERIFY'}
    if preprocessed.get('idat_provenance') is not None:
        binding = preprocessed['idat_provenance']
        if sha(binding['path']) != binding['sha256']:
            raise ValueError('IDAT provenance changed during clock prediction')
        result['idat_provenance'] = binding
    write_new(output / 'predictions.json', result)
    return result


def rna_rle_factors(counts):
    """Independently calculate edgeR's RLE factors for a fixed count matrix."""
    require_matrix(counts, beta=False)
    libraries = counts.sum(axis=0).astype(float)
    positive = counts.loc[(counts > 0).all(axis=1)].astype(float)
    if positive.empty or (libraries <= 0).any():
        raise ValueError('RLE requires positive library sizes and genes observed in every sample')
    geometric = np.exp(np.log(positive).mean(axis=1))
    ratios = positive.div(geometric, axis=0).median(axis=0) / libraries
    factors = ratios / np.exp(np.log(ratios).mean())
    if not np.isfinite(factors).all() or (factors <= 0).any():
        raise ValueError('RLE normalization factors are not finite positive numbers')
    return pd.DataFrame({'library_size': libraries, 'norm_factor': factors,
                         'count_multiplier': 1e7 / (libraries * factors)}).rename_axis('sample_id')


def rna_preprocessing_stratum(packet, sample_ids):
    """Require one explicit source stratum; this is not biological validation."""
    from audit_system.source_packet import structured_characteristics, normalized_source_label
    if not isinstance(packet, dict) or not packet.get('samples'):
        raise ValueError('Published mouse preprocessing requires a pinned source packet')
    samples = packet['samples']
    ids = [sample['sample_id'] for sample in samples]
    if len(set(ids)) != len(ids) or set(ids) != set(sample_ids):
        raise ValueError('RNA preprocessing/source sample membership differs')
    tissues = set()
    for sample in samples:
        if (sample.get('organism_ch1') != 'Mus musculus'
                or normalized_source_label(sample.get('molecule_ch1', '')) != 'total rna'):
            raise ValueError('Initial mouse preprocessing requires explicit mouse total-RNA source')
        fields = structured_characteristics(sample)
        declared = {normalized_source_label(row['value']) for key in ('tissue', 'tissue type')
                    for row in fields.get(key, [])}
        if len(declared) != 1 or not next(iter(declared)):
            raise ValueError('RNA preprocessing tissue stratum is missing or ambiguous')
        tissues.update(declared)
    if len(tissues) != 1:
        raise ValueError('RNA preprocessing cannot pool distinct source tissue strata')
    return {'species': 'Mus musculus', 'molecule': 'total RNA', 'tissue': next(iter(tissues)),
            'sample_ids': list(sample_ids), 'reference_ids': list(sample_ids),
            'reference_rule': 'Per-gene median across every source sample in this stratum',
            'matched_controls_verified': False, 'biological_applicability': 'CANNOT VERIFY'}


def mapped_rna_counts(raw, table, *, raw_filter_first=False):
    """Apply the declared order; raw filtering and mapped summation do not commute."""
    entrez = pd.to_numeric(table['Entrez'].dropna(), errors='raise')
    if (not np.isfinite(entrez).all() or (entrez <= 0).any()
            or not np.equal(entrez, np.floor(entrez)).all()):
        raise ValueError('Entrez mapping identifiers must be positive finite integers')
    mapping = entrez.astype(int).astype(str)
    if mapping.index.has_duplicates: raise ValueError('Ambiguous Ensembl mapping')
    selected = raw.loc[(raw >= 10).sum(axis=1) >= raw.shape[1] * .2] if raw_filter_first else raw
    keep = selected.index.intersection(mapping.index, sort=False)
    expr = selected.loc[keep].copy()
    expr.index = mapping.loc[keep].values
    expr = expr.groupby(level=0).sum()
    if not raw_filter_first:
        expr = expr.loc[(expr >= 10).sum(axis=1) >= expr.shape[1] * .2]
    return expr


def snapshot_rna_implementations(paths, output):
    """Keep executed code bytes independent of later worktree evolution."""
    directory = output / 'implementation'
    directory.mkdir(exist_ok=False)
    snapshots, originals = {}, {}
    for label, path in paths.items():
        path = Path(path).resolve()
        originals[str(path)] = sha(path)
        target = directory / (label + path.suffix)
        with path.open('rb') as reader, target.open('xb') as writer:
            shutil.copyfileobj(reader, writer)
        if sha(target) != originals[str(path)] or sha(path) != originals[str(path)]:
            raise PipelineIntegrityError('RNA implementation changed during snapshot: ' + label)
        snapshots[label] = {'path': str(target.resolve()), 'sha256': originals[str(path)],
                            'executed_path': str(path)}
    return snapshots, originals


def preprocess_rna(counts, species, output, *, normalization='clockbase_legacy', source=None):
    """Preserve legacy defaults; allow a source-bound mouse operation-order check."""
    if normalization not in {'clockbase_legacy', 'rle_sensitivity', 'published_mouse_rle'}:
        raise ValueError('Unknown RNA normalization profile')
    if normalization == 'published_mouse_rle' and species != 'mouse':
        raise ValueError('Published operation-order profile currently supports mouse only')
    declarations = {'raw_counts': {'path': str(Path(counts).resolve()), 'sha256': sha(counts)},
                    'gene_mapping': {'path': str(ASSETS / ('Gene_table_' + species + '.csv')),
                                     'sha256': sha(ASSETS / ('Gene_table_' + species + '.csv'))},
                    'feature_registry': {'path': str(ASSETS / 'Gene_list_all_4.6.txt'),
                                         'sha256': sha(ASSETS / 'Gene_list_all_4.6.txt')}}
    snapshots, executed_implementations = snapshot_rna_implementations({
        'preprocessing_python': Path(__file__).resolve(),
        'preprocessing_r': ASSETS / 'FUN.Preprocessing_for_clock_3.0.R',
        'ortholog_r': ASSETS / 'FUN.Map_to_orthologs.R'}, output)
    declarations.update(snapshots)
    if species == 'human':
        declarations['ortholog_mapping'] = {'path': str(ASSETS / 'Table_of_orthologs.csv'),
                                            'sha256': sha(ASSETS / 'Table_of_orthologs.csv')}
    raw = read_matrix(counts)
    stratum = None
    if normalization == 'published_mouse_rle':
        stratum = rna_preprocessing_stratum(source, list(raw.columns))
        write_new(output / 'source-packet.json', source)
        declarations['source_packet'] = {'path': str((output / 'source-packet.json').resolve()),
                                         'sha256': sha(output / 'source-packet.json')}
        for key in ('sample_source', 'metadata_source'):
            if source.get(key): declarations[key] = source[key]
    table = pd.read_csv(ASSETS / ('Gene_table_' + species + '.csv'), index_col=0)
    expr = mapped_rna_counts(raw, table, raw_filter_first=normalization == 'published_mouse_rle')
    if species == 'human':
        orth = pd.read_csv(ASSETS / 'Table_of_orthologs.csv', index_col=0)
        orth = orth.dropna(subset=['Entrez.Human', 'Entrez.Mouse'])
        orth.index = orth['Entrez.Human'].astype(int).astype(str)
        if orth.index.has_duplicates: raise ValueError('Ambiguous human-to-mouse ortholog mapping')
        keep = expr.index.intersection(orth.index, sort=False)
        expr = expr.loc[keep]
        expr.index = orth.loc[keep, 'Entrez.Mouse'].astype(int).astype(str)
    require_matrix(expr, beta=False)
    if expr.index.has_duplicates: raise ValueError('Ortholog mapping is many-to-one; adapter requires adjudication')
    mapped = output / 'mapped_filtered_counts.csv'
    expr.to_csv(mapped, mode='x')
    # Only the scaled branch is used by the original primary RNA analysis clocks.
    script = output / 'preprocess.R'
    normalization_code = ''
    if normalization in {'rle_sensitivity', 'published_mouse_rle'}:
        normalization_code = (
            'if (!requireNamespace("edgeR",quietly=TRUE)) stop("edgeR is required for RLE sensitivity")\n'
            'd <- edgeR::calcNormFactors(edgeR::DGEList(counts=x),method="RLE")\n'
            'nf <- d$samples$norm.factors\n'
            'ls <- colSums(x)\n'
            'multipliers <- 1e7/(ls*nf)\n'
            'write.csv(data.frame(sample_id=colnames(x),library_size=ls,norm_factor=nf,\n'
            'count_multiplier=multipliers),file.path(args[3],"rle_factors.csv"),row.names=FALSE)\n'
            'writeLines(as.character(utils::packageVersion("edgeR")),file.path(args[3],"edgeR-version.txt"))\n'
            'x <- sweep(x,2,multipliers,"*")\n')
    script.write_text('args <- commandArgs(trailingOnly=TRUE)\n'
        'source(file.path(args[1],"FUN.Preprocessing_for_clock_3.0.R"))\n'
        'x <- as.matrix(read.csv(args[2],row.names=1,check.names=FALSE))\n' + normalization_code +
        'g <- readLines(file.path(args[1],"Gene_list_all_4.6.txt"))\n'
        'y <- preprocessing_for_clock(x,g,log_transform=TRUE,show_plots=FALSE,scaling=TRUE,\n'
        'yugene=FALSE,differential=TRUE,control_group=colnames(x),function_path=paste0(args[1],"/"))\n'
        'write.csv(y$Scaled_data,file.path(args[3],"scaled.csv"))\n'
        'write.csv(y$Scaled_differential_data,file.path(args[3],"scaled_diff.csv"))\n')
    with (output / 'preprocess.log').open('x') as log:
        completed = subprocess.run(['Rscript', str(script), str(ASSETS), str(mapped), str(output)],
                                   stdout=log, stderr=log, timeout=180)
    if completed.returncode: raise RuntimeError('Preserved R preprocessing failed; see preprocess.log')
    transformed = expr
    if normalization in {'rle_sensitivity', 'published_mouse_rle'}:
        factors = rna_rle_factors(expr)
        observed_factors = pd.read_csv(output / 'rle_factors.csv', index_col='sample_id')
        if (list(factors.index) != list(observed_factors.index)
                or list(factors.columns) != list(observed_factors.columns)
                or not np.allclose(factors, observed_factors, rtol=1e-10, atol=1e-12)):
            raise ValueError('Independent RLE calculation differs from edgeR')
        transformed = expr.mul(factors['count_multiplier'], axis=1)
    expected = np.log10(transformed + 1)
    deviations = expected.std(axis=0, ddof=1)
    if not np.isfinite(expected.to_numpy()).all() or not np.isfinite(deviations).all() or (deviations <= 0).any():
        raise ValueError('RNA transformation has nonfinite values or zero sample standard deviation')
    expected = (expected - expected.mean(axis=0)) / expected.std(axis=0, ddof=1)
    genes = (ASSETS / 'Gene_list_all_4.6.txt').read_text().splitlines()
    if len(genes) != len(set(genes)): raise ValueError('Clock feature registry has duplicate identifiers')
    expected = expected.reindex(genes).T
    observed = pd.read_csv(output / 'scaled.csv', index_col=0)
    differential = pd.read_csv(output / 'scaled_diff.csv', index_col=0)
    expected_diff = expected - expected.median(axis=0)
    if list(observed.index) != list(expected.index) or list(observed.columns) != genes:
        raise ValueError('R preprocessing changed sample or feature order')
    if not np.allclose(expected, observed, equal_nan=True, rtol=1e-10, atol=1e-12) or not np.allclose(expected_diff, differential, equal_nan=True, rtol=1e-10, atol=1e-12):
        raise ValueError('Independent preprocessing calculation differs from preserved R function')
    record = {'scaled': str(output / 'scaled.csv'), 'scaled_diff': str(output / 'scaled_diff.csv'),
              'status': 'PASS', 'sample_count': len(observed), 'mapped_filtered_genes': len(expr),
              'clock_feature_columns': len(genes), 'log_transform': 'log10(count+1)',
              'normalization_profile': normalization,
              'filter_order': 'raw_then_mapping' if normalization == 'published_mouse_rle' else 'mapping_then_filter',
              'filter_count_minimum': 10, 'filter_sample_fraction': .2,
              'source_stratum': stratum, 'measurement_sources': declarations,
              'executed_implementation_sha256': executed_implementations,
              'scaling': 'Per sample across mapped/filter-passing genes, sample standard deviation',
              'differential_center': 'Gene median across all study samples, as in the inspected producer',
              'raw_count_preprocessing_recomputed': True, 'independent_preprocessing_verification': 'PASS',
              'historical_producing_version_verified': False,
              'RLE_applied': normalization != 'clockbase_legacy', 'YuGene_applied': False}
    if normalization in {'rle_sensitivity', 'published_mouse_rle'}:
        record.update(rle_factors=str(output / 'rle_factors.csv'),
                      rle_formula='count * 1e7 / (mapped-filtered library size * edgeR RLE factor)',
                      independent_rle_verification='PASS', scientific_status='CANNOT VERIFY',
                      scope=('Development sensitivity with the existing map-before-filter universe; not a full published-preprocessing repair'
                             if normalization == 'rle_sensitivity' else
                             'Source-bound published mouse operation order; independent full conformance and biological application are not established by this receipt'))
    outputs = ['mapped_filtered_counts.csv', 'scaled.csv', 'scaled_diff.csv']
    if normalization != 'clockbase_legacy': outputs.append('rle_factors.csv')
    record['output_sha256'] = {name: sha(output / name) for name in outputs}
    for declaration in declarations.values():
        if sha(declaration['path']) != declaration['sha256']:
            raise PipelineIntegrityError('RNA preprocessing input changed during execution')
    for path, digest in executed_implementations.items():
        if sha(path) != digest: raise PipelineIntegrityError('RNA implementation changed during preprocessing')
    write_new(output / 'preprocessing.json', record)
    return record


def predict_rna(preprocessed, species, output):
    import sklearn
    from audit_system.clock_replay_2026_09_10 import independent_prediction
    if sklearn.__version__ != '1.3.2': raise ValueError('RNA model runtime requires scikit-learn 1.3.2')
    producer_path = ASSETS / '2_run_rna_clock.py'
    snapshots, executed_implementations = snapshot_rna_implementations({
        'producer': producer_path, 'prediction_python': Path(__file__).resolve(),
        'parameter_verifier': Path(__file__).resolve().parents[1] / 'audit_system/clock_replay_2026_09_10.py'}, output)
    implementation_sources = {value['path']: value['sha256'] for value in snapshots.values()}
    spec = importlib.util.spec_from_file_location('preserved_rna_producer', producer_path)
    producer = importlib.util.module_from_spec(spec); spec.loader.exec_module(producer)
    selected = ['EN_Hazard_Mouse_All_All_scaleddiff.pkl', 'EN_Hazard_Rodents_All_All_scaleddiff.pkl'] if species == 'mouse' else ['EN_Hazard2_All_All_All_scaleddiff.pkl']
    declarations = dict(preprocessed.get('measurement_sources', {}))
    processed_path = Path(preprocessed['scaled_diff']).parent / 'preprocessing.json'
    if declarations:
        declarations['preprocessing'] = {'path': str(processed_path.resolve()), 'sha256': sha(processed_path)}
        declarations['transformed_input'] = {'path': str(Path(preprocessed['scaled_diff']).resolve()),
                                            'sha256': sha(preprocessed['scaled_diff'])}
        if (declarations['transformed_input']['sha256'] != preprocessed['output_sha256']['scaled_diff.csv']
                or json.loads(processed_path.read_text()) != preprocessed):
            raise PipelineIntegrityError('RNA preprocessing receipt or transformed input changed before prediction')
        for declaration in declarations.values():
            if sha(declaration['path']) != declaration['sha256']:
                raise PipelineIntegrityError('RNA measurement input changed before prediction')
    frame = pd.read_csv(preprocessed['scaled_diff'], index_col=0)
    predictions = pd.DataFrame(index=frame.index); checks = []
    for name in selected:
        model_digest = sha(ASSETS / name)
        model = producer.load_clock(ASSETS / name)
        clock_name = name.removeprefix('EN_').removesuffix('.pkl') + '_EN'
        computed = np.asarray(producer.predict_tAge(frame, model, clock_name)['tAge'])
        # Independent verifier is parameter based, with no sklearn prediction calls.
        verified = independent_prediction(model, frame, inverse_transform=False)
        if not np.isfinite(computed).all() or not np.allclose(computed, verified, rtol=1e-10, atol=1e-12):
            raise ValueError('Independent parameter calculation failed: ' + clock_name)
        predictions[clock_name + '_tAge'] = computed
        observed = np.isfinite(frame[model.feature_names].to_numpy(dtype=float)).sum(axis=1)
        checks.append({'clock': clock_name, 'model_sha256': model_digest, 'numerical_status': 'PASS',
                       'observed_features_min': int(observed.min()), 'required_features': len(model.feature_names),
                       'max_abs_independent_difference': float(np.max(np.abs(computed - verified))),
                       'measurement_applicability': 'CANNOT VERIFY'})
        if declarations:
            features = [str(value) for value in model.feature_names]
            observed_values = np.isfinite(frame[features].to_numpy(dtype=float))
            imputer, _, selector, estimator = [step for _, step in model.steps]
            training_keep = ~np.isnan(imputer.statistics_)
            scores = np.asarray(selector.scores_).copy()
            scores[np.isnan(scores)] = np.finfo(scores.dtype).min
            mask = np.ones(len(scores), dtype=bool) if selector.k == 'all' else np.zeros(len(scores), dtype=bool)
            if selector.k != 'all' and selector.k > 0:
                mask[np.argsort(scores, kind='mergesort')[-selector.k:]] = True
            selected_features = np.asarray(features)[training_keep][mask]
            used_features = selected_features[np.asarray(estimator.coef_) != 0]
            coverage_path = output / (clock_name + '-feature-coverage.json')
            coverage = {'feature_ids': features, 'training_retained_ids': np.asarray(features)[training_keep].tolist(),
                        'training_dropped_ids': np.asarray(features)[~training_keep].tolist(),
                        'selected_ids': selected_features.tolist(), 'nonzero_coefficient_ids': used_features.tolist(),
                        'imputation_strategy': imputer.strategy,
                        'missing_ids_by_sample': {str(sid): np.asarray(features)[~observed_values[index]].tolist()
                                                  for index, sid in enumerate(frame.index)},
                        'imputed_ids_by_sample': {str(sid): np.asarray(features)[~observed_values[index] & training_keep].tolist()
                                                  for index, sid in enumerate(frame.index)},
                        'observed_selected_by_sample': np.isfinite(frame[selected_features]).sum(axis=1).astype(int).to_dict(),
                        'observed_nonzero_by_sample': np.isfinite(frame[used_features]).sum(axis=1).astype(int).to_dict()}
            write_new(coverage_path, coverage)
            checks[-1].update(measurement_sources={**declarations,
                'model': {'path': str(ASSETS / name), 'sha256': model_digest},
                'feature_coverage': {'path': str(coverage_path.resolve()), 'sha256': sha(coverage_path)}},
                score_units=('log10 mortality hazard-ratio score; not years' if species == 'mouse' else 'CANNOT VERIFY'),
                normalization_profile=preprocessed['normalization_profile'])
        if sha(ASSETS / name) != model_digest:
            raise PipelineIntegrityError('RNA model changed during prediction')
    path = output / 'clock_predictions.csv'
    predictions.rename_axis('sample_id').to_csv(path, mode='x')
    result = {'path': str(path), 'clocks': list(predictions.columns), 'sample_count': len(predictions), 'checks': checks,
              'clock_predictions_recomputed': True, 'scientific_status': 'CANNOT VERIFY',
              'implementation_sha256': implementation_sources,
              'prediction_implementation_sources': snapshots,
              'executed_implementation_sha256': executed_implementations,
              'prediction_sha256': sha(path)}
    write_new(output / 'predictions.json', result)
    for check in checks:
        for declaration in check.get('measurement_sources', {}).values():
            if sha(declaration['path']) != declaration['sha256']:
                raise PipelineIntegrityError('RNA measurement input changed during prediction')
    for path, digest in {**implementation_sources, **executed_implementations}.items():
        if sha(path) != digest:
            raise PipelineIntegrityError('RNA prediction implementation changed during execution')
    return result


def _context_cli(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description='Build source-only GEO/Methods context without model or clock calls')
    parser.add_argument('--study', required=True)
    parser.add_argument('--series', type=Path, required=True, help='Existing GEO Series SOFT source')
    parser.add_argument('--output', type=Path, required=True, help='New context directory')
    parser.add_argument('--no-publication-methods', action='store_true')
    args = parser.parse_args(argv)
    try:
        if args.output.exists():
            raise FileExistsError('Use a new context directory; existing sources are preserved')
        result = acquire_study_context(args.study, args.series, args.output,
                                       publication_methods=not args.no_publication_methods)
        print(json.dumps({'study_context': str(result), 'model_calls': 0, 'clock_runs': 0}))
    except Exception as error:
        raise SystemExit(f'Source context stopped: {type(error).__name__}: {error}')


# Added to the existing input-adapter module; imports are supplied by that module.
DNAM_IDAT_PROFILE = {
    'name': 'minfi_ssnoob_beta100_detection01_v1',
    'array': 'IlluminaHumanMethylationEPIC', 'platform': 'GPL23976',
    'import': {'extended': True, 'force': False},
    'noob': {'offset': 15, 'dyeCorr': True, 'dyeMethod': 'single'},
    'beta': {'type': 'Illumina', 'offset': 100, 'betaThreshold': 0, 'minZero': True},
    'qc': {'detection_type': 'm+u', 'failed_if_greater_than': 0.01,
           'mask_nonfinite': True, 'failed_cells': 'NA', 'sample_exclusion': 'none',
           'bead_counts': 'both channels retained; no bead-count filter'},
    'averaging': False, 'imputation': False, 'author_method_reproduction': False,
    'scientific_status': 'CANNOT VERIFY',
}
DNAM_IDAT_SUPPORTED_PROFILE = {
    **DNAM_IDAT_PROFILE,
    'name': 'minfi_ssnoob_beta100_detection01_supported_v2',
    'probe_membership': {
        'raw': 'manifest probes with every required address observed',
        'corrected': 'address-supported probes in matching installed annotation correction groups',
        'normalization_input': 'all observed RG addresses',
        'controls': 'all manifest control addresses required',
        'exclusions': 'address_absent and annotation_absent recorded separately',
    },
}
DNAM_IDAT_PROFILES = {value['name']: value for value in (DNAM_IDAT_PROFILE, DNAM_IDAT_SUPPORTED_PROFILE)}
DNAM_IDAT_LIMITS = {'compressed_bytes_per_file': 128_000_000,
                    'decompressed_bytes_per_file': 256_000_000,
                    'compressed_bytes_total': 2_000_000_000,
                    'processing_seconds': 3600, 'samples': 64}


def _dnam_idat_limits(limits):
    if (not isinstance(limits, dict) or set(limits) != set(DNAM_IDAT_LIMITS)
            or any(type(limits[key]) is not int or not 0 < limits[key] <= ceiling
                   for key, ceiling in DNAM_IDAT_LIMITS.items())):
        raise ValueError('Invalid IDAT resource limits')
    return dict(limits)


def _idat_output_size(output, maximum=None):
    """Bound retained producer files; failed writes stay available for review."""
    if maximum is not None and (type(maximum) is not int or maximum <= 0):
        raise ValueError('IDAT output byte cap must be a positive integer')
    size = sum(path.stat().st_size for path in Path(output).rglob('*') if path.is_file())
    if maximum is not None and size > maximum:
        raise ValueError('IDAT retained output exceeds the explicit byte cap')
    return size


def _idat_running_group_members(group):
    """Exclude dead/zombie members when Linux retains orphan process records."""
    proc = Path('/proc')
    if not (proc / 'self/stat').is_file():
        # Other POSIX systems can only prove cleanup once the group is absent.
        try: os.killpg(group, 0)
        except ProcessLookupError: return []
        return [group]
    members = []
    for entry in proc.iterdir():
        if not entry.name.isdigit(): continue
        try:
            # comm may contain spaces and parentheses; fields after its final
            # closing parenthesis start at state, ppid, pgrp.
            fields = (entry / 'stat').read_text().rsplit(')', 1)[1].split()
        except (FileNotFoundError, ProcessLookupError):
            continue
        # A zombie leader may retain live threads. Only a single-threaded
        # dead record proves that this process cannot keep writing.
        if int(fields[2]) == group and (fields[0] not in {'Z', 'X', 'x'} or int(fields[17]) > 1):
            members.append(int(entry.name))
    return members


def _run_idat_worker(command, stdout, stderr, timeout, output, max_output_bytes=None):
    if max_output_bytes is None:
        return subprocess.run(command, stdout=stdout, stderr=stderr, timeout=timeout, check=False)
    _idat_output_size(output, max_output_bytes)
    started = time.monotonic()
    process = subprocess.Popen(command, stdout=stdout, stderr=stderr, start_new_session=True)
    try:
        while process.poll() is None:
            _idat_output_size(output, max_output_bytes)
            if time.monotonic() - started >= timeout:
                raise TimeoutError('IDAT worker exceeded its frozen processing time limit')
            time.sleep(.1)
    finally:
        # The owned process group can outlive its leader. Always clean the
        # group, even when poll()/wait() already reaped the direct worker.
        try: os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError: pass
        try: process.wait(timeout=2)
        except subprocess.TimeoutExpired: pass
        finally:
            try: os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError: pass
            process.wait(timeout=2)
            # SIGKILL is asynchronous. Waiting for the leader does not wait for
            # its children; require the owned group to stop before returning.
            deadline = time.monotonic() + 2
            while _idat_running_group_members(process.pid):
                if time.monotonic() >= deadline:
                    raise TimeoutError('IDAT worker process group did not stop after SIGKILL')
                time.sleep(.01)
    _idat_output_size(output, max_output_bytes)
    return subprocess.CompletedProcess(command, process.returncode)


def dnam_idat_source_pairs(study, series, sample_source):
    """Extract literal per-GSM channel ownership, without inferring replicates."""
    if not re.fullmatch(r'GSE\d+', study):
        raise ValueError('IDAT input requires a GEO Series accession')
    series_text, sample_text = Path(series).read_text(), Path(sample_source).read_text()
    if re.findall(r'^\^SERIES = (GSE\d+)\s*$', series_text, re.M) != [study]:
        raise ValueError('IDAT Series source does not identify the requested study')
    declared = re.findall(r'^!Series_sample_id = (GSM\d+)\s*$', series_text, re.M)
    if not declared or len(set(declared)) != len(declared):
        raise ValueError('IDAT Series membership is absent or duplicated')
    records, current = {}, None
    for line in sample_text.splitlines():
        key, sep, value = line.partition(' = ')
        if key == '^SAMPLE' and sep:
            if value in records or not re.fullmatch(r'GSM\d+', value):
                raise ValueError('IDAT sample source has invalid or duplicate identities')
            current = value; records[current] = {'urls': [], 'platforms': []}
        elif current and sep and key.startswith('!Sample_supplementary_file'):
            if re.search(r'_(?:Grn|Red)\.idat\.gz$', value): records[current]['urls'].append(value)
        elif current and key == '!Sample_platform_id' and sep:
            records[current]['platforms'].append(value)
    if set(records) != set(declared):
        raise ValueError('IDAT physical Series/Sample membership differs')
    pairs = []
    for gsm in declared:
        record, channels = records[gsm], {}
        if record['platforms'] != [DNAM_IDAT_PROFILE['platform']]:
            raise ValueError('IDAT profile currently supports source-declared GPL23976 EPIC only')
        for url in record['urls']:
            parsed = urllib.parse.urlsplit(url)
            name = posixpath.basename(parsed.path)
            match = re.fullmatch(r'('+gsm+r'_[A-Za-z0-9_.+-]+)_(Grn|Red)\.idat\.gz', name)
            prefix = f'/geo/samples/{gsm[:3]}{gsm[3:-3]}nnn/{gsm}/suppl/'
            if (parsed.scheme not in {'ftp', 'https'} or parsed.netloc != 'ftp.ncbi.nlm.nih.gov'
                    or parsed.query or parsed.fragment or not parsed.path.startswith(prefix)
                    or parsed.path != prefix + name or not match):
                raise ValueError('IDAT supplementary URL lacks exact source sample ownership')
            channel = match.group(2)
            if channel in channels: raise ValueError('Duplicate source IDAT channel')
            channels[channel] = {'url': url, 'https_url': 'https://'+parsed.netloc+parsed.path,
                                 'filename': name, 'basename': match.group(1)}
        if set(channels) != {'Grn', 'Red'} or channels['Grn']['basename'] != channels['Red']['basename']:
            raise ValueError('Each source sample requires one explicit matching Grn/Red pair')
        pairs.append({'sample_id': gsm, 'channels': channels})
    return pairs


def prepare_dnam_idat_plan(study, series, sample_source, output, *, profile_name=DNAM_IDAT_SUPPORTED_PROFILE['name'], limits=None):
    """Freeze source ownership, method and finite resource bounds before decoding."""
    output = Path(output)
    if output.exists(): raise FileExistsError('IDAT plan already exists: '+str(output))
    if profile_name not in DNAM_IDAT_PROFILES: raise ValueError('Unknown IDAT processing profile')
    limits = _dnam_idat_limits(DNAM_IDAT_LIMITS if limits is None else limits)
    pairs = dnam_idat_source_pairs(study, series, sample_source)
    def identity(path):
        path = Path(path).resolve()
        return {'path': str(path), 'sha256': sha(path), 'size_bytes': path.stat().st_size}
    result = {'schema': 'clockbase-dnam-idat-plan-v1', 'study_id': study,
              'created_unix': time.time(), 'profile': DNAM_IDAT_PROFILES[profile_name], 'pairs': pairs,
              'source_files': {'series': identity(series), 'samples': identity(sample_source)},
              'limits': limits,
              'method_reference': 'https://bioconductor.org/packages/release/bioc/manuals/minfi/man/minfi.pdf',
              'method_reference_scope': 'preprocessNoob and detectionP documentation; detection failures >0.01',
              'new_analysis_profile': True, 'model_calls': 0,
              'biological_units_inferred': False, 'scientific_status': 'CANNOT VERIFY'}
    if len(pairs) > result['limits']['samples']: raise ValueError('IDAT source exceeds bounded sample count')
    write_new(output, result)
    return result


def _checked_idat_plan(plan_path):
    plan = json.loads(Path(plan_path).read_text())
    if (plan.get('schema') != 'clockbase-dnam-idat-plan-v1'
            or plan.get('profile') not in DNAM_IDAT_PROFILES.values()):
        raise ValueError('Unknown or changed IDAT processing profile')
    for record in plan['source_files'].values():
        if not Path(record['path']).is_file():
            raise FileNotFoundError('Cannot proceed - required data file not found: '+record['path'])
        if sha(record['path']) != record['sha256']: raise ValueError('IDAT plan source changed')
    sources = plan['source_files']
    if plan['pairs'] != dnam_idat_source_pairs(plan['study_id'], sources['series']['path'], sources['samples']['path']):
        raise ValueError('IDAT plan pairs differ from literal physical sources')
    _dnam_idat_limits(plan['limits'])
    if not 0 < len(plan['pairs']) <= plan['limits']['samples']:
        raise ValueError('IDAT sample bound is invalid')
    return plan


def acquire_dnam_idat(plan_path, output, *, reuse_acquisition=None, fetch_callback=None, max_output_bytes=None):
    """Resume only completed byte-verified files; preserve failed download attempts."""
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    plan = _checked_idat_plan(plan_path)
    if fetch_callback is not None and not callable(fetch_callback):
        raise ValueError('IDAT fetch callback must be callable')
    download = fetch if fetch_callback is None else fetch_callback
    output.mkdir(parents=True, exist_ok=True)
    _idat_output_size(output, max_output_bytes)
    receipt = output / 'acquisition.json'
    if receipt.exists():
        value = json.loads(receipt.read_text())
        if (value.get('schema')!='clockbase-dnam-idat-acquisition-v1' or value['plan_sha256'] != sha(plan_path)
                or value.get('profile')!=plan['profile'] or value.get('study_id')!=plan['study_id']):
            raise ValueError('IDAT acquisition plan changed')
        required={kind+':'+p['sample_id']+':'+channel for p in plan['pairs']
                  for channel in ('Grn','Red') for kind in ('compressed','source_receipt','raw','decompression')}
        if set(value.get('files',{}))!=required or len(value.get('pairs',[]))!=len(plan['pairs']):
            raise ValueError('IDAT acquisition does not retain all required files')
        for frozen,stored in zip(plan['pairs'],value['pairs']):
            if frozen['sample_id']!=stored.get('sample_id') or set(stored.get('channels',{}))!={'Grn','Red'}:
                raise ValueError('IDAT acquisition sample membership changed')
            for channel,decl in frozen['channels'].items():
                record=stored['channels'][channel]
                if any(record.get(k)!=v for k,v in decl.items()): raise ValueError('IDAT acquisition declaration changed')
                raw=value['files']['raw:'+frozen['sample_id']+':'+channel]
                if any(record.get(k)!=raw.get(k) for k in ('path','sha256')):
                    raise ValueError('IDAT acquisition raw file binding changed')
        total = 0
        for role, record in value['files'].items():
            if sha(record['path']) != record['sha256']: raise ValueError('Completed IDAT input changed')
            size = Path(record['path']).stat().st_size
            if role.startswith('compressed:'):
                if not 0 < size <= plan['limits']['compressed_bytes_per_file']:
                    raise ValueError('Completed compressed IDAT exceeds the plan byte bound')
                total += size
            if role.startswith('raw:') and not 0 < size <= plan['limits']['decompressed_bytes_per_file']:
                raise ValueError('Completed raw IDAT exceeds the plan byte bound')
        if total > plan['limits']['compressed_bytes_total'] or total != value.get('compressed_bytes'):
            raise ValueError('Completed IDAT aggregate bytes differ or exceed the plan bound')
        if reuse_acquisition is not None and value.get('reused_acquisition_sha256') != sha(reuse_acquisition):
            raise ValueError('IDAT resumed acquisition reuse binding changed')
        _idat_output_size(output, max_output_bytes)
        return value
    if reuse_acquisition is not None:
        prior_path = Path(reuse_acquisition).resolve()
        if prior_path.name != 'acquisition.json' or prior_path == receipt:
            raise ValueError('IDAT reuse requires a distinct completed acquisition receipt')
        prior_hash = sha(prior_path)
        prior = json.loads(prior_path.read_text())
        prior_plan_path = Path(prior['plan_path'])
        prior_plan = _checked_idat_plan(prior_plan_path)
        if (prior_plan['study_id'] != plan['study_id'] or prior_plan['pairs'] != plan['pairs']
                or prior_plan['source_files'] != plan['source_files']):
            raise ValueError('IDAT reusable acquisition sources differ from the new plan')
        prior = acquire_dnam_idat(prior_plan_path, prior_path.parent, fetch_callback=fetch_callback)
        total = 0
        for role, record in prior['files'].items():
            size = Path(record['path']).stat().st_size
            if role.startswith('compressed:'):
                if size > plan['limits']['compressed_bytes_per_file']:
                    raise ValueError('Reused compressed IDAT exceeds the new plan bound')
                total += size
            if role.startswith('raw:') and size > plan['limits']['decompressed_bytes_per_file']:
                raise ValueError('Reused raw IDAT exceeds the new plan bound')
        if total > plan['limits']['compressed_bytes_total'] or sha(prior_path) != prior_hash:
            raise ValueError('IDAT reuse aggregate bound or source receipt changed')
        value = {'schema':'clockbase-dnam-idat-acquisition-v1','study_id':plan['study_id'],
                 'plan_path':str(plan_path),'plan_sha256':sha(plan_path),'profile':plan['profile'],
                 'files':prior['files'],'pairs':prior['pairs'],'compressed_bytes':total,
                 'recorded_unix':time.time(),'new_clock_runs':0,'model_calls':0,
                 'scientific_status':'CANNOT VERIFY','new_downloads':0,
                 'reused_acquisition_sha256':prior_hash}
        _checked_idat_plan(plan_path)
        write_new(receipt,value)
        _idat_output_size(output, max_output_bytes)
        return value
    files, pairs, total = {}, [], 0
    for pair in plan['pairs']:
        item = {'sample_id': pair['sample_id'], 'channels': {}}
        for channel, declaration in pair['channels'].items():
            role = pair['sample_id']+':'+channel
            folder = output / pair['sample_id'] / channel
            folder.mkdir(parents=True, exist_ok=True)
            completed = sorted(folder.glob('attempt-*/'+declaration['filename']))
            if len(completed) > 1: raise ValueError('Ambiguous completed raw input attempts')
            if completed:
                compressed = completed[0]
                prior = json.loads(Path(str(compressed)+'.source.json').read_text())
                if prior.get('url') != declaration['https_url'] or prior.get('sha256') != sha(compressed):
                    raise ValueError('Cached IDAT source identity differs')
            else:
                attempt = folder / ('attempt-%03d' % (len(list(folder.glob('attempt-*')))+1))
                attempt.mkdir(exist_ok=False)
                cap = min(plan['limits']['compressed_bytes_per_file'], plan['limits']['compressed_bytes_total']-total)
                if max_output_bytes is not None:
                    cap = min(cap, max_output_bytes - _idat_output_size(output))
                if cap <= 0: raise ValueError('IDAT acquisition aggregate byte limit reached')
                destination = attempt/declaration['filename']
                compressed = Path(download(declaration['https_url'], destination, limit=cap)).resolve()
                if compressed != destination.resolve() or not 0 < compressed.stat().st_size <= cap:
                    raise ValueError('IDAT fetch callback returned an invalid or oversized file')
                _idat_output_size(output, max_output_bytes)
            if not 0 < compressed.stat().st_size <= plan['limits']['compressed_bytes_per_file']:
                raise ValueError('Cached compressed IDAT exceeds the plan byte bound')
            total += compressed.stat().st_size
            if total > plan['limits']['compressed_bytes_total']: raise ValueError('IDAT aggregate byte limit reached')
            target_dir = output / 'uncompressed'
            target_dir.mkdir(exist_ok=True)
            raw = target_dir / declaration['filename'].removesuffix('.gz')
            marker = Path(str(raw)+'.json')
            if marker.exists():
                previous = json.loads(marker.read_text())
                if previous['compressed_sha256'] != sha(compressed) or previous['sha256'] != sha(raw):
                    raise ValueError('Decompressed IDAT bytes changed')
            else:
                part = Path(str(raw)+'.partial')
                if raw.exists() or part.exists(): raise FileExistsError('Unsealed IDAT decompression remains: '+str(raw))
                size = 0
                with gzip.open(compressed,'rb') as src, part.open('xb') as dst:
                    for block in iter(lambda: src.read(1024*1024), b''):
                        size += len(block)
                        if size > plan['limits']['decompressed_bytes_per_file']:
                            raise ValueError('IDAT decompression byte limit reached')
                        dst.write(block)
                        if max_output_bytes is not None:
                            dst.flush()
                            _idat_output_size(output, max_output_bytes)
                with part.open('rb') as handle:
                    if handle.read(4) != b'IDAT': raise ValueError('Source is not an IDAT binary')
                part.rename(raw)
                write_new(marker, {'compressed_sha256':sha(compressed), 'sha256':sha(raw), 'size_bytes':size})
            if not 0 < raw.stat().st_size <= plan['limits']['decompressed_bytes_per_file']:
                raise ValueError('Cached raw IDAT exceeds the plan byte bound')
            _idat_output_size(output, max_output_bytes)
            for kind, path in [('compressed',compressed),('source_receipt',Path(str(compressed)+'.source.json')),
                               ('raw',raw),('decompression',marker)]:
                files[kind+':'+role] = {'path':str(path),'sha256':sha(path),'size_bytes':path.stat().st_size}
            item['channels'][channel] = {'path':str(raw),'sha256':sha(raw),**declaration}
        pairs.append(item)
    value = {'schema':'clockbase-dnam-idat-acquisition-v1','study_id':plan['study_id'],
             'plan_path':str(plan_path),'plan_sha256':sha(plan_path),'profile':plan['profile'],
             'files':files,'pairs':pairs,'compressed_bytes':total,'recorded_unix':time.time(),
             'new_clock_runs':0,'model_calls':0,'scientific_status':'CANNOT VERIFY'}
    _checked_idat_plan(plan_path)
    write_new(receipt,value)
    _idat_output_size(output, max_output_bytes)
    return value


def process_dnam_idat(plan_path, acquisition_path, output, *, max_output_bytes=None):
    """Compute the explicit sample-level candidate, with independent readback."""
    plan_path, acquisition_path, output = map(lambda p: Path(p).resolve(), (plan_path,acquisition_path,output))
    plan = _checked_idat_plan(plan_path)
    if acquisition_path.name != 'acquisition.json': raise ValueError('IDAT acquisition receipt name differs')
    if (output/'provenance.json').is_file():
        from audit_system.report_2026_09_10 import verify_dnam_idat_provenance
        saved=json.loads((output/'provenance.json').read_text())
        if saved['files']['protocol']['sha256']!=sha(plan_path) or saved['files']['acquisition']['sha256']!=sha(acquisition_path):
            raise ValueError('Completed IDAT processing inputs changed')
        verification=verify_dnam_idat_provenance(output/'provenance.json',matrix_path=saved['files']['matrix']['path'])
        _idat_output_size(output, max_output_bytes)
        return {'matrix':saved['files']['matrix'],'dnam_idat_provenance':{'path':str(output/'provenance.json'),
                'sha256':sha(output/'provenance.json')},'verification':verification,'resumed':True}
    if not acquisition_path.is_file():
        raise FileNotFoundError('Processing requires a completed IDAT acquisition: '+str(acquisition_path))
    def no_processing_download(*args, **kwargs):
        raise RuntimeError('IDAT processing cannot download or retry an incomplete acquisition')
    acquisition = acquire_dnam_idat(plan_path,acquisition_path.parent,fetch_callback=no_processing_download)
    if output.exists(): raise FileExistsError('IDAT processing output exists; preserve it and use a new version')
    output.mkdir(parents=True)
    worker = Path(__file__).resolve().parent / 'dnam_idat_preprocess_2026_09_24.R'
    worker_snapshot=output/'implementation-rworker.R'
    inputs_snapshot=output/'implementation-inputs.py'
    shutil.copyfile(worker,worker_snapshot)
    shutil.copyfile(Path(__file__).resolve(),inputs_snapshot)
    files = dict(plan['source_files'])
    files.update(acquisition['files'])
    for role, path in [('protocol',plan_path),('acquisition',acquisition_path),('implementation:rworker',worker_snapshot),
                       ('implementation:inputs',inputs_snapshot)]:
        files[role]={'path':str(path),'sha256':sha(path),'size_bytes':path.stat().st_size}
    config = {'acquisition':str(acquisition_path),'output':str(output),'profile':plan['profile']}
    write_new(output/'worker-config.json',config)
    started=time.time()
    try:
        with (output/'stdout.txt').open('x') as stdout, (output/'stderr.txt').open('x') as stderr:
            result=_run_idat_worker(['Rscript','--vanilla',str(worker_snapshot),str(output/'worker-config.json')],
                stdout,stderr,plan['limits']['processing_seconds'],output,max_output_bytes)
        if result.returncode: raise RuntimeError('IDAT worker failed; inspect preserved stderr, exit '+str(result.returncode))
        processing=json.loads((output/'processing.json').read_text())
        if processing.get('status')!='PASS': raise ValueError('IDAT processing did not pass compatibility checks')
        for role,name in processing['outputs'].items():
            path=(output/name).resolve()
            if path.parent!=output or not path.is_file(): raise ValueError('Invalid IDAT worker output locator')
            files[role]={'path':str(path),'sha256':sha(path),'size_bytes':path.stat().st_size}
        for role,name in [('processing','processing.json'),('worker_config','worker-config.json'),
                          ('worker_stdout','stdout.txt'),('worker_stderr','stderr.txt')]:
            path=output/name; files[role]={'path':str(path),'sha256':sha(path),'size_bytes':path.stat().st_size}
        for record in files.values():
            if sha(record['path'])!=record['sha256']: raise ValueError('IDAT transitive source changed during processing')
        provenance={'schema':'clockbase-dnam-idat-provenance-v1','study_id':plan['study_id'],
            'sample_ids':[p['sample_id'] for p in plan['pairs']],'profile':plan['profile'],'files':files,
            'elapsed_seconds':time.time()-started,'scientific_status':'CANNOT VERIFY','accepted_results':0}
        write_new(output/'provenance.json',provenance)
        from audit_system.report_2026_09_10 import verify_dnam_idat_provenance
        verification=verify_dnam_idat_provenance(output/'provenance.json',matrix_path=files['matrix']['path'])
        write_new(output/'independent-verification.json',verification)
        _idat_output_size(output, max_output_bytes)
        return {'matrix':files['matrix'],'dnam_idat_provenance':{'path':str(output/'provenance.json'),
                'sha256':sha(output/'provenance.json')},'verification':verification}
    except Exception as error:
        write_new(output/'failure.json',{'status':'FAIL','error':type(error).__name__+': '+str(error),
            'elapsed_seconds':time.time()-started,'scientific_status':'CANNOT VERIFY'})
        raise


def _idat_cli(argv=None):
    import argparse
    parser=argparse.ArgumentParser(description='Prepare, acquire and reprocess source-declared EPIC IDATs; no model calls')
    commands=parser.add_subparsers(dest='action',required=True)
    prepare=commands.add_parser('prepare-idat')
    prepare.add_argument('--study',required=True)
    prepare.add_argument('--series',type=Path,required=True)
    prepare.add_argument('--samples',type=Path,required=True)
    prepare.add_argument('--output',type=Path,required=True)
    prepare.add_argument('--profile',choices=sorted(DNAM_IDAT_PROFILES),default=DNAM_IDAT_SUPPORTED_PROFILE['name'])
    acquire=commands.add_parser('acquire-idat')
    acquire.add_argument('--plan',type=Path,required=True)
    acquire.add_argument('--output',type=Path,required=True)
    acquire.add_argument('--reuse-acquisition',type=Path)
    process=commands.add_parser('process-idat')
    process.add_argument('--plan',type=Path,required=True)
    process.add_argument('--acquisition',type=Path,required=True)
    process.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    try:
        if args.action=='prepare-idat':
            value=prepare_dnam_idat_plan(args.study,args.series,args.samples,args.output,profile_name=args.profile)
        elif args.action=='acquire-idat': value=acquire_dnam_idat(args.plan,args.output,reuse_acquisition=args.reuse_acquisition)
        else: value=process_dnam_idat(args.plan,args.acquisition,args.output)
        print(json.dumps({'status':'PASS','action':args.action,'output':str(args.output.resolve()),
            'study_id':value.get('study_id',value.get('verification',{}).get('study_id')),
            'scientific_status':'CANNOT VERIFY','model_calls':0},allow_nan=False))
    except Exception as error:
        raise SystemExit('IDAT operation stopped: '+type(error).__name__+': '+str(error))


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1] in {'prepare-idat','acquire-idat','process-idat'}:
        _idat_cli()
    else:
        _context_cli()
