"""Download an ETag-pinned ARCHS4 release and inventory sample overlap, not dates."""
import argparse
import csv
import contextlib
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.metadata
import json
import os
import platform
from pathlib import Path
import re
import shutil
import subprocess
import sys
import urllib.request


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def download(species, declaration, destination, *, request=None):
    if request is not None:
        return _bounded_download(species, declaration, destination, request)
    url = declaration['url']
    if url != f'https://s3.k8s.maayanlab.cloud/archs4/files/{species}_gene_v2.latest.h5':
        raise ValueError('Unrecognized official download URL')
    target = destination / f'{species}_gene_v2.latest.h5'
    receipt = destination / f'{species}_download.json'
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved['declaration'] != declaration or not target.is_file() or digest(target) != saved['sha256']:
            raise ValueError(f'Completed download changed: {target}')
        return saved
    if target.exists():
        raise FileExistsError(f'Unverified existing target: {target}')
    part = target.with_suffix('.h5.part')
    started = datetime.now(timezone.utc).isoformat()
    etag = declaration['http_headers']['Etag']
    command = (['aria2c', '--continue=true', '--auto-file-renaming=false', '--file-allocation=none',
                '--split=8', '--max-connection-per-server=8', '--min-split-size=16M',
                '--max-tries=3', '--retry-wait=3', '--console-log-level=warn',
                '--summary-interval=0', '--download-result=hide', '--header=If-Match: ' + etag,
                '--dir=' + str(destination), '--out=' + part.name, url]
               if shutil.which('aria2c') else
               ['curl', '--fail', '--show-error', '--silent', '--retry', '2',
                '--retry-delay', '3', '--connect-timeout', '30', '--continue-at', '-',
                '--header', 'If-Match: ' + etag, '--output', str(part), url])
    with (destination / f'{species}_curl.log').open('ab') as log:
        result = subprocess.run(command, stdout=log, stderr=log)
    if result.returncode or not part.exists() or part.stat().st_size != declaration['file_size']:
        raise RuntimeError(f'Incomplete {species} download; partial file retained for resume')
    with urllib.request.urlopen(urllib.request.Request(url, method='HEAD'), timeout=30) as response:
        if response.headers.get('Etag') != etag:
            raise RuntimeError('Remote ETag changed during download')
    import h5py
    with h5py.File(part, 'r') as h5:
        samples = h5['meta/samples/geo_accession'].shape[0]
        shape = list(h5['data/expression'].shape)
    if samples != declaration['samples'] or shape[1] != samples:
        raise RuntimeError('Downloaded HDF5 sample dimension differs from catalog')
    checksum = digest(part)
    part.rename(target)
    saved = {'declaration': declaration, 'started_at_utc': started,
             'completed_at_utc': datetime.now(timezone.utc).isoformat(),
             'path': str(target), 'size_bytes': target.stat().st_size,
             'transport': command[0], 'maximum_connections_per_file': 8 if command[0]=='aria2c' else 1,
             'sha256': checksum, 'expression_shape': shape,
             'first_public_dates_verified': False, 'source_etag_verified': True}
    receipt.write_text(json.dumps(saved, indent=2) + '\n')
    print(json.dumps({'species': species, 'download_complete': True, 'shape': shape}), flush=True)
    return saved


def _bounded_download(species, declaration, destination, request):
    """Reuse source validation with an explicitly budgeted, no-retry transport.

    The caller owns the cache lock, lifetime reservations and each HTTP receipt.
    A partial object is retained and cannot be resumed implicitly.
    """
    url = f'https://s3.k8s.maayanlab.cloud/archs4/files/{species}_gene_v2.latest.h5'
    if species not in {'human', 'mouse'} or declaration.get('url') != url:
        raise ValueError('Unrecognized official ARCHS4 source')
    expected_size, expected_samples = declaration.get('file_size'), declaration.get('samples')
    if any(type(v) is not int or v <= 0 for v in (expected_size, expected_samples)):
        raise ValueError('ARCHS4 catalog size/sample declaration is invalid')
    etag = declaration['http_headers']['Etag']
    if not isinstance(etag, str) or not etag or etag.startswith('W/'):
        raise ValueError('ARCHS4 requires a nonempty strong ETag')
    destination = Path(destination); destination.mkdir(parents=True, exist_ok=True)
    target = destination / f'{species}_gene_v2.latest.h5'; receipt = destination / f'{species}_download.json'
    if receipt.exists():
        if digest(receipt) != Path(str(receipt)+'.sha256').read_text().strip():
            raise ValueError('Cached ARCHS4 receipt changed')
        saved = json.loads(receipt.read_text())
        if (saved.get('schema') != 'clockbase-bounded-archs4-download-v1' or saved['declaration'] != declaration
                or saved['path'] != str(target.resolve()) or not target.is_file()
                or target.stat().st_size != expected_size or digest(target) != saved['sha256']):
            raise ValueError('Cached ARCHS4 object changed')
        return saved
    part = target.with_suffix('.h5.part')
    if target.exists() or part.exists():
        raise FileExistsError('Unverified/partial ARCHS4 object retained; no implicit retry')
    began = datetime.now(timezone.utc).isoformat()
    response = request('GET', url, part, expected_size, {'If-Match': etag})
    if response['headers'].get('etag') != etag or part.stat().st_size != expected_size:
        raise ValueError('Downloaded ARCHS4 ETag or size differs from selection')
    after = request('HEAD', url, None, 0, {})
    if (after['headers'].get('etag') != etag
            or after['headers'].get('content-length') != str(expected_size)):
        raise ValueError('ARCHS4 source changed during bounded download')
    import h5py
    with h5py.File(part, 'r') as handle:
        samples = handle['meta/samples/geo_accession'].shape[0]
        shape = list(handle['data/expression'].shape)
    if len(shape) != 2 or not all(shape) or samples != expected_samples or shape[1] != samples:
        raise ValueError('ARCHS4 HDF5 dimensions differ from catalog')
    checksum = digest(part); part.rename(target)
    saved = {'schema': 'clockbase-bounded-archs4-download-v1', 'declaration': declaration,
             'started_at_utc': began, 'completed_at_utc': datetime.now(timezone.utc).isoformat(),
             'path': str(target.resolve()), 'size_bytes': expected_size, 'sha256': checksum,
             'expression_shape': shape, 'maximum_connections_per_file': 1, 'automatic_retries': 0,
             'source_etag_verified': True, 'first_public_dates_verified': False,
             'response_receipt': response['receipt'], 'post_download_head_receipt': after['receipt']}
    with receipt.open('x') as handle: json.dump(saved, handle, sort_keys=True, indent=2); handle.write('\n')
    Path(str(receipt)+'.sha256').write_text(digest(receipt)+'\n')
    return saved


def sample_inventory(path):
    import h5py
    with h5py.File(path, 'r') as h5:
        decode = lambda ds: [v.decode() if isinstance(v, bytes) else str(v) for v in ds[:]]
        samples = decode(h5['meta/samples/geo_accession'])
        series = decode(h5['meta/samples/series_id'])
    if len(samples) != len(series) or len(set(samples)) != len(samples):
        raise ValueError(f'Duplicate sample identifiers or mismatched metadata: {path}')
    # Preserve all study memberships; the historical extractor kept only the first.
    return {sample: sorted({g.strip() for g in memberships.split(',') if g.strip()})
            for sample, memberships in zip(samples, series)}


HARDENED_INTAKE_SCHEMA = 'clockbase-new-rna-intake-plan-v2'
HARDENED_INTAKE_SCHEMA_V21 = 'clockbase-new-rna-intake-plan-v2.1'
BATCH_EXTRACTION_SCHEMA = 'clockbase-archs4-batch-extraction-plan-v1'
BATCH_SOURCE_PROFILE_SCHEMA = 'clockbase-rna-source-profile-census-v1'
BATCH_STUDY_EXTRACTION_SCHEMA = 'clockbase-archs4-batch-study-extraction-v1'
BATCH_SUMMARY_SCHEMA = 'clockbase-archs4-batch-extraction-summary-v1'
BATCH_IMPLEMENTATIONS = {
    'batch_coordinator': Path(__file__).resolve(),
    'agent_workflow_helpers': Path(__file__).resolve().parent / 'agent_workflow_2026_09_10.py',
    'extension_eligibility': Path(__file__).resolve().parent / 'extension_eligibility_2026_09_10.py',
    'source_packet_dependency': Path(__file__).resolve().parent / 'source_packet.py',
    'pipeline_inputs': Path(__file__).resolve().parents[1] / 'src/pipeline_inputs_2026_09_10.py',
    'pipeline_runtime': Path(__file__).resolve().parents[1] / 'src/pipeline_runtime_2026_09_10.py',
    'metadata_adapter': Path(__file__).resolve().parents[1] / 'src/update_pipeline_2026_09_20.py',
}
HARDENED_V21_NUMERIC_CHARS = 64
HARDENED_V21_IDENTIFIER_CHARS = 128
HARDENED_V21_IMPLEMENTATIONS = (
    'archs4_update_2026_09_10.py', 'pipeline_inputs_2026_09_10.py',
    'pipeline_runtime_2026_09_10.py', 'agent_workflow_2026_09_10.py',
    'source_packet.py', 'clock_replay_2026_09_10.py', '2_run_rna_clock.py')
HARDENED_V21_DEPENDENCIES = (
    'FUN.Preprocessing_for_clock_3.0.R', 'FUN.Map_to_orthologs.R',
    'Gene_list_all_4.6.txt', 'Table_of_orthologs.csv')
HARDENED_V21_MODELS = {
    'human': ('EN_Hazard2_All_All_All_scaleddiff.pkl',),
    'mouse': ('EN_Hazard_Mouse_All_All_scaleddiff.pkl',
              'EN_Hazard_Rodents_All_All_scaleddiff.pkl'),
}
RNA_MODALITY_BY_SPECIES = {'human': 'human_rna', 'mouse': 'mouse_rna'}
RNA_ORGANISM_BY_SPECIES = {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}


def _intake_write_new(path, value):
    """Create one immutable JSON receipt without rewriting a prior stage."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def _intake_json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _intake_file_declaration(declaration, label):
    if not isinstance(declaration, dict) or set(declaration) - {'path', 'sha256', 'size_bytes', 'kind', 'expression_shape', 'sample_count'}:
        raise ValueError(f'Invalid file declaration: {label}')
    path = Path(declaration.get('path', '')).expanduser().resolve()
    expected = declaration.get('sha256')
    if not path.is_file() or not isinstance(expected, str) or len(expected) != 64:
        raise ValueError(f'Missing or malformed file declaration: {label}')
    size = path.stat().st_size
    if declaration.get('size_bytes') is not None and declaration['size_bytes'] != size:
        raise ValueError(f'File size changed before hardened intake: {label}')
    actual = digest(path)
    if actual != expected:
        raise ValueError(f'File SHA256 changed before hardened intake: {label}')
    return {'label': label, 'path': str(path), 'size_bytes': size, 'sha256': actual}


def _intake_runtime_versions(plan):
    runtime = plan.get('runtime')
    if not isinstance(runtime, dict):
        raise ValueError('Hardened intake requires a runtime declaration')
    actual = {'python': platform.python_version(), 'packages': {}}
    packages = runtime.get('packages')
    if not isinstance(packages, dict) or not packages:
        raise ValueError('Hardened intake requires pinned Python packages')
    for name, expected in packages.items():
        if not isinstance(name, str) or not isinstance(expected, str):
            raise ValueError('Runtime package declarations must be strings')
        try:
            actual['packages'][name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError as error:
            raise ValueError(f'Pinned runtime package is unavailable: {name}') from error
        if actual['packages'][name] != expected:
            raise ValueError(f'Runtime package version differs: {name}')
    if runtime.get('python') != actual['python']:
        raise ValueError('Python runtime version differs from hardened plan')
    if runtime.get('Rscript'):
        result = subprocess.run(['Rscript', '--version'], capture_output=True, text=True, check=False)
        observed = (result.stdout + result.stderr).strip()
        if observed != runtime['Rscript']:
            raise ValueError('Rscript runtime version differs from hardened plan')
        actual['Rscript'] = observed
    if runtime.get('edgeR'):
        if not runtime.get('Rscript'):
            raise ValueError('edgeR pin requires an Rscript pin')
        result = subprocess.run(['Rscript', '-e', "cat(as.character(packageVersion('edgeR')))"],
                                capture_output=True, text=True, check=False)
        observed = result.stdout.strip()
        if result.returncode or observed != runtime['edgeR']:
            raise ValueError('edgeR runtime version differs from hardened plan')
        actual['edgeR'] = observed
    return actual


def _intake_read_metadata(path, study, species):
    try:
        payload = Path(path).read_bytes()
        records = list(csv.reader(payload.decode('utf-8-sig').splitlines()))
    except (OSError, UnicodeError, csv.Error) as error:
        raise ValueError(f'Cannot read hardened metadata: {study}') from error
    if len(records) < 2 or not records[0] or any(not key.strip() for key in records[0]):
        raise ValueError(f'Empty or malformed hardened metadata: {study}')
    header = records[0]
    if len(set(header)) != len(header):
        raise ValueError(f'Duplicate hardened metadata columns: {study}')
    if any(len(row) != len(header) for row in records[1:]):
        raise ValueError(f'Ragged hardened metadata rows: {study}')
    rows = [dict(zip(header, row)) for row in records[1:]]
    id_key = next((key for key in ('GSM', 'gsm_id', 'sample_id', 'geo_accession') if key in header), None)
    if not id_key:
        raise ValueError(f'No GEO sample identity column: {study}')
    ids = [row[id_key] for row in rows]
    if not ids or len(set(ids)) != len(ids) or any(not value.strip() for value in ids):
        raise ValueError(f'Missing or duplicate GEO sample identities: {study}')
    expected_organism = {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}.get(species)
    if expected_organism is None:
        raise ValueError(f'Unknown hardened species: {species}')
    if {row.get('organism_ch1', '') for row in rows} != {expected_organism}:
        raise ValueError(f'Metadata organism differs from plan: {study}')
    if {row.get('molecule_ch1', '').casefold() for row in rows} != {'total rna'}:
        raise ValueError(f'Metadata molecule differs from plan: {study}')
    return rows, ids, id_key


def _intake_release_preflight(release, study, species, source_ids, limits):
    """Verify source membership and bound CSV size before expression extraction."""
    import h5py
    import numpy as np
    path = Path(release['path'])
    expected_shape = release.get('expression_shape')
    expected_samples = release.get('sample_count')
    with h5py.File(path, 'r') as h5:
        shape = list(h5['data/expression'].shape)
        if expected_shape is not None and shape != expected_shape:
            raise ValueError(f'ARCHS4 expression shape differs: {species}')
        if expected_samples is not None and int(h5['meta/samples/geo_accession'].shape[0]) != expected_samples:
            raise ValueError(f'ARCHS4 sample count differs: {species}')
        decode = lambda values: [v.decode() if isinstance(v, bytes) else str(v) for v in values]
        geo_ids = decode(h5['meta/samples/geo_accession'][:])
        memberships = decode(h5['meta/samples/series_id'][:])
        if len(geo_ids) != len(set(geo_ids)) or len(geo_ids) != len(memberships):
            raise ValueError(f'ARCHS4 sample identities are not unique/aligned: {species}')
        indices = [i for i, line in enumerate(memberships)
                   if study in {part.strip() for part in line.split(',')}]
        selected = [geo_ids[i] for i in indices]
        if set(selected) != set(source_ids) or len(selected) != len(source_ids):
            raise ValueError(f'ARCHS4/GEO full roster differs; no partial extraction: {study}')
        expected_organism = {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}[species]
        organisms = decode(h5['meta/samples/organism_ch1'][indices])
        if set(organisms) != {expected_organism}:
            raise ValueError(f'ARCHS4 organism differs from plan: {study}')
        gene_count = int(h5['data/expression'].shape[0])
        dtype = np.dtype(h5['data/expression'].dtype)
        if dtype.kind not in 'biufc' or dtype.itemsize > 16:
            raise ValueError(f'Unsupported ARCHS4 expression dtype for conservative bound: {dtype}')
        gene_ids = decode(h5['meta/genes/ensembl_gene'][:])
        if len(gene_ids) != gene_count or len(set(gene_ids)) != gene_count:
            raise ValueError(f'ARCHS4 gene identities are not unique/aligned: {species}')
    numeric_chars = int(limits.get('max_numeric_chars', 64))
    identifier_chars = int(limits.get('max_identifier_chars', 128))
    safe_bound = (gene_count * len(source_ids) * (numeric_chars + 2)
                  + gene_count * (identifier_chars + 1)
                  + len(source_ids) * (identifier_chars + 1)
                  + 1024)
    cap = int(limits['max_extracted_bytes'])
    if safe_bound > cap:
        raise ValueError(f'Conservative extraction bound exceeds cap before extraction: {study}')
    return {'study_id': study, 'species': species, 'sample_count': len(source_ids),
            'gene_count': gene_count, 'expression_dtype': str(dtype),
            'expression_shape': shape, 'selected_sample_ids': selected,
            'safe_upper_bound_bytes': safe_bound, 'max_extracted_bytes': cap,
            'preflight_status': 'PASS'}


def _intake_terminal_stage(store, name, function):
    """Run one immutable terminal stage; failures become non-retryable receipts."""
    def execute(directory):
        try:
            value = function(directory)
            envelope = {'status': 'PASS', 'phase': name, 'value': value, 'automatic_retry': False}
        except Exception as error:  # noqa: BLE001 - receipt must preserve every failure class
            envelope = {'status': 'CANNOT VERIFY', 'phase': name,
                        'error': f'{type(error).__name__}: {error}', 'automatic_retry': False}
        _intake_write_new(directory / 'terminal-attempt.json', envelope)
        return envelope
    return store.run(name, execute)


def _intake_verify_frozen_files(plan):
    frozen = plan.get('frozen_files')
    if not isinstance(frozen, list) or not frozen:
        raise ValueError('Hardened intake requires implementation, dependency and model asset declarations')
    records = []
    kinds = set()
    for index, declaration in enumerate(frozen):
        kind = declaration.get('kind') if isinstance(declaration, dict) else None
        if kind not in {'implementation', 'dependency_asset', 'model_asset'}:
            raise ValueError(f'Unsupported frozen file kind at index {index}')
        kinds.add(kind)
        records.append({**_intake_file_declaration(declaration, f'frozen_files[{index}]'), 'kind': kind})
    if kinds != {'implementation', 'dependency_asset', 'model_asset'}:
        raise ValueError('Hardened intake requires implementation, dependency and model hashes')
    return records


def _intake_v21_expected_assets(species):
    expected = {(name, 'implementation') for name in HARDENED_V21_IMPLEMENTATIONS}
    expected.update((name, 'dependency_asset') for name in HARDENED_V21_DEPENDENCIES)
    expected.add((f'Gene_table_{species}.csv', 'dependency_asset'))
    expected.update((name, 'model_asset') for name in HARDENED_V21_MODELS[species])
    return expected


def _intake_v21_verify_frozen_files(plan):
    records = _intake_verify_frozen_files(plan)
    studies = plan.get('studies') or []
    species = {study.get('species') for study in studies if isinstance(study, dict)}
    if not species.issubset(HARDENED_V21_MODELS):
        raise ValueError('V2.1 requires human or mouse study species')
    expected = set().union(*(_intake_v21_expected_assets(item) for item in species))
    observed = {}
    for record in records:
        name = Path(record['path']).name
        if name in observed:
            raise ValueError(f'V2.1 frozen asset is declared more than once: {name}')
        observed[name] = (record['kind'], record['sha256'])
    missing = sorted(name for name, kind in expected
                     if observed.get(name, (None, None))[0] != kind)
    if missing:
        raise ValueError('V2.1 frozen execution assets are incomplete: ' + ', '.join(missing))
    return records


def _intake_v21_recheck(declarations):
    current = []
    for label, declaration in declarations:
        current.append((label, _intake_file_declaration(declaration, label)))
    return current


def _intake_v21_release_preflight(release, study, species, source_ids, limits):
    """Apply fixed conservative limits and reject unsupported complex matrices."""
    import h5py
    import numpy as np

    preflight = _intake_release_preflight(
        release, study, species, source_ids,
        {**limits, 'max_numeric_chars': HARDENED_V21_NUMERIC_CHARS,
         'max_identifier_chars': HARDENED_V21_IDENTIFIER_CHARS})
    if np.dtype(preflight['expression_dtype']).kind not in 'biuf':
        raise ValueError('V2.1 rejects complex or unsupported ARCHS4 expression dtype')
    with h5py.File(release['path'], 'r') as h5:
        decode = lambda values: [v.decode() if isinstance(v, bytes) else str(v) for v in values]
        gene_ids = decode(h5['meta/genes/ensembl_gene'][:])
        selected_ids = decode(h5['meta/samples/geo_accession'][:])
    max_gene_identifier = max((len(value) for value in gene_ids), default=0)
    max_sample_identifier = max((len(value) for value in selected_ids), default=0)
    if max(max_gene_identifier, max_sample_identifier, len(study)) > HARDENED_V21_IDENTIFIER_CHARS:
        raise ValueError(f'V2.1 source identifier exceeds {HARDENED_V21_IDENTIFIER_CHARS} characters: {study}')
    safe_bound = (preflight['gene_count'] * len(source_ids) *
                  (HARDENED_V21_NUMERIC_CHARS + 2)
                  + preflight['gene_count'] * (HARDENED_V21_IDENTIFIER_CHARS + 1)
                  + len(source_ids) * (HARDENED_V21_IDENTIFIER_CHARS + 1) + 1024)
    if safe_bound > int(limits['max_extracted_bytes']):
        raise ValueError(f'V2.1 conservative extraction bound exceeds cap: {study}')
    return {**preflight, 'bound_profile': 'fixed-v2.1',
            'max_numeric_chars': HARDENED_V21_NUMERIC_CHARS,
            'max_identifier_chars': HARDENED_V21_IDENTIFIER_CHARS,
            'max_gene_identifier_length': max_gene_identifier,
            'max_sample_identifier_length': max_sample_identifier,
            'safe_upper_bound_bytes': safe_bound}


def _intake_v21_lock_path(destination):
    return destination.parent / (destination.name + '.v2.1.lock')


@contextlib.contextmanager
def _intake_v21_run_lock(destination):
    """Serialize V2.1 callers across processes with an advisory lock."""
    lock_path = _intake_v21_lock_path(Path(destination))
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    handle = lock_path.open('a+')
    try:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f'V2.1 intake is already running: {lock_path}') from error
        yield lock_path
    finally:
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        handle.close()


def validate_hardened_intake_plan_v21(plan):
    """Read-only strict preflight for plans that have not used V2 output paths."""
    if not isinstance(plan, dict) or plan.get('schema') != HARDENED_INTAKE_SCHEMA_V21:
        raise ValueError('Unsupported hardened intake V2.1 plan schema')
    limits = plan.get('limits')
    if not isinstance(limits, dict):
        raise ValueError('V2.1 intake limits are required')
    if int(limits.get('max_numeric_chars', 0)) < HARDENED_V21_NUMERIC_CHARS:
        raise ValueError('V2.1 numeric bound cannot be weakened by the caller')
    if int(limits.get('max_identifier_chars', 0)) < HARDENED_V21_IDENTIFIER_CHARS:
        raise ValueError('V2.1 identifier bound cannot be weakened by the caller')
    if limits.get('concurrency') != 1:
        raise ValueError('V2.1 requires serial execution')
    studies = plan.get('studies')
    if not isinstance(studies, list) or not studies:
        raise ValueError('V2.1 requires at least one study')
    if len(studies) > int(limits.get('max_studies', 0)):
        raise ValueError('V2.1 study count exceeds its bound')
    study_ids = [item.get('study_id') for item in studies if isinstance(item, dict)]
    if len(study_ids) != len(studies) or any(not re.fullmatch(r'GSE[0-9]+', value or '') for value in study_ids):
        raise ValueError('V2.1 study IDs must match GSE[0-9]+')
    if len(set(study_ids)) != len(study_ids):
        raise ValueError('V2.1 study IDs must be unique')
    _intake_runtime_versions(plan)
    frozen = _intake_v21_verify_frozen_files(plan)
    declarations = [('frozen_files[%d]' % index, item)
                    for index, item in enumerate(plan['frozen_files'])]
    releases = {}
    for species, declaration in (plan.get('releases') or {}).items():
        releases[species] = _intake_file_declaration(declaration, f'releases.{species}')
        releases[species].update({key: declaration[key] for key in ('expression_shape', 'sample_count')
                                  if key in declaration})
        declarations.append((f'releases.{species}', declaration))
    if not releases:
        raise ValueError('V2.1 requires at least one release declaration')
    prepared = []
    preflights = []
    for study in studies:
        study_id, species = study['study_id'], study.get('species')
        if species not in releases:
            raise ValueError(f'V2.1 study release is missing: {study_id}')
        source_files = study.get('source_files') or {}
        metadata_decl, context_decl = source_files.get('metadata'), source_files.get('context')
        metadata = _intake_file_declaration(metadata_decl, f'{study_id}.metadata')
        context = _intake_file_declaration(context_decl, f'{study_id}.context')
        declarations.extend(((f'{study_id}.metadata', metadata_decl),
                             (f'{study_id}.context', context_decl)))
        reuse_decl = study.get('reuse_extraction')
        if reuse_decl:
            _intake_file_declaration(reuse_decl, f'{study_id}.reuse_extraction')
            declarations.append((f'{study_id}.reuse_extraction', reuse_decl))
        rows, source_ids, _ = _intake_read_metadata(metadata['path'], study_id, species)
        if any(len(value) > HARDENED_V21_IDENTIFIER_CHARS for value in source_ids):
            raise ValueError(f'V2.1 metadata identifier exceeds {HARDENED_V21_IDENTIFIER_CHARS}: {study_id}')
        if len(rows) > int(limits.get('max_samples_per_study', 0)):
            raise ValueError(f'V2.1 sample cap exceeded before extraction: {study_id}')
        preflights.append(_intake_v21_release_preflight(
            releases[species], study_id, species, source_ids, limits))
        prepared.append({'study_id': study_id, 'metadata': metadata, 'context': context,
                         'sample_count': len(source_ids)})
    final_frozen = _intake_v21_verify_frozen_files(plan)
    if [(x['path'], x['sha256'], x['size_bytes']) for x in frozen] != \
            [(x['path'], x['sha256'], x['size_bytes']) for x in final_frozen]:
        raise ValueError('V2.1 frozen execution hashes changed during preflight')
    final_inputs = _intake_v21_recheck(declarations)
    return {'schema': 'clockbase-hardened-intake-preflight-v2.1',
            'status': 'PASS', 'study_count': len(prepared), 'studies': prepared,
            'preflights': preflights, 'frozen_files': final_frozen,
            'final_input_hash_recheck': 'PASS', 'concurrency_lock': 'REQUIRED',
            'fixed_bound_profile': {'max_numeric_chars': HARDENED_V21_NUMERIC_CHARS,
                                    'max_identifier_chars': HARDENED_V21_IDENTIFIER_CHARS},
            'verified_declarations': len(final_inputs)}


def run_hardened_intake_v21(plan, output=None):
    """Run a strict V2.1 plan under an exclusive lock; V2 outputs are untouched."""
    if not isinstance(plan, dict) or plan.get('schema') != HARDENED_INTAKE_SCHEMA_V21:
        raise ValueError('run_hardened_intake_v21 requires the V2.1 schema')
    destination = Path(output or plan.get('output', '')).expanduser().resolve()
    with _intake_v21_run_lock(destination):
        validate_hardened_intake_plan_v21(plan)
        compatible = dict(plan)
        compatible['schema'] = HARDENED_INTAKE_SCHEMA
        result = run_hardened_intake(compatible, output=destination)
        # Recheck all frozen/source declarations after stages, before reporting success.
        validate_hardened_intake_plan_v21(plan)
        return result


def run_hardened_intake(plan, output=None):
    """Execute a plan-driven, bounded ARCHS4 intake in a new v2 output directory.

    The plan is the caller's portable resource declaration.  Existing v1 output
    directories cannot be reused; a v2 directory may only resume its own frozen
    configuration, and a failed terminal stage is never retried automatically.
    """
    if not isinstance(plan, dict) or plan.get('schema') != HARDENED_INTAKE_SCHEMA:
        raise ValueError('Unsupported hardened intake plan schema')
    limits = plan.get('limits')
    if not isinstance(limits, dict):
        raise ValueError('Hardened intake limits are required')
    required_limits = ('max_studies', 'max_samples_per_study', 'max_extracted_bytes', 'concurrency')
    if any(key not in limits for key in required_limits):
        raise ValueError('Hardened intake limits are incomplete')
    if limits['concurrency'] != 1 or int(limits['max_studies']) <= 0 or int(limits['max_samples_per_study']) <= 0:
        raise ValueError('Hardened intake requires positive serial limits')
    studies = plan.get('studies')
    if not isinstance(studies, list) or not studies or len(studies) > int(limits['max_studies']):
        raise ValueError('Study count exceeds hardened intake bound')
    destination = Path(output or plan.get('output', '')).expanduser().resolve()
    if not str(destination) or destination == Path('/').resolve():
        raise ValueError('A dedicated hardened intake output directory is required')
    if destination.exists() and any(destination.iterdir()) and not (destination / 'config.json').is_file():
        raise ValueError('Hardened intake requires a new empty output directory; v1 output cannot be reused')

    runtime = _intake_runtime_versions(plan)
    frozen = _intake_verify_frozen_files(plan)
    releases = {}
    for species, declaration in (plan.get('releases') or {}).items():
        releases[species] = _intake_file_declaration(declaration, f'releases.{species}')
        releases[species].update({key: declaration[key] for key in ('expression_shape', 'sample_count') if key in declaration})
    if not releases:
        raise ValueError('Hardened intake requires at least one release declaration')

    prepared = []
    for study in studies:
        if not isinstance(study, dict) or not study.get('study_id') or study.get('species') not in releases:
            raise ValueError('Malformed hardened study declaration')
        study_id, species = study['study_id'], study['species']
        metadata_decl = (study.get('source_files') or {}).get('metadata')
        context_decl = (study.get('source_files') or {}).get('context')
        metadata = _intake_file_declaration(metadata_decl, f'{study_id}.metadata')
        context = _intake_file_declaration(context_decl, f'{study_id}.context')
        reuse_decl = None
        if study.get('reuse_extraction'):
            reuse_decl = _intake_file_declaration(study['reuse_extraction'], f'{study_id}.reuse_extraction')
            if reuse_decl['size_bytes'] > int(limits['max_extracted_bytes']):
                raise ValueError(f'Reused extraction exceeds byte cap before stage: {study_id}')
        rows, source_ids, id_key = _intake_read_metadata(metadata['path'], study_id, species)
        if len(rows) > int(limits['max_samples_per_study']):
            raise ValueError(f'Sample cap exceeded before extraction: {study_id}')
        prepared.append({'plan': study, 'study_id': study_id, 'species': species,
                         'metadata': metadata, 'context': context, 'rows': rows,
                         'source_ids': source_ids, 'id_key': id_key, 'reuse_decl': reuse_decl})

    plan_identity = _intake_json_hash(plan)
    config = {'schema': 'clockbase-new-rna-intake-v2', 'plan_sha256': plan_identity,
              'limits': limits, 'runtime': runtime, 'frozen_files': frozen,
              'releases': releases,
              'studies': [{key: item[key] for key in ('study_id', 'species', 'metadata', 'context', 'source_ids')}
                          for item in prepared],
              'scientific_status': 'CANNOT VERIFY', 'atlas_status': 'QUARANTINED'}
    destination.mkdir(parents=True, exist_ok=True)
    from src.pipeline_runtime_2026_09_10 import StageStore
    store = StageStore(destination, config)
    os.environ['OMP_NUM_THREADS'] = str(plan.get('omp_threads', 1))
    os.environ['OPENBLAS_NUM_THREADS'] = str(plan.get('omp_threads', 1))
    os.environ['MKL_NUM_THREADS'] = str(plan.get('omp_threads', 1))

    preflighted = []
    safe_total = 0
    for item in prepared:
        study = item['plan']; study_id = item['study_id']; species = item['species']
        release = releases[species]
        preflight = _intake_release_preflight(release, study_id, species, item['source_ids'], limits)
        safe_total += preflight['safe_upper_bound_bytes']
        preflighted.append((item, preflight))
    if safe_total > int(limits['max_extracted_bytes']):
        raise ValueError('Aggregate conservative extraction bound exceeds cap before extraction')

    results = []
    extracted_total = 0
    for item, preflight in preflighted:
        study = item['plan']; study_id = item['study_id']; species = item['species']
        release = releases[species]
        target = study.get('target', 'Bounded source design audit')
        context_path = item['context']['path']
        metadata_path = item['metadata']['path']
        study_store = StageStore(destination / 'studies' / study_id,
                                  {**config, 'study_id': study_id, 'preflight': preflight})

        def make_packet(directory):
            from audit_system.agent_workflow_2026_09_10 import source_packet as build_source_packet
            from audit_system.source_packet import compact_source_packet
            raw_packet = build_source_packet(metadata_path, study_id, target, context=context_path)
            model_packet = compact_source_packet(raw_packet)
            path = directory / 'source-packet.json'
            _intake_write_new(path, model_packet)
            return {'path': str(path), 'sha256': digest(path), 'sample_count': len(item['source_ids']),
                    'representation': 'compact_source_packet', 'raw_packet_written': False,
                    'metadata_sha256': item['metadata']['sha256'], 'context_sha256': item['context']['sha256'],
                    'model_view_keys': sorted(model_packet)}

        packet_result = _intake_terminal_stage(study_store, 'source_packet', make_packet)
        if packet_result['status'] != 'PASS':
            results.append({'study_id': study_id, 'preflight': preflight, 'source_packet': packet_result})
            continue

        reuse = study.get('reuse_extraction')
        if reuse:
            reuse_decl = item['reuse_decl']

            def extract(directory):
                target_counts = directory / 'counts.csv'
                shutil.copyfile(reuse_decl['path'], target_counts)
                columns, _ = _intake_csv_ids(target_counts)
                if set(columns) != set(item['source_ids']) or len(columns) != len(item['source_ids']):
                    raise ValueError('Reused count matrix roster differs from metadata')
                receipt = {'path': str(target_counts), 'sha256': digest(target_counts),
                           'sample_count': len(columns), 'input_mode': 'REUSED_VERIFIED_INGESTED_MATRIX',
                           'source_counts': reuse_decl, 'source_hdf5_sha256': release['sha256'],
                           'partial_study': False}
                _intake_write_new(directory / 'hardened-extraction.json', receipt)
                return receipt
        else:
            def extract(directory):
                from src.pipeline_inputs_2026_09_10 import extract_rna
                result = extract_rna(Path(release['path']), study_id, species, item['source_ids'], directory)
                counts = Path(result['path'])
                size = counts.stat().st_size
                if size > int(limits['max_extracted_bytes']):
                    raise ValueError('Extracted matrix exceeds hardened byte cap')
                receipt = {**result, 'input_mode': 'FRESH_ARCHS4_EXTRACTION',
                           'source_hdf5_sha256': release['sha256'], 'source_hdf5_size_bytes': release['size_bytes'],
                           'safe_upper_bound_bytes': preflight['safe_upper_bound_bytes'],
                           'post_extraction_bytes': size}
                _intake_write_new(directory / 'hardened-extraction.json', receipt)
                return receipt
        extraction_result = _intake_terminal_stage(study_store, 'extract', extract)
        if extraction_result['status'] != 'PASS':
            results.append({'study_id': study_id, 'preflight': preflight,
                            'source_packet': packet_result, 'extract': extraction_result})
            continue
        extraction = extraction_result['value']
        extracted_total += int(extraction.get('post_extraction_bytes', Path(extraction['path']).stat().st_size))
        if extracted_total > int(limits['max_extracted_bytes']):
            raise ValueError('Aggregate extracted byte cap exceeded')
        counts = Path(extraction['path'])
        normalization = study.get('normalization', 'clockbase_legacy')

        def preprocess(directory):
            from src.pipeline_inputs_2026_09_10 import preprocess_rna
            raw_packet = None
            if normalization == 'published_mouse_rle':
                from audit_system.agent_workflow_2026_09_10 import source_packet as build_source_packet
                raw_packet = build_source_packet(metadata_path, study_id, target, context=context_path)
            return preprocess_rna(counts, species, directory, normalization=normalization,
                                  source=raw_packet if normalization == 'published_mouse_rle' else None)
        preprocessing_result = _intake_terminal_stage(study_store, 'preprocess', preprocess)
        current = {'study_id': study_id, 'preflight': preflight, 'source_packet': packet_result,
                   'extract': extraction_result, 'preprocess': preprocessing_result}
        if preprocessing_result['status'] == 'PASS':
            preprocessed = preprocessing_result['value']

            def predict(directory):
                from src.pipeline_inputs_2026_09_10 import predict_rna
                return predict_rna(preprocessed, species, directory)
            current['predict'] = _intake_terminal_stage(study_store, 'predict', predict)
        results.append(current)
    summary = {'schema': 'clockbase-new-rna-intake-summary-v2', 'plan_sha256': plan_identity,
               'config_sha256': _intake_json_hash(config), 'runtime': runtime, 'results': results,
               'aggregate_extracted_bytes': extracted_total, 'model_api_calls': 0,
               'atlas_status': 'QUARANTINED', 'scientific_status': 'CANNOT VERIFY'}
    summary_path = destination / 'summary.json'
    if summary_path.exists():
        if json.loads(summary_path.read_text()) != summary:
            raise ValueError('Hardened summary changed; use a new output directory')
    else:
        _intake_write_new(summary_path, summary)
    return summary


def _intake_csv_ids(path):
    with Path(path).open(newline='') as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = list(reader)
    if len(header) < 2 or not rows:
        raise ValueError('Count matrix has no sample columns')
    ids = header[1:]
    if len(ids) != len(set(ids)) or any(not value for value in ids):
        raise ValueError('Count matrix sample identities are invalid')
    return ids, [row[0] for row in rows]


def _batch_seal(path, value):
    """Write one immutable JSON artifact and its detached SHA256."""
    path = Path(path)
    _intake_write_new(path, value)
    with path.with_suffix(path.suffix + '.sha256').open('x') as handle:
        handle.write(digest(path) + '\n')


def _batch_read_sealed(path):
    path = Path(path)
    sidecar = path.with_suffix(path.suffix + '.sha256')
    if not path.is_file() or not sidecar.is_file() or digest(path) != sidecar.read_text().strip():
        raise ValueError('Batch artifact changed or lacks its SHA256: ' + str(path))
    return json.loads(path.read_text())


def _batch_stat(path):
    stat = Path(path).stat()
    return {'device': stat.st_dev, 'inode': stat.st_ino, 'size_bytes': stat.st_size,
            'mtime_ns': stat.st_mtime_ns, 'ctime_ns': stat.st_ctime_ns}


def _batch_file_declaration(declaration, label):
    if (not isinstance(declaration, dict)
            or set(declaration) - {'path', 'sha256', 'size_bytes'}):
        raise ValueError('Invalid batch file declaration: ' + label)
    path = Path(declaration.get('path', '')).expanduser().resolve()
    expected = declaration.get('sha256')
    if not path.is_file() or not isinstance(expected, str) or not re.fullmatch(r'[0-9a-f]{64}', expected):
        raise ValueError('Missing or malformed batch file declaration: ' + label)
    size = path.stat().st_size
    if declaration.get('size_bytes') is not None and declaration['size_bytes'] != size:
        raise ValueError('Batch file size changed: ' + label)
    observed = digest(path)
    if observed != expected:
        raise ValueError('Batch file SHA256 changed: ' + label)
    return {'path': str(path), 'sha256': observed, 'size_bytes': size}


def _batch_release_declaration(declaration):
    allowed = {'path', 'sha256', 'size_bytes', 'expression_shape', 'sample_count', 'species'}
    if not isinstance(declaration, dict) or set(declaration) - allowed:
        raise ValueError('Invalid batch ARCHS4 release declaration')
    path = Path(declaration.get('path', '')).expanduser().resolve()
    checksum = declaration.get('sha256')
    if (not path.is_file() or not re.fullmatch(r'[0-9a-f]{64}', checksum or '')
            or declaration.get('species') not in {'human', 'mouse'}):
        raise ValueError('Missing or malformed batch ARCHS4 release declaration')
    if declaration.get('size_bytes') != path.stat().st_size:
        raise ValueError('Batch ARCHS4 release size changed before execution')
    shape = declaration.get('expression_shape')
    if (not isinstance(shape, list) or len(shape) != 2
            or any(type(value) is not int or value <= 0 for value in shape)
            or declaration.get('sample_count') != shape[1]):
        raise ValueError('Batch ARCHS4 release shape/sample declaration is invalid')
    return {**declaration, 'path': str(path)}


def _batch_source_profile(declaration, release):
    descriptor_decl = _batch_file_declaration(declaration, 'source_profile_descriptor')
    try:
        descriptor_bytes = Path(descriptor_decl['path']).read_bytes()
        if hashlib.sha256(descriptor_bytes).hexdigest() != descriptor_decl['sha256']:
            raise ValueError('Source-profile descriptor changed while it was read')
        descriptor = json.loads(descriptor_bytes)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError('Cannot read source-profile descriptor') from error
    required = {'schema', 'profile', 'modality', 'studies', 'supported_study_ids',
                'stopped_study_ids', 'counts', 'source_files'}
    if (not required <= set(descriptor)
            or descriptor.get('schema') != BATCH_SOURCE_PROFILE_SCHEMA
            or descriptor.get('profile') != 'existing_total_rna_intake_v21'
            or descriptor.get('modality') != RNA_MODALITY_BY_SPECIES[release['species']]):
        raise ValueError('Unsupported batch source-profile descriptor')
    studies = descriptor['studies']
    supported = descriptor['supported_study_ids']
    stopped = descriptor['stopped_study_ids']
    if (not isinstance(studies, list) or not studies
            or not isinstance(supported, list) or not isinstance(stopped, list)
            or len(supported) != len(set(supported)) or len(stopped) != len(set(stopped))
            or set(supported) & set(stopped)):
        raise ValueError('Source-profile partition is malformed')
    rows = {}
    all_samples = set()
    molecule_counts = {}
    resolved_sources = {}
    for row in studies:
        if not isinstance(row, dict):
            raise ValueError('Source-profile study record is malformed')
        study = row.get('study_id')
        if (not re.fullmatch(r'GSE[0-9]+', study or '') or study in rows
                or row.get('modality') != descriptor['modality']):
            raise ValueError('Source-profile study identity is invalid or duplicated')
        sample_ids = row.get('sample_ids')
        samples = row.get('samples')
        if (not isinstance(sample_ids, list) or not sample_ids
                or len(sample_ids) != len(set(sample_ids))
                or any(not re.fullmatch(r'GSM[0-9]+', value or '') for value in sample_ids)
                or row.get('sample_count') != len(sample_ids)
                or not isinstance(samples, list) or len(samples) != len(sample_ids)):
            raise ValueError('Source-profile sample inventory is invalid: ' + study)
        records = {}
        for sample in samples:
            if not isinstance(sample, dict) or sample.get('sample_id') in records:
                raise ValueError('Source-profile sample record is malformed: ' + study)
            sample_id = sample.get('sample_id')
            organism = sample.get('organism_ch1')
            molecule = sample.get('molecule_ch1')
            if (not re.fullmatch(r'GSM[0-9]+', sample_id or '')
                    or not isinstance(organism, list) or not organism
                    or not isinstance(molecule, list) or not molecule
                    or any(not isinstance(value, str) for value in organism + molecule)):
                raise ValueError('Source-profile sample annotations are invalid: ' + study)
            records[sample_id] = sample
            for value in molecule:
                molecule_counts[value] = molecule_counts.get(value, 0) + 1
        if set(records) != set(sample_ids):
            raise ValueError('Source-profile sample records differ from sample IDs: ' + study)
        overlap = all_samples & set(sample_ids)
        if overlap:
            raise ValueError('Source-profile studies reuse sample identities: ' + ', '.join(sorted(overlap)))
        all_samples.update(sample_ids)
        sources = row.get('sources')
        if not isinstance(sources, dict) or set(sources) != {'series', 'samples'}:
            raise ValueError('Source-profile study lacks exact Series/Samples sources: ' + study)
        resolved = {}
        for role in ('series', 'samples'):
            resolved[role] = _batch_file_declaration(sources[role], f'{study}.{role}')
        status = row.get('profile_status')
        if status not in {'SUPPORTED_SOURCE_PROFILE', 'UNSUPPORTED_SOURCE_PROFILE'}:
            raise ValueError('Source-profile status is invalid: ' + study)
        if status == 'SUPPORTED_SOURCE_PROFILE' and (
                {value for sample in samples for value in sample['organism_ch1']}
                    != {RNA_ORGANISM_BY_SPECIES[release['species']]}
                or {value.casefold() for sample in samples for value in sample['molecule_ch1']}
                    != {'total rna'}):
            raise ValueError('Supported source profile conflicts with organism or total-RNA assay: ' + study)
        if row.get('scientific_status') not in {'CANNOT', 'CANNOT VERIFY'}:
            raise ValueError('Source-profile scientific status must remain unresolved: ' + study)
        rows[study] = row
        resolved_sources[study] = resolved
    if set(rows) != set(supported) | set(stopped):
        raise ValueError('Source-profile partition does not cover every enrolled study')
    if any(rows[study].get('profile_status') != 'SUPPORTED_SOURCE_PROFILE' for study in supported):
        raise ValueError('Supported partition contains a non-supported source profile')
    if any(rows[study].get('profile_status') != 'UNSUPPORTED_SOURCE_PROFILE' for study in stopped):
        raise ValueError('Stopped partition contains a supported source profile')
    counts = descriptor['counts']
    expected_counts = {'enrolled_studies': len(rows), 'enrolled_samples': len(all_samples),
                       'supported_studies': len(supported), 'stopped_studies': len(stopped),
                       'molecule_label_occurrences': molecule_counts}
    if any(counts.get(key) != value for key, value in expected_counts.items()):
        raise ValueError('Source-profile counts do not reproduce its records')
    source_files = descriptor['source_files']
    if not isinstance(source_files, dict) or not source_files:
        raise ValueError('Source-profile source-file ledger is absent')
    verified_files = {}
    release_path = Path(release['path']).resolve()
    for name, checksum in source_files.items():
        path = Path(name).expanduser().resolve()
        if not re.fullmatch(r'[0-9a-f]{64}', checksum or '') or not path.is_file():
            raise ValueError('Source-profile source-file declaration is malformed')
        if path == release_path:
            if checksum != release['sha256']:
                raise ValueError('Source-profile release identity differs from batch plan')
            verified_files[str(path)] = {'path': str(path), 'sha256': checksum,
                                         'size_bytes': path.stat().st_size, 'release_hash_deferred': True}
        else:
            verified_files[str(path)] = _batch_file_declaration(
                {'path': str(path), 'sha256': checksum}, 'source_profile.source_files')
    for study, sources in resolved_sources.items():
        for role, record in sources.items():
            if verified_files.get(record['path'], {}).get('sha256') != record['sha256']:
                raise ValueError(f'{study} {role} is absent from the source-profile source-file ledger')
    return descriptor_decl, descriptor, rows, resolved_sources, verified_files


def _batch_implementation(plan, destination):
    declarations = plan.get('implementation')
    if not isinstance(declarations, dict) or set(declarations) != set(BATCH_IMPLEMENTATIONS):
        raise ValueError('Batch extraction requires the exact executed implementation inventory')
    verified = {}
    frozen = {}
    frozen_root = destination / 'implementation'
    frozen_root.mkdir(parents=True, exist_ok=True)
    for role, expected_path in BATCH_IMPLEMENTATIONS.items():
        record = _batch_file_declaration(declarations[role], 'implementation.' + role)
        if Path(record['path']) != expected_path:
            raise ValueError('Batch implementation path differs for role: ' + role)
        target = frozen_root / (role + expected_path.suffix)
        if target.exists():
            if digest(target) != record['sha256']:
                raise ValueError('Frozen batch implementation changed: ' + role)
        else:
            with expected_path.open('rb') as source, target.open('xb') as output:
                shutil.copyfileobj(source, output)
            if digest(target) != record['sha256']:
                raise ValueError('Batch implementation changed during snapshot: ' + role)
        verified[role] = record
        frozen[role] = {'path': str(target), 'sha256': record['sha256']}
    return verified, frozen


def _batch_verify_small_inputs(implementation, source_files, descriptor_decl,
                               plan_path, plan_sha256):
    if digest(plan_path) != plan_sha256:
        raise ValueError('Batch extraction plan changed during execution')
    for role, record in implementation.items():
        if digest(Path(record['path'])) != record['sha256']:
            raise ValueError('Batch implementation changed during execution: ' + role)
    if digest(Path(descriptor_decl['path'])) != descriptor_decl['sha256']:
        raise ValueError('Source-profile descriptor changed during execution')
    for record in source_files.values():
        if record.get('release_hash_deferred'):
            continue
        if digest(Path(record['path'])) != record['sha256']:
            raise ValueError('Source-profile source file changed during execution: ' + record['path'])


def _batch_stage_receipt(destination, study_store, stage):
    path = study_store.root / stage / 'complete.json'
    return {'path': str(path.relative_to(destination)), 'sha256': digest(path)}


def _batch_relative_file(destination, path):
    path = Path(path).resolve()
    if not path.is_relative_to(destination):
        raise ValueError('Batch output escaped its destination: ' + str(path))
    return {'path': str(path.relative_to(destination)), 'sha256': digest(path),
            'size_bytes': path.stat().st_size}


def _batch_extract_from_index(destination, directory, *, study, release, h5, genes,
                              geo_ids, memberships, indices, descriptor_record,
                              metadata_result, source_records, descriptor_sha256,
                              source_profile_record_sha256, preflight):
    """Use extract_rna's matrix construction and require_matrix semantics.

    The batch coordinator supplies indices from one shared release metadata scan;
    the historical extract_rna entry point remains unchanged.
    """
    import pandas as pd
    from src.pipeline_inputs_2026_09_10 import require_matrix
    before = _batch_stat(release['path'])
    selected = [geo_ids[index] for index in indices]
    frame = pd.DataFrame(h5['data/expression'][:, indices], index=genes, columns=selected)
    after = _batch_stat(release['path'])
    if before != after:
        from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError
        raise PipelineIntegrityError('ARCHS4 source changed during batch study extraction')
    require_matrix(frame, beta=False)
    counts = directory / 'counts.csv'
    frame.rename_axis('ensembl_gene').to_csv(counts, mode='x')
    if counts.stat().st_size > preflight['max_extracted_bytes_per_study']:
        raise ValueError('Extracted matrix exceeds per-study byte cap')
    matrix = _batch_relative_file(destination, counts)
    matrix.update(feature_count=len(frame), sample_count=len(frame.columns),
                  cell_count=len(frame) * len(frame.columns),
                  feature_id_sha256=_intake_json_hash(list(frame.index)),
                  sample_id_sha256=_intake_json_hash(list(frame.columns)))
    receipt = {
        'schema': BATCH_STUDY_EXTRACTION_SCHEMA, 'status': 'PASS',
        'path_base': 'batch_root',
        'study_id': study, 'modality': RNA_MODALITY_BY_SPECIES[release['species']], 'species': release['species'],
        'matrix': matrix, 'sample_ids_in_matrix_order': list(frame.columns),
        'all_study_memberships': {geo_ids[index]: [value.strip() for value in memberships[index].split(',')]
                                  for index in indices},
        'input_mode': 'FRESH_ARCHS4_BATCH_EXTRACTION',
        'source': {'release_sha256': release['sha256'],
                   'release_size_bytes': release['size_bytes'],
                   'expression_shape': release['expression_shape'],
                   'expression_dtype': str(h5['data/expression'].dtype),
                   'series_sha256': source_records['series']['sha256'],
                   'samples_sha256': source_records['samples']['sha256'],
                   'metadata_sha256': digest(Path(metadata_result['metadata'])),
                   'context_sha256': digest(Path(metadata_result['study_context'])),
                   'source_profile_descriptor_sha256': descriptor_sha256,
                   'source_profile_record_sha256': source_profile_record_sha256},
        'preflight': preflight,
        'scientific_status': 'CANNOT VERIFY',
        'original_paper_2025_eligibility': 'QUARANTINED'}
    receipt_path = directory / 'extraction.json'
    _intake_write_new(receipt_path, receipt)
    return {'receipt_path': str(receipt_path.relative_to(destination)),
            'receipt_sha256': digest(receipt_path), 'matrix': matrix,
            'sample_ids': list(frame.columns)}


def _batch_verify_completed(destination, summary, plan_sha256, release, descriptor_decl,
                            implementation, source_files, plan_path):
    if summary.get('schema') != BATCH_SUMMARY_SCHEMA or summary.get('plan_sha256') != plan_sha256:
        raise ValueError('Completed batch summary belongs to another plan')
    for relative, checksum in summary.get('artifacts', {}).items():
        path = (destination / relative).resolve()
        if not path.is_relative_to(destination) or not path.is_file() or digest(path) != checksum:
            raise ValueError('Completed batch artifact changed: ' + relative)
    _batch_verify_small_inputs(implementation, source_files, descriptor_decl,
                               plan_path, plan_sha256)
    observed = digest(Path(release['path']))
    if observed != release['sha256']:
        raise ValueError('ARCHS4 release changed after completed batch')
    return summary


def run_archs4_batch_extraction(plan_path, output=None):
    """Extract a source-profile-gated ARCHS4 cohort without clocks or models."""
    plan_path = Path(plan_path).expanduser().resolve()
    try:
        plan_bytes = plan_path.read_bytes()
        plan = json.loads(plan_bytes)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError('Cannot read ARCHS4 batch extraction plan') from error
    plan_sha256 = hashlib.sha256(plan_bytes).hexdigest()
    if plan.get('schema') != BATCH_EXTRACTION_SCHEMA:
        raise ValueError('Unsupported ARCHS4 batch extraction schema')
    declared_output = output if output is not None else plan.get('output')
    if not isinstance(declared_output, (str, os.PathLike)) or not str(declared_output):
        raise ValueError('A dedicated batch extraction output is required')
    destination = Path(declared_output).expanduser().resolve()
    if destination == Path('/').resolve():
        raise ValueError('A dedicated batch extraction output is required')
    limits = plan.get('limits')
    required_limits = {'max_studies', 'max_samples_per_study', 'max_extracted_bytes_per_study',
                       'max_total_extracted_bytes', 'concurrency'}
    if (not isinstance(limits, dict) or set(limits) != required_limits
            or limits.get('concurrency') != 1
            or any(type(limits.get(key)) is not int or limits[key] <= 0
                   for key in required_limits - {'concurrency'})):
        raise ValueError('Batch extraction limits are incomplete or invalid')
    selected = plan.get('selected_study_ids')
    if (not isinstance(selected, list) or not selected or len(selected) > limits['max_studies']
            or len(selected) != len(set(selected))
            or any(not re.fullmatch(r'GSE[0-9]+', value or '') for value in selected)):
        raise ValueError('Selected batch study inventory is invalid')
    release = _batch_release_declaration(plan.get('release'))
    modality = RNA_MODALITY_BY_SPECIES[release['species']]
    organism = RNA_ORGANISM_BY_SPECIES[release['species']]
    if plan.get('modality', 'human_rna') != modality:
        raise ValueError('Batch plan modality differs from its declared ARCHS4 species')
    runtime = _intake_runtime_versions(plan)
    for package in ('h5py', 'numpy', 'pandas'):
        if package not in runtime['packages']:
            raise ValueError('Batch runtime lacks required package pin: ' + package)
    if destination.exists() and any(destination.iterdir()) and not (destination / 'config.json').is_file():
        raise ValueError('Batch extraction requires a new output or its own resumable configuration')
    destination.mkdir(parents=True, exist_ok=True)
    with _intake_v21_run_lock(destination):
        descriptor_decl, descriptor, profile_rows, profile_sources, source_files = _batch_source_profile(
            plan.get('source_profile_descriptor'), release)
        if selected != descriptor['supported_study_ids']:
            raise ValueError('Selected studies must exactly preserve the supported source-profile order')
        if len(profile_rows) > limits['max_studies']:
            raise ValueError('Enrolled source-profile study count exceeds the batch cap')
        if any(len(profile_rows[study]['sample_ids']) > limits['max_samples_per_study']
               for study in profile_rows):
            raise ValueError('Source-profile sample count exceeds the batch cap')
        implementation, frozen_implementation = _batch_implementation(plan, destination)
        config = {'schema': 'clockbase-archs4-batch-extraction-config-v1',
                  'modality': modality,
                  'plan_path': str(plan_path), 'plan_sha256': plan_sha256,
                  'source_profile_descriptor': descriptor_decl, 'release': release,
                  'selected_study_ids': selected, 'stopped_study_ids': descriptor['stopped_study_ids'],
                  'limits': limits, 'runtime': runtime, 'implementation': implementation,
                  'scientific_status': 'CANNOT VERIFY',
                  'original_paper_2025_eligibility': 'QUARANTINED'}
        config_path = destination / 'config.json'
        if config_path.exists():
            if _batch_read_sealed(config_path) != config:
                raise ValueError('Batch extraction configuration changed; use a new output')
        else:
            _batch_seal(config_path, config)
        summary_path = destination / 'batch-summary.json'
        if summary_path.exists():
            summary = _batch_read_sealed(summary_path)
            return _batch_verify_completed(destination, summary, plan_sha256, release,
                                           descriptor_decl, implementation, source_files, plan_path)

        attempt_root = destination / 'execution-attempts'
        attempt_root.mkdir(exist_ok=True)
        attempt = attempt_root / ('attempt-%03d' % len(list(attempt_root.glob('attempt-*'))))
        attempt.mkdir()
        release_stat_before = _batch_stat(release['path'])
        release_sha_before = digest(Path(release['path']))
        before = {'schema': 'clockbase-archs4-batch-release-check-v1', 'phase': 'before',
                  'path': release['path'], 'sha256': release_sha_before,
                  'size_bytes': release_stat_before['size_bytes'], 'hash_scope': 'FULL_HDF5',
                  'stat': release_stat_before,
                  'expected_sha256': release['sha256'], 'status': 'PASS' if release_sha_before == release['sha256'] else 'INTEGRITY_FAILURE'}
        before_path = attempt / 'release-before.json'
        _batch_seal(before_path, before)
        if release_sha_before != release['sha256']:
            _batch_seal(attempt / 'attempt.json', {
                'schema': 'clockbase-archs4-batch-execution-attempt-v1',
                'status': 'INTEGRITY_FAILURE', 'plan_sha256': plan_sha256,
                'release_before': {'path': str(before_path.relative_to(destination)),
                                   'sha256': digest(before_path)},
                'release_after': None, 'result_count': 0,
                'error': 'ValueError: ARCHS4 release SHA256 differs before batch extraction',
                'model_api_calls': 0, 'preprocessing_runs': 0, 'clock_runs': 0})
            raise ValueError('ARCHS4 release SHA256 differs before batch extraction')

        results = []
        study_stores = {}
        metadata = {}
        failure = None
        try:
            from src.pipeline_runtime_2026_09_10 import StageStore
            from src.update_pipeline_2026_09_20 import frozen_metadata

            for study in descriptor['stopped_study_ids']:
                row = profile_rows[study]
                store = StageStore(destination / 'studies' / study,
                    {'schema': BATCH_EXTRACTION_SCHEMA, 'study_id': study,
                     'plan_sha256': plan_sha256,
                     'source_profile_descriptor_sha256': descriptor_decl['sha256'],
                     'source_profile_record_sha256': _intake_json_hash(row),
                     'profile_status': row['profile_status']})
                study_stores[study] = store

                def stop(directory, study=study, row=row):
                    record = {'schema': 'clockbase-archs4-source-profile-stop-v1',
                              'study_id': study, 'modality': modality, 'species': release['species'],
                              'status': 'SOURCE_PROFILE_UNRESOLVED',
                              'profile_status': row['profile_status'],
                              'limitations': row.get('limitations', []),
                              'sample_ids': row['sample_ids'],
                              'source_profile_descriptor_sha256': descriptor_decl['sha256'],
                              'source_profile_record_sha256': _intake_json_hash(row),
                              'scientific_status': 'CANNOT VERIFY',
                              'original_paper_2025_eligibility': 'QUARANTINED'}
                    path = directory / 'source-profile-stop.json'
                    _intake_write_new(path, record)
                    return {'path': str(path.relative_to(destination)), 'sha256': digest(path)}
                stop_result = store.run('source_profile_stop', stop)
                results.append({'study_id': study, 'modality': modality,
                    'source_profile_status': row['profile_status'],
                    'source_profile_record_sha256': _intake_json_hash(row),
                    'status': 'SOURCE_PROFILE_UNRESOLVED',
                    'series': profile_sources[study]['series'], 'samples': profile_sources[study]['samples'],
                    'matrix': None, 'extraction_receipt': None, 'stop_receipt': stop_result,
                    'stage_receipt': _batch_stage_receipt(destination, store, 'source_profile_stop')})

            for study in selected:
                row = profile_rows[study]
                store = StageStore(destination / 'studies' / study,
                    {'schema': BATCH_EXTRACTION_SCHEMA, 'study_id': study,
                     'plan_sha256': plan_sha256,
                     'source_profile_descriptor_sha256': descriptor_decl['sha256'],
                     'source_profile_record_sha256': _intake_json_hash(row),
                     'profile_status': row['profile_status']})
                study_stores[study] = store
                source_records = profile_sources[study]

                def prepare_metadata(directory, study=study, source_records=source_records):
                    value = frozen_metadata(study, directory, source_records)
                    rows, sample_ids, _ = _intake_read_metadata(value['metadata'], study, release['species'])
                    if set(sample_ids) != set(profile_rows[study]['sample_ids']):
                        raise ValueError('Prepared metadata differs from the source-profile GSM membership')
                    return {**value, 'sample_ids': sample_ids,
                            'metadata_reindexed_from_profile': sample_ids != profile_rows[study]['sample_ids']}
                metadata_result = _intake_terminal_stage(store, 'metadata', prepare_metadata)
                metadata[study] = metadata_result

            import h5py
            import numpy as np
            with h5py.File(release['path'], 'r') as h5:
                expression = h5['data/expression']
                shape = list(expression.shape)
                dtype = np.dtype(expression.dtype)
                if shape != release['expression_shape'] or shape[1] != release['sample_count']:
                    raise ValueError('ARCHS4 expression shape differs from the batch plan')
                if dtype.kind not in 'biuf' or dtype.itemsize > 16:
                    raise ValueError('Unsupported ARCHS4 expression dtype for batch extraction')
                decode = lambda values: [value.decode() if isinstance(value, bytes) else str(value) for value in values]
                geo_ids = decode(h5['meta/samples/geo_accession'][:])
                memberships = decode(h5['meta/samples/series_id'][:])
                organisms = decode(h5['meta/samples/organism_ch1'][:])
                genes = decode(h5['meta/genes/ensembl_gene'][:])
                if (len(geo_ids) != shape[1] or len(geo_ids) != len(set(geo_ids))
                        or len(memberships) != len(geo_ids) or len(organisms) != len(geo_ids)):
                    raise ValueError('ARCHS4 sample identities are duplicated or unaligned')
                if (len(genes) != shape[0] or len(genes) != len(set(genes))
                        or any(not value for value in genes)):
                    raise ValueError('ARCHS4 gene identities are duplicated or unaligned')
                if max([len(value) for value in genes + geo_ids + selected], default=0) > HARDENED_V21_IDENTIFIER_CHARS:
                    raise ValueError(f'ARCHS4 batch identifier exceeds {HARDENED_V21_IDENTIFIER_CHARS} characters')
                selected_set = set(selected)
                indices_by_study = {study: [] for study in selected}
                for index, membership in enumerate(memberships):
                    for study in {value.strip() for value in membership.split(',') if value.strip()} & selected_set:
                        indices_by_study[study].append(index)
                safe_total = 0
                preflights = {}
                for study in selected:
                    indices = indices_by_study[study]
                    observed_ids = [geo_ids[index] for index in indices]
                    sample_count = len(profile_rows[study]['sample_ids'])
                    safe_bound = (len(genes) * sample_count * (HARDENED_V21_NUMERIC_CHARS + 2)
                                  + len(genes) * (HARDENED_V21_IDENTIFIER_CHARS + 1)
                                  + sample_count * (HARDENED_V21_IDENTIFIER_CHARS + 1) + 1024)
                    supported = (set(observed_ids) == set(profile_rows[study]['sample_ids'])
                                 and len(observed_ids) == sample_count
                                 and {organisms[index] for index in indices} == {organism})
                    within_cap = safe_bound <= limits['max_extracted_bytes_per_study']
                    if within_cap:
                        safe_total += safe_bound
                    preflights[study] = {'study_id': study, 'species': release['species'],
                        'sample_count': sample_count, 'gene_count': len(genes),
                        'expression_dtype': str(dtype), 'expression_shape': shape,
                        'selected_sample_ids': observed_ids, 'membership_and_organism_match': supported,
                        'safe_upper_bound_bytes': safe_bound,
                        'max_extracted_bytes_per_study': limits['max_extracted_bytes_per_study'],
                        'within_per_study_cap': within_cap,
                        'bound_profile': 'fixed-v2.1'}
                if safe_total > limits['max_total_extracted_bytes']:
                    raise ValueError('Aggregate conservative extraction bound exceeds the batch cap')

                extracted_bytes = 0
                for study in selected:
                    row = profile_rows[study]
                    store = study_stores[study]
                    metadata_result = metadata[study]
                    preflight = preflights[study]

                    def extract(directory, study=study, row=row, metadata_result=metadata_result,
                                preflight=preflight):
                        if metadata_result['status'] != 'PASS':
                            raise ValueError('Prepared metadata stage is unresolved')
                        if not preflight['membership_and_organism_match']:
                            raise ValueError('ARCHS4/GEO full roster or organism differs; no partial extraction')
                        if not preflight['within_per_study_cap']:
                            raise ValueError('Conservative extraction bound exceeds the per-study cap')
                        return _batch_extract_from_index(
                            destination, directory, study=study, release=release, h5=h5, genes=genes,
                            geo_ids=geo_ids, memberships=memberships, indices=indices_by_study[study],
                            descriptor_record=row, metadata_result=metadata_result['value'],
                            source_records=profile_sources[study], descriptor_sha256=descriptor_decl['sha256'],
                            source_profile_record_sha256=_intake_json_hash(row), preflight=preflight)
                    extraction_result = _intake_terminal_stage(store, 'extract', extract)
                    result = {'study_id': study, 'modality': modality,
                        'source_profile_status': row['profile_status'],
                        'source_profile_record_sha256': _intake_json_hash(row),
                        'series': profile_sources[study]['series'], 'samples': profile_sources[study]['samples'],
                        'stage_receipt': _batch_stage_receipt(destination, store, 'extract'),
                        'stop_receipt': None}
                    if extraction_result['status'] == 'PASS':
                        value = extraction_result['value']
                        extracted_bytes += value['matrix']['size_bytes']
                        if extracted_bytes > limits['max_total_extracted_bytes']:
                            raise ValueError('Aggregate extracted bytes exceed the batch cap')
                        result.update(status='EXTRACTED', matrix=value['matrix'],
                                      extraction_receipt={'path': value['receipt_path'],
                                                          'sha256': value['receipt_sha256']})
                    else:
                        result.update(status='EXTRACTION_UNRESOLVED', matrix=None,
                                      extraction_receipt=None,
                                      error=extraction_result.get('error'))
                    results.append(result)
        except Exception as error:  # preserve a batch-level failure before returning it
            failure = error

        release_stat_after = _batch_stat(release['path'])
        release_sha_after = digest(Path(release['path']))
        after_status = 'PASS' if release_sha_after == release['sha256'] else 'INTEGRITY_FAILURE'
        after = {'schema': 'clockbase-archs4-batch-release-check-v1', 'phase': 'after',
                 'path': release['path'], 'sha256': release_sha_after,
                 'size_bytes': release_stat_after['size_bytes'], 'hash_scope': 'FULL_HDF5',
                 'stat': release_stat_after,
                 'expected_sha256': release['sha256'], 'status': after_status}
        after_path = attempt / 'release-after.json'
        _batch_seal(after_path, after)
        if release_sha_after != release_sha_before:
            failure = ValueError('ARCHS4 release changed between full-file batch hashes')
        try:
            _batch_verify_small_inputs(implementation, source_files, descriptor_decl,
                                       plan_path, plan_sha256)
        except Exception as error:
            failure = failure or error
        ordered_results = None
        counts = None
        if failure is None:
            try:
                by_study = {row['study_id']: row for row in results}
                if set(by_study) != set(profile_rows) or len(by_study) != len(results):
                    raise ValueError('Batch results do not account for every enrolled study exactly once')
                ordered_results = [by_study[row['study_id']] for row in descriptor['studies']]
                extracted = sum(row['status'] == 'EXTRACTED' for row in ordered_results)
                source_stops = sum(row['status'] == 'SOURCE_PROFILE_UNRESOLVED'
                                   for row in ordered_results)
                extraction_unresolved = sum(row['status'] == 'EXTRACTION_UNRESOLVED'
                                            for row in ordered_results)
                counts = {'enrolled_studies': len(ordered_results), 'selected_studies': len(selected),
                          'source_profile_stops': source_stops, 'extracted_studies': extracted,
                          'extraction_unresolved_studies': extraction_unresolved,
                          'enrolled_samples': descriptor['counts']['enrolled_samples'],
                          'selected_samples': sum(len(profile_rows[study]['sample_ids'])
                                                  for study in selected),
                          'extracted_samples': sum(row['matrix']['sample_count'] for row in ordered_results
                                                   if row['matrix'] is not None),
                          'extracted_bytes': sum(row['matrix']['size_bytes'] for row in ordered_results
                                                 if row['matrix'] is not None)}
            except Exception as error:
                failure = error
        attempt_record = {'schema': 'clockbase-archs4-batch-execution-attempt-v1',
                          'status': 'INTEGRITY_FAILURE' if release_sha_after != release_sha_before else
                                    ('FAIL' if failure else 'COMPLETED'),
                          'plan_sha256': plan_sha256,
                          'release_before': {'path': str(before_path.relative_to(destination)),
                                             'sha256': digest(before_path)},
                          'release_after': {'path': str(after_path.relative_to(destination)),
                                            'sha256': digest(after_path)},
                          'result_count': len(results),
                          'error': None if failure is None else type(failure).__name__ + ': ' + str(failure),
                          'model_api_calls': 0, 'preprocessing_runs': 0, 'clock_runs': 0}
        attempt_path = attempt / 'attempt.json'
        _batch_seal(attempt_path, attempt_record)
        if failure is not None:
            raise failure
        artifacts = {str(path.relative_to(destination)): digest(path)
                     for path in destination.rglob('*')
                     if path.is_file() and path.name not in {'batch-summary.json', 'batch-summary.json.sha256'}
                     and path != _intake_v21_lock_path(destination)}
        summary = {'schema': BATCH_SUMMARY_SCHEMA, 'modality': modality,
                   'path_base': 'batch_root',
                   'status': 'COMPLETED_WITH_UNRESOLVED' if source_stops or extraction_unresolved else 'COMPLETED',
                   'execution_status': 'COMPLETED', 'plan_sha256': plan_sha256,
                   'config_sha256': digest(config_path),
                   'source_profile_descriptor': descriptor_decl,
                   'source_hdf5': {'sha256': release['sha256'], 'size_bytes': release['size_bytes'],
                                   'species': release['species'],
                                   'expression_shape': release['expression_shape'],
                                   'sample_count': release['sample_count'],
                                   'expression_dtype': str(dtype),
                                   'hash_scope': 'full_file_before_and_after_single_batch'},
                   'release_before': attempt_record['release_before'],
                   'release_after': attempt_record['release_after'],
                   'implementation': frozen_implementation,
                   'counts': counts, 'results': ordered_results,
                   'execution_attempt': {'path': str(attempt_path.relative_to(destination)),
                                         'sha256': digest(attempt_path)},
                   'artifacts': artifacts, 'model_api_calls': 0,
                   'preprocessing_runs': 0, 'clock_runs': 0,
                   'scientific_status': 'CANNOT VERIFY',
                   'original_paper_2025_eligibility': 'QUARANTINED'}
        _batch_seal(summary_path, summary)
        return summary


def remaining_storage(declarations, destination):
    remaining=0
    for species,declaration in declarations.items():
        target=destination/f'{species}_gene_v2.latest.h5'
        part=target.with_suffix('.h5.part')
        path=target if target.exists() else part
        # aria2 parts are sparse: logical file length is not downloaded bytes.
        allocated=min(path.stat().st_blocks*512,declaration['file_size']) if path.exists() else 0
        remaining+=max(0,declaration['file_size']-allocated)
    return remaining


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog', type=Path)
    parser.add_argument('--destination', type=Path)
    parser.add_argument('--old-root', type=Path)
    parser.add_argument('--batch-extraction-plan', type=Path)
    parser.add_argument('--batch-output', type=Path)
    args = parser.parse_args()
    if args.batch_extraction_plan is not None:
        if any(value is not None for value in (args.catalog, args.destination, args.old_root)):
            parser.error('--batch-extraction-plan cannot be combined with release-download arguments')
        result = run_archs4_batch_extraction(args.batch_extraction_plan, output=args.batch_output)
        print(json.dumps({'status': result['status'], 'counts': result['counts']}, sort_keys=True))
        return
    if args.batch_output is not None:
        parser.error('--batch-output requires --batch-extraction-plan')
    if any(value is None for value in (args.catalog, args.destination, args.old_root)):
        parser.error('--catalog, --destination and --old-root are required for release download mode')
    declarations = json.loads(args.catalog.read_text())
    args.destination.mkdir(parents=True, exist_ok=True)
    lock_path = args.destination / 'release-selection.json'
    if lock_path.exists() and json.loads(lock_path.read_text()) != declarations:
        raise ValueError('Resume rejected: release declaration changed')
    if not lock_path.exists():
        lock_path.write_text(json.dumps(declarations, indent=2) + '\n')
    remaining = remaining_storage(declarations,args.destination)
    if shutil.disk_usage(args.destination).free < remaining + 20 * 1024**3:
        raise RuntimeError('Insufficient free space for selected files plus a 20 GiB reserve')
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(download, s, v, args.destination) for s, v in declarations.items()]
        for f in futures:
            f.result()
    for species in declarations:
        target = args.destination / f'{species}_overlap.json'
        if target.exists():
            continue
        old = sample_inventory(args.old_root / f'{species}_gene_v2.4.h5')
        new = sample_inventory(args.destination / f'{species}_gene_v2.latest.h5')
        old_studies = {g for gs in old.values() for g in gs}
        added = set(new) - set(old)
        new_studies = {g for gs in new.values() for g in gs} - old_studies
        shared_study_samples = {g for gsm in set(new) & set(old) for g in new[gsm]}
        # This is an accession/sample screen. Related families and first-public dates remain unchecked.
        summary = {'old_samples': len(old), 'new_samples': len(new), 'added_samples': len(added),
                   'removed_samples': len(set(old)-set(new)), 'new_accessions': len(new_studies),
                   'new_accessions_with_shared_old_samples': len(new_studies & shared_study_samples),
                   'eligible_extension_studies': None, 'first_public_dates_verified': False}
        with (args.destination / f'{species}_new_sample_memberships.jsonl').open('x') as handle:
            for gsm in sorted(added):
                handle.write(json.dumps({'gsm_id': gsm, 'gse_ids': new[gsm]}) + '\n')
        target.write_text(json.dumps(summary, indent=2) + '\n')
        print(json.dumps({'species': species, **summary}), flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        raise SystemExit(f'ARCHS4 update stopped: {type(error).__name__}: {error}')
