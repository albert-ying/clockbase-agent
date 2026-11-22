"""Verify extracted counts or discovered study membership against pinned HDF5."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
import time

import h5py
import numpy as np


def sha(path):
    digest = hashlib.sha256()
    if not Path(path).is_file():
        raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def bound_json(path, suffix):
    if not path.is_file():
        raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
    if sha(path) != path.with_suffix(suffix).read_text().strip():
        raise ValueError('Receipt hash differs: ' + str(path))
    return json.loads(path.read_text())


def verify(execution_path):
    started = time.monotonic()
    execution = json.loads(execution_path.read_text())
    if len(execution['studies']) != 1 or execution['until'] != 'input':
        raise ValueError('This readback expects one input-only study execution')
    study = execution['studies'][0]
    if study['status'] != 'PIPELINE_COMPLETED_THROUGH_INPUT':
        raise ValueError('Input extraction did not complete')
    version_root = Path(study['version_root'])
    version = bound_json(version_root / 'version.json', '.json.sha256')
    stage_path = version_root / 'pipeline/02_molecular_input/complete.json'
    stage = bound_json(stage_path, '.sha256')
    for name, expected in stage['outputs'].items():
        path = (stage_path.parent / name).resolve()
        if not path.is_relative_to(stage_path.parent.resolve()) or sha(path) != expected:
            raise ValueError('Extraction output differs from its completed stage')
    counts = Path(stage['result']['path'])
    counts_relative = str(counts.relative_to(stage_path.parent))
    if counts_relative not in stage['outputs']:
        raise ValueError('Counts file is not a declared completed-stage output')
    release = Path(version['sources']['archs4']['path'])
    before = release.stat()
    release_hash = sha(release)
    if release_hash != version['sources']['archs4']['sha256']:
        raise ValueError('Full HDF5 source digest changed')
    with counts.open(newline='') as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = list(reader)
    columns = header[1:]
    gene_ids = [row[0] for row in rows]
    values = np.asarray([row[1:] for row in rows], dtype=np.float64)
    if len(columns) != len(set(columns)) or len(gene_ids) != len(set(gene_ids)):
        raise ValueError('Duplicate identity in extracted matrix')
    decode = lambda data: [v.decode() if isinstance(v, bytes) else str(v) for v in data]
    with h5py.File(release, 'r') as handle:
        source_ids = decode(handle['meta/samples/geo_accession'][:])
        if len(source_ids) != len(set(source_ids)):
            raise ValueError('Duplicate identity in HDF5')
        membership = decode(handle['meta/samples/series_id'][:])
        positions = [index for index, labels in enumerate(membership)
                     if study['study_id'] in {v.strip() for v in labels.split(',')}]
        if [source_ids[index] for index in positions] != columns:
            raise ValueError('Study HDF5 GSM membership/order differs from extraction')
        if decode(handle['meta/genes/ensembl_gene'][:]) != gene_ids:
            raise ValueError('HDF5 gene identities/order differ from extraction')
        if set(decode(handle['meta/samples/organism_ch1'][positions])) != {'Homo sapiens'}:
            raise ValueError('Source organism differs from this human RNA verification')
        expected_values = handle['data/expression'][:, positions]
    if values.shape != expected_values.shape or not np.array_equal(values, expected_values):
        raise ValueError('Extracted matrix differs from HDF5 source cells')
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('RNA count range is invalid')
    if (before.st_size, before.st_mtime_ns) != (release.stat().st_size, release.stat().st_mtime_ns):
        raise ValueError('HDF5 changed during independent readback')
    series_path = Path(version['sources']['series']['path'])
    sample_path = Path(version['sources']['samples']['path'])
    for role, path in [('series', series_path), ('samples', sample_path)]:
        if sha(path) != version['sources'][role]['sha256']:
            raise ValueError('GEO source changed: ' + role)
    series = series_path.read_text()
    series_samples = re.findall(r'^!Series_sample_id = (GSM[0-9]+)$', series, re.M)
    samples = re.findall(r'^\^SAMPLE = (GSM[0-9]+)$', sample_path.read_text(), re.M)
    if set(series_samples) != set(columns) or set(samples) != set(columns):
        raise ValueError('GEO source identity differs from HDF5 extraction')
    dates = re.findall(r'^!Series_status = Public on (.+)$', series, re.M)
    if len(dates) != 1:
        raise ValueError('Unique source first-public date absent')
    first_public = datetime.strptime(dates[0], '%b %d %Y').date().isoformat()
    if first_public != study['first_public_date']:
        raise ValueError('First-public source parsing differs')
    return {
        'status': 'PASS', 'scope': 'Source identity and every extracted count; no scientific acceptance',
        'study_id': study['study_id'], 'source_features': len(gene_ids), 'source_samples': len(columns),
        'source_sample_ids': columns, 'count_cells_compared': int(values.size),
        'all_count_cells_exactly_equal': True, 'full_hdf5_hash_recomputed': release_hash,
        'hdf5_source': str(release), 'hdf5_size_bytes': before.st_size,
        'first_public_date': first_public, 'first_public_source': str(series_path),
        'counts_sha256': sha(counts), 'counts_path': str(counts),
        'execution_path': str(execution_path), 'execution_sha256': sha(execution_path),
        'version_sha256': sha(version_root / 'version.json'),
        'verification_script_sha256': sha(Path(__file__)),
        'reserved_model_calls': execution['reserved_model_calls'],
        'reserved_cost_usd': execution['reserved_cost_usd'],
        'clock_rerun_performed': False, 'scientific_status': 'CANNOT VERIFY',
        'original_paper_2025_eligibility': 'QUARANTINED',
        'elapsed_seconds': time.monotonic() - started,
    }


def discovery_availability(discovery_path, release_path, expected_release_sha256):
    """Read metadata only; matching membership does not validate an assay or clock."""
    started = time.monotonic()
    discovery_path, release_path = discovery_path.resolve(), release_path.resolve()
    discovery_hash = sha(discovery_path)
    verifier_hash = sha(Path(__file__))
    summary = bound_json(discovery_path, '.json.sha256')
    if (summary.get('schema') != 'clockbase-update-discovery-result-v1'
            or summary.get('status') != 'COMPLETED'):
        raise ValueError('Expected a completed discovery summary')
    for name, expected in summary['output_sha256'].items():
        path = Path(name).resolve()
        if not path.is_relative_to(discovery_path.parent) or sha(path) != expected:
            raise ValueError('Discovery output changed or escapes its directory: ' + name)
    for name, expected in summary['input_sha256'].items():
        if sha(name) != expected:
            raise ValueError('Discovery input changed: ' + name)
    ledger_path = discovery_path.parent / 'candidate-ledger.jsonl'
    if summary['output_sha256'].get(str(ledger_path)) != sha(ledger_path):
        raise ValueError('Candidate ledger is not bound to the completed discovery')
    config_path = discovery_path.parent / 'config.json'
    config_hash = sha(config_path)
    if summary['output_sha256'].get(str(config_path)) != config_hash:
        raise ValueError('Discovery configuration is not bound to its completion')
    config = bound_json(config_path, '.json.sha256')
    plan = config['plan']
    if plan['modality'] not in {'human_rna', 'mouse_rna'}:
        raise ValueError('ARCHS4 availability requires an RNA discovery')
    species = 'Homo sapiens' if plan['modality'] == 'human_rna' else 'Mus musculus'
    with ledger_path.open() as handle:
        ledger = [json.loads(line) for line in handle if line.strip()]
    if (len(ledger) != summary['enumerated_uids']
            or len({r['geo_uid'] for r in ledger}) != len(ledger)):
        raise ValueError('Discovery UID inventory is incomplete or duplicated')
    accessions = [r['accession'] for r in ledger if r.get('accession') is not None]
    if (len(accessions) != len(set(accessions))
            or any(not isinstance(value, str) or not re.fullmatch(r'GSE[0-9]+', value) for value in accessions)):
        raise ValueError('Discovery Series accessions are noncanonical or duplicated')
    seen_studies, seen_samples, seen_hashes = set(), set(), {}
    for name in ('seen_worklist', 'seen_sample_sources'):
        declaration = plan.get(name)
        if declaration is None:
            continue
        path = Path(declaration['path']).expanduser().resolve()
        digest = sha(path)
        if (digest != declaration['sha256'] or summary['input_sha256'].get(str(path)) != digest):
            raise ValueError('Operational seen source is not bound to the completed discovery: ' + name)
        seen_hashes[str(path)] = digest
        with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as handle:
            if name == 'seen_worklist':
                for line in handle:
                    value = json.loads(line).get('study_id', '')
                    if re.fullmatch(r'GSE[0-9]+', value):
                        seen_studies.add(value)
            else:
                for item in json.load(handle):
                    if item.get('accession') in seen_studies:
                        seen_samples.update(v for v in item.get('sample_ids', [])
                                            if re.fullmatch(r'GSM[0-9]+', v))
    if not re.fullmatch(r'[0-9a-f]{64}', expected_release_sha256 or ''):
        raise ValueError('A pinned release SHA256 is required')
    before = release_path.stat()
    observed_hash = sha(release_path)
    if observed_hash != expected_release_sha256:
        raise ValueError('Pinned ARCHS4 release SHA256 differs')
    target_studies = {r['accession'] for r in ledger if r.get('accession')}
    source_sets = {study: set() for study in target_studies}
    organisms = {}
    decode = lambda data: [v.decode() if isinstance(v, bytes) else str(v) for v in data]
    with h5py.File(release_path, 'r') as handle:
        samples = decode(handle['meta/samples/geo_accession'][:])
        memberships = decode(handle['meta/samples/series_id'][:])
        species_values = decode(handle['meta/samples/organism_ch1'][:])
        expression_shape = list(handle['data/expression'].shape)
    if (len(samples) != len(memberships) or len(samples) != len(species_values)
            or len(set(samples)) != len(samples)
            or any(not re.fullmatch(r'GSM[0-9]+', sample) for sample in samples)
            or len(expression_shape) != 2 or expression_shape[1] != len(samples)):
        raise ValueError('ARCHS4 sample identities, metadata or expression dimensions disagree')
    for sample, labels, organism in zip(samples, memberships, species_values):
        targets = target_studies.intersection(v.strip() for v in labels.split(',') if v.strip())
        for study in targets:
            source_sets[study].add(sample)
            organisms[sample] = organism
    rows, molecular_sources, source_hashes = [], {}, {}
    eligible_states = {'UNRESOLVED_NO_ASSAY_INPUT', 'READY_FOR_INGESTION_ATTEMPT'}
    for record in ledger:
        study = record.get('accession')
        observed = source_sets.get(study, set())
        row = {'geo_uid': record['geo_uid'], 'study_id': study,
               'source_candidate_status': record['status'],
               'first_public_date': record.get('first_public_date'),
               'archs4_sample_ids': sorted(observed), 'geo_sample_ids': [],
               'availability': 'SOURCE_MEMBERSHIP_UNVERIFIED',
               'declared_as_available_input': False}
        receipt = record.get('source_receipt')
        if receipt:
            source_path = Path(receipt['response_path']).resolve()
            source_digest = sha(source_path)
            if (receipt.get('status') != 'PASS' or source_digest != receipt.get('sha256')
                    or source_digest != record.get('series_source_sha256')):
                raise ValueError('GEO source bytes differ for ' + str(study))
            source_hashes[str(source_path)] = source_digest
            text = source_path.read_text()
            source_study = re.findall(r'^\^SERIES = (GSE[0-9]+)\s*$', text, re.M)
            source_samples = re.findall(r'^!Series_sample_id = (GSM[0-9]+)\s*$', text, re.M)
            public_dates = re.findall(r'^!Series_status = Public on (.+)$', text, re.M)
            ledger_samples = record.get('sample_ids')
            if (source_study != [study] or not source_samples
                    or len(source_samples) != len(set(source_samples))
                    or not isinstance(ledger_samples, list)
                    or any(not isinstance(value, str) for value in ledger_samples)
                    or len(ledger_samples) != len(set(ledger_samples))
                    or set(source_samples) != set(ledger_samples)
                    or len(public_dates) != 1):
                raise ValueError('GEO source membership differs from discovery for ' + str(study))
            public_date = datetime.strptime(public_dates[0].strip(), '%b %d %Y').date().isoformat()
            if public_date != record['first_public_date']:
                raise ValueError('Source first-public date differs for ' + str(study))
            expected = set(source_samples)
            row.update(geo_sample_ids=sorted(expected), missing_sample_ids=sorted(expected-observed),
                       extra_sample_ids=sorted(observed-expected), series_source_sha256=source_digest,
                       first_public_date=public_date, organism_values=sorted({organisms[s] for s in observed}))
            if observed == expected:
                row['availability'] = ('EXACT_GEO_MEMBERSHIP' if row['organism_values'] == [species]
                                       else 'ORGANISM_MISMATCH')
            elif not observed:
                row['availability'] = 'STUDY_ABSENT_FROM_RELEASE'
            else:
                row['availability'] = 'SAMPLE_MEMBERSHIP_MISMATCH'
        rows.append(row)
    owners = {}
    for row in rows:
        for sample in row['geo_sample_ids']:
            owners.setdefault(sample, set()).add(row['study_id'])
    for row in rows:
        exclusions = []
        study = row['study_id']
        if study in seen_studies:
            exclusions.append('SEEN_ACCESSION')
        row['shared_seen_sample_ids'] = sorted(set(row['geo_sample_ids']) & seen_samples)
        row['shared_candidate_sample_ids'] = sorted(s for s in row['geo_sample_ids'] if len(owners[s]) > 1)
        if row['shared_seen_sample_ids']:
            exclusions.append('SHARED_SEEN_SAMPLES')
        if row['shared_candidate_sample_ids']:
            exclusions.append('SHARED_CANDIDATE_SAMPLES')
        row['operational_exclusions'] = exclusions
        if (row['availability'] == 'EXACT_GEO_MEMBERSHIP' and not exclusions
                and row['source_candidate_status'] in eligible_states
                and plan['date_from'] <= row['first_public_date'] <= plan['date_to']):
            molecular_sources[study] = {'archs4': {'path': str(release_path), 'sha256': observed_hash}}
            row['declared_as_available_input'] = True
    after = release_path.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
        raise ValueError('ARCHS4 release changed during membership inspection')
    for name, expected in source_hashes.items():
        if sha(name) != expected:
            raise ValueError('GEO source changed during membership inspection')
    for name, expected in seen_hashes.items():
        if sha(name) != expected:
            raise ValueError('Operational overlap source changed during membership inspection')
    if sha(ledger_path) != summary['output_sha256'][str(ledger_path)]:
        raise ValueError('Discovery ledger changed during membership inspection')
    if (sha(discovery_path) != discovery_hash or bound_json(discovery_path, '.json.sha256') != summary
            or sha(config_path) != config_hash or bound_json(config_path, '.json.sha256') != config
            or sha(Path(__file__)) != verifier_hash):
        raise ValueError('Discovery summary, configuration or verifier changed during membership inspection')
    return {'schema': 'clockbase-discovered-rna-input-availability-v1', 'status': 'PASS',
            'recorded_at_utc': datetime.now(timezone.utc).isoformat(),
            'scope': 'Pinned ARCHS4 metadata and source GEO membership; no expression-cell readback, preprocessing, clock or scientific validation',
            'discovery': {'path': str(discovery_path), 'sha256': discovery_hash},
            'discovery_config': {'path': str(config_path), 'sha256': config_hash},
            'release': {'path': str(release_path), 'sha256': observed_hash, 'bytes': before.st_size,
                        'mtime_ns': before.st_mtime_ns, 'samples': len(samples), 'expression_shape': expression_shape},
            'source_files': source_hashes, 'operational_overlap_source_files': seen_hashes, 'rows': rows,
            'counts': {'enumerated_uids': len(rows), 'available_input_studies': len(molecular_sources),
                       'by_availability': dict(sorted(Counter(r['availability'] for r in rows).items()))},
            'molecular_sources': molecular_sources, 'paid_model_calls': 0, 'molecular_downloads': 0,
            'scientific_status': 'CANNOT VERIFY', 'original_paper_2025_eligibility': 'QUARANTINED',
            'verification_script_sha256': verifier_hash, 'elapsed_seconds': time.monotonic()-started}


def discovery_source_profiles(discovery_path):
    """Census literal GEO assay annotations; no model, expression or clock read."""
    started = time.monotonic()
    discovery_path = Path(discovery_path).resolve()
    discovery_hash, implementation_hash = sha(discovery_path), sha(Path(__file__))
    summary = bound_json(discovery_path, '.json.sha256')
    if (summary.get('schema') != 'clockbase-update-discovery-result-v1'
            or summary.get('status') != 'COMPLETED'):
        raise ValueError('Source profiles require a completed discovery summary')
    files = {str(discovery_path): discovery_hash}
    for name, expected in summary['output_sha256'].items():
        path = Path(name).resolve()
        if not path.is_relative_to(discovery_path.parent) or sha(path) != expected:
            raise ValueError('Discovery output changed or escapes its directory: ' + name)
        files[str(path)] = expected
    config_path = discovery_path.parent / 'config.json'
    candidate_path = discovery_path.parent / 'candidates.json'
    ledger_path = discovery_path.parent / 'candidate-ledger.jsonl'
    if any(str(path) not in files for path in (config_path, candidate_path, ledger_path)):
        raise ValueError('Source-profile configuration, candidates and ledger must be bound outputs')
    config = bound_json(config_path, '.json.sha256')
    modality = config['plan']['modality']
    if modality not in {'human_rna', 'mouse_rna'}:
        raise ValueError('The existing total-RNA source profile requires an RNA modality')
    expected_species = 'Homo sapiens' if modality == 'human_rna' else 'Mus musculus'
    candidates = json.loads(candidate_path.read_text())
    if candidates.get('schema') != 'clockbase-update-candidates-v1':
        raise ValueError('Unexpected source-profile candidate schema')
    studies = candidates.get('studies')
    if not isinstance(studies, list) or len(studies) != summary['ready_for_ingestion_attempt']:
        raise ValueError('Candidate count differs from the completed discovery')
    candidate_ids = [r.get('study_id') for r in studies]
    if (len(set(candidate_ids)) != len(candidate_ids)
            or any(not re.fullmatch(r'GSE[0-9]+', value or '') for value in candidate_ids)):
        raise ValueError('Source-profile candidate accessions are invalid or duplicated')
    ledger = [json.loads(line) for line in ledger_path.read_text().splitlines() if line.strip()]
    if len(ledger) != summary['enumerated_uids']:
        raise ValueError('Source-profile ledger does not account for the query inventory')
    ready = [r for r in ledger if r['status'] == 'READY_FOR_INGESTION_ATTEMPT']
    if (len(ready) != len(studies)
            or {r['accession'] for r in ready} != set(candidate_ids)):
        raise ValueError('Source-profile candidates differ from the ready ledger entries')
    ready_by_id = {r['accession']: r for r in ready}
    records, molecule_counts = [], Counter()
    for entry in studies:
        study = entry['study_id']
        if entry['modality'] != modality:
            raise ValueError('Mixed source-profile candidate modalities')
        sources = {}
        for role in ('series', 'samples'):
            declaration = entry['sources'][role]
            if set(declaration) != {'path', 'sha256'}:
                raise ValueError('Source profiles require locally bound GEO ' + role)
            path = Path(declaration['path']).resolve()
            if sha(path) != declaration['sha256']:
                raise ValueError('Source-profile GEO hash differs: ' + study + ' ' + role)
            sources[role] = {'path': str(path), 'sha256': declaration['sha256']}
            files[str(path)] = declaration['sha256']
        series_text = Path(sources['series']['path']).read_text()
        series_ids = re.findall(r'^\^SERIES = (GSE[0-9]+)\s*$', series_text, re.M)
        ids = re.findall(r'^!Series_sample_id = (GSM[0-9]+)\s*$', series_text, re.M)
        dates = re.findall(r'^!Series_status = Public on (.+)$', series_text, re.M)
        source_record = ready_by_id[study]
        if (series_ids != [study] or not ids or len(ids) != len(set(ids))
                or source_record['sample_ids'] != ids or len(dates) != 1
                or source_record['series_source_sha256'] != sources['series']['sha256']):
            raise ValueError('Source-profile Series identity or membership differs: ' + study)
        first_public = datetime.strptime(dates[0], '%b %d %Y').date().isoformat()
        if first_public != source_record['first_public_date']:
            raise ValueError('Source-profile first-public date differs: ' + study)
        sample_text = Path(sources['samples']['path']).read_text()
        blocks = re.split(r'^\^SAMPLE = ', sample_text, flags=re.M)[1:]
        samples, issues = [], []
        for block in blocks:
            lines = block.splitlines()
            sample = lines[0].strip()
            if not re.fullmatch(r'GSM[0-9]+', sample):
                raise ValueError('Noncanonical SOFT sample identity: ' + study)
            # A later platform/series section cannot supply sample annotations.
            sample_body = '\n'.join(re.split(r'^\^', '\n'.join(lines[1:]), maxsplit=1, flags=re.M)[:1])
            fields = {}
            for key in ('organism_ch1', 'molecule_ch1'):
                fields[key] = re.findall(r'^!Sample_' + key + r' = (.*)$', sample_body, re.M)
                if len(fields[key]) != 1 or not fields[key][0].strip():
                    issues.append(sample + ': missing or ambiguous ' + key)
            molecule_counts.update(fields['molecule_ch1'])
            samples.append({'sample_id': sample, **fields})
        observed_ids = [r['sample_id'] for r in samples]
        if len(observed_ids) != len(set(observed_ids)) or set(observed_ids) != set(ids):
            raise ValueError('Source-profile sample SOFT membership differs: ' + study)
        organisms = {v for r in samples for v in r['organism_ch1']}
        molecules = {v for r in samples for v in r['molecule_ch1']}
        if organisms != {expected_species}:
            issues.append('Recorded organism does not match the declared modality')
        if {v.casefold() for v in molecules} != {'total rna'}:
            issues.append('Recorded molecule labels are outside the existing total-RNA intake profile')
        records.append({'study_id': study, 'modality': modality, 'sample_ids': ids,
            'sample_count': len(ids), 'first_public_date': first_public,
            'source_status_quote': '!Series_status = Public on ' + dates[0],
            'sources': sources, 'samples': samples,
            'molecule_label_counts': dict(sorted(Counter(v for r in samples for v in r['molecule_ch1']).items())),
            'profile_status': 'UNSUPPORTED_SOURCE_PROFILE' if issues else 'SUPPORTED_SOURCE_PROFILE',
            'profile': 'existing_total_rna_intake_v21', 'limitations': issues,
            'scientific_status': 'CANNOT VERIFY'})
    for name, expected in files.items():
        if sha(name) != expected:
            raise ValueError('Source-profile input changed during inspection: ' + name)
    if (bound_json(discovery_path, '.json.sha256') != summary
            or sha(Path(__file__)) != implementation_hash):
        raise ValueError('Source-profile summary or implementation changed during inspection')
    supported = [r['study_id'] for r in records if r['profile_status'] == 'SUPPORTED_SOURCE_PROFILE']
    stopped = [r['study_id'] for r in records if r['profile_status'] != 'SUPPORTED_SOURCE_PROFILE']
    return {'schema': 'clockbase-rna-source-profile-census-v1', 'status': 'PASS',
        'recorded_at_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Literal GEO assay annotations against the existing total-RNA intake profile; unsupported does not mean erroneous or scientifically unsuitable',
        'discovery': {'path': str(discovery_path), 'sha256': discovery_hash},
        'candidates': {'path': str(candidate_path), 'sha256': files[str(candidate_path)]},
        'profile': 'existing_total_rna_intake_v21', 'modality': modality, 'studies': records,
        'supported_study_ids': supported, 'stopped_study_ids': stopped,
        'counts': {'enrolled_studies': len(records), 'enrolled_samples': sum(r['sample_count'] for r in records),
                   'supported_studies': len(supported), 'stopped_studies': len(stopped),
                   'molecule_label_occurrences': dict(sorted(molecule_counts.items()))},
        'source_files': files, 'verification_script_sha256': implementation_hash,
        'model_calls': 0, 'expression_cells_read': 0, 'clock_runs': 0,
        'scientific_status': 'CANNOT VERIFY', 'original_paper_2025_eligibility': 'QUARANTINED',
        'elapsed_seconds': time.monotonic()-started}


def verify_batch_extraction(batch_path, release_path=None):
    """Read every derived CSV cell directly against HDF5, without producer imports."""
    started = time.monotonic()
    batch_path = Path(batch_path).resolve()
    batch_root = batch_path.parent
    batch_hash, code_hash = sha(batch_path), sha(Path(__file__))
    batch = bound_json(batch_path, '.json.sha256')
    if (batch.get('schema') != 'clockbase-archs4-batch-extraction-summary-v1'
            or batch.get('path_base') != 'batch_root'
            or batch.get('execution_status') != 'COMPLETED'
            or batch.get('scientific_status') != 'CANNOT VERIFY'
            or batch.get('original_paper_2025_eligibility') != 'QUARANTINED'
            or any(batch.get(key) != 0 for key in ('model_api_calls', 'preprocessing_runs', 'clock_runs'))):
        raise ValueError('Every-cell readback requires a completed scoped extraction batch')
    bound_files = {str(batch_path): batch_hash}

    def local(declaration, *, member=True):
        path = Path(declaration['path'])
        path = (batch_root / path).resolve() if not path.is_absolute() else path.resolve()
        digest = declaration['sha256']
        if member and (not path.is_relative_to(batch_root)
                       or batch['artifacts'].get(str(path.relative_to(batch_root))) != digest):
            raise ValueError('Batch declaration is not an exact archived member: ' + str(path))
        if sha(path) != digest:
            raise ValueError('Batch source or artifact hash differs: ' + str(path))
        bound_files[str(path)] = digest
        return path

    def stage_member(stage_path, stage, declaration):
        path = local(declaration)
        if (not path.is_relative_to(stage_path.parent)
                or stage['outputs'].get(str(path.relative_to(stage_path.parent))) != declaration['sha256']):
            raise ValueError('Declared artifact is not a sealed stage member: ' + str(path))
        return path

    def terminal_stage(study, phase, declaration=None):
        expected_path = batch_root / 'studies' / study / phase / 'complete.json'
        if declaration is None:
            name = str(expected_path.relative_to(batch_root))
            declaration = {'path': name, 'sha256': batch['artifacts'].get(name)}
        path = local(declaration)
        if path != expected_path:
            raise ValueError('Terminal stage path belongs to a different study or phase')
        stage = bound_json(path, '.sha256')
        members = [name for name in stage['outputs'] if Path(name).name == 'terminal-attempt.json']
        if len(members) != 1:
            raise ValueError('Terminal stage does not bind exactly one terminal attempt')
        member = path.parent / members[0]
        terminal_path = stage_member(path, stage, {'path': str(member), 'sha256': stage['outputs'][members[0]]})
        terminal = json.loads(terminal_path.read_text())
        if (terminal != stage['result'] or terminal.get('phase') != phase
                or terminal.get('automatic_retry') is not False):
            raise ValueError('Terminal stage result and retained attempt disagree')
        return path, stage, terminal

    for name, digest in batch['artifacts'].items():
        local({'path': name, 'sha256': digest})
    profile_path = local(batch['source_profile_descriptor'], member=False)
    profile = json.loads(profile_path.read_text())
    if (profile.get('schema') != 'clockbase-rna-source-profile-census-v1'
            or profile.get('status') != 'PASS'):
        raise ValueError('Batch source profile descriptor is unavailable')
    species_by_modality = {'human_rna': ('human', 'Homo sapiens'),
                           'mouse_rna': ('mouse', 'Mus musculus')}
    modality = profile.get('modality')
    if modality not in species_by_modality:
        raise ValueError('Batch source profile has no supported RNA modality')
    species, organism = species_by_modality[modality]
    # Archived human v1 summaries lack these additive fields. A mouse batch
    # must declare both; an absent field can never silently mean mouse.
    if (batch.get('modality', 'human_rna') != modality
            or batch['source_hdf5'].get('species', 'human') != species):
        raise ValueError('Batch/profile/release modality or species differs')
    config_path = local({'path': 'config.json', 'sha256': batch['config_sha256']}, member=False)
    config = bound_json(config_path, '.json.sha256')
    plan_path = local({'path': config['plan_path'], 'sha256': batch['plan_sha256']}, member=False)
    # Older batch plans were SHA-bound by config/summary but not sidecar-sealed.
    # The source hash above remains mandatory for both versions.
    plan = json.loads(plan_path.read_text())
    if (config.get('schema') != 'clockbase-archs4-batch-extraction-config-v1'
            or plan.get('schema') != 'clockbase-archs4-batch-extraction-plan-v1'
            or config.get('plan_sha256') != batch['plan_sha256']
            or not (config.get('modality', 'human_rna') == plan.get('modality', 'human_rna') == modality)):
        raise ValueError('Batch config/plan modality or identity differs')
    if (not (config.get('release', {}).get('species', 'human')
             == plan.get('release', {}).get('species', 'human') == species)
            or config.get('release', {}).get('sha256') != plan.get('release', {}).get('sha256')
            or config.get('release', {}).get('sha256') != batch['source_hdf5']['sha256']):
        raise ValueError('Batch config/plan release species or hash differs')
    profile_rows = profile['studies']
    profile_ids = [row['study_id'] for row in profile_rows]
    rows = batch['results']
    if (len(profile_ids) != len(set(profile_ids))
            or [row['study_id'] for row in rows] != profile_ids
            or len(rows) != batch['counts']['enrolled_studies']):
        raise ValueError('Batch result inventory differs from the complete source cohort')
    by_study = {row['study_id']: row for row in profile_rows}
    supported, stopped = profile['supported_study_ids'], profile['stopped_study_ids']
    if (len(set(supported)) != len(supported) or len(set(stopped)) != len(stopped)
            or set(supported) & set(stopped) or set(supported) | set(stopped) != set(profile_ids)):
        raise ValueError('Batch source-profile partition is incomplete or overlapping')
    status_counts = Counter(row['status'] for row in rows)
    if set(status_counts) - {'EXTRACTED', 'SOURCE_PROFILE_UNRESOLVED', 'EXTRACTION_UNRESOLVED'}:
        raise ValueError('Batch inventory has an unsupported terminal status')
    expected_counts = {
        'enrolled_studies': len(rows), 'selected_studies': len(supported),
        'source_profile_stops': status_counts['SOURCE_PROFILE_UNRESOLVED'],
        'extracted_studies': status_counts['EXTRACTED'],
        'extraction_unresolved_studies': status_counts['EXTRACTION_UNRESOLVED'],
        'enrolled_samples': sum(len(row['sample_ids']) for row in profile_rows),
        'selected_samples': sum(len(by_study[study]['sample_ids']) for study in supported),
        'extracted_samples': sum(row['matrix']['sample_count'] for row in rows if row['matrix'] is not None),
        'extracted_bytes': sum(row['matrix']['size_bytes'] for row in rows if row['matrix'] is not None)}
    expected_status = ('COMPLETED_WITH_UNRESOLVED' if any(
        row['status'] != 'EXTRACTED' for row in rows) else 'COMPLETED')
    if batch['counts'] != expected_counts or batch.get('status') != expected_status:
        raise ValueError('Batch counts or terminal summary differ from the retained inventory')
    checks = {}
    for role, phase in [('release_before', 'before'), ('release_after', 'after')]:
        path = local(batch[role])
        value = bound_json(path, '.json.sha256')
        if (value.get('phase') != phase or value.get('status') != 'PASS'
                or value.get('hash_scope') != 'FULL_HDF5'
                or value.get('sha256') != batch['source_hdf5']['sha256']
                or value.get('size_bytes') != batch['source_hdf5']['size_bytes']):
            raise ValueError('Batch full-file release check does not support its origin')
        checks[role] = value
    if checks['release_before']['path'] != checks['release_after']['path']:
        raise ValueError('Batch release locator differs between full-file checks')
    release = Path(release_path or checks['release_before']['path']).resolve()
    release_hash = sha(release)
    before = release.stat()
    if release_hash != batch['source_hdf5']['sha256'] or before.st_size != batch['source_hdf5']['size_bytes']:
        raise ValueError('Independent release hash differs from the batch origin')
    decode = lambda values: [value.decode() if isinstance(value, bytes) else str(value) for value in values]
    results = []
    with h5py.File(release, 'r') as handle:
        expression = handle['data/expression']
        # Current release counts are unsigned integers. Other source types need
        # an explicit verifier profile; no rounding or tolerance is introduced.
        if expression.dtype.kind != 'u' or expression.dtype.itemsize > 8:
            raise ValueError('Exact every-cell readback requires unsigned integer ARCHS4 counts')
        ids = decode(handle['meta/samples/geo_accession'][:])
        labels = decode(handle['meta/samples/series_id'][:])
        organisms = decode(handle['meta/samples/organism_ch1'][:])
        genes = decode(handle['meta/genes/ensembl_gene'][:])
        if (list(expression.shape) != batch['source_hdf5']['expression_shape']
                or str(expression.dtype) != batch['source_hdf5']['expression_dtype']
                or len(ids) != expression.shape[1] or len(genes) != expression.shape[0]
                or len(ids) != len(set(ids)) or len(genes) != len(set(genes))
                or len(labels) != len(ids) or len(organisms) != len(ids)):
            raise ValueError('Independent HDF5 metadata identities or dimensions differ')
        positions = {study: [] for study in profile_ids}
        wanted = set(profile_ids)
        for index, membership in enumerate(labels):
            for study in wanted.intersection(part.strip() for part in membership.split(',')):
                positions[study].append(index)
        for row in rows:
            study = row['study_id']
            source_profile = by_study[study]
            result = {'study_id': study, 'modality': row['modality'], 'status': 'CANNOT VERIFY',
                      'scientific_status': 'CANNOT VERIFY', 'cells_compared': 0}
            try:
                record_hash = hashlib.sha256(json.dumps(source_profile, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
                if record_hash != row['source_profile_record_sha256']:
                    raise ValueError('Per-study source profile hash differs')
                if (row['modality'] != modality or source_profile['modality'] != modality
                        or row['source_profile_status'] != source_profile['profile_status']):
                    raise ValueError('Per-study modality or source-profile status differs')
                for role in ('series', 'samples'):
                    path = local(row[role], member=False)
                    declared = source_profile['sources'][role]
                    if path != Path(declared['path']).resolve() or sha(path) != declared['sha256']:
                        raise ValueError('Per-study GEO source differs from the profile')
                if study in stopped:
                    if (row['status'] != 'SOURCE_PROFILE_UNRESOLVED' or row['matrix'] is not None
                            or row['extraction_receipt'] is not None
                            or source_profile['profile_status'] != 'UNSUPPORTED_SOURCE_PROFILE'):
                        raise ValueError('Unsupported source profile was omitted or converted to extraction')
                    stop = json.loads(local(row['stop_receipt']).read_text())
                    if (stop.get('sample_ids') != source_profile['sample_ids']
                            or stop.get('study_id') != study
                            or stop.get('modality') != modality
                            or stop.get('species', 'human') != species
                            or stop.get('status') != 'SOURCE_PROFILE_UNRESOLVED'
                            or stop.get('source_profile_descriptor_sha256') != sha(profile_path)
                            or stop.get('source_profile_record_sha256') != record_hash
                            or stop.get('scientific_status') != 'CANNOT VERIFY'):
                        raise ValueError('Source-profile stop changed its samples or scientific scope')
                    stage_path = local(row['stage_receipt'])
                    stage = bound_json(stage_path, '.sha256')
                    stop_path = local(row['stop_receipt'])
                    if (not stop_path.is_relative_to(stage_path.parent)
                            or stage['outputs'].get(str(stop_path.relative_to(stage_path.parent)))
                            != row['stop_receipt']['sha256']):
                        raise ValueError('Source-profile stop is not a sealed stage member')
                    result.update(reason='Unsupported existing source profile; no expression read for this study',
                                  source_profile_status=source_profile['profile_status'],
                                  sample_ids=source_profile['sample_ids'])
                    results.append(result)
                    continue
                if row['status'] == 'EXTRACTION_UNRESOLVED':
                    if row['matrix'] is not None or row['extraction_receipt'] is not None:
                        raise ValueError('Unresolved extraction must not claim a derived matrix')
                    _, _, terminal = terminal_stage(study, 'extract', row['stage_receipt'])
                    if (terminal.get('status') != 'CANNOT VERIFY' or 'value' in terminal
                            or not isinstance(terminal.get('error'), str)
                            or terminal['error'] != row.get('error')):
                        raise ValueError('Unresolved extraction disagrees with its terminal result')
                    result.update(reason=row.get('error', 'Extraction did not produce a verified matrix'))
                    results.append(result)
                    continue
                if row['status'] != 'EXTRACTED' or source_profile['profile_status'] != 'SUPPORTED_SOURCE_PROFILE':
                    raise ValueError('Unexpected extracted-study or source-profile status')
                counts = local(row['matrix'])
                extraction_path = local(row['extraction_receipt'])
                extraction = json.loads(extraction_path.read_text())
                if (extraction.get('schema') != 'clockbase-archs4-batch-study-extraction-v1'
                        or extraction.get('path_base') != 'batch_root'
                        or extraction.get('study_id') != study or extraction.get('status') != 'PASS'
                        or extraction.get('modality') != modality
                        or extraction.get('species') != species
                        or extraction.get('scientific_status') != 'CANNOT VERIFY'
                        or extraction['matrix'] != row['matrix']
                        or extraction['source']['release_sha256'] != release_hash
                        or extraction['source']['release_size_bytes'] != before.st_size
                        or extraction['source']['source_profile_record_sha256'] != record_hash
                        or any(extraction['source'][role + '_sha256'] != row[role]['sha256']
                               for role in ('series', 'samples'))
                        or extraction['source']['source_profile_descriptor_sha256'] != sha(profile_path)):
                    raise ValueError('Extraction receipt does not bind this matrix, profile and origin')
                stage_path, stage, terminal = terminal_stage(study, 'extract', row['stage_receipt'])
                value = terminal.get('value') or {}
                if (terminal.get('status') != 'PASS' or value.get('matrix') != row['matrix']
                        or value.get('receipt_path') != row['extraction_receipt']['path']
                        or value.get('receipt_sha256') != row['extraction_receipt']['sha256']
                        or value.get('sample_ids') != extraction['sample_ids_in_matrix_order']):
                    raise ValueError('Extracted matrix disagrees with its successful terminal result')
                for path, digest in [(counts, row['matrix']['sha256']),
                                     (extraction_path, row['extraction_receipt']['sha256'])]:
                    stage_member(stage_path, stage, {'path': str(path), 'sha256': digest})
                metadata_path, metadata_stage, metadata_terminal = terminal_stage(study, 'metadata')
                metadata_value = metadata_terminal.get('value') or {}
                if metadata_terminal.get('status') != 'PASS':
                    raise ValueError('Extracted matrix lacks a successful metadata stage')
                for role, field in [('metadata', 'metadata'), ('context', 'study_context'),
                                    ('series', 'series_source'), ('samples', 'source_sample_records')]:
                    stage_member(metadata_path, metadata_stage, {
                        'path': metadata_value[field], 'sha256': extraction['source'][role + '_sha256']})
                indices = positions[study]
                source_ids = [ids[index] for index in indices]
                full_memberships = {ids[index]: [part.strip() for part in labels[index].split(',')]
                                    for index in indices}
                if (set(source_ids) != set(source_profile['sample_ids'])
                        or len(source_ids) != len(source_profile['sample_ids'])
                        or {organisms[index] for index in indices} != {organism}
                        or extraction['sample_ids_in_matrix_order'] != source_ids
                        or extraction.get('all_study_memberships') != full_memberships):
                    raise ValueError('Full HDF5 study membership, order or organism differs')
                cells = 0
                with counts.open(newline='') as csv_handle:
                    reader = csv.reader(csv_handle, strict=True)
                    header = next(reader)
                    if header != ['ensembl_gene', *source_ids]:
                        raise ValueError('Literal CSV sample order or header differs')
                    for start in range(0, len(genes), 2048):
                        end = min(len(genes), start + 2048)
                        block = [next(reader) for _ in range(end-start)]
                        if (any(len(record) != len(header) for record in block)
                                or [record[0] for record in block] != genes[start:end]):
                            raise ValueError('Literal CSV feature identity, order or width differs')
                        if any(re.fullmatch(r'[0-9]+', value) is None
                               for record in block for value in record[1:]):
                            raise ValueError('Derived CSV count is not an unsigned integer literal')
                        observed = np.asarray([record[1:] for record in block], dtype=np.uint64)
                        expected = expression[start:end, indices]
                        if observed.shape != expected.shape or not np.array_equal(observed, expected):
                            raise ValueError('Derived CSV count cell differs from the HDF5 source')
                        cells += observed.size
                    if next(reader, None) is not None:
                        raise ValueError('Derived CSV contains additional feature records')
                if (row['matrix']['feature_count'] != len(genes)
                        or row['matrix']['sample_count'] != len(source_ids)
                        or row['matrix']['cell_count'] != cells
                        or row['matrix']['size_bytes'] != counts.stat().st_size):
                    raise ValueError('Derived matrix dimensions or byte count differ')
                result.update(status='PASS', matrix_sha256=sha(counts),
                    extraction_receipt_sha256=sha(extraction_path), feature_count=len(genes),
                    sample_count=len(source_ids), sample_ids=source_ids, cells_compared=int(cells),
                    feature_ids_equal=True, sample_ids_equal=True, all_cells_equal=True)
            except (OSError, ValueError, TypeError, KeyError, StopIteration, OverflowError, csv.Error) as error:
                result.update(status='FAIL', reason=type(error).__name__ + ': ' + str(error))
            results.append(result)
    final_hash = sha(release)
    after = release.stat()
    if (final_hash != release_hash
            or (before.st_size, before.st_ino, before.st_mtime_ns, before.st_ctime_ns)
            != (after.st_size, after.st_ino, after.st_mtime_ns, after.st_ctime_ns)):
        raise ValueError('HDF5 source changed during independent every-cell verification')
    for path, expected in bound_files.items():
        if sha(path) != expected:
            raise ValueError('Bound input changed during independent every-cell verification: ' + path)
    if bound_json(batch_path, '.json.sha256') != batch or sha(Path(__file__)) != code_hash:
        raise ValueError('Batch summary or independent verifier changed during readback')
    status_counts = dict(sorted(Counter(row['status'] for row in results).items()))
    return {'schema': 'clockbase-archs4-independent-every-cell-v1',
        'modality': modality, 'species': species,
        'status': 'FAIL' if status_counts.get('FAIL') else 'PASS',
        'scope': 'Exact derived-count source-cell readback and complete unresolved-case accounting; no preprocessing or clock validation',
        'recorded_at_utc': datetime.now(timezone.utc).isoformat(),
        'batch_summary': {'path': str(batch_path), 'sha256': batch_hash}, 'batch_summary_sha256': batch_hash,
        'source_profile_sha256': sha(profile_path), 'release_sha256': release_hash,
        'release_sha256_before': release_hash, 'release_sha256_after': final_hash,
        'release_path': str(release), 'release_origin_path': checks['release_before']['path'],
        'verification_script_sha256': code_hash,
        'verification_script': {'path': str(Path(__file__).resolve()), 'sha256': code_hash},
        'results': results, 'study_status_counts': status_counts,
        'enrolled_studies': len(results), 'cells_compared': sum(row['cells_compared'] for row in results),
        'source_files': bound_files, 'model_calls': 0, 'preprocessing_runs': 0, 'clock_runs': 0,
        'scientific_status': 'CANNOT VERIFY', 'original_paper_2025_eligibility': 'QUARANTINED',
        'elapsed_seconds': time.monotonic()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--execution', type=Path)
    mode.add_argument('--discovery', type=Path)
    mode.add_argument('--profile-discovery', type=Path)
    mode.add_argument('--batch-extraction', type=Path)
    parser.add_argument('--release', type=Path)
    parser.add_argument('--release-sha256')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.output.exists():
            raise FileExistsError('Existing verification output is preserved: ' + str(args.output))
        if args.discovery:
            if not args.release or not args.release_sha256:
                raise ValueError('--discovery requires --release and --release-sha256')
            result = discovery_availability(args.discovery, args.release, args.release_sha256)
        elif args.batch_extraction:
            if args.release_sha256:
                raise ValueError('Batch origin SHA256 comes from its sealed source receipts')
            result = verify_batch_extraction(args.batch_extraction, args.release)
        else:
            if args.release or args.release_sha256:
                raise ValueError('Release options apply only to --discovery')
            result = (discovery_source_profiles(args.profile_discovery)
                      if args.profile_discovery else verify(args.execution.resolve()))
        with args.output.open('x') as handle:
            json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write('\n')
    except Exception as error:
        print(type(error).__name__ + ': ' + str(error), file=sys.stderr)
        return 1
    print(json.dumps(result))
    return 0


if __name__ == '__main__':
    sys.exit(main())
