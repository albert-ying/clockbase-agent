"""Prospective comparison with source references frozen before numerical work.

The controller reuses the original CodingAgent and four existing comparison
methods. Enrollment includes failed attempts. Reference judgments never enter
the producer or method prompts. A committed joint release is required even for
the development exercise of this route.
"""
import argparse
import fcntl
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

from .agent_workflow_2026_09_10 import sha, source_packet, write_new
from src.pipeline_runtime_2026_09_10 import BoundedOpenAI, StageStore, PipelineIntegrityError

ROOT = Path(__file__).resolve().parents[1]
PLAN_SCHEMA = 'clockbase-prospective-evaluation-plan-v1'
CONFIG_SCHEMA = 'clockbase-prospective-comparison-v1'
RELEASE_SCHEMA = 'clockbase-prospective-release-v1'
UNAVAILABLE = 'UNAVAILABLE_PROSPECTIVE_PREPARATION'
ANALYSIS = {'primary_endpoint': 'accepted_result_error_rate', 'study_weighting': 'equal',
            'uncertainty': 'paired_family_cluster_percentile_bootstrap', 'resamples': 2000, 'seed': 20260910}
IMPLEMENTATIONS = tuple('audit_system/' + name for name in (
    'heldout_2026_09_17.py', 'evaluation_2026_09_10.py', 'agent_workflow_2026_09_10.py',
    'report_2026_09_10.py', 'source_packet.py', 'adapters/measurement_pilot.py',
    'claim_gate.py', 'policy.py', 'models.py', 'ledger.py', 'hashing.py',
    'dependency_graph.py', 'extension_eligibility_2026_09_10.py')) + tuple(
    'src/' + name for name in ('pipeline_2026_09_10.py', 'pipeline_inputs_2026_09_10.py',
    'pipeline_runtime_2026_09_10.py', 'coding_agent.py', 'detailed_trace_inspection.py'))


IntegrityError = PipelineIntegrityError


class ProspectiveStore(StageStore):
    def run(self, name, function):
        folder = self.root / name
        if not (folder / 'complete.json').exists() and any(folder.glob('attempt-*')):
            raise IntegrityError('Interrupted stage cannot retry automatically; preserve it and use an explicit recovery run: ' + name)
        def with_api_evidence(directory):
            result = function(directory)
            if _api_patterns(name) is not None:
                write_new(directory / 'api-ledger-identities.json', _api_inventory(self.root, name))
            return result
        result = super().run(name, with_api_evidence)
        if _api_patterns(name) is not None:
            _verify_stage_api(self.root, name)
        return result


def _api_patterns(name):
    if re.fullmatch(r'GSE\d+/coding_agent', name):
        return ['baseline-api/' + name.split('/')[0] + '/**/*']
    if re.fullmatch(r'GSE\d+/(independent_audit|self_review)/review', name):
        study, method, _ = name.split('/')
        return ['api/' + study + '-' + method + '-*/**/*']
    if name == '99_summary': return ['api/**/*', 'baseline-api/**/*']
    return None


def _api_inventory(archive, name):
    return {str(path.relative_to(archive)): sha(path)
            for pattern in _api_patterns(name) or []
            for path in sorted(Path(archive).glob(pattern)) if path.is_file()}


def _verify_stage_api(archive, name):
    from .evaluation_2026_09_10 import reference_stage
    identities, _ = reference_stage(archive, name, 'api-ledger-identities.json')
    if identities != _api_inventory(Path(archive), name):
        raise IntegrityError('API request, response, reservation or usage evidence changed: ' + name)


def _verify_completed_api(archive):
    archive = Path(archive)
    for pattern in ('GSE*/coding_agent/complete.json', 'GSE*/*/review/complete.json', '99_summary/complete.json'):
        for path in archive.glob(pattern):
            name = str(path.parent.relative_to(archive))
            if _api_patterns(name) is not None: _verify_stage_api(archive, name)


def read_json(path):
    path = Path(path)
    try:
        with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as handle:
            def nonfinite(value): raise ValueError('Nonfinite JSON constant: ' + value)
            return json.load(handle, parse_constant=nonfinite)
    except (OSError, ValueError) as error:
        raise IntegrityError(f'Cannot read declared JSON artifact {path}: {error}') from error


def checked_path(item):
    if (not isinstance(item, dict) or not isinstance(item.get('path'), str)
            or not re.fullmatch(r'[0-9a-f]{64}', str(item.get('sha256', '')))):
        raise IntegrityError('Artifact requires a path and SHA256')
    path = Path(item['path']).resolve()
    if not path.is_file():
        raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
    if sha(path) != item['sha256']:
        raise IntegrityError('Frozen artifact changed: ' + str(path))
    return path


def git_payload(path, revision):
    path = Path(path).resolve()
    if not path.is_relative_to(ROOT):
        raise IntegrityError('Joint freeze artifact must be inside the repository: ' + str(path))
    process = subprocess.run(['git', 'show', revision + ':' + str(path.relative_to(ROOT))],
                             cwd=ROOT, capture_output=True)
    if process.returncode:
        raise IntegrityError('Joint freeze artifact is not committed: ' + str(path))
    return process.stdout


def require_committed(path, revision):
    if hashlib.sha256(git_payload(path, revision)).hexdigest() != sha(path):
        raise IntegrityError('Working artifact differs from its committed freeze: ' + str(path))


def _limits(value, calls):
    keys = {'max_calls', 'max_cost', 'output_tokens', 'prompt_bytes',
            'input_per_million', 'output_per_million'}
    if not isinstance(value, dict) or set(value) != keys or value['max_calls'] != calls:
        raise IntegrityError('Explicit per-actor call/cost/token/byte/tariff limits required')
    if any(type(v) not in (int, float) or not math.isfinite(v) or v <= 0 for v in value.values()):
        raise IntegrityError('Resource limits must be positive and finite')
    if any(type(value[k]) is not int for k in ('max_calls', 'output_tokens', 'prompt_bytes')):
        raise IntegrityError('Call, token and byte ceilings must be integers')
    return value


def validate_release(release_path):
    """Check the joint freeze before reading scientific sources or API access."""
    from .evaluation_2026_09_10 import COMPARISON_METHODS, verify_family_enrollment
    release_path = Path(release_path).resolve(); release = read_json(release_path)
    if (release.get('schema') != RELEASE_SCHEMA or release.get('state') != 'OPEN'
            or not re.fullmatch(r'[0-9a-f]{40}', str(release.get('freeze_commit', '')))):
        raise IntegrityError('A committed OPEN joint release is required')
    require_committed(release_path, 'HEAD')
    commit = release['freeze_commit']
    check = subprocess.run(['git', 'merge-base', '--is-ancestor', commit, 'HEAD'], cwd=ROOT)
    if check.returncode:
        raise IntegrityError('Joint freeze is not an ancestor of the current checkout')
    plan_path = checked_path(release['plan']); require_committed(plan_path, commit)
    plan = read_json(plan_path)
    if plan.get('schema') != PLAN_SCHEMA or plan.get('split') not in {'development', 'held_out'}:
        raise IntegrityError('Unrecognized prospective comparison plan')
    if (plan.get('methods') != list(COMPARISON_METHODS)
            or plan.get('baseline_rule') != 'retain_producer_RESULT_rows_without_audit'
            or plan.get('preprocessing') != 'clockbase_legacy'
            or plan.get('producer_rule') != 'one_fresh_original_CodingAgent_per_ready_study'
            or plan.get('comparison_scope') != 'source_defined_interventions'
            or not isinstance(plan.get('target'), str) or not plan['target'].strip()):
        raise IntegrityError('Primary methods, fresh producer and unchanged legacy recipe must be frozen')
    if (plan.get('provider') != 'openai' or plan.get('model') != 'gpt-5.6-sol'
            or plan.get('concurrency') != 1 or plan.get('automatic_retries') != 0):
        raise IntegrityError('No provider/model substitution or automatic retry is allowed')
    dependencies = {}
    for key in ('protocol', 'family_manifest', 'exposure_manifest', 'cohort_manifest', 'asset_manifest'):
        path = checked_path(plan[key]); require_committed(path, commit); dependencies[key] = path
    if dependencies['asset_manifest'] != (ROOT / 'pipeline_assets/rna_4_6_2026_09_10/manifest.json').resolve():
        raise IntegrityError('Asset freeze must name the executed pipeline manifest')
    if set(plan.get('implementation_sha256', {})) != set(IMPLEMENTATIONS):
        raise IntegrityError('All executed producer/audit/controller implementations must be frozen')
    for relative, digest in plan['implementation_sha256'].items():
        path = checked_path({'path': str(ROOT / relative), 'sha256': digest})
        require_committed(path, commit)
    families = read_json(dependencies['family_manifest'])
    cohort = read_json(dependencies['cohort_manifest'])
    if (cohort.get('schema') != 'clockbase-prospective-cohort-v1'
            or cohort.get('split') != plan['split'] or not isinstance(cohort.get('cases'), list)
            or not cohort['cases']):
        raise IntegrityError('A nonempty prior cohort is required')
    cases = cohort['cases']; ids = [row.get('study_id') for row in cases]
    if (any(not isinstance(x, str) or not re.fullmatch(r'GSE\d+', x) for x in ids)
            or len(ids) != len(set(ids))):
        raise IntegrityError('Cohort contains missing or duplicate study identifiers')
    molecular_identities = {}
    for row in cases:
        if row.get('modality') not in {'human_rna', 'mouse_rna'}:
            raise IntegrityError('This prospective cohort uses the supported RNA modalities only')
        if row.get('metadata_source', {}).get('kind') not in {'geo', 'pinned'}:
            raise IntegrityError('Source acquisition strategy must be frozen')
        if row.get('molecular_source', {}).get('kind') not in {'archs4', 'raw_counts'}:
            raise IntegrityError('Molecular input identity must be frozen')
        molecular = row['molecular_source']
        if (not isinstance(molecular.get('path'), str) or not Path(molecular['path']).is_absolute()
                or not re.fullmatch(r'[0-9a-f]{64}', str(molecular.get('sha256', '')))):
            raise IntegrityError('Molecular source requires an absolute path and SHA256 before source review')
        if molecular['kind'] == 'archs4' and any(type(molecular.get(k)) is not int or molecular[k] < 1
                                               for k in ('size_bytes', 'mtime_ns')):
            raise IntegrityError('ARCHS4 byte size and modification identity must be frozen')
        key = str(Path(molecular['path']).resolve())
        if key in molecular_identities and molecular_identities[key] != molecular:
            raise IntegrityError('The same molecular source has conflicting frozen identities')
        molecular_identities[key] = molecular
        metadata = row['metadata_source']
        if metadata['kind'] == 'pinned':
            record, identities = metadata.get('record', {}), metadata.get('source_sha256', {})
            if (not isinstance(record, dict) or not isinstance(identities, dict)
                    or any(not isinstance(record.get(k), str) or record[k] not in identities
                           for k in ('metadata','source_sample_records','series_source','study_context'))
                    or any(not isinstance(p, str) or not Path(p).is_absolute()
                           or not re.fullmatch(r'[0-9a-f]{64}', str(d)) for p, d in identities.items())):
                raise IntegrityError('Pinned metadata requires all primary/context identities before source access')
        linkage = row.get('expected_linkage', {})
        if (set(linkage) != {'sample_ids', 'pubmed_ids', 'relations'}
                or any(not isinstance(values, list) or any(not isinstance(x, str) for x in values)
                       or len(values) != len(set(values)) for values in linkage.values())
                or not linkage['sample_ids']
                or any(not re.fullmatch(r'GSM\d+', x) for x in linkage['sample_ids'])):
            raise IntegrityError('Frozen primary sample/publication/relation linkage is required')
    enrollment = {'schema': 'clockbase-comparison-enrollment-v1', 'split': plan['split'],
                  'methods': list(COMPARISON_METHODS), 'studies': [
                      {'study_id': row['study_id'], 'family_id': row['family_id']} for row in cases]}
    verify_family_enrollment(enrollment['studies'], families, plan['split'])
    exposed = read_json(dependencies['exposure_manifest'])
    if not isinstance(exposed.get('exposed_studies'), list):
        raise IntegrityError('Explicit development exposure registry is required')
    if not set(exposed['exposed_studies']).issubset(set(families['exposed_studies'])):
        raise IntegrityError('Source-family freeze omits a known exposed study')
    primary_limits = _limits(plan['primary_limits'], 6)
    review_limits = _limits(plan['review_limits'], 1)
    budget = plan.get('budget', {})
    allocation = len(cases) * (primary_limits['max_cost'] + 2 * review_limits['max_cost'])
    if (budget.get('initial_authorized_cap_usd') != 5000
            or type(budget.get('previously_allocated_usd')) not in (int, float)
            or not math.isfinite(budget['previously_allocated_usd'])
            or type(budget.get('this_allocation_usd')) not in (int, float)
            or not math.isfinite(budget['this_allocation_usd']) or not math.isfinite(allocation)
            or budget['previously_allocated_usd'] < 2
            or budget.get('this_allocation_usd') != allocation
            or budget['previously_allocated_usd'] + allocation > 5000):
        raise IntegrityError('Per-study allocations must fit the remaining initial USD5,000 cap')
    if plan.get('reference_procedure', {}).get('kind') != 'separate_source_first_agent':
        raise IntegrityError('Independent source-first reference procedure must be declared')
    if plan.get('analysis') != ANALYSIS:
        raise IntegrityError('Primary endpoint, study weighting and paired-family uncertainty must match the frozen analysis')
    if type(plan.get('max_samples_per_study')) is not int or plan['max_samples_per_study'] < 1:
        raise IntegrityError('A positive common preparation resource ceiling is required')
    return {'release': release, 'release_path': str(release_path), 'release_sha256': sha(release_path),
            'plan': plan, 'plan_path': str(plan_path), 'plan_sha256': sha(plan_path),
            'split': plan['split'], 'enrollment': enrollment, 'cases': cases,
            'family_manifest': families, 'family_manifest_path': str(dependencies['family_manifest'])}


def archive_config(context):
    return {'schema': CONFIG_SCHEMA, 'split': context['split'],
            'development_only': context['split'] == 'development',
            'release_path': context['release_path'], 'release_sha256': context['release_sha256'],
            'plan_path': context['plan_path'], 'plan_sha256': context['plan_sha256'],
            'implementation_sha256': context['plan']['implementation_sha256'],
            'provider': context['plan']['provider'], 'model': context['plan']['model'],
            'reference_labels_supplied_to_methods': False}


def validate_archive_context(archive):
    from .evaluation_2026_09_10 import reference_stage
    archive = Path(archive).resolve(); config = read_json(archive / 'config.json')
    if config.get('schema') != CONFIG_SCHEMA:
        raise IntegrityError('Explicit prospective archive schema required')
    context = validate_release(config['release_path'])
    if config != archive_config(context):
        raise IntegrityError('Prospective archive differs from its joint release')
    enrollment, result = reference_stage(archive, '00_enrollment', 'enrollment.json')
    if (enrollment != context['enrollment'] or result.get('enrollment') != enrollment
            or result.get('sha256') != sha(Path(result['path']))):
        raise IntegrityError('Enrolled denominator differs from its joint freeze')
    _verify_completed_api(archive)
    return context


def validate_reference_packet_context(packet):
    from .evaluation_2026_09_10 import reference_stage, reference_source_view
    if packet.get('schema') != 'clockbase-reference-design-packet-v2' or packet.get('phase') != 'SOURCE_ONLY':
        raise IntegrityError('Source-first reference packet v2 required')
    declared = packet.get('evaluation_context', {})
    archive = Path(declared.get('archive', '')).resolve()
    if sha(archive / 'config.json') != declared.get('archive_config_sha256'):
        raise IntegrityError('Source-only packet names a changed evaluation archive')
    context = validate_archive_context(archive)
    if packet.get('split') != context['split']:
        raise IntegrityError('Source-only packet names another cohort split')
    selected = [x for x in context['cases'] if x['study_id'] == packet.get('study_id')]
    if len(selected) != 1 or selected[0]['family_id'] != packet.get('family_id'):
        raise IntegrityError('Source-only packet is outside frozen enrollment')
    name = packet['study_id'] + '/source_packet'
    if declared.get('source_stage') != name:
        raise IntegrityError('Source-only packet names another source stage')
    envelope, receipt = reference_stage(archive, name, 'source-envelope.json')
    if (envelope != receipt or envelope.get('status') != 'PASS'
            or packet['source_packet'] != reference_source_view(envelope['packet'])):
        raise IntegrityError('Reference packet differs from its frozen source-only evidence')
    _verify_sources(envelope['source_sha256'])
    return context


def _verify_sources(identities):
    for path, digest in identities.items():
        checked_path({'path': path, 'sha256': digest})


def _stage(store, name, function):
    (store.root / name).parent.mkdir(parents=True, exist_ok=True)
    return store.run(name, function)


def terminal_stage(store, name, function):
    """A failed attempt becomes immutable evidence, never an automatic retry."""
    def execute(directory):
        try:
            value = function(directory)
            envelope = {'status': 'PASS', 'phase': name, 'value': value}
        except IntegrityError:
            raise
        except Exception as error:
            # Transport messages from BoundedOpenAI suppress response bodies/keys.
            envelope = {'status': 'CANNOT VERIFY', 'phase': name,
                        'error': f'{type(error).__name__}: {error}', 'automatic_retry': False}
        path = directory / 'terminal-attempt.json'; write_new(path, envelope)
        return {**envelope, 'terminal_receipt': str(path), 'terminal_sha256': sha(path)}
    return _stage(store, name, execute)


def _pinned_metadata(source, directory, study):
    metadata = dict(source['record'])
    _verify_sources(source['source_sha256'])
    for key in ('metadata', 'source_sample_records', 'series_source', 'study_context'):
        if metadata.get(key) not in source['source_sha256']:
            raise IntegrityError('Pinned metadata omits a source identity: ' + key)
    if metadata.get('accession', metadata.get('study_id', study)) != study:
        raise IntegrityError('Pinned metadata names another study')
    write_new(directory / 'source-membership.json', metadata)
    return metadata


def prepare_sources(store, context):
    from .evaluation_2026_09_10 import reference_source_view, REFERENCE_DESIGN_DOMAINS
    from src.pipeline_inputs_2026_09_10 import acquire_metadata
    plan = context['plan']; results = []
    for case in context['cases']:
        study = case['study_id']; source = case['metadata_source']
        def acquire(directory):
            return (_pinned_metadata(source, directory, study) if source['kind'] == 'pinned' else
                    acquire_metadata(study, directory, publication_methods=True))
        metadata = terminal_stage(store, study + '/metadata', acquire)
        def packet_stage(directory):
            if metadata['status'] != 'PASS':
                value = {'status': 'CANNOT VERIFY', 'study_id': study, 'family_id': case['family_id'],
                         'reason': 'Primary metadata unavailable', 'metadata_attempt': metadata}
            else:
                data = metadata['value']
                try:
                    if any(set(data.get(key, [])) != set(expected)
                           for key, expected in case['expected_linkage'].items()):
                        raise ValueError('Current primary linkage differs from the frozen family evidence')
                    packet = source_packet(Path(data['metadata']), study, plan['target'],
                        context=data['study_context'], comparison_scope=plan['comparison_scope'],
                        sample_source=data['source_sample_records'])
                    identities = {str(Path(data[key]).resolve()): sha(data[key]) for key in
                                  ('metadata', 'source_sample_records', 'series_source', 'study_context')}
                    for item in read_json(data['study_context'])['sources']:
                        path = (Path(data['study_context']).parent / item['path']).resolve()
                        if sha(path) != item['sha256']:
                            raise IntegrityError('Study-context source changed: ' + str(path))
                        identities[str(path)] = item['sha256']
                    value = {'status': 'PASS', 'study_id': study, 'family_id': case['family_id'],
                             'packet': packet, 'metadata': data, 'source_sha256': identities}
                except IntegrityError:
                    raise
                except Exception as error:
                    value = {'status': 'CANNOT VERIFY', 'study_id': study, 'family_id': case['family_id'],
                             'reason': f'Source packet failed: {type(error).__name__}: {error}',
                             'metadata_attempt': metadata}
            write_new(directory / 'source-envelope.json', value)
            return value
        envelope = _stage(store, study + '/source_packet', packet_stage)
        if envelope['status'] != 'PASS':
            results.append(envelope); continue
        _verify_sources(envelope['source_sha256'])
        def forms(directory):
            packet = {'schema': 'clockbase-reference-design-packet-v2', 'phase': 'SOURCE_ONLY',
                'study_id': study, 'family_id': case['family_id'], 'split': context['split'],
                'evaluation_context': {'archive': str(store.root),
                    'archive_config_sha256': sha(store.root / 'config.json'),
                    'source_stage': study + '/source_packet'},
                'source_packet': reference_source_view(envelope['packet'])}
            path = directory / 'source-design.json'; write_new(path, packet)
            form = {'schema': 'clockbase-reference-design-judgment-v1', 'completed': False,
                'packet_sha256': sha(path), 'study_id': study, 'family_id': case['family_id'],
                'reviewer': {'name': '', 'role': 'independent_source_adjudication',
                             'independent_from_methods': None},
                'domains': {name: {'status': 'CANNOT VERIFY', 'material': None,
                                  'rationale': '', 'evidence': []} for name in REFERENCE_DESIGN_DOMAINS},
                'comparisons': [], 'additional_sources_needed': []}
            write_new(directory / 'design-form.json', form)
            return {'study_id': study, 'family_id': case['family_id'], 'status': 'SOURCE_ONLY',
                    'source_packet': str(path), 'packet_sha256': sha(path),
                    'form': str(directory / 'design-form.json')}
        results.append(_stage(store, study + '/reference_packet', forms))
    return results


def freeze_source_states(store, context, state_manifest):
    """Require a completed or unavailable source reference before all numbers."""
    from .evaluation_2026_09_10 import read_reference_design, reference_stage
    states = read_json(state_manifest)
    if (states.get('schema') != 'clockbase-source-reference-states-v1'
            or states.get('archive_config_sha256') != sha(store.root / 'config.json')):
        raise IntegrityError('Source reference states must name the enrolled archive')
    supplied = states.get('studies', [])
    by_study = {row['study_id']: row for row in supplied}
    if len(by_study) != len(supplied) or set(by_study) != {x['study_id'] for x in context['cases']}:
        raise IntegrityError('Every enrolled study requires a source reference state before numbers')
    for case in context['cases']:
        study = case['study_id']; declared = by_study[study]
        references, missingness = [], []
        if declared.get('kind') == 'FROZEN':
            path = checked_path(declared['artifact']); frozen = read_reference_design(path)
            packet = frozen['packet']
            if (packet['study_id'] != study or packet['family_id'] != case['family_id']
                    or packet.get('evaluation_context', {}).get('archive') != str(store.root)):
                raise IntegrityError('Source reference belongs to another study or archive')
            references.append({'phase': 'SOURCE_DESIGN', 'path': str(path), 'sha256': sha(path)})
        elif declared.get('kind') == 'MISSING':
            path = checked_path(declared['artifact']); receipt = read_json(path)
            if (receipt.get('schema') != 'clockbase-reference-missingness-v1'
                    or receipt.get('phase') != 'SOURCE_DESIGN' or receipt.get('study_id') != study
                    or receipt.get('family_id') != case['family_id']
                    or receipt.get('archive_config_sha256') != sha(store.root / 'config.json')
                    or not isinstance(receipt.get('reason'), str) or not receipt['reason'].strip()
                    or not receipt.get('evidence')):
                raise IntegrityError('Unavailable source reference requires a pinned administrative receipt')
            for item in receipt['evidence']: checked_path(item)
            missingness.append({'phase': 'SOURCE_DESIGN', 'path': str(path), 'sha256': sha(path)})
        else:
            raise IntegrityError('Source reference is neither frozen nor explicitly unavailable')
        expected = {'study_id': study, 'references': references, 'missingness': missingness}
        def save(directory):
            write_new(directory / 'reference-state.json', expected)
            return expected
        observed = _stage(store, study + '/reference_state', save)
        if observed != expected:
            raise IntegrityError('Source reference changed after its pre-numerical freeze')
    return {'source_reference_states_complete': True,
            'manifest': str(Path(state_manifest).resolve()), 'sha256': sha(state_manifest)}


class ActorClient:
    """Use one shared review ledger, with additional strict actor ceilings."""
    def __init__(self, directory, actor, *, limits, total_calls, total_cost, model):
        self.directory, self.actor, self.limit = Path(directory), actor, limits['max_calls']
        self.actor_cost = limits['max_cost']; self.calls = 0
        self.client = BoundedOpenAI(directory, actor, model=model,
            **{**limits, 'max_calls': total_calls, 'max_cost': total_cost})

    def invoke(self, messages):
        if self.calls >= self.limit:
            raise RuntimeError('Per-actor model turn limit reached before transport')
        # Existing slots are rechecked by BoundedOpenAI; no failed request is retried.
        prior = list(self.directory.glob(self.actor + '-*/reservation.json'))
        charged = sum(read_json(p.parent / 'receipt.json')['estimated_cost_usd']
                      if (p.parent / 'receipt.json').exists() else read_json(p)['reserved_usd'] for p in prior)
        wire_bytes = len(json.dumps([{'role': m.type, 'content': m.content} for m in messages],
                                    ensure_ascii=False).encode())
        # Include conservative framing and transport keys; BoundedOpenAI also
        # reserves its exact request under the common ledger before transport.
        reserve = ((wire_bytes + 2048 + 100 * len(messages)) * self.client.input_rate
                   + self.client.output_tokens * self.client.output_rate) / 1e6
        slot = self.directory / (self.actor + '-%03d' % self.calls)
        if not slot.exists() and (len(prior) >= self.limit or charged + reserve > self.actor_cost):
            raise RuntimeError('Per-actor allocation exhausted before transport')
        self.calls += 1
        return self.client.invoke(messages)


def _input_path(source):
    path = Path(source['path']).resolve()
    if not path.is_file():
        raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
    if source['kind'] == 'raw_counts':
        return checked_path(source)
    # Full HDF5 hashes are verified once in the shared input stage. Every
    # extraction then checks size/mtime and hashes its selected count matrix.
    if (path.stat().st_size != source['size_bytes'] or path.stat().st_mtime_ns != source['mtime_ns']):
        raise IntegrityError('Pinned ARCHS4 file identity changed: ' + str(path))
    return path


def prepare_case(store, context, case):
    from .evaluation_2026_09_10 import reference_stage
    from src.pipeline_inputs_2026_09_10 import read_matrix, extract_rna, preprocess_rna, predict_rna
    from src.pipeline_2026_09_10 import primary_analysis
    study = case['study_id']; plan = context['plan']
    envelope, _ = reference_stage(store.root, study + '/source_packet', 'source-envelope.json')
    enrollment, er = reference_stage(store.root, '00_enrollment', 'enrollment.json')
    def failure(attempt):
        receipt = Path(attempt['terminal_receipt'])
        config = store.root / 'config.json'
        record = {'schema': 'clockbase-prospective-unavailable-attempt-v1',
            'enrollment_sha256': er['sha256'], 'phase': attempt['phase'],
            'failure_receipt': str(receipt), 'failure_sha256': sha(receipt),
            'recorded_failure': read_json(receipt), 'upstream_config': str(config),
            'upstream_config_sha256': sha(config)}
        return {'study_id': study, 'family_id': case['family_id'], 'preparation_status': UNAVAILABLE,
                'attempt': record, 'source_sha256': {str(receipt): sha(receipt), str(config): sha(config)}}
    if envelope['status'] != 'PASS':
        def unavailable(_directory): raise RuntimeError(envelope.get('reason', 'Source packet unavailable'))
        return failure(terminal_stage(store, study + '/preparation', unavailable))
    _verify_sources(envelope['source_sha256'])
    metadata, packet = envelope['metadata'], envelope['packet']; species = case['modality'].split('_')[0]
    source = case['molecular_source']
    def acquire(directory):
        if len(metadata['sample_ids']) > plan['max_samples_per_study']:
            raise RuntimeError('Study exceeds the frozen common preparation sample ceiling')
        path = _input_path(source)
        if source['kind'] == 'archs4':
            return extract_rna(path, study, species, metadata['sample_ids'], directory)
        frame = read_matrix(path)
        if set(frame.columns) != set(metadata['sample_ids']):
            raise ValueError('Raw count/GEO membership differs; no silent sample subset')
        destination = directory / 'counts.csv'; frame.to_csv(destination, mode='x')
        return {'path': str(destination), 'sha256': sha(destination), 'source': str(path)}
    raw = terminal_stage(store, study + '/molecular_input', acquire)
    if raw['status'] != 'PASS': return failure(raw)
    processed = terminal_stage(store, study + '/preprocessing', lambda d: preprocess_rna(
        Path(raw['value']['path']), species, d, normalization='clockbase_legacy', source=packet))
    if processed['status'] != 'PASS': return failure(processed)
    clocks = terminal_stage(store, study + '/clock_predictions', lambda d: predict_rna(processed['value'], species, d))
    if clocks['status'] != 'PASS': return failure(clocks)
    limits = plan['primary_limits']
    client = ActorClient(store.root / 'baseline-api' / study, 'primary', limits=limits,
                        total_calls=limits['max_calls'], total_cost=limits['max_cost'], model=plan['model'])
    primary = terminal_stage(store, study + '/coding_agent', lambda d: primary_analysis(
        packet, clocks['value'], client, d))
    if primary['status'] != 'PASS': return failure(primary)
    identities = dict(envelope['source_sha256'])
    for item in (raw, processed, clocks, primary):
        stage = Path(item['terminal_receipt']).parent
        identities.update({str(path): sha(path) for path in stage.rglob('*') if path.is_file()})
    return {'study_id': study, 'family_id': case['family_id'], 'preparation_status': 'PASS',
            'packet': packet, 'primary': primary['value'], 'clocks': clocks['value'],
            'source_sha256': identities, 'historical_release_state': 'CANNOT VERIFY',
            'producer_rule': plan['producer_rule'], 'preprocessing': 'clockbase_legacy'}


def run(release_path, output, *, until='sources', reference_states=None):
    from .evaluation_2026_09_10 import (compare_development_case, compare_unavailable_case,
                                      comparison_resource_record)
    from src.pipeline_inputs_2026_09_10 import checked_assets
    context = validate_release(release_path); output = Path(output).resolve(); plan = context['plan']
    output.mkdir(parents=True, exist_ok=True)
    with (output / '.run.lock').open('a') as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error: raise RuntimeError('Another process owns the prospective run') from error
        store = ProspectiveStore(output, archive_config(context))
        def enrollment(directory):
            path = directory / 'enrollment.json'; write_new(path, context['enrollment'])
            return {'enrollment': context['enrollment'], 'path': str(path), 'sha256': sha(path)}
        _stage(store, '00_enrollment', enrollment)
        _verify_completed_api(output)
        sources = prepare_sources(store, context)
        if until == 'sources':
            return {'phase': 'SOURCE_ONLY', 'sources': sources, 'generation_calls': 0,
                    'clock_execution_performed': False, 'reference_labels_completed': False}
        if reference_states is None:
            raise IntegrityError('A complete source-reference state manifest is required before numerical preparation')
        freeze_source_states(store, context, reference_states)
        if until == 'references': return {'phase': 'SOURCE_REFERENCE_STATES_FROZEN'}
        def data_identity(directory):
            checked_assets(); identities = {}
            for case in context['cases']:
                source = case['molecular_source']; path = Path(source['path']).resolve()
                if str(path) in identities: continue
                try:
                    _input_path(source)
                    if sha(path) != source['sha256']:
                        raise IntegrityError('Frozen molecular source checksum changed: ' + str(path))
                    identities[str(path)] = {'status': 'PASS', 'sha256': source['sha256']}
                except FileNotFoundError as error:
                    identities[str(path)] = {'status': 'CANNOT VERIFY', 'error': str(error)}
            write_new(directory / 'molecular-source-identities.json', identities)
            return identities
        _stage(store, '01_molecular_source_identities', data_identity)
        def inputs(directory):
            prepared = [prepare_case(store, context, case) for case in context['cases']]
            write_new(directory / 'prepared-cases.json', prepared)
            return {'cases': prepared, 'reference_labels_supplied': False,
                    'fresh_producer_attempted': True, 'preprocessing': 'clockbase_legacy'}
        prepared = _stage(store, '00_inputs', inputs)['cases']
        for case in prepared: _verify_sources(case['source_sha256'])
        if until == 'producer': return {'phase': 'FRESH_PRODUCER_COMPLETE', 'studies': len(prepared)}
        limits = plan['review_limits']
        def client_factory(study, method):
            return ActorClient(output / 'api', study + '-' + method, limits=limits,
                total_calls=len(prepared) * 2, total_cost=len(prepared) * 2 * limits['max_cost'], model=plan['model'])
        rows = []
        for case in prepared:
            rows.extend(compare_unavailable_case(store, case) if case.get('preparation_status') == UNAVAILABLE
                        else compare_development_case(store, case, client_factory))
        def summary(directory):
            records = [comparison_resource_record(output, row) for row in rows]
            baseline = {}
            for case in context['cases']:
                study = case['study_id']; slots = list((output / 'baseline-api' / study).glob('primary-*'))
                receipts = [read_json(p / 'receipt.json') for p in slots if (p / 'receipt.json').exists()]
                baseline[study] = {'generation_attempts': sum((p / 'reservation.json').is_file() for p in slots),
                    'responses_with_usage': len(receipts),
                    'estimated_known_cost_usd': sum(x['estimated_cost_usd'] for x in receipts),
                    'unknown_billing_attempts': sum(not (p / 'receipt.json').exists() for p in slots),
                    'billed_cost_usd': None if slots else 0}
            result = {'schema': 'clockbase-prospective-method-results-v1', 'split': context['split'],
                'attempted_studies': len(prepared), 'records': records, 'baseline_resources': baseline,
                'accepted_result_error_rate': None, 'reference_evaluation_complete': False,
                'scope': 'Method decisions only; independently frozen references are required for error rates'}
            write_new(directory / 'comparison-summary.json', result)
            write_new(directory / 'comparison-records.json', records)
            return result
        result = _stage(store, '99_summary', summary)
        validate_archive_context(output)
        for case in prepared: _verify_sources(case['source_sha256'])
        return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--until', choices=('sources', 'references', 'producer', 'methods'), default='sources')
    parser.add_argument('--reference-states', type=Path)
    args = parser.parse_args(argv)
    return run(args.release, args.output, until=args.until, reference_states=args.reference_states)


if __name__ == '__main__':
    try:
        print(json.dumps(main(), indent=2, allow_nan=False))
    except Exception as error:
        raise SystemExit(f'Prospective comparison stopped: {type(error).__name__}: {error}') from error
