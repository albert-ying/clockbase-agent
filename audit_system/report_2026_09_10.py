"""Build source-linked development evidence and explicitly budgeted source reviews.

The explicit measurement-profile mode recomputes clocks from source CpGs.
The source-design plan and batch modes can make their declared model calls. Other modes
recheck or summarize preserved outputs. No mode estimates audit
error rates without independent reference judgments.
"""
import argparse
import ast
from collections import Counter
import csv
from decimal import Decimal
import gzip
import hashlib
from datetime import datetime, timezone
from importlib.metadata import version
import json
import fcntl
from itertools import zip_longest
from pathlib import Path
import re
import shutil
import time
import sys

import numpy as np
import pandas as pd
from .agent_workflow_2026_09_10 import sha, write_new


def check_rna_measurement_profile(plan_path, output):
    """Run a pinned exposed-source measurement check without design/model calls."""
    from .agent_workflow_2026_09_10 import source_packet
    from .adapters.measurement_pilot import verify_clock_measurement, MOUSE_RNA_MODELS
    from src.pipeline_inputs_2026_09_10 import checked_assets, preprocess_rna, predict_rna, ASSETS
    from src.pipeline_runtime_2026_09_10 import StageStore
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    plan = json.loads(plan_path.read_text())
    if (plan.get('schema') != 'clockbase-rna-measurement-development-v1' or plan.get('split') != 'development'
            or plan.get('species') != 'mouse' or plan.get('normalization') != 'published_mouse_rle'
            or not re.fullmatch(r'GSE\d+', plan.get('study_id', ''))):
        raise ValueError('RNA measurement check requires a pinned exposed mouse development plan')
    paths = {key: Path(plan[key]).resolve() for key in ('counts', 'metadata', 'sample_source')}
    sources = plan['source_sha256']
    if any(str(path) not in sources for path in paths.values()):
        raise ValueError('RNA measurement input is absent from the frozen source manifest')
    for path, digest in sources.items():
        if not Path(path).is_file(): raise FileNotFoundError('Cannot proceed - required data file not found: ' + path)
        if sha(path) != digest: raise ValueError('RNA measurement source changed: ' + path)
    checked_assets()
    packet = source_packet(paths['metadata'], plan['study_id'], 'Declared RNA numerical application check',
                           sample_source=paths['sample_source'])
    source_ids = {row['sample_id'] for row in packet['samples']}
    for name in ('processing_sample_ids', 'reference_sample_ids', 'verification_sample_ids'):
        if (not isinstance(plan.get(name), list) or set(plan[name]) != source_ids
                or len(plan[name]) != len(source_ids)):
            raise ValueError('RNA measurement plan must freeze the complete ' + name)
    if plan.get('clock_model_sha256') != MOUSE_RNA_MODELS:
        raise ValueError('RNA measurement plan must freeze both public model identities')
    root = Path(__file__).resolve().parents[1]
    code = ['audit_system/report_2026_09_10.py', 'audit_system/agent_workflow_2026_09_10.py',
            'audit_system/source_packet.py', 'audit_system/adapters/measurement_pilot.py', 'audit_system/hashing.py',
            'audit_system/clock_replay_2026_09_10.py', 'src/pipeline_inputs_2026_09_10.py', 'src/pipeline_runtime_2026_09_10.py']
    config = {'plan': plan, 'plan_sha256': sha(plan_path), 'source_packet': packet,
              'implementation_sha256': {name: sha(root/name) for name in code},
              'asset_manifest_sha256': sha(ASSETS/'manifest.json'),
              'dependencies': {name: version(name) for name in ('numpy', 'pandas', 'scikit-learn', 'joblib')},
              'development_only': True, 'model_calls': 0, 'scientific_status': 'CANNOT VERIFY'}
    store = StageStore(output, config)
    processed = store.run('01_preprocessing', lambda d: preprocess_rna(paths['counts'], 'mouse', d,
                         normalization='published_mouse_rle', source=packet))
    if (processed['source_stratum']['sample_ids'] != plan['processing_sample_ids']
            or processed['source_stratum']['reference_ids'] != plan['reference_sample_ids']):
        raise ValueError('RNA preprocessing membership/order differs from the frozen plan')
    clocks = store.run('02_clock_predictions', lambda d: predict_rna(processed, 'mouse', d))
    if {row['clock']: row['model_sha256'] for row in clocks['checks']} != plan['clock_model_sha256']:
        raise ValueError('RNA executed model identities differ from the frozen plan')
    def verify(directory):
        ids = plan['verification_sample_ids']
        rows = [verify_clock_measurement(packet, clocks, clock, ids) for clock in clocks['clocks']]
        result = {'study_id': plan['study_id'], 'rows': rows, 'source_packet': packet, 'clocks': clocks,
                  'clock_predictions_recomputed': True, 'raw_count_preprocessing_recomputed': True,
                  'scientific_status': 'CANNOT VERIFY', 'accepted_scientific_repairs': 0,
                  'scope': 'Source-supported model-score application only; no intervention/design adjudication or audit error-rate estimate',
                  'development_only': True, 'model_calls': 0}
        write_new(directory/'measurement-verification.json', result)
        compact = [{'clock': row['clock'], 'measurement_status': row['status'],
                    'processing_samples': len(row.get('processing_sample_ids', [])),
                    'required_features': row.get('required_features'),
                    'observed_min': min(row['observed_features_by_sample'].values()) if row.get('observed_features_by_sample') else None,
                    'imputed_max': max(row['imputed_features_by_sample'].values()) if row.get('imputed_features_by_sample') else None,
                    'scientific_status': 'CANNOT VERIFY'} for row in rows]
        with (directory/'measurement-table.csv').open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(compact[0])); writer.writeheader(); writer.writerows(compact)
        text = ['| Clock | Measurement | Samples | Required features | Observed minimum | Imputed maximum | Scientific status |',
                '|---|---|---:|---:|---:|---:|---|']
        text += ['| ' + ' | '.join(str(v) for v in row.values()) + ' |' for row in compact]
        (directory/'measurement-table.md').write_text('\n'.join(text) + '\n')
        return result
    result = store.run('03_measurement_verification', verify)
    for path, digest in sources.items():
        if sha(path) != digest: raise ValueError('RNA input changed during the measurement run')
    for path, digest in config['implementation_sha256'].items():
        if sha(root/path) != digest: raise ValueError('RNA implementation changed during the measurement run')
    return result


def compare_rna_normalization(plan_path, output):
    """Quantify one declared preprocessing difference on exposed source data.

    This preserves the legacy mapping/filtering universe. It is a sensitivity
    analysis, not an automatically accepted preprocessing or biological repair.
    """
    from .agent_workflow_2026_09_10 import source_packet
    from src.pipeline_inputs_2026_09_10 import checked_assets, preprocess_rna, predict_rna, read_matrix, ASSETS
    from src.pipeline_runtime_2026_09_10 import StageStore
    from src.pipeline_2026_09_10 import recompute_declared_contrast
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    plan = json.loads(plan_path.read_text())
    conformance = plan.get('schema') == 'clockbase-rna-preprocessing-conformance-v1'
    if (plan.get('schema') not in {'clockbase-rna-normalization-sensitivity-v1', 'clockbase-rna-preprocessing-conformance-v1'}
            or plan.get('split') != 'development' or plan.get('species') not in {'mouse', 'human'}
            or not re.fullmatch(r'GSE\d+', plan.get('study_id', ''))):
        raise ValueError('RNA normalization sensitivity requires an explicitly exposed development study')
    if conformance and (plan['species'] != 'mouse' or plan.get('profiles') !=
                        ['clockbase_legacy', 'rle_sensitivity', 'published_mouse_rle']):
        raise ValueError('Mouse conformance plan must freeze the three preprocessing profiles')
    inputs = plan['source_sha256']
    required = [Path(plan[key]).resolve() for key in ('counts', 'metadata', 'archived_predictions')]
    if any(str(path) not in inputs for path in required):
        raise ValueError('Normalization input is absent from the pinned source manifest')
    for path, digest in inputs.items():
        if not Path(path).is_file(): raise FileNotFoundError('Cannot proceed - required data file not found: ' + path)
        if sha(path) != digest: raise ValueError('Normalization source changed: ' + path)
    checked_assets()
    source = source_packet(required[1], plan['study_id'], 'Development preprocessing sensitivity')
    source_ids = {row['sample_id'] for row in source['samples']}
    counts = read_matrix(required[0])
    if set(counts.columns) != source_ids:
        raise ValueError('RNA counts and source sample membership differ')
    comparisons = plan.get('comparisons')
    if not isinstance(comparisons, list) or not comparisons:
        raise ValueError('Sensitivity contrasts must be declared before computation')
    for row in comparisons:
        members = row['case_ids'] + row['control_ids']
        if (len(set(members)) != len(members) or set(members) - source_ids
                or min(len(row['case_ids']), len(row['control_ids'])) < 2 or row.get('test') != 'Welch'):
            raise ValueError('Sensitivity requires declared source groups and the supported conditional Welch calculation')
    root = Path(__file__).resolve().parents[1]
    code = [root/name for name in ['audit_system/report_2026_09_10.py', 'audit_system/agent_workflow_2026_09_10.py',
            'audit_system/source_packet.py', 'audit_system/clock_replay_2026_09_10.py',
            'src/pipeline_inputs_2026_09_10.py', 'src/pipeline_runtime_2026_09_10.py', 'src/pipeline_2026_09_10.py']]
    config = {'plan': plan, 'plan_sha256': sha(plan_path), 'implementation_sha256': {str(p): sha(p) for p in code},
              'asset_manifest_sha256': sha(ASSETS/'manifest.json'), 'model_calls': 0,
              'scientific_status': 'CANNOT VERIFY', 'development_only': True}
    store = StageStore(output, config)
    variants = {}
    profiles = plan['profiles'] if conformance else ['clockbase_legacy', 'rle_sensitivity']
    for index, normalization in enumerate(profiles):
        processed = store.run(f'{index*2+1:02d}_{normalization}_preprocessing',
            lambda directory, profile=normalization: preprocess_rna(required[0], plan['species'], directory,
                                                                   normalization=profile, source=source))
        variants[normalization] = store.run(f'{index*2+2:02d}_{normalization}_clocks',
            lambda directory, data=processed: predict_rna(data, plan['species'], directory))

    def compare(directory):
        frames = {profile: pd.read_csv(value['path']).set_index('sample_id') for profile, value in variants.items()}
        legacy, sensitivity = frames['clockbase_legacy'], frames['rle_sensitivity']
        if (set(legacy.index) != source_ids or not legacy.index.equals(sensitivity.index)
                or not legacy.columns.equals(sensitivity.columns)):
            raise ValueError('Normalization variants changed sample or clock identity')
        if any(not frame.index.equals(legacy.index) or not frame.columns.equals(legacy.columns)
               for frame in frames.values()):
            raise ValueError('Preprocessing profiles changed sample or clock identity')
        old = pd.read_csv(required[2]).set_index('sample_id')
        if set(old.index) != source_ids or set(old.columns) != set(legacy.columns):
            raise ValueError('Archived prediction identities differ from the fresh legacy run')
        old = old.loc[legacy.index, legacy.columns]
        legacy_delta = float(np.max(np.abs(old.to_numpy() - legacy.to_numpy())))
        if not np.allclose(old, legacy, rtol=1e-10, atol=1e-12):
            raise ValueError('Fresh legacy predictions do not reproduce the archived input')
        delta = sensitivity - legacy
        delta.rename_axis('sample_id').to_csv(directory/'per-sample-clock-differences.csv', mode='x')
        if conformance:
            (frames['published_mouse_rle'] - legacy).rename_axis('sample_id').to_csv(
                directory/'published-order-minus-legacy.csv', mode='x')
            (frames['published_mouse_rle'] - sensitivity).rename_axis('sample_id').to_csv(
                directory/'published-order-minus-rle-sensitivity.csv', mode='x')
        rows = []
        for index, comparison in enumerate(comparisons):
            for profile, frame in frames.items():
                for clock in frame.columns:
                    result = recompute_declared_contrast(frame, comparison['case_ids'], comparison['control_ids'], clock, 'Welch')
                    rows.append({'study_id': plan['study_id'], 'comparison_index': index,
                                 'target': comparison['target'], 'normalization': profile, **result,
                                 'scientific_status': 'CANNOT VERIFY',
                                 'scope': 'Conditional arithmetic; source design and biological interpretation remain unresolved'})
        write_new(directory/'conditional-contrasts.json', rows)
        with (directory/'conditional-contrasts.csv').open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader()
            writer.writerows({key: json.dumps(value) if isinstance(value, list) else value
                              for key, value in row.items()} for row in rows)
        result = {'study_id': plan['study_id'], 'sample_count': len(legacy), 'clock_count': len(legacy.columns),
                  'conditional_contrast_rows': len(rows), 'archived_legacy_max_abs_difference': legacy_delta,
                  'normalization_max_abs_difference_by_clock': {clock: float(delta[clock].abs().max()) for clock in delta},
                  'clock_execution_performed': True, 'model_calls': 0, 'accepted_scientific_repairs': 0,
                  'held_out_studies_evaluated': 0, 'scientific_status': 'CANNOT VERIFY',
                  'scope': 'Observed preprocessing sensitivity under fixed legacy mapping/filtering; no biological-validity or corrected-preprocessing claim'}
        if conformance:
            result.update(profiles=profiles,
                published_order_max_abs_difference_by_clock={clock: float((frames['published_mouse_rle'][clock] - legacy[clock]).abs().max()) for clock in legacy},
                filter_order_max_abs_difference_by_clock={clock: float((frames['published_mouse_rle'][clock] - sensitivity[clock]).abs().max()) for clock in legacy},
                scope='Development preprocessing comparison; source design and biological application remain unresolved; independent conformance assessment is separate')
        write_new(directory/'summary.json', result)
        return result
    result = store.run(f'{len(profiles)*2+1:02d}_normalization_comparison', compare)
    for path, digest in {**inputs, **config['implementation_sha256']}.items():
        if sha(path) != digest: raise ValueError('Normalization input or implementation changed during execution: '+path)
    return result


def export_source_audit_trace(run, output):
    """Export a completed development batch child without API envelopes or paths.

    A parent artifact manifest is required even for failed model responses, so
    negative results receive the same integrity check as completed reviews.
    Source-only checks never become a scientific acceptance badge.
    """
    from .agent_workflow_2026_09_10 import AUDIT_DOMAINS, verify_domain_review
    from .source_packet import source_json_bytes
    from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError
    from urllib.parse import urlparse
    run, output = Path(run).resolve(), Path(output).resolve()

    def read_stage(directory):
        complete = directory/'complete.json'
        if sha(complete) != (directory/'complete.sha256').read_text().strip():
            raise PipelineIntegrityError('Export rejected a changed stage receipt')
        record = json.loads(complete.read_text())
        for name, checksum in record['outputs'].items():
            path = (directory/name).resolve()
            if not path.is_relative_to(directory) or sha(path) != checksum:
                raise PipelineIntegrityError('Export rejected a changed or escaped stage artifact')
        return record

    if run.parent.name != 'studies':
        raise ValueError('Public source trace requires a completed batch child')
    parent = run.parent.parent
    completed = read_stage(parent/'01_batch_execution')
    names = [name for name in completed['outputs'] if name.endswith('/child-artifact-hashes.json')]
    if len(names) != 1:
        raise ValueError('Public export requires the completed child artifact manifest')
    batch_manifest = parent/'01_batch_execution'/names[0]
    expected = json.loads(batch_manifest.read_text())
    observed = {str(path.relative_to(parent)): sha(path) for path in (parent/'studies').rglob('*') if path.is_file()}
    if expected != observed:
        raise PipelineIntegrityError('Public export rejected changed batch children')
    config = json.loads((run/'config.json').read_text())
    if output.exists():
        receipt = json.loads((output/'export-receipt.json').read_text())
        if (receipt['source_run_config_sha256'] != sha(run/'config.json')
                or receipt['source_batch_manifest_sha256'] != sha(batch_manifest)):
            raise PipelineIntegrityError('Existing public export belongs to a different source run')
        for name, checksum in receipt['output_sha256'].items():
            path = (output/name).resolve()
            if not path.is_relative_to(output) or sha(path) != checksum:
                raise PipelineIntegrityError('Existing public export changed')
        return {'study_id': config['plan']['study_id'], 'trace':str(output/'audit_trace.json'),
                'scientific_status':'CANNOT VERIFY', 'export_sha256':sha(output/'audit_trace.json')}
    plan = config['plan']
    if plan.get('schema') not in {'clockbase-source-design-run-v3', 'clockbase-source-design-run-v4'} or plan.get('split') != 'development':
        raise ValueError('Public domain export requires the opt-in version3 development plan')
    study = plan['study_id']
    if not re.fullmatch(r'GSE\d+', study):
        raise ValueError('Public empirical source trace requires a GEO study identity')
    preflight = read_stage(run/'00_source_preflight')
    names = [name for name in preflight['outputs'] if name.endswith('/source-packet.json')]
    if len(names) != 1:
        raise ValueError('Source packet is unavailable for public trace')
    packet = json.loads((run/'00_source_preflight'/names[0]).read_text())
    if (packet['study_id'] != study or hashlib.sha256(source_json_bytes(packet)).hexdigest() != config['source_packet_sha256']):
        raise PipelineIntegrityError('Public trace source packet differs from preflight')
    design, verification = None, None
    if (run/'01_source_design/complete.json').exists():
        stage = read_stage(run/'01_source_design')
        design = stage['result']
        verification = verify_domain_review(design, packet)
        names = [name for name in stage['outputs'] if name.endswith('/domain-verification.json')]
        if len(names) != 1 or json.loads((run/'01_source_design'/names[0]).read_text()) != verification:
            raise PipelineIntegrityError('Domain verification does not reproduce from the source packet')
    else:
        matches = [row for row in completed['result']['rows'] if row['study_id'] == study]
        if len(matches) != 1 or matches[0]['response_contract_status'] != 'FAIL':
            raise ValueError('Neither a completed review nor an attested failed review exists')

    def public_url(url):
        parsed = urlparse(url)
        if (parsed.scheme != 'https' or parsed.hostname not in {'www.ncbi.nlm.nih.gov', 'pmc.ncbi.nlm.nih.gov', 'pubmed.ncbi.nlm.nih.gov', 'www.ebi.ac.uk'}
                or parsed.username or parsed.password):
            raise ValueError('Source URL is outside the declared public-source hosts')
        if parsed.hostname == 'www.ebi.ac.uk' and not re.fullmatch(r'/europepmc/webservices/rest/PMC\d+/fullTextXML', parsed.path):
            raise ValueError('Unrecognized Europe PMC public-source path')
        return url

    document_urls = {doc['source_id']: public_url(doc['source_url']) for doc in packet.get('study_documents', [])}
    def evidence(entries):
        rows = []
        for item in entries:
            row = {key: item[key] for key in ('sample_id', 'source_id', 'field', 'quote') if key in item}
            if 'sample_id' in row:
                if not re.fullmatch(r'GSM\d+', row['sample_id']):
                    raise ValueError('Public sample evidence requires a GEO sample identity')
                row['source_url'] = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + row['sample_id']
            else:
                row['source_url'] = document_urls[row['source_id']]
            rows.append(row)
        return rows

    public_sources = {'schema': 'clockbase-public-source-records-v1', 'study_id': study,
        'source_metadata_sha256': packet['source_sha256'], 'source_url': public_url(packet['source_url']),
        'samples': packet['samples'],
        'study_documents': [{'source_id': doc['source_id'], 'source_url': document_urls[doc['source_id']],
            'source_sha256': doc['source_sha256'], 'fields': doc['fields']} for doc in packet.get('study_documents', [])]}
    domains = []
    for name in AUDIT_DOMAINS:
        assessment = (dict(design['audit_domains'][name]) if design else {
            'status':'CANNOT VERIFY', 'summary':'Model review failed its response contract; no domain assessment accepted',
            'evidence':[], 'unresolved':['Complete validated model review unavailable'],
            'proposed_repair':None, 'verification_needed':['A new explicitly budgeted review, retaining this failed attempt']})
        assessment['evidence'] = evidence(assessment['evidence'])
        domains.append({'id':name, 'title':name.replace('_',' ').capitalize(),
            'model_assessment':assessment, 'verification':verification['domains'][name] if verification else {
                'status':'CANNOT VERIFY', 'scope':'No validated model assessment to check', 'checks':[]}})
    comparisons = []
    for item in (design or {}).get('comparisons', []):
        row = {key:item[key] for key in ('comparison_id','comparison_name','target','case_ids','control_ids',
                                       'excluded_ids','biological_unit','pairing','unresolved')}
        row['status'] = 'CANNOT VERIFY'
        comparisons.append(row)
    covariates = [{**item, 'evidence':evidence(item['evidence'])} for item in (design or {}).get('covariate_plan', [])]
    receipts = [json.loads(path.read_text()) for path in sorted((run/'api').glob('*/receipt.json'))]
    reservation_paths = list((run/'api').glob('*/reservation.json'))
    for path in sorted((run/'api').glob('*/receipt.json')):
        receipt = json.loads(path.read_text())
        if sha(path.parent/'response.json') != receipt['response_sha256']:
            raise PipelineIntegrityError('API response hash differs from its retained receipt')
    trace = {'schema':'clockbase-audit-trace-v1', 'study_id':study, 'run_id':'sha256:'+sha(run/'config.json'),
        'recorded_at_utc':datetime.now(timezone.utc).isoformat().replace('+00:00','Z'), 'evaluation_split':'development',
        'scope':'Source-only development audit; numerical outputs and molecular matrices were not supplied to the reviewer',
        'models':[{'provider':plan['provider'], 'model':plan['model'], 'role':'source-design reviewer'}],
        'overall':{'status':'CANNOT VERIFY', 'reason':'Source interpretation, molecular/clock applicability and independent scientific reference judgments remain unresolved'},
        'domains':domains, 'comparisons':comparisons, 'covariates':covariates,
        'artifacts':[], 'events':[
            {'stage':'source preparation','status':'PASS','summary':'Frozen source packet and complete batch artifact manifest verified'},
            {'stage':'model response contract','status':'PASS' if design else 'FAIL',
             'summary':'Exact citations, sample membership and required review domains validated' if design else 'All failed responses are retained; no repair or scientific result accepted'},
            {'stage':'repair application','status':'CANNOT VERIFY','summary':'No repair applied by this source-only review'}],
        'resource_usage':{'generation_attempts':len(reservation_paths),'responses_received':len(receipts),
            'estimated_cost_usd':sum(row['estimated_cost_usd'] for row in receipts),'billed_cost_usd':None},
        'limitations':['Model assessments are advisory; exact quotations do not establish semantic support.',
            'Published source prose may contain findings; this development run is not certified outcome-blinded.',
            'No molecular QC, clock execution, adjusted analysis or independent reference adjudication was performed by this source-only run.',
            'Costs are API estimates, not invoices; collaboration-runtime charges are not included.',
            'No scientific error rate can be estimated from these development cases.']}
    for path in sorted((run/'01_source_design').glob('attempt-*/response-validation/*-validation.json')):
        validation = json.loads(path.read_text())
        trace['events'].append({'stage':'response validation '+path.name.removesuffix('-validation.json'),
            'status':validation['validation_status'],
            'summary':validation.get('error', 'Response schema, sample assignment and quotation contract passed; scientific status remains unresolved')})
    def safe_content(value):
        if isinstance(value, dict):
            for key, item in value.items(): safe_content(key); safe_content(item)
        elif isinstance(value, list):
            for item in value: safe_content(item)
        elif isinstance(value, str) and re.search(r'/Entropy/|/Lasagna/|/home/kying|/Users/|Bearer\s+|sk-[A-Za-z0-9_-]{20,}', value):
            raise ValueError('Export refused potentially private text; source quotations were not rewritten')
    safe_content(public_sources); safe_content(trace)
    if verification: safe_content(verification)
    output.mkdir(parents=True, exist_ok=False)
    (output/'public').mkdir()
    for label, value in [('source-records',public_sources), ('source-verification',verification)]:
        if value is None: continue
        path = output/'public'/(label+'.json'); write_new(path, value)
        trace['artifacts'].append({'id':label, 'kind':'public source evidence', 'label':label.replace('-',' ').capitalize(),
                                   'sha256':sha(path), 'download_path':'public/'+path.name})
    trace['artifacts'].append({'id':'source-metadata', 'kind':'source hash', 'label':'Original source metadata digest',
                              'sha256':packet['source_sha256'], 'url':public_url(packet['source_url'])})
    write_new(output/'audit_trace.json', trace)
    write_new(output/'export-receipt.json', {'schema':'clockbase-public-trace-export-v1',
        'source_run_config_sha256':sha(run/'config.json'), 'source_batch_manifest_sha256':sha(batch_manifest),
        'output_sha256':{str(path.relative_to(output)):sha(path) for path in output.rglob('*') if path.is_file()},
        'scientific_status':'CANNOT VERIFY'})
    return {'study_id':study, 'trace':str(output/'audit_trace.json'), 'scientific_status':'CANNOT VERIFY',
            'export_sha256':sha(output/'audit_trace.json')}


def run_source_design_plan(plan_path, output, *, prepare_only=False):
    """Budgeted source review, with opt-in validation feedback and frozen inputs."""
    from .agent_workflow_2026_09_10 import source_packet, review_messages
    from .source_packet import source_json_bytes
    from src.pipeline_2026_09_10 import assess_design
    from src.pipeline_runtime_2026_09_10 import StageStore, BoundedOpenAI
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    plan = json.loads(plan_path.read_text())
    schema = plan.get('schema')
    if schema == 'clockbase-source-design-run-v1':
        review_profile, validation_retries = 'legacy', 0
        if 'review_profile' in plan or 'validation_retries' in plan:
            raise ValueError('Source-review options require a version2 plan')
    elif schema == 'clockbase-source-design-run-v2':
        review_profile, validation_retries = plan.get('review_profile'), plan.get('validation_retries')
        if review_profile != 'source_fields_v2' or type(validation_retries) is not int or validation_retries not in (0, 1):
            raise ValueError('Version2 requires source_fields_v2 and zero or one validation retry')
    elif schema == 'clockbase-source-design-run-v3':
        review_profile, validation_retries = plan.get('review_profile'), plan.get('validation_retries')
        if review_profile != 'source_domains_v3' or type(validation_retries) is not int or validation_retries not in (0, 1):
            raise ValueError('Version3 requires source_domains_v3 and zero or one validation retry')
    elif schema == 'clockbase-source-design-run-v4':
        review_profile, validation_retries = plan.get('review_profile'), plan.get('validation_retries')
        if review_profile != 'source_domains_v4' or type(validation_retries) is not int or validation_retries not in (0, 1):
            raise ValueError('Version4 requires source_domains_v4 and zero or one validation retry')
    else:
        raise ValueError('Unsupported source design plan schema')
    source_binding_profile = plan.get('source_binding_profile', 'legacy')
    if source_binding_profile != 'legacy' and schema != 'clockbase-source-design-run-v4':
        raise ValueError('Opt-in source binding requires a version4 development plan')
    if (plan.get('split') != 'development' or plan.get('provider') != 'openai'
            or type(plan.get('max_calls')) is not int or plan['max_calls'] != 1 + validation_retries):
        raise ValueError('Source design plan requires bounded OpenAI development calls')
    for path, digest in plan['source_sha256'].items():
        if sha(path) != digest: raise ValueError('Frozen design source changed: '+path)
    paths = [Path(plan['metadata']).resolve()]
    if plan.get('context'): paths.append(Path(plan['context']).resolve())
    if any(str(path) not in plan['source_sha256'] for path in paths):
        raise ValueError('Design input is absent from source manifest')
    packet_args = {'comparison_scope':plan.get('comparison_scope', 'exact_target')} if schema == 'clockbase-source-design-run-v4' else {}
    packet = source_packet(paths[0], plan['study_id'], plan['target'], context=plan.get('context'), **packet_args)
    messages = review_messages(packet, review_profile=review_profile, source_binding_profile=source_binding_profile)
    request = {'model':plan['model'], 'messages':[
        {'role':'system' if m.type=='system' else 'user','content':m.content} for m in messages],
        'max_completion_tokens':plan['output_tokens'], 'reasoning_effort':'low', 'store':False}
    size = len(json.dumps(request, ensure_ascii=False).encode())
    limits = {key:plan[key] for key in ('max_calls','max_cost','output_tokens','prompt_bytes','input_per_million','output_per_million')}
    if any(type(limits[key]) is not int or limits[key] <= 0 for key in ('max_calls','output_tokens','prompt_bytes')):
        raise ValueError('Design token, byte and call limits must be positive integers')
    if any(type(limits[key]) not in (int,float) or not np.isfinite(limits[key]) or limits[key] <= 0
           for key in ('max_cost','input_per_million','output_per_million')):
        raise ValueError('Design cost limit and tariffs must be positive finite numbers')
    reserve = ((size+100*len(messages))*limits['input_per_million']
               +limits['output_tokens']*limits['output_per_million'])/1e6
    if size > limits['prompt_bytes'] or reserve > limits['max_cost']:
        raise ValueError('Design request exceeds frozen byte/cost limits before transport')
    root = Path(__file__).resolve().parents[1]
    code = [root/name for name in ['audit_system/report_2026_09_10.py','audit_system/agent_workflow_2026_09_10.py',
        'audit_system/source_packet.py','src/pipeline_2026_09_10.py','src/pipeline_runtime_2026_09_10.py']]
    config = {'plan':plan,'plan_sha256':sha(plan_path),'source_packet_sha256':hashlib.sha256(source_json_bytes(packet)).hexdigest(),
              'implementation_sha256':{str(path):sha(path) for path in code},'python':sys.version,
              'versions':{name:version(name) for name in ('numpy','pandas','langchain-core')},
              'request_bytes':size,'reserved_cost_upper_estimate_usd':reserve,
              'request_sha256':hashlib.sha256(source_json_bytes(request)).hexdigest()}
    output.mkdir(parents=True, exist_ok=True)
    with (output/'.run.lock').open('a') as lock:
        try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError('Another process owns this source design run')
        store = StageStore(output,config)
        def prepare(directory):
            write_new(directory/'source-packet.json',packet)
            write_new(directory/'prepared-request.json',request)
            return {'status':'PASS','samples':len(packet['samples']),'request_bytes':size,
                    'reserved_cost_upper_estimate_usd':reserve,'generation_calls':0}
        prepared = store.run('00_source_preflight',prepare)
        if prepare_only: return prepared
        def review(directory):
            client = BoundedOpenAI(output/'api','source_design',model=plan['model'],**limits,
                                   require_model_match=schema in {'clockbase-source-design-run-v2', 'clockbase-source-design-run-v3', 'clockbase-source-design-run-v4'})
            result = assess_design(packet,client,directory, review_profile=review_profile,
                                   validation_retries=validation_retries, source_binding_profile=source_binding_profile)
            if schema in {'clockbase-source-design-run-v2', 'clockbase-source-design-run-v3', 'clockbase-source-design-run-v4'}:
                write_new(directory/'api-artifact-hashes.json', {
                    str(path.relative_to(output/'api')): sha(path)
                    for path in (output/'api').rglob('*') if path.is_file()})
            return result
        design = store.run('01_source_design',review)
        if schema in {'clockbase-source-design-run-v2', 'clockbase-source-design-run-v3', 'clockbase-source-design-run-v4'}:
            completed = json.loads((output/'01_source_design/complete.json').read_text())
            manifests = [name for name in completed['outputs'] if name.endswith('/api-artifact-hashes.json')]
            if len(manifests) != 1:
                raise ValueError('Completed source review lacks its API artifact manifest')
            expected = json.loads((output/'01_source_design'/manifests[0]).read_text())
            observed = {str(path.relative_to(output/'api')): sha(path)
                        for path in (output/'api').rglob('*') if path.is_file()}
            if not expected or expected != observed:
                raise ValueError('Completed source-review API artifacts changed; no new call')
        saved_request = json.loads((output/'api/source_design-000/request.json').read_text())
        if saved_request != request: raise ValueError('Executed model request differs from preflight')
        for path,digest in {**plan['source_sha256'],**config['implementation_sha256']}.items():
            if sha(path) != digest: raise ValueError('Design input or implementation changed during execution')
        return {'study_id':plan['study_id'],'design_status':design['status'],'scientific_status':'CANNOT VERIFY',
                'samples':len(packet['samples']),'design':design,
                'scope':'Source-only design proposal; no clocks, numerical repair or reference adjudication'}


def run_source_design_batch(plan_path, output, *, prepare_only=False):
    """Preflight every child before transport; freeze their aggregate allocation.

    Each child has its own persistent allowance, including validation feedback.
    Their sum cannot exceed the parent cap. Failed/uncertain requests remain
    charged to that child's allowance; completed batches resume without calls.
    """
    from src.pipeline_runtime_2026_09_10 import StageStore, PipelineIntegrityError
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    batch = json.loads(plan_path.read_text())
    if (batch.get('schema') != 'clockbase-source-design-batch-v1'
            or batch.get('split') != 'development' or batch.get('provider') != 'openai'
            or batch.get('concurrency') != 1 or type(batch.get('concurrency')) is not int):
        raise ValueError('Source design batch requires sequential OpenAI development runs')
    cap = batch.get('max_total_cost_usd')
    if type(cap) not in (int, float) or not np.isfinite(cap) or cap <= 0:
        raise ValueError('Batch cost cap must be positive and finite')
    if type(batch.get('max_total_calls')) is not int or batch['max_total_calls'] <= 0:
        raise ValueError('Batch call cap must be a positive integer')
    if not isinstance(batch.get('plans'), list) or not batch['plans']:
        raise ValueError('Batch requires frozen child plans')
    children, seen, allocation, calls = [], set(), Decimal('0'), 0
    for entry in batch['plans']:
        path = Path(entry['path']).resolve()
        if sha(path) != entry['sha256']:
            raise PipelineIntegrityError('Frozen child plan changed: ' + str(path))
        child = json.loads(path.read_text())
        study = child.get('study_id')
        if not isinstance(study, str) or not re.fullmatch(r'[A-Za-z0-9_-]+', study) or study in seen:
            raise ValueError('Batch study IDs must be unique safe directory names')
        if (child.get('schema') not in {'clockbase-source-design-run-v2', 'clockbase-source-design-run-v3', 'clockbase-source-design-run-v4'}
                or child.get('provider') != batch['provider'] or child.get('model') != batch.get('model')):
            raise ValueError('Batch child must use the declared provider/model and a version2 or version3 plan')
        value = child.get('max_cost')
        if type(value) not in (int, float) or not np.isfinite(value) or value <= 0:
            raise ValueError('Child cost cap must be positive and finite')
        if type(child.get('max_calls')) is not int or child['max_calls'] <= 0:
            raise ValueError('Child call cap must be a positive integer')
        allocation += Decimal(str(value))
        calls += child['max_calls']
        seen.add(study)
        children.append((path, child, output/'studies'/study))
    if allocation > Decimal(str(cap)) or calls > batch['max_total_calls']:
        raise ValueError('Child allocations exceed the aggregate batch cost/call cap; nothing sent')
    output.mkdir(parents=True, exist_ok=True)
    with (output/'.batch.lock').open('a') as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError('Another process owns this source design batch')
        # Every source, schema, implementation and initial request is checked
        # before the first child can make a generation call.
        for path, _, directory in children:
            run_source_design_plan(path, directory, prepare_only=True)
        frozen = {str(path): sha(path) for path, _, _ in children}
        frozen[str(plan_path)] = sha(plan_path)
        for _, child, directory in children:
            frozen.update(child['source_sha256'])
            frozen.update(json.loads((directory/'config.json').read_text())['implementation_sha256'])
        config = {'batch': batch, 'frozen_sha256': frozen,
                  'allocated_upper_limit_usd': float(allocation), 'maximum_calls': calls}
        store = StageStore(output, config)
        if prepare_only:
            return {'status': 'PASS', 'studies': len(children), 'generation_calls': 0,
                    'allocated_upper_limit_usd': float(allocation), 'maximum_calls': calls}
        def check_frozen():
            for path, digest in frozen.items():
                if sha(path) != digest:
                    raise PipelineIntegrityError('Batch input or implementation changed: ' + path)
        def child_hashes():
            return {str(path.relative_to(output)): sha(path)
                    for path in (output/'studies').rglob('*') if path.is_file()}
        def execute(directory):
            rows = []
            for path, child, child_directory in children:
                check_frozen()
                try:
                    result = run_source_design_plan(path, child_directory)
                    rows.append({key: value for key, value in result.items() if key != 'design'})
                    rows[-1]['response_contract_status'] = 'PASS'
                except PipelineIntegrityError:
                    raise
                except Exception as error:
                    rows.append({'study_id': child['study_id'], 'response_contract_status': 'FAIL',
                                 'scientific_status': 'CANNOT VERIFY',
                                 'error': f'{type(error).__name__}: {error}'})
                check_frozen()
            reservations = list((output/'studies').glob('*/api/*/reservation.json'))
            receipts = list((output/'studies').glob('*/api/*/receipt.json'))
            estimated = sum(json.loads(path.read_text())['estimated_cost_usd'] for path in receipts)
            uncertain = sum(json.loads(path.read_text())['reserved_usd'] for path in reservations
                            if not (path.parent/'receipt.json').exists())
            write_new(directory/'child-artifact-hashes.json', child_hashes())
            return {'rows': rows, 'generation_attempts': len(reservations), 'responses_received': len(receipts),
                    'estimated_cost_usd': estimated, 'uncertain_billing_reserved_usd': uncertain,
                    'billed_cost_usd': None, 'allocated_upper_limit_usd': float(allocation),
                    'scientific_status': 'CANNOT VERIFY', 'accepted_scientific_repairs': 0,
                    'scope': 'Development response-contract check; no scientific error-rate estimate'}
        result = store.run('01_batch_execution', execute)
        complete = json.loads((output/'01_batch_execution/complete.json').read_text())
        manifests = [path for path in complete['outputs'] if path.endswith('/child-artifact-hashes.json')]
        if len(manifests) != 1 or json.loads((output/'01_batch_execution'/manifests[0]).read_text()) != child_hashes():
            raise PipelineIntegrityError('Completed batch child artifacts changed; no new call')
        check_frozen()
        return result


def _xlsx_recovery_chunks(source, layout, columns, chunk_size=4096):
    """Independent openpyxl reader; never use cached formula results."""
    from openpyxl import load_workbook
    workbook = load_workbook(source, read_only=True, data_only=False, keep_links=False)
    try:
        if workbook.sheetnames != [layout['worksheet']]:
            raise ValueError('Independent XLSX reader found a different worksheet selection')
        sheet = workbook[layout['worksheet']]
        # XLSX dimensions are author-supplied metadata, not an observed row count.
        sheet.reset_dimensions()
        rows = sheet.iter_rows()
        header = next(rows, ())
        if ([cell.value for cell in header] != layout['header']
                or any(cell.data_type not in {'s', 'inlineStr'} for cell in header)):
            raise ValueError('Independent XLSX header differs from the source layout')
        width = len(header)
        features, values = [], []
        for cells in rows:
            if len(cells) > width:
                raise ValueError('Independent XLSX row extends beyond the header')
            if not cells or not isinstance(cells[0].value, str) or cells[0].data_type not in {'s', 'inlineStr'}:
                raise ValueError('Independent XLSX reader requires a text feature identity')
            row = []
            for cell in cells[1:]:
                value = cell.value
                if cell.data_type in {'f', 'b', 'e', 'd'} or isinstance(value, bool):
                    raise ValueError('Independent XLSX reader rejects formula, boolean, error and date cells')
                if value is None or (isinstance(value, str) and value.strip().lower() in {'', 'na', 'nan', 'null'}):
                    row.append(np.nan)
                elif type(value) in (str, int, float):
                    row.append(float(value))
                else:
                    raise ValueError('Independent XLSX cell is not numeric or explicitly missing')
            row.extend([np.nan] * (width-1-len(row)))
            features.append(cells[0].value); values.append(row)
            if len(features) == chunk_size:
                yield pd.DataFrame(values, index=features, columns=columns)
                features, values = [], []
        if features:
            yield pd.DataFrame(values, index=features, columns=columns)
        # The read-only openpyxl iterator omits merged-range metadata. Check
        # the exact worksheet part selected by that independent reader before
        # accepting its values; never interpret a merge as ordinary missingness.
        import xml.etree.ElementTree as ET
        import zipfile
        namespace = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
        with zipfile.ZipFile(source) as archive, archive.open(sheet._worksheet_path) as handle:
            sheet_data = None
            for event, element in ET.iterparse(handle, events=('start', 'end')):
                if event == 'start' and element.tag == namespace+'mergeCell':
                    raise ValueError('Independent XLSX reader rejects merged cells across sample or feature identities')
                if event == 'start' and element.tag == namespace+'sheetData':
                    sheet_data = element
                if event == 'end' and element.tag == namespace+'row':
                    element.clear()
                    if sheet_data is not None: sheet_data.clear()
    finally:
        workbook.close()


def _idat_verification_profile(name='minfi_ssnoob_beta100_detection01_v1'):
    """Independent declarations of strict v1 and frozen supported-subset v2."""
    profile = {'name': 'minfi_ssnoob_beta100_detection01_v1',
        'array': 'IlluminaHumanMethylationEPIC', 'platform': 'GPL23976',
        'import': {'extended': True, 'force': False},
        'noob': {'offset': 15, 'dyeCorr': True, 'dyeMethod': 'single'},
        'beta': {'type': 'Illumina', 'offset': 100, 'betaThreshold': 0, 'minZero': True},
        'qc': {'detection_type': 'm+u', 'failed_if_greater_than': 0.01,
            'mask_nonfinite': True, 'failed_cells': 'NA', 'sample_exclusion': 'none',
            'bead_counts': 'both channels retained; no bead-count filter'},
        'averaging': False, 'imputation': False, 'author_method_reproduction': False,
        'scientific_status': 'CANNOT VERIFY'}
    if name == 'minfi_ssnoob_beta100_detection01_supported_v2':
        profile.update(name=name, probe_membership={
            'raw': 'manifest probes with every required address observed',
            'corrected': 'address-supported probes in matching installed annotation correction groups',
            'normalization_input': 'all observed RG addresses',
            'controls': 'all manifest control addresses required',
            'exclusions': 'address_absent and annotation_absent recorded separately'})
    elif name != profile['name']:
        raise ValueError('Unsupported IDAT numerical candidate profile')
    return profile


def _idat_literal_pairs(study, series_path, samples_path):
    """Read GEO ownership without the acquisition adapter's source parser."""
    from urllib.parse import urlsplit
    series_ids, members, samples, active = [], [], {}, None
    for line in Path(series_path).read_text().splitlines():
        key, separator, value = line.partition(' = ')
        if separator and key == '^SERIES': series_ids.append(value.strip())
        if separator and key == '!Series_sample_id': members.append(value.strip())
    if (series_ids != [study] or not members or len(members) != len(set(members))
            or any(not re.fullmatch(r'GSM\d+', value) for value in members)):
        raise ValueError('IDAT independent Series identity or membership differs')
    for line in Path(samples_path).read_text().splitlines():
        key, separator, value = line.partition(' = ')
        if not separator: continue
        if key == '^SAMPLE':
            if not re.fullmatch(r'GSM\d+', value) or value in samples:
                raise ValueError('IDAT independent Sample identity is invalid or duplicated')
            active = value; samples[active] = {'platforms': [], 'urls': []}
        elif active and key == '!Sample_platform_id': samples[active]['platforms'].append(value)
        elif active and key.startswith('!Sample_supplementary_file'):
            if re.search(r'_(Grn|Red)\.idat\.gz$', value): samples[active]['urls'].append(value)
    if set(samples) != set(members):
        raise ValueError('IDAT independent Series and physical Sample records disagree')
    pairs = []
    for gsm in members:
        channels = {}
        if samples[gsm]['platforms'] != ['GPL23976']:
            raise ValueError('IDAT independent source platform is not the declared EPIC profile')
        for url in samples[gsm]['urls']:
            parsed = urlsplit(url); filename = parsed.path.rsplit('/', 1)[-1]
            match = re.fullmatch(re.escape(gsm)+r'_([A-Za-z0-9_.+-]+)_(Grn|Red)\.idat\.gz', filename)
            expected_parent = '/geo/samples/'+gsm[:3]+gsm[3:-3]+'nnn/'+gsm+'/suppl/'
            if (parsed.scheme not in {'ftp', 'https'} or parsed.netloc != 'ftp.ncbi.nlm.nih.gov'
                    or parsed.query or parsed.fragment or not match
                    or parsed.path != expected_parent+filename):
                raise ValueError('IDAT source URL does not declare exact GSM/channel ownership')
            channel = match[2]
            if channel in channels: raise ValueError('IDAT source repeats a channel')
            channels[channel] = {'url': url, 'https_url': 'https://ftp.ncbi.nlm.nih.gov'+parsed.path,
                'filename': filename, 'basename': filename.rsplit('_', 1)[0]}
        if (set(channels) != {'Grn', 'Red'}
                or channels['Grn']['basename'] != channels['Red']['basename']):
            raise ValueError('IDAT source does not declare one matching channel pair')
        pairs.append({'sample_id': gsm, 'channels': channels})
    return pairs


def _idat_binary_readback(path):
    """Separate IDAT v3 decoder: field directory, identity, addresses and quants.

    The format was checked against installed illuminaio's readIDAT_nonenc source;
    this routine does not call that library or the producer's decoded tables.
    """
    import struct
    try:
        data = Path(path).read_bytes()
    except OSError as error:
        raise OSError('Cannot read IDAT binary: '+str(path)) from error
    if len(data) < 16 or data[:4] != b'IDAT' or struct.unpack_from('<Q', data, 4)[0] != 3:
        raise ValueError('Independent reader requires an unencrypted IDAT v3 binary')
    count = struct.unpack_from('<I', data, 12)[0]
    directory_end = 16+10*count
    if not 1 <= count <= 256 or directory_end > len(data):
        raise ValueError('Invalid IDAT field directory length')
    fields = {}
    for index in range(count):
        code, offset = struct.unpack_from('<HQ', data, 16+10*index)
        if code in fields or not directory_end <= offset < len(data):
            raise ValueError('Invalid or repeated IDAT field offset')
        fields[code] = offset
    if len(set(fields.values())) != count:
        raise ValueError('IDAT fields overlap at their starting offsets')
    required = {1000, 102, 103, 104, 107, 400, 402, 403, 404}
    if not required <= fields.keys(): raise ValueError('IDAT lacks required identity or intensity fields')
    offsets = sorted(fields.values())+[len(data)]
    ends = dict(zip(offsets[:-1], offsets[1:]))
    def field(code, size):
        start = fields[code]
        if start+size > ends[start]: raise ValueError('IDAT typed field is truncated or overlapping')
        return memoryview(data)[start:start+size]
    n = struct.unpack('<I', field(1000, 4))[0]
    if not 0 < n <= 5_000_000: raise ValueError('IDAT address count exceeds the bounded decoder')
    address = np.frombuffer(field(102, n*4), dtype='<i4').copy()
    if (address <= 0).any() or len(np.unique(address)) != n:
        raise ValueError('IDAT address identities are invalid or duplicated')
    quants = {'address': address,
        'Mean': np.frombuffer(field(104, n*2), dtype='<u2').copy(),
        'SD': np.frombuffer(field(103, n*2), dtype='<u2').copy(),
        'NBeads': np.frombuffer(field(107, n), dtype='u1').copy()}
    def string(code):
        start, end = fields[code], ends[fields[code]]
        cursor, length, shift = start, 0, 0
        while True:
            if cursor >= end or shift > 28: raise ValueError('Malformed IDAT string length')
            value = data[cursor]; cursor += 1
            length |= (value & 127) << shift
            if value < 128: break
            shift += 7
        if not length or cursor+length > end: raise ValueError('Truncated or empty IDAT identity string')
        try: result = data[cursor:cursor+length].decode('ascii')
        except UnicodeError as error: raise ValueError('Non-ASCII IDAT identity') from error
        if any(ord(char) < 32 for char in result): raise ValueError('Control character in IDAT identity')
        return result
    header = {'version': 3, 'nSNPsRead': n, 'RedGreen': struct.unpack('<i', field(400, 4))[0],
              'Barcode': string(402), 'ChipType': string(403), 'position': string(404)}
    return header, quants


def _idat_csv_chunks(path, columns, index_name, chunksize=16384):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        header = next(csv.reader(handle), None)
    if header != [index_name]+list(columns):
        raise ValueError('IDAT table header/order differs: '+str(path))
    yield from pd.read_csv(path, index_col=0, converters={0: str}, chunksize=chunksize,
        float_precision='round_trip', na_values=['NA'], keep_default_na=False, on_bad_lines='error')


def _idat_verify_supported_outputs(paths, ids, type_i, type_ii, controls,
                                   common_addresses, channels, physical, processing):
    """Derive v2 membership before reading raw/corrected output rows.

    Annotation groups are a sealed producer export, not an independent annotation
    database. The joins, address rules, arithmetic and comparisons below are
    independent; the installed annotation and normalization remain unvalidated.
    """
    manifest = pd.concat([
        type_i[['probe_id', 'Color']].assign(Type='I')[['probe_id', 'Type', 'Color']],
        type_ii[['probe_id']].assign(Type='II', Color='')], ignore_index=True)
    full_ids = manifest['probe_id'].tolist()
    required = np.unique(np.concatenate([type_i['AddressA'], type_i['AddressB'],
                                        type_ii['AddressA'], controls['Address']]))
    missing = np.setdiff1d(required, common_addresses)
    if not np.isin(controls['Address'], common_addresses).all():
        raise ValueError('IDAT supported v2 still requires every control address')
    address_ok = np.concatenate([
        np.isin(type_i['AddressA'], common_addresses) & np.isin(type_i['AddressB'], common_addresses),
        np.isin(type_ii['AddressA'], common_addresses)])
    raw_ids = manifest.loc[address_ok, 'probe_id'].tolist()
    if not raw_ids:
        raise ValueError('IDAT supported v2 has no address-supported probes')
    # Both OOB pools are needed by this declared Noob method; do not replace a
    # missing pool or silently change the normalization input.
    if not {'Grn', 'Red'} <= set(manifest.loc[address_ok & (manifest['Type'] == 'I'), 'Color']):
        raise ValueError('IDAT supported v2 lacks an observed Type I out-of-band pool')

    def table(role, columns, allow_empty=False):
        frame = pd.read_csv(paths[role], dtype=str, keep_default_na=False)
        if list(frame.columns) != columns or (frame.empty and not allow_empty):
            raise ValueError('IDAT v2 table header or empty roster differs: '+role)
        return frame
    annotation = table('annotation_groups', ['probe_id', 'Type', 'Color'])
    if (annotation['probe_id'].duplicated().any()
            or not annotation['probe_id'].str.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.:+-]*').all()
            or not (((annotation['Type'] == 'I') & annotation['Color'].isin(['Grn', 'Red']))
                    | ((annotation['Type'] == 'II') & (annotation['Color'] == ''))).all()):
        raise ValueError('IDAT annotation groups contain duplicate identities or invalid type/color')
    joined = annotation.set_index('probe_id').reindex(full_ids)
    annotation_ok = joined['Type'].notna().to_numpy()
    matching = ((joined['Type'].to_numpy() == manifest['Type'].to_numpy())
                & (joined['Color'].to_numpy() == manifest['Color'].to_numpy()))
    if np.any(annotation_ok & ~matching):
        raise ValueError('IDAT annotation and manifest type/color groups conflict')
    retained = address_ok & annotation_ok
    corrected_ids = manifest.loc[retained, 'probe_id'].tolist()
    if not corrected_ids:
        raise ValueError('IDAT supported v2 has no annotation-supported corrected probes')
    disposition = np.where(~address_ok, 'address_absent',
                           np.where(~annotation_ok, 'annotation_absent', 'retained'))
    expected_membership = manifest.assign(
        address_supported=np.where(address_ok, 'TRUE', 'FALSE'),
        annotation_supported=np.where(annotation_ok, 'TRUE', 'FALSE'), disposition=disposition)
    membership = table('probe_membership', list(expected_membership.columns))
    if not membership.equals(expected_membership):
        raise ValueError('IDAT recorded probe exclusions or manifest order differ from independent derivation')
    missing_table = table('missing_addresses', ['address'], allow_empty=True)
    if list(missing_table['address']) != missing.astype(str).tolist():
        raise ValueError('IDAT recorded missing address roster/order differs from the binary')
    expected_compatibility = {'channel_headers': 'PASS', 'address_sets': 'PASS',
        'address_count': len(common_addresses), 'importer_values': 'PASS', 'manifest_mapping': 'PASS',
        'required_manifest_addresses': len(required), 'missing_manifest_addresses': len(missing),
        'probe_count': len(corrected_ids), 'control_count': len(controls),
        'typeI_count': len(type_i), 'typeII_count': len(type_ii), 'raw_array_ids': physical,
        'raw_probe_count': len(raw_ids), 'annotation_excluded_probe_count': int(np.sum(address_ok & ~annotation_ok)),
        'address_excluded_probe_count': int(np.sum(~address_ok)), 'corrected_probe_count': len(corrected_ids)}
    compatibility = processing.get('compatibility')
    if (not isinstance(compatibility, dict) or compatibility != expected_compatibility
            or any(type(compatibility.get(key)) is not int
                   for key, value in expected_compatibility.items() if type(value) is int)):
        raise ValueError('IDAT v2 worker compatibility counts differ from independent readback')
    type_i = type_i.set_index('probe_id'); type_ii = type_ii.set_index('probe_id')

    def matrices(roles, roster):
        readers = {role: _idat_csv_chunks(paths[role], ids, 'probe_id') for role in roles}
        offset = 0
        for chunks in zip_longest(*readers.values()):
            if any(chunk is None for chunk in chunks):
                raise ValueError('IDAT v2 matrices have unequal row counts')
            stop = offset+len(chunks[0]); expected = roster[offset:stop]
            if not expected or any(list(frame.index) != expected or list(frame.columns) != ids for frame in chunks):
                raise ValueError('IDAT v2 matrix probe/sample identities or roster order differ')
            yield offset, stop, expected, {role: frame.to_numpy(float) for role, frame in zip(readers, chunks)}
            offset = stop
        if offset != len(roster):
            raise ValueError('IDAT v2 matrix omits required roster rows')

    raw_detection = np.empty((len(raw_ids), len(ids)), dtype=float)
    raw_values = {role: [] for role in ('raw_meth', 'raw_unmeth')}
    for start, stop, names, values in matrices(('raw_meth', 'raw_unmeth', 'raw_detection_p'), raw_ids):
        expected_m, expected_u = np.empty_like(values['raw_meth']), np.empty_like(values['raw_unmeth'])
        i_mask = np.array([name in type_i.index for name in names]); ii_mask = ~i_mask
        i_rows = type_i.reindex(np.array(names)[i_mask]); ii_rows = type_ii.reindex(np.array(names)[ii_mask])
        for column, gsm in enumerate(ids):
            for color in ('Grn', 'Red'):
                chosen = i_rows['Color'].to_numpy() == color; positions = np.flatnonzero(i_mask)[chosen]
                index, means = channels[(gsm, color)]
                expected_m[positions, column] = means[index.get_indexer(i_rows.loc[chosen, 'AddressB'])]
                expected_u[positions, column] = means[index.get_indexer(i_rows.loc[chosen, 'AddressA'])]
            for color, target in [('Grn', expected_m), ('Red', expected_u)]:
                index, means = channels[(gsm, color)]
                target[ii_mask, column] = means[index.get_indexer(ii_rows['AddressA'])]
        if not np.array_equal(values['raw_meth'], expected_m) or not np.array_equal(values['raw_unmeth'], expected_u):
            raise ValueError('IDAT v2 raw methylated/unmethylated cells differ from binary manifest mappings')
        detection = values['raw_detection_p']
        if np.any((detection[np.isfinite(detection)] < 0) | (detection[np.isfinite(detection)] > 1)):
            raise ValueError('IDAT full-raw detection probabilities exceed the unit interval')
        raw_detection[start:stop] = detection
        for role in raw_values: raw_values[role].append(values[role])
    raw_medians = {role: np.median(np.concatenate(chunks), axis=0) for role, chunks in raw_values.items()}
    del raw_values
    selected_raw = pd.Index(raw_ids).get_indexer(corrected_ids)
    tallies = {name: np.zeros(len(ids), dtype=np.int64) for name in
        ('detection_failed_count', 'detection_nonfinite_count', 'beta_nonfinite_count', 'masked_count', 'observed_count')}
    max_delta = np.zeros(len(ids)); paired_finite = np.zeros(len(ids), dtype=np.int64)
    corrected_roles = ('matrix', 'beta_unmasked', 'meth_corrected', 'unmeth_corrected', 'detection_p', 'beta_single_array')
    for start, stop, names, values in matrices(corrected_roles, corrected_ids):
        detection = values['detection_p']; beta = values['beta_unmasked']
        if not np.array_equal(detection, raw_detection[selected_raw[start:stop]], equal_nan=True):
            raise ValueError('IDAT corrected detection values differ from the full-raw subset')
        with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
            m, u = np.maximum(values['meth_corrected'], 0), np.maximum(values['unmeth_corrected'], 0)
            expected_beta = m/(m+u+100)
        if (not np.allclose(beta, expected_beta, rtol=0, atol=1e-12, equal_nan=True)
                or np.any((beta[np.isfinite(beta)] < 0) | (beta[np.isfinite(beta)] > 1))):
            raise ValueError('IDAT v2 corrected beta arithmetic differs')
        single = values['beta_single_array']; finite = np.isfinite(beta)
        if (not np.array_equal(np.isfinite(single), finite)
                or not np.array_equal(single[~finite], beta[~finite], equal_nan=True)):
            raise ValueError('IDAT single-array versus joint nonfinite identities differ')
        difference = np.zeros_like(beta)
        np.subtract(single, beta, out=difference, where=finite)
        max_delta = np.maximum(max_delta, np.max(np.abs(difference), axis=0))
        paired_finite += finite.sum(axis=0)
        failed = np.isfinite(detection) & (detection > 0.01)
        det_na, beta_na = ~np.isfinite(detection), ~finite
        mask = failed | det_na | beta_na; expected = beta.copy(); expected[mask] = np.nan
        if not np.allclose(values['matrix'], expected, rtol=0, atol=1e-12, equal_nan=True):
            raise ValueError('IDAT v2 detection/nonfinite masking differs')
        for name, flags in [('detection_failed_count', failed), ('detection_nonfinite_count', det_na),
                            ('beta_nonfinite_count', beta_na), ('masked_count', mask), ('observed_count', ~mask)]:
            tallies[name] += flags.sum(axis=0)
    recorded_delta = np.array([processing['sample_vs_joint_max_abs_delta'][gsm] for gsm in ids])
    # R write.table emits 15 significant decimal digits even with digits=17.
    # For two beta values in [0,1], their subtraction can shift by <=1e-15
    # from that rounding, plus one binary spacing at 1 for parsing/subtraction.
    serialization_bound = 1e-15+np.spacing(1.0)
    if ((paired_finite == 0).any() or (max_delta > 1e-12+serialization_bound).any()
            or not np.allclose(recorded_delta, max_delta, rtol=0, atol=serialization_bound)):
        raise ValueError('IDAT saved single-array versus joint beta delta differs from its record')
    if (tallies['observed_count'] == 0).any():
        raise ValueError('IDAT candidate loses all measurements for a sample')
    qc = pd.read_csv(paths['sample_qc'], keep_default_na=False)
    if list(qc.columns) != ['sample_id', 'probe_count', *tallies, 'raw_median_M', 'raw_median_U'] or list(qc['sample_id']) != ids:
        raise ValueError('IDAT v2 sample QC identities/order differ')
    expected_qc = {'probe_count': np.full(len(ids), len(corrected_ids)), **tallies,
                  'raw_median_M': raw_medians['raw_meth'], 'raw_median_U': raw_medians['raw_unmeth']}
    if any(not np.array_equal(qc[name].to_numpy(float), values) for name, values in expected_qc.items()):
        raise ValueError('IDAT v2 sample QC counts or full-raw medians differ')
    summary = {'full_manifest_probe_count': len(full_ids), 'raw_probe_count': len(raw_ids),
        'corrected_probe_count': len(corrected_ids),
        'address_excluded_probe_count': int(np.sum(~address_ok)),
        'annotation_excluded_probe_count': int(np.sum(address_ok & ~annotation_ok)),
        'missing_manifest_addresses': len(missing), 'probe_membership_status': 'PASS',
        'detection_subset_status': 'PASS', 'single_array_beta_cells_compared': len(corrected_ids)*len(ids),
        'sample_vs_joint_max_abs_delta': dict(zip(ids, max_delta.tolist())),
        'single_array_delta_csv_rounding_bound': float(serialization_bound),
        'sample_vs_joint_scope': 'Saved beta matrices compared; same-minfi repeatability only, not independent normalization'}
    return len(corrected_ids), tallies, qc, summary


def verify_dnam_idat_provenance(provenance_path, matrix_path=None):
    """Verify the closed raw-IDAT chain and candidate arithmetic, not ssNoob validity.

    Independent binary/source parsers check raw identities and intensities.
    Manifest annotation, detection-P generation and normalization are recorded
    minfi evidence; arithmetic or same-library repeatability cannot validate them.
    """
    provenance_path = Path(provenance_path).resolve()
    def read_json(path):
        def invalid(value): raise ValueError('Nonfinite JSON token in IDAT evidence: '+value)
        try:
            value = json.loads(Path(path).read_text(), parse_constant=invalid)
        except OSError as error:
            raise FileNotFoundError('Cannot proceed - required data file not found: '+str(path)) from error
        if not isinstance(value, dict): raise ValueError('IDAT evidence must be a JSON object')
        return value
    def canonical(value): return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)
    def finite(value, label):
        if type(value) not in (int, float) or not np.isfinite(value) or value < 0:
            raise ValueError('Invalid nonnegative IDAT '+label)
    provenance = read_json(provenance_path)
    declared_profile = provenance.get('profile')
    if not isinstance(declared_profile, dict): raise ValueError('IDAT provenance profile must be an object')
    profile = _idat_verification_profile(declared_profile.get('name'))
    supported_v2 = profile['name'] == 'minfi_ssnoob_beta100_detection01_supported_v2'
    study = provenance.get('study_id', '')
    if (set(provenance) != {'schema','study_id','sample_ids','profile','files','elapsed_seconds','scientific_status','accepted_results'}
            or provenance['schema'] != 'clockbase-dnam-idat-provenance-v1'
            or not isinstance(study, str) or not re.fullmatch(r'GSE\d+', study)
            or canonical(provenance['profile']) != canonical(profile)
            or provenance['scientific_status'] != 'CANNOT VERIFY'
            or type(provenance['accepted_results']) is not int or provenance['accepted_results'] != 0):
        raise ValueError('Unsupported IDAT provenance identity, profile or scientific scope')
    finite(provenance['elapsed_seconds'], 'elapsed time')
    files = provenance['files']
    if not isinstance(files, dict) or not files: raise ValueError('IDAT provenance lacks file closure')
    checked = {}
    def declaration(record):
        if (not isinstance(record, dict) or set(record) != {'path','sha256','size_bytes'}
                or not isinstance(record['path'], str) or not Path(record['path']).is_absolute()
                or str(Path(record['path']).resolve()) != record['path']
                or not isinstance(record['sha256'], str) or not re.fullmatch(r'[0-9a-f]{64}', record['sha256'])
                or type(record['size_bytes']) is not int or record['size_bytes'] < 0):
            raise ValueError('IDAT file declaration is not canonical and complete')
        path = Path(record['path'])
        if not path.is_file():
            raise FileNotFoundError('Cannot proceed - required data file not found: '+str(path))
        if path.stat().st_size != record['size_bytes'] or sha(path) != record['sha256']:
            raise ValueError('IDAT declared source bytes changed: '+str(path))
        checked[str(path)] = record['sha256']
        return path
    paths = {role: declaration(record) for role, record in files.items()}
    if len(set(paths.values())) != len(paths): raise ValueError('IDAT roles alias the same physical file')
    base_roles = {'series','samples','protocol','acquisition','processing','implementation:rworker',
                  'implementation:inputs','worker_config','worker_stdout','worker_stderr'}
    if not base_roles <= paths.keys(): raise ValueError('IDAT provenance omits a required source role')
    pairs = _idat_literal_pairs(study, paths['series'], paths['samples'])
    ids = [pair['sample_id'] for pair in pairs]
    if provenance['sample_ids'] != ids: raise ValueError('IDAT sample roster/order differs from physical source')
    plan, acquisition, processing = (read_json(paths[name]) for name in ('protocol','acquisition','processing'))
    if (plan.get('schema') != 'clockbase-dnam-idat-plan-v1' or plan.get('study_id') != study
            or canonical(plan.get('profile')) != canonical(profile) or plan.get('pairs') != pairs
            or plan.get('source_files') != {name: files[name] for name in ('series','samples')}
            or plan.get('new_analysis_profile') is not True or plan.get('model_calls') != 0
            or plan.get('biological_units_inferred') is not False or plan.get('scientific_status') != 'CANNOT VERIFY'):
        raise ValueError('IDAT plan differs from independent source ownership or profile')
    ceilings = {'compressed_bytes_per_file':128_000_000, 'decompressed_bytes_per_file':256_000_000,
                'compressed_bytes_total':2_000_000_000, 'processing_seconds':3600, 'samples':64}
    limits = plan.get('limits', {})
    if set(limits) != set(ceilings) or any(type(limits[key]) is not int or not 0 < limits[key] <= cap for key, cap in ceilings.items()):
        raise ValueError('IDAT resource limits are absent or exceed the supported bounds')
    if len(ids) > limits['samples']: raise ValueError('IDAT sample bound exceeded')
    acquisition_roles = {kind+':'+gsm+':'+channel for gsm in ids for channel in ('Grn','Red')
                         for kind in ('compressed','raw','source_receipt','decompression')}
    if (acquisition.get('schema') != 'clockbase-dnam-idat-acquisition-v1' or acquisition.get('study_id') != study
            or acquisition.get('plan_path') != str(paths['protocol']) or acquisition.get('plan_sha256') != files['protocol']['sha256']
            or canonical(acquisition.get('profile')) != canonical(profile)
            or acquisition.get('files') != {role: files.get(role) for role in acquisition_roles}
            or set(acquisition.get('files', {})) != acquisition_roles
            or acquisition.get('model_calls') != 0 or acquisition.get('new_clock_runs') != 0
            or acquisition.get('scientific_status') != 'CANNOT VERIFY'):
        raise ValueError('IDAT acquisition is not closed over the frozen plan and channel files')
    reuse_fields = {'reused_acquisition_sha256', 'new_downloads'}
    if supported_v2 and reuse_fields.intersection(acquisition) and (
            not reuse_fields <= acquisition.keys()
            or not isinstance(acquisition.get('reused_acquisition_sha256'), str)
            or not re.fullmatch(r'[0-9a-f]{64}', acquisition['reused_acquisition_sha256'])
            or type(acquisition.get('new_downloads')) is not int or acquisition['new_downloads'] != 0):
        raise ValueError('IDAT v2 acquisition reuse requires a prior receipt hash and zero new downloads')
    matrix_names = {'matrix':'beta-qc.csv.gz', 'beta_unmasked':'beta-unmasked.csv.gz',
        'meth_corrected':'meth-corrected.csv.gz', 'unmeth_corrected':'unmeth-corrected.csv.gz',
        'raw_meth':'raw-meth.csv.gz', 'raw_unmeth':'raw-unmeth.csv.gz', 'detection_p':'detection-p.csv.gz'}
    outputs = {**matrix_names, 'manifest_I':'manifest-I.csv','manifest_II':'manifest-II.csv',
        'manifest_controls':'manifest-controls.csv','sample_qc':'sample-qc.csv',
        'session_info':'session-info.txt','function_source':'function-source.R.txt'}
    if supported_v2:
        outputs.update(annotation_groups='annotation-groups.csv', probe_membership='probe-membership.csv',
            missing_addresses='missing-addresses.csv', raw_detection_p='raw-detection-p.csv.gz',
            beta_single_array='beta-single-array.csv.gz')
    for gsm in ids:
        for channel in ('Grn','Red'):
            outputs['channel_header:'+gsm+':'+channel] = gsm+'-'+channel+'-header.json'
            outputs['channel_decoded:'+gsm+':'+channel] = gsm+'-'+channel+'-decoded.csv.gz'
    if set(files) != base_roles | acquisition_roles | set(outputs):
        raise ValueError('IDAT provenance has omitted, unknown or extra file roles')
    if (processing.get('schema') != 'clockbase-dnam-idat-processing-v1' or processing.get('status') != 'PASS'
            or processing.get('sample_ids') != ids or processing.get('array') != profile['array']
            or processing.get('annotation', {}).get('array') != profile['array']
            or canonical(processing.get('profile')) != canonical(profile) or processing.get('outputs') != outputs
            or processing.get('scientific_status') != 'CANNOT VERIFY' or processing.get('accepted_results') != 0
            or processing.get('author_method_reproduction') is not False or processing.get('model_calls') != 0
            or processing.get('new_clock_runs') != 0):
        raise ValueError('IDAT worker profile, outputs or scientific scope differs')
    root = provenance_path.parent
    for role, filename in {**outputs,'processing':'processing.json','worker_config':'worker-config.json',
                           'worker_stdout':'stdout.txt','worker_stderr':'stderr.txt'}.items():
        if paths[role] != root/filename: raise ValueError('IDAT worker output escapes its closed processing directory')
    config = read_json(paths['worker_config'])
    if config != {'acquisition':str(paths['acquisition']), 'output':str(root), 'profile':profile}:
        raise ValueError('IDAT worker configuration is not bound to the acquisition/profile')
    if matrix_path is not None and Path(matrix_path).resolve() != paths['matrix']:
        raise ValueError('IDAT caller matrix is not the provenance-bound matrix')
    packages = processing.get('packages', {})
    if not isinstance(packages, dict) or any(not re.fullmatch(r'\d+(?:\.\d+)+(?:[-+].+)?', str(packages.get(name, '')))
            for name in ('minfi','illuminaio','IlluminaHumanMethylationEPICmanifest')):
        raise ValueError('IDAT processing lacks its recorded package versions')
    deltas = processing.get('sample_vs_joint_max_abs_delta', {})
    if not isinstance(deltas, dict) or set(deltas) != set(ids): raise ValueError('IDAT repeatability record loses sample identities')
    for value in deltas.values():
        finite(value, 'recorded single/joint delta')
        if value > 1e-12: raise ValueError('Recorded single/joint ssNoob check did not pass')
    expected_acquisition_pairs, total_compressed, channels, physical, common_addresses = [], 0, {}, [], None
    for pair in pairs:
        gsm = pair['sample_id']; joined = {'sample_id':gsm,'channels':{}}; headers = {}
        for channel in ('Grn','Red'):
            suffix = gsm+':'+channel; source = pair['channels'][channel]
            compressed, raw = paths['compressed:'+suffix], paths['raw:'+suffix]
            if (compressed.name != source['filename'] or raw.name != source['filename'][:-3]
                    or paths['source_receipt:'+suffix] != Path(str(compressed)+'.source.json')
                    or paths['decompression:'+suffix] != Path(str(raw)+'.json')):
                raise ValueError('IDAT channel file locator differs from its literal source')
            receipt, marker = read_json(paths['source_receipt:'+suffix]), read_json(paths['decompression:'+suffix])
            if (receipt.get('url') != source['https_url'] or receipt.get('sha256') != files['compressed:'+suffix]['sha256']
                    or receipt.get('size_bytes') != compressed.stat().st_size
                    or marker != {'compressed_sha256':files['compressed:'+suffix]['sha256'],
                                  'sha256':files['raw:'+suffix]['sha256'],'size_bytes':raw.stat().st_size}):
                raise ValueError('IDAT acquisition receipt or decompression binding differs')
            finite(receipt.get('retrieved_unix'), 'source retrieval time')
            total_compressed += compressed.stat().st_size
            if not 0 < compressed.stat().st_size <= limits['compressed_bytes_per_file']:
                raise ValueError('IDAT compressed input violates its frozen resource bound')
            expanded = hashlib.sha256(); size = 0
            try:
                with gzip.open(compressed, 'rb') as handle:
                    for block in iter(lambda: handle.read(1024*1024), b''):
                        size += len(block)
                        if size > limits['decompressed_bytes_per_file']: raise ValueError('IDAT decompression exceeds its bound')
                        expanded.update(block)
            except (OSError, EOFError) as error:
                raise ValueError('Invalid IDAT compressed stream: '+str(compressed)) from error
            if size != raw.stat().st_size or expanded.hexdigest() != files['raw:'+suffix]['sha256']:
                raise ValueError('Compressed IDAT content differs from its declared uncompressed bytes')
            header, binary = _idat_binary_readback(raw)
            # The acquired EPIC files record RedGreen=0 for both named channels.
            # Retain this integer without assigning unsupported channel semantics;
            # channel ownership comes from the exact physical GEO declaration.
            if read_json(paths['channel_header:'+suffix]) != {'sample_id':gsm,'channel':channel,**header}:
                raise ValueError('IDAT exported header differs from the independent binary decoder')
            headers[channel] = header
            addresses = np.sort(binary['address'])
            if common_addresses is None: common_addresses = addresses
            if not np.array_equal(addresses, common_addresses): raise ValueError('IDAT channel address sets disagree')
            offset = 0
            for chunk in _idat_csv_chunks(paths['channel_decoded:'+suffix], ['Mean','SD','NBeads'], 'address'):
                stop = offset+len(chunk)
                if (stop > len(binary['address']) or list(chunk.index) != binary['address'][offset:stop].astype(str).tolist()
                        or any(not np.array_equal(chunk[name].to_numpy(float), binary[name][offset:stop]) for name in ('Mean','SD','NBeads'))):
                    raise ValueError('IDAT decoded channel rows/intensities differ from its binary')
                offset = stop
            if offset != header['nSNPsRead']: raise ValueError('IDAT decoded table omits binary addresses')
            channels[(gsm,channel)] = (pd.Index(binary['address']), binary['Mean'])
            joined['channels'][channel] = {'path':str(raw),'sha256':files['raw:'+suffix]['sha256'],**source}
        if any(headers['Grn'][name] != headers['Red'][name] for name in ('Barcode','ChipType','position','nSNPsRead')):
            raise ValueError('IDAT paired binary headers refer to different physical arrays')
        physical.append(headers['Grn']['Barcode']+':'+headers['Grn']['position'])
        expected_acquisition_pairs.append(joined)
    if (len(set(physical)) != len(ids) or acquisition.get('pairs') != expected_acquisition_pairs
            or acquisition.get('compressed_bytes') != total_compressed or total_compressed > limits['compressed_bytes_total']):
        raise ValueError('IDAT acquisition duplicates physical arrays or changes its source/channel/byte accounting')
    def manifest(role, expected_columns):
        frame = pd.read_csv(paths[role], dtype=str, keep_default_na=False)
        if list(frame.columns) != expected_columns or frame.empty:
            raise ValueError('IDAT manifest columns are unsupported: '+role)
        for name in [name for name in expected_columns if name.startswith('Address')]:
            if not frame[name].str.fullmatch(r'[1-9]\d*').all(): raise ValueError('Invalid IDAT manifest address')
            frame[name] = frame[name].astype(np.int64)
        return frame
    type_i = manifest('manifest_I',['probe_id','AddressA','AddressB','Color'])
    type_ii = manifest('manifest_II',['probe_id','AddressA'])
    controls = manifest('manifest_controls',['Address','Type','Color','ExtendedType'])
    if (not type_i['Color'].isin(['Grn','Red']).all() or (type_i['AddressA'] == type_i['AddressB']).any()
            or not {'NEGATIVE','NORM_A','NORM_C','NORM_G','NORM_T'} <= set(controls['Type'])):
        raise ValueError('IDAT manifest color, Type I addresses or control families differ')
    probe_ids = list(type_i['probe_id'])+list(type_ii['probe_id'])
    probe_set = set(probe_ids)
    if len(probe_ids) != len(probe_set) or any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.:+-]*', probe) for probe in probe_ids):
        raise ValueError('IDAT manifest probe identities are invalid or duplicated')
    subset_summary = {}
    if supported_v2:
        count, tallies, qc, subset_summary = _idat_verify_supported_outputs(
            paths, ids, type_i, type_ii, controls, common_addresses, channels, physical, processing)
        reused = reuse_fields <= acquisition.keys()
        subset_summary['acquisition_mode'] = 'reused' if reused else 'fresh'
        subset_summary['acquisition_reuse_scope'] = (
            'Source bytes independently reread; prior acquisition receipt identified by hash only' if reused else
            'Fresh acquisition; source bytes independently reread, no prior receipt declared')
    else:
        required_addresses = np.unique(np.concatenate([type_i['AddressA'],type_i['AddressB'],type_ii['AddressA'],controls['Address']]))
        if not np.isin(required_addresses, common_addresses).all(): raise ValueError('IDAT manifest/control address absent from binary inputs')
        expected_compatibility = {'channel_headers':'PASS','address_sets':'PASS','address_count':len(common_addresses),
            'importer_values':'PASS','manifest_mapping':'PASS','required_manifest_addresses':len(required_addresses),
            'missing_manifest_addresses':0,'probe_count':len(probe_ids),'control_count':len(controls),
            'typeI_count':len(type_i),'typeII_count':len(type_ii),'raw_array_ids':physical}
        if processing.get('compatibility') != expected_compatibility:
            raise ValueError('IDAT worker compatibility counts differ from independent readback')
        type_i = type_i.set_index('probe_id'); type_ii = type_ii.set_index('probe_id')
        readers = {role: _idat_csv_chunks(paths[role], ids, 'probe_id') for role in matrix_names}
        seen = set(); count = 0
        tallies = {name: np.zeros(len(ids),dtype=np.int64) for name in
                   ('detection_failed_count','detection_nonfinite_count','beta_nonfinite_count','masked_count','observed_count')}
        raw_medians = {'raw_meth':[], 'raw_unmeth':[]}
        for chunks in zip_longest(*readers.values()):
            if any(chunk is None for chunk in chunks): raise ValueError('IDAT matrices have unequal row counts')
            tables = dict(zip(readers, chunks)); names = list(chunks[0].index)
            if (not names or len(names) != len(set(names)) or seen.intersection(names)
                    or any(not chunk.index.equals(chunks[0].index) or list(chunk.columns) != ids for chunk in chunks)):
                raise ValueError('IDAT matrix probe/sample identities or order differ')
            if not set(names) <= probe_set: raise ValueError('IDAT matrix contains a probe absent from its manifest')
            seen.update(names); count += len(names)
            values = {role: frame.to_numpy(float) for role,frame in tables.items()}
            expected_m, expected_u = np.empty_like(values['raw_meth']), np.empty_like(values['raw_unmeth'])
            i_mask = np.array([name in type_i.index for name in names]); ii_mask = ~i_mask
            i_rows, ii_rows = type_i.reindex(np.array(names)[i_mask]), type_ii.reindex(np.array(names)[ii_mask])
            for column,gsm in enumerate(ids):
                for color in ('Grn','Red'):
                    chosen = i_rows['Color'].to_numpy() == color; positions = np.flatnonzero(i_mask)[chosen]
                    index, means = channels[(gsm,color)]
                    expected_m[positions,column] = means[index.get_indexer(i_rows.loc[chosen,'AddressB'])]
                    expected_u[positions,column] = means[index.get_indexer(i_rows.loc[chosen,'AddressA'])]
                for color,target in [('Grn',expected_m),('Red',expected_u)]:
                    index, means = channels[(gsm,color)]
                    target[ii_mask,column] = means[index.get_indexer(ii_rows['AddressA'])]
            if not np.array_equal(values['raw_meth'],expected_m) or not np.array_equal(values['raw_unmeth'],expected_u):
                raise ValueError('IDAT raw methylated/unmethylated cells do not match binary manifest mappings')
            for role in raw_medians: raw_medians[role].append(values[role])
            with np.errstate(invalid='ignore',divide='ignore',over='ignore'):
                m, u = np.maximum(values['meth_corrected'],0), np.maximum(values['unmeth_corrected'],0)
                expected_beta = m/(m+u+100)
            beta, detection = values['beta_unmasked'], values['detection_p']
            if (not np.allclose(beta,expected_beta,rtol=0,atol=1e-12,equal_nan=True)
                    or np.any((beta[np.isfinite(beta)] < 0) | (beta[np.isfinite(beta)] > 1))
                    or np.any((detection[np.isfinite(detection)] < 0) | (detection[np.isfinite(detection)] > 1))):
                raise ValueError('IDAT corrected beta arithmetic or probability range differs')
            failed = np.isfinite(detection) & (detection > 0.01)
            det_na, beta_na = ~np.isfinite(detection), ~np.isfinite(beta)
            mask = failed | det_na | beta_na; target = beta.copy(); target[mask] = np.nan
            if not np.allclose(values['matrix'],target,rtol=0,atol=1e-12,equal_nan=True):
                raise ValueError('IDAT detection/nonfinite masking differs from the frozen profile')
            for name,flags in [('detection_failed_count',failed),('detection_nonfinite_count',det_na),
                               ('beta_nonfinite_count',beta_na),('masked_count',mask),('observed_count',~mask)]:
                tallies[name] += flags.sum(axis=0)
        if seen != probe_set or count != len(probe_ids) or (tallies['observed_count'] == 0).any():
            raise ValueError('IDAT candidate loses manifest probes or all measurements for a sample')
        qc_columns = ['sample_id','probe_count',*tallies,'raw_median_M','raw_median_U']
        qc = pd.read_csv(paths['sample_qc'], keep_default_na=False)
        if list(qc.columns) != qc_columns or list(qc['sample_id']) != ids:
            raise ValueError('IDAT sample QC identities/order differ')
        expected_qc = {'probe_count': np.full(len(ids),count), **tallies,
            'raw_median_M':np.median(np.concatenate(raw_medians['raw_meth']),axis=0),
            'raw_median_U':np.median(np.concatenate(raw_medians['raw_unmeth']),axis=0)}
        if any(not np.array_equal(qc[name].to_numpy(float),values) for name,values in expected_qc.items()):
            raise ValueError('IDAT sample QC counts or raw-intensity medians differ')
    for path,digest in checked.items():
        if sha(path) != digest: raise ValueError('IDAT source changed during independent verification: '+path)
    if read_json(provenance_path) != provenance: raise ValueError('IDAT provenance changed during verification')
    result_files = {role:{'path':record['path'],'sha256':record['sha256']} for role,record in files.items()}
    result_files['provenance'] = {'path':str(provenance_path),'sha256':sha(provenance_path)}
    measurement = {'schema':'clockbase-dnam-idat-independent-readback-v1','status':'PASS',
        'profile':profile['name'],'array':profile['array'],'sample_count':len(ids),'probe_count':count,
        'channel_count':2*len(ids),'addresses_per_channel':len(common_addresses),
        'raw_channel_values_compared':2*len(ids)*len(common_addresses)*3,
        'raw_methylation_cells_compared':2*subset_summary.get('raw_probe_count',count)*len(ids),'beta_cells_compared':count*len(ids),
        'masked_cells':int(tallies['masked_count'].sum()),'observed_cells':int(tallies['observed_count'].sum()),
        'sample_qc':qc.to_dict(orient='records'),'binary_decoder':'independent Python IDAT v3 field reader',
        'channel_ownership':'Exact GEO supplementary URLs and filenames; raw RedGreen header retained uninterpreted',
        'raw_manifest_mapping_status':'PASS','beta_arithmetic_status':'PASS','qc_mask_status':'PASS',
        'normalization_validation':'CANNOT VERIFY','detection_probability_validation':'CANNOT VERIFY',
        'manifest_annotation_validation':'CANNOT VERIFY','specimen_identity_beyond_GEO':'CANNOT VERIFY',
        'normalization_scope':'Recorded minfi ssNoob; no independent normalization implementation or scientific validation',
        'sample_vs_joint_scope':'Recorded same-minfi repeatability only; separate matrices are not available to this verifier',
        'scientific_status':'CANNOT VERIFY','accepted_results':0,'model_calls':0,'clock_execution_performed':False}
    measurement.update(subset_summary)
    return {'study_id':study,'sample_ids':ids,'files':result_files,'measurement':measurement,'scientific_status':'CANNOT VERIFY'}


def verify_dnam_recovery(proposal, reconstruction):
    """Independently parse both tables and compare every row/value, not a sample."""
    from urllib.parse import urlparse
    from src.pipeline_inputs_2026_09_10 import beta_table_header, reconcile_recovery_metadata, recovery_alias_owners
    for record in proposal['source_files'].values():
        if sha(record['path']) != record['sha256']: raise ValueError('Recovery verification source changed')
    original = beta_table_header(proposal['source_files']['failed_matrix']['path'])
    if not original['empty']: raise ValueError('Original is not an empty molecular table')
    by_gsm, primary = reconcile_recovery_metadata(proposal['source_files']['metadata']['path'],
        proposal['source_files']['series']['path'], proposal['study_id'],
        proposal['source_files'].get('sample_records', {}).get('path'))
    crosswalk = proposal['crosswalk']
    primary_owners = recovery_alias_owners(primary)
    ids = [row['gsm_id'] for row in crosswalk]
    original_ids = original['header'][1:] if original['feature_column'] else original['header']
    if len(set(ids)) != len(ids) or set(ids) != set(by_gsm) or set(ids) != set(original_ids):
        raise ValueError('Recovery verification found inconsistent sample membership')
    for binding in crosswalk:
        row, token = by_gsm[binding['gsm_id']], binding['matrix_sample_id']
        if primary_owners.get(token) != {binding['gsm_id']}:
            raise ValueError('Recovery mapping alias is ambiguous in primary Series')
        if binding['mapping_rule'] == 'explicit_GSM':
            valid = token == binding['gsm_id']
        elif binding['mapping_rule'] == 'exact_sample_supplementary_idat_coordinate':
            value = binding['quote']
            match = re.fullmatch(r'(GSM\d+)_([^_]+(?:_[^_]+)*)_(?:Grn|Red)\.idat\.gz',
                                 Path(urlparse(value).path).name)
            valid = (binding['field'] == 'supplementary_file'
                     and value in primary[binding['gsm_id']].get('supplementary_file', [])
                     and match is not None and match[1] == binding['gsm_id'] and match[2] == token)
        else:
            value = row.get(binding['field'], '')
            if value != binding['quote']: raise ValueError('Recovery mapping quote differs from source')
            if binding['mapping_rule'] == 'exact_bracketed_source_title_token':
                valid = (binding['field'] == 'title' and primary[binding['gsm_id']].get('title') == [value]
                         and value.rstrip().endswith(']') and value.rsplit('[', 1)[-1].rstrip()[:-1] == token)
            elif binding['mapping_rule'] == 'explicit_source_sample_label':
                label, separator, value = value.partition(':')
                valid = (binding['field'].startswith('characteristics_ch1')
                         and binding['quote'] in primary[binding['gsm_id']].get('characteristics_ch1', [])
                         and bool(separator) and label.strip().casefold() in {'sample id','sample identifier','sample code'}
                         and value.strip() == token)
            else: valid = False
        if not valid: raise ValueError('Recovery mapping lacks an exact source binding')
    source = proposal['source_files']['supplement']['path']
    source_layout = beta_table_header(source)
    if source_layout != proposal['source_layout']:
        raise ValueError('Supplement layout changed before recovery verification')
    target = Path(reconstruction['path'])
    if sha(target) != reconstruction['sha256']: raise ValueError('Reconstructed file changed before verification')
    options = {'index_col': 0, 'chunksize': 4096, 'float_precision': 'round_trip', 'converters': {0: str},
               'na_values': ['', 'na', 'nan', 'null']}
    matrix_ids = [row['matrix_sample_id'] for row in crosswalk]
    if source_layout.get('format') == 'xlsx':
        left = _xlsx_recovery_chunks(source, source_layout, matrix_ids)
    elif source_layout['delimiter'] == 'space':
        if source_layout['feature_column']:
            raise ValueError('Whitespace recovery layout unexpectedly declares a feature header')
        left = pd.read_csv(source, sep=r'\s+', header=None, skiprows=1,
                           names=['ID_REF', *matrix_ids], quotechar='"', **options)
    else:
        left = pd.read_csv(source, sep=source_layout['delimiter'], **options)
    right = pd.read_csv(target, sep='\t', **options)
    count, missing, seen = 0, 0, set()
    categories = Counter()
    observed = np.zeros(len(ids), dtype=np.int64)
    for source_chunk, repaired_chunk in zip_longest(left, right):
        if source_chunk is None or repaired_chunk is None: raise ValueError('Recovered table silently loses/adds rows')
        if list(source_chunk.columns) != matrix_ids or list(repaired_chunk.columns) != ids:
            raise ValueError('Recovered columns differ from source-bound order')
        if not source_chunk.index.equals(repaired_chunk.index): raise ValueError('Recovered probe identities/order changed')
        probes = list(source_chunk.index)
        if len(probes) != len(set(probes)) or seen.intersection(probes) or any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.:+-]*', str(p)) for p in probes):
            raise ValueError('Recovery verifier found duplicate/invalid feature identities')
        for probe in probes:
            categories['cg_digits' if re.fullmatch(r'cg\d+', str(probe)) else 'ch_prefix' if str(probe).startswith('ch.') else 'other'] += 1
        a, b = source_chunk.to_numpy(float), repaired_chunk.to_numpy(float)
        if a.shape != b.shape or not np.array_equal(a, b, equal_nan=True):
            raise ValueError('Recovered molecular values differ from source')
        finite = np.isfinite(a)
        if np.isinf(a).any() or np.any((a[finite] < 0) | (a[finite] > 1)):
            raise ValueError('Source is not a finite/missing beta matrix within [0,1]')
        seen.update(probes); count += len(probes)
        observed += finite.sum(axis=0); missing += int((~finite).sum())
    if not count or (observed == 0).any(): raise ValueError('Source has no measured values for a sample')
    if count != reconstruction['feature_count'] or len(ids) != reconstruction['sample_count'] or missing != reconstruction['missing_values']:
        raise ValueError('Reconstruction counts differ from independent parser')
    if dict(zip(ids, observed.tolist())) != reconstruction['finite_values_per_sample']:
        raise ValueError('Per-sample measurement coverage differs')
    if dict(categories) != reconstruction['feature_id_categories']:
        raise ValueError('Source feature identity categories differ')
    for record in proposal['source_files'].values():
        if sha(record['path']) != record['sha256']: raise ValueError('Source changed during recovery verification')
    if sha(target) != reconstruction['sha256']: raise ValueError('Recovered table changed during verification')
    return {'status': 'PASS', 'scope': 'Lossless measured-beta input reconstruction with source-bound assay sample IDs',
            'source_rows': count, 'source_samples': len(ids), 'values_compared': count*len(ids),
            'missing_values': missing, 'value_changes': 0, 'sample_drops': 0, 'probe_drops': 0,
            'feature_id_categories': dict(categories), 'feature_annotation_status': 'CANNOT VERIFY',
            'parser': ('openpyxl cell reader and pandas float_precision=round_trip versus worksheet XML token reconstruction'
                       if source_layout.get('format') == 'xlsx' else
                       'pandas float_precision=round_trip versus csv/NumPy reconstruction'),
            'biological_unit_status': 'CANNOT VERIFY', 'clock_applicability': 'CANNOT VERIFY',
            'scientific_status': 'CANNOT VERIFY', 'clock_execution_performed': False,
            'accepted_result_error_rate': None}


def recover_dnam_input(study, failed_matrix, supplement, metadata, series, output, *, prepare_only=False, sample_source=None):
    """Resumable detection/proposal, reconstruction and separate verification."""
    from src.pipeline_inputs_2026_09_10 import propose_dnam_recovery, reconstruct_dnam_input
    from src.pipeline_runtime_2026_09_10 import StageStore
    output = Path(output).resolve()
    proposal = propose_dnam_recovery(study, failed_matrix, supplement, metadata, series, sample_source=sample_source)
    root = Path(__file__).resolve().parents[1]
    code = [Path(__file__), root/'src/pipeline_inputs_2026_09_10.py', root/'src/pipeline_runtime_2026_09_10.py',
            root/'audit_system/agent_workflow_2026_09_10.py', root/'audit_system/source_packet.py']
    config = {'proposal': proposal, 'implementation_sha256': {str(path): sha(path) for path in code},
              'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__, 'generation_calls': 0}
    if proposal['source_layout'].get('format') == 'xlsx':
        import openpyxl
        config['openpyxl'] = openpyxl.__version__
    output.mkdir(parents=True, exist_ok=True)
    with (output/'.run.lock').open('a') as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError('Another process owns this recovery run')
        store = StageStore(output, config)
        def propose(directory):
            write_new(directory/'recovery-proposal.json', proposal)
            return proposal
        frozen = store.run('00_detection_and_proposal', propose)
        if prepare_only: return {'prepared': True, 'study_id': study, 'source_samples': len(frozen['crosswalk']), 'generation_calls': 0}
        reconstructed = store.run('01_reconstruction', lambda directory: reconstruct_dnam_input(frozen, directory))
        def verify(directory):
            result = verify_dnam_recovery(frozen, reconstructed)
            write_new(directory/'verification.json', result)
            return {**result, 'recovered_input': reconstructed['path'], 'recovered_sha256': reconstructed['sha256'],
                    'generation_calls': 0, 'original_preserved': True}
        result = store.run('02_independent_verification', verify)
        for path, digest in config['implementation_sha256'].items():
            if sha(path) != digest: raise ValueError('Recovery implementation changed during execution')
        return result


def load_unavailable_comparison(case, expected_inputs, enrollment_sha256):
    """Pin one recorded upstream failure without inventing an analysis result.

    Call only while freezing a new comparison. Later availability does not
    change the saved attempt; its existing evidence is still rehashed on resume.
    """
    from datetime import datetime, timezone
    directory = Path(case['archive']).resolve()
    config_path = directory / 'config.json'
    failure_path = Path(case['failure_receipt']).resolve()
    if not failure_path.is_relative_to(directory):
        raise ValueError('Unavailable receipt escapes its upstream archive')
    relative = failure_path.relative_to(directory)
    allowed = {'01_metadata', '02_molecular_input', '03_blinded_design', '04_preprocessing',
               '05_clock_predictions', '06_coding_agent'}
    if (len(relative.parts) != 3 or relative.parts[0] not in allowed
            or not re.fullmatch(r'attempt-\d{3}', relative.parts[1]) or relative.parts[2] != 'failure.json'):
        raise ValueError('Unavailable case requires a recorded upstream prerequisite attempt')
    used = {}
    def pin(path):
        key = str(path.resolve())
        if key not in expected_inputs or sha(path) != expected_inputs[key]:
            raise ValueError('Unavailable-attempt evidence differs from its pinned identity: ' + key)
        used[key] = expected_inputs[key]
    pin(config_path); pin(failure_path)
    config, failure = json.loads(config_path.read_text()), json.loads(failure_path.read_text())
    if config.get('study') != case['study_id']:
        raise ValueError('Unavailable upstream study differs from enrollment')
    elapsed = failure.get('elapsed_seconds')
    if (failure.get('status') != 'CANNOT VERIFY' or not isinstance(failure.get('error'), str)
            or not failure['error'].strip() or type(elapsed) not in (int, float)
            or not np.isfinite(elapsed) or elapsed < 0):
        raise ValueError('Unavailable attempt lacks a valid recorded failure')
    stage = directory / relative.parts[0]
    attempts = sorted(stage.glob('attempt-*'))
    if not attempts or failure_path.parent != attempts[-1]:
        raise ValueError('Unavailable attempt is not the latest upstream attempt at enrollment')
    completed = [directory / name / 'complete.json' for name in allowed
                 if name >= relative.parts[0]]
    if any(path.exists() for path in completed):
        raise ValueError('Unavailable attempt conflicts with a completed upstream prerequisite')
    partial = {}
    for path in sorted(failure_path.parent.rglob('*')):
        if path.is_file():
            if not path.resolve().is_relative_to(failure_path.parent):
                raise ValueError('Unavailable partial artifact escapes the selected attempt')
            pin(path); partial[str(path)] = used[str(path.resolve())]
    return {'study_id': case['study_id'], 'family_id': case['family_id'],
            'preparation_status': 'UNAVAILABLE_SHARED_PREREQUISITE',
            'source_sha256': used,
            'attempt': {'schema': 'clockbase-unavailable-attempt-v1',
                'enrollment_sha256': enrollment_sha256, 'upstream_config': str(config_path),
                'upstream_config_sha256': used[str(config_path)], 'failure_receipt': str(failure_path),
                'failure_sha256': used[str(failure_path)], 'stage': relative.parts[0],
                'selected_attempt': relative.parts[1], 'recorded_failure': failure,
                'partial_artifact_sha256': partial,
                'absent_completed_prerequisites_at_freeze': [str(path) for path in sorted(completed)],
                'observation_utc': datetime.now(timezone.utc).isoformat(),
                'scope': 'Recorded shared prerequisite failure; partial files are not validated numerical outputs',
                'later_recovery_policy': 'Retain this frozen attempt; a recovered analysis requires a separate comparison'}}


def load_archived_comparison(case, expected_inputs):
    """Load one pinned original analysis without using its previous audit verdict."""
    from .agent_workflow_2026_09_10 import source_packet
    directory = Path(case['archive']).resolve()
    used = {}

    def read(path):
        path = path.resolve()
        expected = expected_inputs.get(str(path))
        if expected is None or sha(path) != expected:
            raise ValueError('Archived input is absent from or differs from its identity manifest: ' + str(path))
        used[str(path)] = expected
        return json.loads(path.read_text())

    config = read(directory/'config.json')
    if config['study'] != case['study_id']:
        raise ValueError('Archived case identity differs from enrollment')

    def source(stage, field):
        receipt = read(directory/stage/'complete.json')
        filename = Path(receipt['result'][field]).name
        names = [name for name in receipt['outputs'] if Path(name).name == filename]
        if len(names) != 1:raise ValueError('Missing/ambiguous archived source: ' + filename)
        path = (directory/stage/names[0]).resolve()
        if not path.is_relative_to(directory/stage):raise ValueError('Archived input escapes its stage')
        if sha(path) != receipt['outputs'][names[0]] or expected_inputs.get(str(path)) != sha(path):
            raise ValueError('Archived data differ from saved source identity: ' + str(path))
        used[str(path)] = sha(path)
        return receipt['result'], path

    metadata, metadata_path = source('01_metadata','metadata')
    clocks, clock_path = source('05_clock_predictions','path')
    primary, primary_path = source('06_coding_agent','path')
    packet = source_packet(metadata_path,case['study_id'],config['target'],context=case.get('context'),
                           comparison_scope=config.get('comparison_scope', 'exact_target'))
    if case.get('context'):
        context_path = Path(case['context']).resolve()
        used[str(context_path)] = sha(context_path)
        context = json.loads(context_path.read_text())
        for entry in context['sources']:
            used[str((context_path.parent/entry['path']).resolve())] = entry['sha256']
    from src.pipeline_2026_09_10 import read_primary_outputs
    manifest_path = source('06_coding_agent', 'manifest_path')[1] if primary.get('manifest_path') else None
    checked_primary = read_primary_outputs(primary_path, clocks, packet, manifest_path)
    rows = checked_primary['rows']
    if rows != primary['rows']:raise ValueError('Original analysis differs from preserved CSV')
    if manifest_path:
        if any(checked_primary[key] != primary[key] for key in ('schema', 'comparisons', 'unassigned_ids')):
            raise ValueError('Original comparison manifest differs from preserved receipt')
        primary['manifest_path'] = str(manifest_path)
    predictions = pd.read_csv(clock_path)
    if predictions.sample_id.duplicated().any() or set(predictions.sample_id) != set(metadata['sample_ids']):
        raise ValueError('Prediction/source sample membership differs')
    if set(predictions.sample_id) != {row['sample_id'] for row in packet['samples']}:
        raise ValueError('Design/source sample membership differs')
    primary['path'], clocks['path'] = str(primary_path), str(clock_path)
    return {'study_id':case['study_id'],'packet':packet,'primary':primary,'clocks':clocks,
            'source_sha256':used,'modality':config['modality'],
            'historical_release_state':'CANNOT VERIFY'}


def verify_study_families(evidence, output, extension_receipts=None):
    """Check preserved link witnesses and rebuild components with SciPy.

    This is source-identity/graph verification, not audit accuracy or a complete
    historical-corpus determination. No numerical study outcomes are loaded.
    """
    from collections import defaultdict
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from .evaluation_2026_09_10 import verify_family_enrollment
    evidence, output = Path(evidence).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((evidence/'manifest.json').read_text())
    for name, expected in manifest['output_sha256'].items():
        if sha(evidence/name) != expected:
            raise ValueError('Family manifest output changed: ' + name)
    def read(name):
        with gzip.open(evidence/name,'rt') as handle:return json.load(handle)
    families = read('families.json.gz'); records = read('source-records.json.gz')
    sources = read('source-catalog.json.gz'); exposures = read('development-exposure.json.gz')
    by_source = {record['source_id']: record for record in records}
    if len(by_source) != len(records):raise ValueError('Duplicate family source identity')
    if not set(by_source) <= {source['source_id'] for source in sources}:
        raise ValueError('Family graph cites unavailable source records')
    # Recheck the source files after the completed graph build, including input
    # files whose parsing failed. Failure is retained; it is not clean evidence.
    for index, source in enumerate(sources,1):
        if sha(source['path']) != source['sha256']:
            raise ValueError('Family input changed: ' + source['path'])
        if index % 10000 == 0:print(f'Family verification: {index}/{len(sources)} input hashes',file=sys.stderr,flush=True)
    exposure_ids = set()
    for source in exposures:
        if sha(source['path']) != source['sha256']:raise ValueError('Exposure source changed')
        exposure_ids.update(source['accessions'])
    if exposure_ids != set(families['exposed_studies']):raise ValueError('Exposure registry differs from families')
    verify_family_enrollment([],families,'development')
    studies = sorted(row['study_id'] for row in families['studies'])
    indices = {study:i for i,study in enumerate(studies)}
    owner, left, right = {}, [], []
    for record in records:
        index = indices[record['accession']]
        for field in ('sample_ids','pubmed_ids','bioproject_ids'):
            for value in record.get(field,[]):
                key = (field,value)
                if key in owner:
                    left.append(index);right.append(owner[key])
                else:owner[key]=index
        for accession in record.get('linked_accessions',[]):
            left.append(index);right.append(indices[accession])
    graph = coo_matrix((np.ones(len(left),dtype=bool),(left,right)),shape=(len(studies),len(studies))).tocsr()
    component_count, labels = connected_components(graph,directed=False)
    independently_grouped = defaultdict(list)
    for study,label in zip(studies,labels):independently_grouped[int(label)].append(study)
    expected = {tuple(sorted(family['members'])) for family in families['families']}
    observed = {tuple(values) for values in independently_grouped.values()}
    if expected != observed:raise ValueError('Independent source-graph components differ')
    for edge in families['edges']:
        cited = [by_source[source] for source in edge['source_ids']]
        if edge['kind']=='GEO_relation':
            if not any({record['accession'], linked}=={edge['study_a'],edge['study_b']}
                       for record in cited for linked in record.get('linked_accessions',[])):
                raise ValueError('GEO relation witness is not present in its cited source')
        else:
            field={'GSM':'sample_ids','PubMed':'pubmed_ids','BioProject':'bioproject_ids'}[edge['kind']]
            supported={record['accession'] for record in cited if edge['identifier'] in record.get(field,[])}
            if supported != {edge['study_a'],edge['study_b']}:
                raise ValueError('Shared-identifier witness is not present in both cited sources')
    newly_exposed=[row for row in families['studies'] if row['family_exposed'] and not row['directly_exposed']]
    write_new(output/'additional-family-exclusions.json',newly_exposed)
    failures=[source for source in sources if source['status']!='PASS']
    write_new(output/'source-parse-failures.json',failures)
    extension_rows=[]; receipt_sources={}
    if extension_receipts is not None:
        by_study={study:family for family in families['families'] for study in family['members']}
        local={record['accession'] for record in records if record['source_kind']=='sample_metadata'}
        for path in sorted(Path(extension_receipts).glob('GSE*.json')):
            receipt_sources[str(path.resolve())]=sha(path)
            record=json.loads(path.read_text());study=record['accession']
            if study not in by_study:raise ValueError('Extension study missing from family screening ledger')
            family=by_study[study];overlap=sorted(set(family['members'])&local)
            date=record.get('series_first_public_date','')
            extension_rows.append({'study':study,'series_public_date':date,
                'archived_date_agreement':record.get('release_date_agreement',False),
                'date_supported_2025':record.get('release_date_agreement',False) and '2025-01-01'<=date<='2025-12-31',
                'family_id':family['family_id'],'accession_in_available_local_metadata':study in local,
                'family_overlaps_available_local_metadata':bool(overlap),'local_family_members':overlap,
                'family_exposed':bool(family['exposed_members']),'source_linkage_status':family['source_linkage_status'],
                'final_extension_eligibility':'CANNOT VERIFY',
                'scope':'Local source overlap screen; original manuscript corpus boundary remains unverified'})
        write_new(output/'extension-family-screen.json',extension_rows)
    result={'status':'PASS','source_files_verified':len(sources),'exposure_files_verified':len(exposures),
            'source_parse_failures':len(failures),'studies':len(studies),'independent_graph_components':int(component_count),
            'source_link_witnesses_verified':len(families['edges']),'additional_family_exclusions':len(newly_exposed),
            'directly_exposed_studies':len(exposure_ids),'family_exposed_studies':sum(s['family_exposed'] for s in families['studies']),
            'extension_candidates_screened':len(extension_rows),'extension_receipt_sha256':receipt_sources,
            'extension_counts':dict(Counter(('2025' if row['date_supported_2025'] else 'other')+':'+
                ('local_family_overlap' if row['family_overlaps_available_local_metadata'] else 'no_observed_local_family_overlap')
                for row in extension_rows)),
            'family_manifest_sha256':sha(evidence/'families.json.gz'),'implementation_sha256':sha(__file__),
            'held_out_studies_evaluated':0,'audited_extension_studies':0,'accepted_result_error_rate':None,
            'clock_execution_performed':False,'model_calls':0,'evaluation_release':'CLOSED'}
    write_new(output/'verification.json',result)
    return result


def extract_clock_sites(series_matrix, metadata, sites, output):
    """Stream a GEO matrix, retaining measured clock sites and every source GSM.

    Checks row widths, unique probes, GEO identity, complete sample membership,
    table closure and source hashes. Unselected probe values are not validated
    or processed. Missing selected measurements remain missing.
    """
    from .agent_workflow_2026_09_10 import source_packet
    from src.pipeline_inputs_2026_09_10 import require_matrix
    series_matrix, metadata = Path(series_matrix).resolve(), Path(metadata).resolve()
    for path in (series_matrix, metadata):
        if not path.is_file():
            raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
    sources = {str(path): {'sha256': sha(path), 'bytes': path.stat().st_size,
               'mtime_ns': path.stat().st_mtime_ns} for path in (series_matrix, metadata)}
    wanted, retained, seen = set(sites), {}, set()
    accessions, source_samples, header, ended = [], [], None, False
    opener = gzip.open if series_matrix.suffix == '.gz' else open
    print('Streaming GEO table; retaining only requested clock CpGs', file=sys.stderr, flush=True)
    with opener(series_matrix, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('!Series_geo_accession\t'):
                accessions += next(csv.reader([line], delimiter='\t'))[1:]
            elif line.startswith('!Sample_geo_accession\t'):
                source_samples += next(csv.reader([line], delimiter='\t'))[1:]
            elif line.strip() == '!series_matrix_table_begin':
                header = next(csv.reader([next(handle)], delimiter='\t'))
                break
        if len(accessions) != 1 or not re.fullmatch(r'GSE\d+', accessions[0]):
            raise ValueError('Missing or ambiguous GEO series identity')
        if not header or header[0] != 'ID_REF' or len(header) < 2 or len(set(header)) != len(header):
            raise ValueError('Missing, empty or duplicated GEO matrix header')
        ids = header[1:]
        if source_samples != ids or any(not re.fullmatch(r'GSM\d+', sid) for sid in ids):
            raise ValueError('GEO source sample declarations differ from matrix header')
        for line in handle:
            if line.strip() == '!series_matrix_table_end':
                ended = True
                if any(line.strip() for line in handle):
                    raise ValueError('Unexpected content after GEO matrix end')
                break
            if line.count('\t') != len(ids):
                raise ValueError('GEO matrix row width differs from its header')
            probe = line.split('\t', 1)[0].strip().strip('"')
            if not probe or probe in seen:
                raise ValueError('Blank or duplicate GEO probe identity: ' + probe)
            seen.add(probe)
            if probe in wanted:
                values = next(csv.reader([line], delimiter='\t'))[1:]
                retained[probe] = [np.nan if value.strip().lower() in {'', 'null', 'na', 'nan'}
                                   else float(value) for value in values]
            if len(seen) % 100000 == 0:
                print(f'Read {len(seen)} source rows; retained {len(retained)} clock sites', file=sys.stderr, flush=True)
    if not ended or not seen:
        raise ValueError('GEO table is empty or its end marker is missing')
    packet = source_packet(metadata, accessions[0], 'Measured-feature clock implementation check; no intervention contrast')
    if set(ids) != {sample['sample_id'] for sample in packet['samples']}:
        raise ValueError('GEO matrix membership differs from source metadata')
    raw = pd.DataFrame.from_dict(retained, orient='index', columns=ids)
    require_matrix(raw, beta=True)
    for path, identity in sources.items():
        if sha(path) != identity['sha256']:
            raise ValueError('Source changed during measurement extraction: ' + path)
    raw_path = output / 'measured-clock-cpgs.csv'
    raw.rename_axis('ID_REF').to_csv(raw_path, mode='x')
    with metadata.open('rb') as source, (output / 'source-metadata.csv').open('xb') as dest:
        shutil.copyfileobj(source, dest)
    if sha(output / 'source-metadata.csv') != sources[str(metadata)]['sha256']:
        raise ValueError('Source metadata changed during archival copy')
    write_new(output / 'source-packet.json', packet)
    receipt = {'study': accessions[0], 'source_files': sources, 'source_url': packet['source_url'],
               'source_rows': len(seen), 'source_samples': len(ids), 'selected_sites': len(raw),
               'requested_sites': sorted(wanted), 'missing_sites': sorted(wanted-set(raw.index)),
               'path': str(raw_path), 'sha256': sha(raw_path), 'imputation': 'none',
               'sample_selection': 'All source samples; no outcome-dependent selection',
               'validation_scope': 'All source row identities and widths; numerical values only at retained clock sites'}
    write_new(output / 'extraction.json', receipt)
    return raw, packet, receipt


def check_measurement_profiles(series_matrix, metadata, output):
    """Recompute two blood clocks from source CpGs; no design or error-rate claim."""
    from biolearn.data_library import GeoData
    from biolearn.model_gallery import ModelGallery
    from .adapters.measurement_pilot import BLOOD_CLOCK_PROFILES, verify_clock_measurement
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    root = Path(__file__).resolve().parents[1]
    code = [Path(__file__), root/'audit_system/source_packet.py', root/'audit_system/agent_workflow_2026_09_10.py',
            root/'audit_system/adapters/measurement_pilot.py', root/'src/pipeline_inputs_2026_09_10.py']
    code_hashes = {str(path): sha(path) for path in code}
    models = {name: ModelGallery().get(name, imputation_method='none') for name in BLOOD_CLOCK_PROFILES}
    sites = {site for model in models.values() for site in model.coefficients.index if site != 'intercept'}
    raw, packet, extraction = extract_clock_sites(series_matrix, metadata, sites, output)
    predictions, checks = pd.DataFrame(index=raw.columns), []
    for name, model in models.items():
        coefficients = model.coefficients.iloc[:, 0]
        required = [site for site in coefficients.index if site != 'intercept']
        measured = raw.reindex(required)
        observed = measured.notna().sum(axis=0)
        complete = observed.index[observed == len(required)]
        predictions[name] = np.nan
        max_difference = None
        if len(complete):
            values = measured.loc[:, complete]
            calculated = model.predict(GeoData(metadata=pd.DataFrame(index=complete), dnam=values.copy()))['Predicted'].reindex(complete)
            independent = coefficients.loc[required].to_numpy(float) @ values.to_numpy(float) + float(coefficients.get('intercept', 0))
            if not np.isfinite(calculated).all() or not np.allclose(calculated, independent, rtol=1e-10, atol=1e-10):
                raise ValueError('Independent measured-feature prediction failed: ' + name)
            predictions.loc[complete, name] = calculated
            max_difference = float(np.max(np.abs(calculated.to_numpy(float)-independent)))
        coefficients_path = output / (name + '-coefficients.csv')
        coefficients.to_csv(coefficients_path, mode='x')
        checks.append({'clock': name, 'numerical_status': 'PASS' if len(complete) else 'CANNOT VERIFY',
                       'numerical_scope': 'Only samples with all required CpGs measured; all others retain missing predictions',
                       'required_features': len(required), 'observed_features_min': int(observed.min()),
                       'complete_measurement_samples': len(complete), 'incomplete_measurement_samples': len(raw.columns)-len(complete),
                       'max_abs_independent_difference': max_difference, 'measurement_sources': {
                           'raw_beta': {'path': extraction['path'], 'sha256': extraction['sha256']},
                           'coefficients': {'path': str(coefficients_path), 'sha256': sha(coefficients_path)}}})
    predictions_path = output / 'clock_predictions.csv'
    predictions.rename_axis('sample_id').to_csv(predictions_path, mode='x')
    clocks = {'path': str(predictions_path), 'clocks': list(models), 'sample_count': len(raw.columns),
              'checks': checks, 'clock_predictions_recomputed': True}
    rows, verification = [], []
    for check in checks:
        result = verify_clock_measurement(packet, clocks, check['clock'], list(raw.columns))
        verification.append(result)
        rows.append({'study': packet['study_id'], 'clock': check['clock'], 'source_samples': len(raw.columns),
                     **{key: check[key] for key in ('required_features', 'observed_features_min',
                        'complete_measurement_samples', 'incomplete_measurement_samples', 'max_abs_independent_difference')},
                     'measurement_status': result['status'], 'unresolved_checks': len(result['unresolved']),
                     'failed_checks': len(result['failures']), 'accepted_intervention_results': 0})
    write_new(output / 'clock-receipts.json', clocks)
    write_new(output / 'measurement-verification.json', verification)
    with (output / 'measurement-validation.csv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    table = ['| Study | Clock | Source samples | Complete measurements | Measured CpGs (min/required) | Measurement status |',
             '|---|---|---:|---:|---:|---|']
    for row in rows:
        table.append(f"| {row['study']} | {row['clock']} | {row['source_samples']} | {row['complete_measurement_samples']} | "
                     f"{row['observed_features_min']}/{row['required_features']} | {row['measurement_status']} |")
    (output / 'measurement-validation.md').write_text('\n'.join(table) + '\n\n'
        'Clocks were recalculated from measured source CpGs, without imputation. All source samples are in the denominator. '
        'This is development evidence for clock execution and the declared measurement profiles; no intervention design, '
        'aging-reversal claim, out-of-sample accuracy or audit error rate is established. '
        + ('GSE40279 is the original Hannum training cohort and was already development-exposed as a ClockBase imputation reference. '
           if packet['study_id'] == 'GSE40279' else 'Source cohort training-set overlap has not been adjudicated. ')
        + 'Measurement PASS is a method decision, not an independently adjudicated reference label.\n')
    for path in code:
        if sha(path) != code_hashes[str(path)]:
            raise ValueError('Implementation changed during measurement validation: ' + str(path))
        destination = output / 'executed-source' / path.relative_to(root)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with path.open('rb') as source, destination.open('xb') as dest:
            shutil.copyfileobj(source, dest)
    result = {'created_at_utc': datetime.now(timezone.utc).isoformat(), 'study': packet['study_id'], 'rows': rows,
              'clock_execution_performed': True, 'model_calls': 0, 'estimated_model_cost_usd': 0,
              'elapsed_seconds': time.monotonic()-started, 'development_only': True,
              'accepted_intervention_results': 0, 'accepted_result_error_rate': None, 'held_out_studies_evaluated': 0,
              'source_sha256': {path: identity['sha256'] for path, identity in extraction['source_files'].items()},
              'implementation_sha256': code_hashes,
              'versions': {name: version(name) for name in ('biolearn', 'numpy', 'pandas')},
              'python_version': sys.version, 'parameters': {'imputation': 'none', 'rtol': 1e-10, 'atol': 1e-10},
              'output_sha256': {str(path.relative_to(output)): sha(path) for path in output.rglob('*') if path.is_file()}}
    write_new(output / 'manifest.json', result)
    return result


def enrich_pipeline_sources(evidence, output, receipt_manifest=None):
    """Retrieve public Methods for archived development inputs, without outcomes."""
    from .agent_workflow_2026_09_10 import source_packet
    from src.pipeline_inputs_2026_09_10 import acquire_study_context
    evidence, output = Path(evidence).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    rows, provenance = [], {}
    started = time.monotonic()
    prior_hashes = {}
    if receipt_manifest is not None:
        receipt_manifest = Path(receipt_manifest).resolve()
        prior_hashes = json.loads(receipt_manifest.read_text())['source_sha256']
        provenance[str(receipt_manifest)] = sha(receipt_manifest)
    for directory in sorted(evidence.iterdir()):
        stage = directory / '01_metadata'
        if not (stage / 'complete.json').is_file() or not (directory / '08_report/complete.json').is_file():
            continue
        config_path = directory / 'config.json'
        config = json.loads(config_path.read_text())
        provenance[str(config_path)] = sha(config_path)
        receipt_path = stage / 'complete.json'
        if (stage / 'complete.sha256').is_file():
            checksum = (stage / 'complete.sha256').read_text().strip()
        else:
            candidates = [value for key, value in prior_hashes.items()
                          if Path(key).parts[-3:] == (directory.name, '01_metadata', 'complete.json')]
            if len(candidates) != 1:
                raise ValueError('Legacy receipt requires one matching --receipt-manifest hash: ' + str(stage))
            checksum = candidates[0]
        if sha(receipt_path) != checksum:
            raise ValueError('Archived metadata receipt changed: ' + str(stage))
        receipt = json.loads(receipt_path.read_text())
        provenance[str(receipt_path)] = sha(receipt_path)

        def source(field):
            filename = Path(receipt['result'][field]).name
            names = [p for p in receipt['outputs'] if Path(p).name == filename]
            if len(names) != 1:
                raise ValueError('Archived metadata source missing/ambiguous: ' + filename)
            path = (stage / names[0]).resolve()
            if not path.is_relative_to(stage) or sha(path) != receipt['outputs'][names[0]]:
                raise ValueError('Archived metadata source changed: ' + str(path))
            provenance[str(path)] = sha(path)
            return path

        metadata, series = source('metadata'), source('series_source')
        destination = output / directory.name
        context_path = acquire_study_context(config['study'], series, destination)
        packet = source_packet(metadata, config['study'], config['target'], context=context_path)
        write_new(destination / 'source-packet.json', packet)
        methods = [doc for doc in packet['study_documents'] if doc['source_kind'] == 'publication_methods']
        rows.append({'study': config['study'], 'samples': len(packet['samples']),
                     'methods_publications': len({doc['source_sha256'] for doc in methods}),
                     'methods_sections': len(methods),
                     'source_acquisition_gaps': sum(row['status'] == 'CANNOT VERIFY' for row in packet['source_acquisition']),
                     'packet_bytes': len(json.dumps(packet, ensure_ascii=False).encode()),
                     'design_reassessed': False, 'scientific_status': 'CANNOT VERIFY'})
    if not rows:
        raise ValueError('No archived source-metadata stages found')
    with (output / 'source-context.csv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    table = ['| Study | Samples | Methods publications | Methods sections | Acquisition gaps | Packet bytes |',
             '|---|---:|---:|---:|---:|---:|']
    for row in rows:
        table.append('| ' + ' | '.join(str(row[k]) for k in ('study','samples','methods_publications',
                      'methods_sections','source_acquisition_gaps','packet_bytes')) + ' |')
    (output / 'source-context.md').write_text('\n'.join(table) + '\n\n'
        'Archived development metadata with newly retrieved public Methods. Source identity, sample membership and '
        'quotation extraction were checked. Acquisition gaps count recorded lookup/extraction failures, not '
        'scientific errors or unique missing papers. No model/clock calls, design reassessment or held-out evaluation.\n')
    root = Path(__file__).resolve().parents[1]
    code = [Path(__file__), root/'src/pipeline_inputs_2026_09_10.py',
            root/'audit_system/source_packet.py', root/'audit_system/agent_workflow_2026_09_10.py']
    for path in code:
        destination = output / 'executed-source' / path.resolve().relative_to(root)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
    result = {'created_at_utc': datetime.now(timezone.utc).isoformat(), 'study_count': len(rows),
              'studies': rows, 'model_calls': 0, 'estimated_model_cost_usd': 0,
              'clock_execution_performed': False, 'held_out_studies_evaluated': 0,
              'design_reassessed': False, 'accepted_result_error_rate': None,
              'elapsed_seconds': time.monotonic()-started, 'source_sha256': provenance,
              'implementation_sha256': {str(path): sha(path) for path in code},
              'output_sha256': {str(path.relative_to(output)): sha(path) for path in output.rglob('*') if path.is_file()}}
    write_new(output / 'manifest.json', result)
    return result


def recheck_pipeline_audits(evidence, output):
    """Re-audit archived predictions and analyses; do not rerun clocks/models."""
    from .agent_workflow_2026_09_10 import parse_design, source_packet
    from src.pipeline_2026_09_10 import numerical_audit, read_primary_outputs
    evidence, output = Path(evidence).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    rows, provenance = [], {}
    started = time.monotonic()
    for directory in sorted(evidence.iterdir()):
        if not (directory / '08_report/complete.json').is_file():
            continue

        def read(path):
            provenance[str(path)] = sha(path)
            return json.loads(path.read_text())

        def source_file(stage, field):
            receipt = read(directory / stage / 'complete.json')
            original = Path(receipt['result'][field])
            names = [name for name in receipt['outputs'] if Path(name).name == original.name]
            if len(names) != 1:
                raise ValueError('Ambiguous archived stage source: ' + str(original))
            path = (directory / stage / names[0]).resolve()
            if not path.is_relative_to(directory / stage):
                raise ValueError('Archived source escapes stage directory')
            if not path.is_file():
                raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
            if sha(path) != receipt['outputs'][names[0]]:
                raise ValueError('Archived source differs from execution receipt: ' + str(path))
            provenance[str(path)] = sha(path)
            return receipt['result'], path

        config = read(directory / 'config.json')
        metadata, metadata_path = source_file('01_metadata', 'metadata')
        clocks, clock_path = source_file('05_clock_predictions', 'path')
        primary, primary_path = source_file('06_coding_agent', 'path')
        design_receipt = read(directory / '03_blinded_design/complete.json')
        names = [name for name in design_receipt['outputs'] if Path(name).name == 'design.json']
        if len(names) != 1:
            raise ValueError('Missing or ambiguous archived design')
        design_path = (directory / '03_blinded_design' / names[0]).resolve()
        if not design_path.is_relative_to(directory / '03_blinded_design'):
            raise ValueError('Archived design escapes stage directory')
        if sha(design_path) != design_receipt['outputs'][names[0]]:
            raise ValueError('Archived design differs from execution receipt')
        design = read(design_path)
        if design != design_receipt['result']:
            raise ValueError('Design result differs from saved response')
        context_path = source_file('01_metadata', 'study_context')[1] if metadata.get('study_context') else None
        packet = source_packet(metadata_path, config['study'], config['target'], context=context_path)
        parse_design(json.dumps(design), packet)
        manifest_path = source_file('06_coding_agent', 'manifest_path')[1] if primary.get('manifest_path') else None
        checked_primary = read_primary_outputs(primary_path, clocks, packet, manifest_path)
        original_rows = checked_primary['rows']
        if original_rows != primary['rows']:
            raise ValueError('Primary result differs from the saved analysis CSV')
        if manifest_path:
            if any(checked_primary[key] != primary[key] for key in ('schema', 'comparisons', 'unassigned_ids')):
                raise ValueError('Primary comparison manifest differs from receipt')
            primary['manifest_path'] = str(manifest_path)
        predictions = pd.read_csv(clock_path)
        if predictions.sample_id.duplicated().any() or set(predictions.sample_id) != set(metadata['sample_ids']):
            raise ValueError('Archived prediction samples differ from source metadata')
        if len(predictions) != clocks['sample_count']:
            raise ValueError('Archived prediction sample count differs from receipt')
        target = output / directory.name
        target.mkdir()
        primary['path'], clocks['path'] = str(primary_path), str(clock_path)
        audit = numerical_audit(design, primary, clocks, target, packet=packet)
        write_new(target / 'result.json', audit)
        states = Counter(check['calculation_status'] for check in audit['checks'])
        primary_states = Counter(check['primary_analysis_state']['state'] for check in audit['checks'])
        target_states = Counter(check['primary_target_binding']['status'] for check in audit['checks'])
        rows.append({'study': config['study'], 'modality': config['modality'],
                     'clock_results': len(audit['checks']), 'calculation_pass': states['PASS'],
                     'calculation_fail': states['FAIL'], 'calculation_unresolved': states['CANNOT VERIFY'],
                     'computed_rows': primary_states['COMPUTED'], 'abstained_rows': primary_states['ABSTAINED'],
                     'inconsistent_rows': primary_states['INCONSISTENT'], 'target_binding_pass': target_states['PASS'],
                     'target_binding_fail': target_states['FAIL'], 'target_binding_unresolved': target_states['CANNOT VERIFY'],
                     'accepted_results': audit['accepted_results'], 'scientific_status': audit['scientific_status'],
                     'source': str(directory), 'result': str(target / 'result.json')})
    if not rows:
        raise ValueError('No complete archived runs found')
    with (output / 'audit-recheck.csv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    table = ['| Study | Clock rows | Computed | Abstained | Inconsistent | Target bound | Calculation PASS | Accepted | Scientific status |',
             '|---|---:|---:|---:|---:|---:|---:|---:|---|']
    for row in rows:
        table.append('| ' + ' | '.join(str(row[k]) for k in ('study', 'clock_results', 'computed_rows',
            'abstained_rows', 'inconsistent_rows', 'target_binding_pass', 'calculation_pass',
            'accepted_results', 'scientific_status')) + ' |')
    (output / 'audit-recheck.md').write_text('\n'.join(table) +
        '\n\nDevelopment replay of archived predictions and primary analyses. Conditional arithmetic checks do not validate study design or aging interpretation. No new clock/model execution or held-out error-rate estimate.\n')
    root = Path(__file__).resolve().parents[1]
    implementations = [Path(__file__), root/'src/pipeline_2026_09_10.py',
                      *[root/'audit_system'/name for name in ('adapters/measurement_pilot.py','claim_gate.py','policy.py','models.py',
                         'ledger.py','hashing.py','dependency_graph.py','agent_workflow_2026_09_10.py','source_packet.py')]]
    result = {'created_at_utc': datetime.now(timezone.utc).isoformat(), 'studies': rows,
              'study_count': len(rows), 'clock_execution_performed': False, 'model_calls': 0,
              'estimated_model_cost_usd': 0, 'elapsed_seconds': time.monotonic()-started,
              'held_out_studies_evaluated': 0, 'accepted_result_error_rate': None,
              'source_sha256': provenance, 'implementation_sha256': {str(p): sha(p) for p in implementations},
              'output_sha256': {str(p.relative_to(output)): sha(p) for p in output.rglob('*') if p.is_file()}}
    write_new(output / 'manifest.json', result)
    return result


def build_pipeline_report(evidence, output):
    """Summarize preserved development runs; no additional model/clock calls."""
    output.mkdir(parents=True, exist_ok=False)
    rows, features, provenance = [], [], {}
    for directory in sorted(evidence.iterdir()):
        if not (directory / '08_report/complete.json').exists(): continue
        def read(name):
            path = directory / name
            provenance[str(path.resolve())] = sha(path)
            return json.loads(path.read_text())
        config = read('config.json')
        result = read('08_report/complete.json')['result']
        clocks = read('05_clock_predictions/complete.json')['result']
        source = read('01_metadata/complete.json')['result']
        stages = [read(str(p.relative_to(directory))) for p in sorted(directory.glob('0*/complete.json'))]
        calls = [read(str(p.relative_to(directory))) for p in sorted(directory.glob('api/*/receipt.json'))]
        if len(stages) != 8 or len(calls) != result['generation_calls']:
            raise ValueError('Incomplete source-to-audit receipts: ' + str(directory))
        clock_file = next(directory.glob('05_clock_predictions/attempt-*/clock_predictions.csv'))
        table = pd.read_csv(clock_file)
        if len(table) != result['sample_count'] or table.sample_id.duplicated().any():
            raise ValueError('Prediction sample count/identity differs from result')
        if set(table.sample_id) != set(source['sample_ids']):
            raise ValueError('Prediction samples differ from source GEO membership')
        provenance[str(clock_file.resolve())] = sha(clock_file)
        original_location = Path(clocks['path'])
        prediction_receipt = read('05_clock_predictions/complete.json')
        expected_key = str(original_location.relative_to(original_location.parent.parent))
        if prediction_receipt['outputs'].get(expected_key) != sha(clock_file):
            raise ValueError('Prediction content differs from execution receipt')
        estimate = sum(c['estimated_cost_usd'] for c in calls)
        if not np.isclose(estimate, result['estimated_cost_usd'], rtol=0, atol=1e-12):
            raise ValueError('Usage cost differs from reported cost')
        rows.append({'study':config['study'], 'modality':config['modality'], 'samples':len(table),
                     'clocks':len(clocks['clocks']), 'completed_stages':len(stages),
                     'execution_status':result['execution_status'], 'scientific_status':result['scientific_status'],
                     'model_calls':len(calls), 'estimated_cost_usd':round(estimate,8),
                     'stage_elapsed_seconds':sum(s['elapsed_seconds'] for s in stages),
                     'accepted_results':result['accepted_results'], 'accepted_result_error_rate':None,
                     'source':str(directory.relative_to(evidence))})
        for check in clocks['checks']:
            numerator, denominator = check['observed_features_min'], check['required_features']
            if not 0 <= numerator <= denominator or denominator <= 0: raise ValueError('Invalid feature coverage')
            features.append({'study':config['study'], 'modality':config['modality'], 'clock':check['clock'],
                             'samples':len(table), 'observed_features_min':numerator, 'required_features':denominator,
                             'min_observed_percent':100*numerator/denominator,
                             'max_abs_independent_prediction_difference':check['max_abs_independent_difference'],
                             'scientific_status':check['measurement_applicability']})
    if not rows: raise ValueError('No complete source-linked runs found')
    for name, records in [('pipeline-validation.csv',rows),('feature-coverage.csv',features)]:
        with (output/name).open('x') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    lines=['| Study | Modality | Samples | Clocks | Stages | Calls | Estimated USD | Scientific status |',
           '|---|---|---:|---:|---:|---:|---:|---|']
    for row in rows:
        lines.append('| '+' | '.join(str(row[k]) for k in ['study','modality','samples','clocks','completed_stages','model_calls','estimated_cost_usd','scientific_status'])+' |')
    (output/'pipeline-validation.md').write_text('\n'.join(lines)+'\n\nDevelopment execution evidence; no held-out error rate. Stage time excludes source-file hashing before stage 1 and package installation. Cost estimates are not invoices.\n')
    style = Path(__file__).resolve().parents[1]/'pipeline_assets/figure_style_2026_09_10'
    style_manifest = json.loads((style/'manifest.json').read_text())
    for name, record in style_manifest.items():
        if sha(style/name) != record['sha256']: raise ValueError('Plotting utility snapshot changed')
    sys.path.insert(0,str(style))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from nature_figure_utils import retro_style, trim_axes, figure_preflight, INK
    retro_style()
    fig,ax=plt.subplots(figsize=(5.3,2.15),layout='constrained')
    labels=[]
    for row in features:
        clock=row['clock'].replace('_All_All_scaleddiff_EN','').replace('Hazard_Mouse','Hazard Mouse').replace('Hazard_Rodents','Hazard Rodents').replace('Hazard2_All','Hazard2 All')
        labels.append(f"{row['study']} · {clock} ({row['observed_features_min']:,}/{row['required_features']:,})")
    ax.scatter([r['min_observed_percent'] for r in features], np.arange(len(features)),
               s=22, facecolors='white', edgecolors=INK, linewidths=.8, clip_on=False)
    ax.set_yticks(np.arange(len(features)), labels)
    ax.set_xticks([0,25,50,75,100]);ax.set_xlim(-2,102);ax.set_ylim(-.55,len(features)-.45)
    ax.set_xlabel('Observed required features (minimum across samples, %)')
    ax.set_title(f'Source coverage ($n = {len(rows)}$ development studies)',loc='left',fontsize=9,fontweight='bold')
    trim_axes(fig)
    ax.set_ylim(-.45,len(features)-.55)
    ax.invert_yaxis()
    qc=figure_preflight(fig,{'coverage':ax},verbose=False)
    write_new(output/'figure-preflight.json',qc)
    if not qc['passed']:
        raise ValueError('Figure preflight failed: '+qc['text'])
    fig.savefig(output/'feature-coverage.pdf',bbox_inches='tight')
    fig.savefig(output/'feature-coverage.png',dpi=200,bbox_inches='tight')
    plt.close(fig)
    (output/'feature-coverage-caption.md').write_text('Observed measurement coverage in development runs. Each point is the smallest fraction of required clock features observed across the source samples before imputation; labels show numerator/denominator. The counts come from the preserved prediction receipts. This is neither clock applicability validation nor audit accuracy. The table reports the number of samples for each study. Missing measurements are retained as an applicability question; no biological finding was accepted.\n')
    write_new(output/'manifest.json',{'source_sha256':provenance,'script_sha256':sha(__file__),
              'outputs':{p.name:sha(p) for p in output.iterdir() if p.is_file()},
              'development_runs':len(rows),'generation_calls':sum(r['model_calls'] for r in rows),
              'estimated_cost_usd':sum(r['estimated_cost_usd'] for r in rows),
              'held_out_cases_evaluated':0,'accepted_result_error_rate':None})


def build(evidence,geo_root,output):
    output.mkdir(parents=True,exist_ok=False)
    previous=evidence.parent/'systematic-audit-2026-09-09/verified-source-census'
    inventory=previous/'empty-matrix-verification.jsonl'
    manifest=json.loads((previous/'manifest.json').read_text())
    if sha(inventory)!=manifest['output_sha256'][inventory.name]:
        raise ValueError('Independent source inventory changed')
    producer=geo_root/'src/4_compute_age.py'
    tree=ast.parse(producer.read_text())
    functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='require_observed_beta_matrix']
    if len(functions)!=1:raise ValueError('Producer guard not found unambiguously')
    namespace={'pd':pd,'np':np}
    exec(compile(ast.Module(body=functions,type_ignores=[]),str(producer),'exec'),namespace)
    guard=namespace['require_observed_beta_matrix']
    started=time.monotonic();checks=[]
    for line in inventory.read_text().splitlines():
        row=json.loads(line);path=Path(row['matrix_path'])
        if sha(path)!=row['matrix_sha256']:raise ValueError('Empirical matrix changed: '+str(path))
        frame=pd.read_csv(path,index_col=0)
        error=None
        try:guard(frame)
        except ValueError as exception:error=str(exception)
        if sha(path)!=row['matrix_sha256']:raise ValueError('Empirical matrix changed during guard check')
        checks.append({'study_id':row['study_id'],'matrix_path':str(path),'matrix_sha256':row['matrix_sha256'],
            'probe_rows':frame.shape[0],'sample_columns':frame.shape[1],
            'guard_rejected':bool(error and error.startswith('AUDIT_FAIL: no observed feature rows')),
            'reason':error,'clock_execution_performed':False})
    write_new(output/'producer-source-rechecks.json',{'checks':checks,'elapsed_seconds':time.monotonic()-started,
        'producer_sha256':sha(producer),'inventory_sha256':sha(inventory),'scientific_status':'CANNOT VERIFY'})
    calls=[]
    for marker in sorted(evidence.rglob('request-started.json')):
        directory=marker.parent
        config=json.loads((directory/'config.json').read_text())
        result=json.loads((directory/'result.json').read_text())
        calls.append({'run':str(directory.relative_to(evidence)),'study_id':config['study'],
            'provider':config['provider'],'model':config['model'],'response_received':(directory/'response.json').exists(),
            'design_proposed':result.get('design_status')=='DESIGN_PROPOSED',
            'elapsed_seconds':result['elapsed_seconds'],'usage':result.get('usage'),
            'estimated_cost_usd':result.get('estimated_cost_usd'),'billed_cost_usd':result.get('billed_cost_usd'),
            'error':result.get('error'),'source_sha256':config['inputs']['metadata']['sha256'],
            'scientific_status':result['scientific_status']})
    write_new(output/'model-development-calls.json',calls)
    extension=json.loads((evidence/'geo-series-validation/summary.json').read_text())
    rows=[]
    def add(measure,value,denominator,scope,source):
        rows.append(dict(measure=measure,value=value,denominator=denominator,scope=scope,source=str(source)))
    scan=manifest['source_summary']
    add('DNAm files with zero feature rows',len(checks),scan['attempted_files'],'Local source census; not scientific error rate',previous/'manifest.json')
    add('Distinct DNAm source accessions',scan['study_count'],'','Local directory; original corpus boundary unverified',previous/'manifest.json')
    add('Observed empty inputs rejected by producer guard',sum(c['guard_rejected'] for c in checks),len(checks),'No clock execution; selected observed defects only',output/'producer-source-rechecks.json')
    for provider in sorted({c['provider'] for c in calls}):
        subset=[c for c in calls if c['provider']==provider]
        add(provider+' generation responses received',sum(c['response_received'] for c in subset),len(subset),'Transport development; not model accuracy',output/'model-development-calls.json')
    add('OpenAI conservative estimated cost USD',round(sum(c['estimated_cost_usd'] or 0 for c in calls if c['provider']=='openai'),6),'','List-rate estimate; billing separate, NeoLink cost unknown',output/'model-development-calls.json')
    add('DNAm first-public dates concordant',extension['date_agreement_studies'],extension['candidate_studies'],'Series status versus GEO PDAT',evidence/'geo-series-validation/summary.json')
    for year,count in sorted(extension['by_verified_release_year'].items()):
        add('DNAm candidates first public in '+year,count,'','Family/measurement eligibility incomplete',evidence/'geo-series-validation/summary.json')
    add('Old DNAm GSM identifiers screened',extension['old_sample_ids_inspected'],'','Available metadata; not complete historical family map',evidence/'geo-series-validation/summary.json')
    for species in ('human','mouse'):
        record=json.loads((evidence/f'archs4-{species}_overlap.json').read_text())
        for field in ('old_samples','new_samples','added_samples','removed_samples','new_accessions','new_accessions_with_shared_old_samples'):
            add('ARCHS4 '+species+' '+field,record[field],'','Sample/accession inventory; study dates unverified',evidence/f'archs4-{species}_overlap.json')
    with (output/'validation-table.csv').open('x',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    table=['| Measurement | Value | Denominator | Scope |','|---|---:|---:|---|']
    table += ['| '+ ' | '.join(str(r[k]) for k in ('measure','value','denominator','scope'))+' |' for r in rows]
    (output/'validation-table.md').write_text('\n'.join(table)+'\n')
    sources=[inventory,previous/'manifest.json',producer,evidence/'geo-series-validation/summary.json']
    sources+=list(evidence.glob('archs4-*_overlap.json'))
    for marker in evidence.rglob('request-started.json'):
        sources.extend([marker.parent/'result.json',marker.parent/'config.json'])
    write_new(output/'manifest.json',{'recorded_at_utc':datetime.now(timezone.utc).isoformat(),
        'script_sha256':sha(__file__),'source_sha256':{str(p.resolve()):sha(p) for p in sources},
        'output_sha256':{p.name:sha(p) for p in output.iterdir() if p.is_file()},
        'software_versions':{'numpy':np.__version__,'pandas':pd.__version__},
        'held_out_studies_evaluated':0,'independent_scientific_error_rate':None,
        'clock_execution_performed':False})
    print(json.dumps({'empty_inputs_checked':len(checks),'guard_rejected':sum(c['guard_rejected'] for c in checks),
        'generation_attempts':len(calls),'responses_received':sum(c['response_received'] for c in calls),
        'output':str(output)}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence',type=Path,required=True)
    parser.add_argument('--geo-root',type=Path)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--pipeline-runs',action='store_true',help='Summarize unified development runs and plot source coverage')
    mode.add_argument('--recheck-audits',action='store_true',help='Verify archived predictions/analyses and rerun only the audit; no clock or model calls')
    mode.add_argument('--enrich-sources',action='store_true',help='Retrieve public study Methods for archived development inputs; no model/clock calls')
    mode.add_argument('--check-measurement-profiles',action='store_true',help='Recompute blood clocks from a source GEO series matrix; no model calls or audit error-rate estimate')
    mode.add_argument('--verify-families',action='store_true',help='Verify family source identities and rebuild connected components independently; no study outcomes')
    mode.add_argument('--recover-dnam-input',action='store_true',help='Reconstruct an empty input from a source-linked supplementary beta table; no model calls')
    mode.add_argument('--source-design-plan',action='store_true',help='Run one frozen, budgeted source-only development review; --evidence supplies the plan')
    mode.add_argument('--source-design-batch',action='store_true',help='Run frozen source-design development plans under an aggregate budget; --evidence supplies the batch manifest')
    mode.add_argument('--export-audit-trace',action='store_true',help='Export one completed version3 batch child for the website; no model calls')
    mode.add_argument('--rna-normalization-sensitivity',action='store_true',help='Rerun exposed RNA counts under legacy and RLE sensitivity transforms; --evidence supplies the pinned plan')
    mode.add_argument('--rna-measurement-profile',action='store_true',help='Execute the bounded source-backed mouse RNA score profile on a frozen development plan; no model API calls')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--receipt-manifest',type=Path,help='Prior source_sha256 manifest for legacy archived stages without receipt hashes')
    parser.add_argument('--source-metadata',type=Path,help='Source GEO metadata CSV for measurement-profile verification')
    parser.add_argument('--extension-receipts',type=Path,help='Archived GEO date/source receipts to screen against the family map')
    parser.add_argument('--study',help='Explicit GEO study identity for input recovery')
    parser.add_argument('--failed-matrix',type=Path,help='Preserved empty molecular input for source recovery')
    parser.add_argument('--series-source',type=Path,help='GEO SOFT or Series Matrix declaring the supplement')
    parser.add_argument('--sample-source',type=Path,help='Primary GEO sample SOFT when stored separately from the series declaration')
    parser.add_argument('--prepare-recovery',action='store_true')
    parser.add_argument('--prepare-design',action='store_true',help='Freeze source-design request and verify limits without a model call')
    args=parser.parse_args()
    try:
        if args.rna_measurement_profile:
            result=check_rna_measurement_profile(args.evidence,args.output)
            print(json.dumps({'study_id': result['study_id'], 'rows': [{k: row.get(k) for k in
                ('clock','status','required_features','failures','unresolved')} for row in result['rows']],
                'scientific_status': result['scientific_status'], 'model_calls': 0}))
        elif args.rna_normalization_sensitivity:
            print(json.dumps(compare_rna_normalization(args.evidence,args.output)))
        elif args.export_audit_trace:
            print(json.dumps(export_source_audit_trace(args.evidence,args.output)))
        elif args.source_design_batch:
            print(json.dumps(run_source_design_batch(args.evidence,args.output,prepare_only=args.prepare_design)))
        elif args.source_design_plan:
            result=run_source_design_plan(args.evidence,args.output,prepare_only=args.prepare_design)
            print(json.dumps({key:value for key,value in result.items() if key!='design'}))
        elif args.recover_dnam_input:
            if not all((args.study,args.failed_matrix,args.source_metadata,args.series_source)):
                raise ValueError('Recovery requires --study, --failed-matrix, --source-metadata and --series-source')
            result=recover_dnam_input(args.study,args.failed_matrix,args.evidence,args.source_metadata,
                                      args.series_source,args.output,prepare_only=args.prepare_recovery,
                                      sample_source=args.sample_source)
            print(json.dumps(result))
        elif args.verify_families:
            result=verify_study_families(args.evidence,args.output,args.extension_receipts)
            print(json.dumps({key:value for key,value in result.items() if key!='extension_receipt_sha256'}))
        elif args.check_measurement_profiles:
            if args.source_metadata is None:raise ValueError('--source-metadata is required for measurement profiles')
            result=check_measurement_profiles(args.evidence,args.source_metadata,args.output)
            print(json.dumps({'study':result['study'],'rows':result['rows'],'output':str(args.output),'model_calls':0,'clock_execution_performed':True}))
        elif args.enrich_sources:
            result=enrich_pipeline_sources(args.evidence,args.output,args.receipt_manifest)
            print(json.dumps({'study_count':result['study_count'],'output':str(args.output),'model_calls':0,'clock_execution_performed':False}))
        elif args.recheck_audits:
            result=recheck_pipeline_audits(args.evidence,args.output)
            print(json.dumps({'study_count':result['study_count'],'output':str(args.output),'model_calls':0,'clock_execution_performed':False}))
        elif args.pipeline_runs:build_pipeline_report(args.evidence,args.output)
        else:
            if args.geo_root is None:raise ValueError('--geo-root is required for producer-census reporting')
            build(args.evidence,args.geo_root,args.output)
    except Exception as error:raise SystemExit(f'Validation report stopped: {type(error).__name__}: {error}')
