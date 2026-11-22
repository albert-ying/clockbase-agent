"""Recover the orientation-validation frame and redraw its recorded sample, without model calls.

  python -m audit_system.validation_frame_2026_10_05 --out NEW_DIRECTORY

The original summaries saved the membership hash, but not a repair hash or directory. Select the
latest surviving pre-sample chain with that membership hash, independently of sample agreement.
The September 30d chain survives even though pass-2 was partly overwritten. Export every row,
the redrawn sample, differences and input provenance; never modify a verifier run.
"""
import argparse
import ast
import hashlib
import inspect
import json
import math
import re
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from audit_system import label_adjudication_2026_09_27 as labels
from audit_system import paths
from audit_system import row_verdicts_2026_09_27 as verdicts
from audit_system.adjudication_files import decision_files
from audit_system.membership_recovery_2026_09_27 import MODALITIES, sha256

SCHEMA = 'clockbase-validation-frame-v1'
EXPECTED = {'dnam': 6188, 'hs_rna': 13425, 'mm_rna': 15107}
CUTOFF = '20261001T162200'
IDENTITY = ['row_ordinal', 'gse_id', 'comparison_name', 'condition_name', 'clock']
# ASTs of the selection and heuristic in a19b79204, the code recorded by the sample summaries.
SELECTION_SHA256 = '63c29f71357b7b85d7ae4705f3393c77e373d485d3e7a29b027427cca67f39d5'
HEURISTIC_SHA256 = 'f0d5df1ca62b230efc036c39328ffbcbee12af8a8c6405bf0f98459719d1828b'
HEURISTIC_BUNDLE_SHA256 = '6ac53a4e4a5d9013fe4ef0ccf29096f4f7ebf9c01a9278165ba4df3d885e4e76'


def read_jsonl(path):
    with open(path) as handle:
        return [json.loads(line) for line in handle if line.strip()]


def provenance(path):
    stat = Path(path).stat()
    return {'path': paths.portable(path), 'sha256': sha256(path), 'bytes': stat.st_size,
            'modified_utc': datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat()}


def before_sample(folder, cutoff=CUTOFF):
    """Use filename timestamps; an undated decision file cannot establish historical availability."""
    selected = []
    for path in decision_files(folder):
        stamp = re.search(r'-(\d{8}T\d{6})\.jsonl$', path.name)
        if not stamp:
            raise ValueError(f'Cannot date adjudication decisions: {paths.portable(path)}')
        if stamp.group(1) < cutoff:
            selected.append(path)
    return selected


def selection_function(source):
    """Execute only label_adjudication's original selection prefix, ending before the ledger/API stage."""
    function = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == 'run')
    stop = next(i for i, n in enumerate(function.body) if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'ledger' for t in n.targets))
    function.body = function.body[:stop]
    digest = hashlib.sha256(ast.dump(function, include_attributes=False).encode()).hexdigest()
    if digest != SELECTION_SHA256:
        raise ValueError('label_adjudication selection differs from the recorded October 1 implementation')
    heuristic = ast.parse(inspect.getsource(labels.label_orientation)).body[0]
    if hashlib.sha256(ast.dump(heuristic, include_attributes=False).encode()).hexdigest() != HEURISTIC_SHA256:
        raise ValueError('label_orientation differs from the recorded October 1 implementation')
    nodes = [n for n in ast.parse(inspect.getsource(verdicts)).body
             if (isinstance(n, ast.FunctionDef) and n.name in ('tokens', 'label_orientation'))
             or (isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'STOP' for t in n.targets))]
    bundle = ast.Module(body=nodes, type_ignores=[])
    if hashlib.sha256(ast.dump(bundle, include_attributes=False).encode()).hexdigest() != HEURISTIC_BUNDLE_SHA256:
        raise ValueError('label_orientation tokenization differs from the recorded October 1 implementation')
    function.body.append(ast.Return(value=ast.Tuple(elts=[ast.Name(id='pool', ctx=ast.Load()),
                                                         ast.Name(id='rows', ctx=ast.Load())], ctx=ast.Load())))
    module = ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[]))
    namespace = dict(vars(labels))
    exec(compile(module, '<label_adjudication selection only>', 'exec'), namespace)
    return namespace['run'], digest


def redraw(membership, repairs, historical_decisions, limit=300):
    select, digest = selection_function(inspect.getsource(labels))
    frame, sample = select('dnam', membership, historical_decisions, select='supported_sample', limit=limit,
                           verdicts_dir=repairs)
    return frame[IDENTITY], sample[IDENTITY], digest


def compare_sample(sample, recorded):
    drawn = sample[IDENTITY].to_dict('records')
    expected = [{k: r.get(k) for k in IDENTITY} for r in recorded]
    keys = lambda records: [tuple(r.get(k) for k in IDENTITY) for r in records]
    a, b = keys(drawn), keys(expected)
    ca, cb = Counter(a), Counter(b)
    return {'ordered_equal': a == b, 'rows_equal': ca == cb,
            'redrawn_only': [dict(zip(IDENTITY, k)) for k in (ca - cb).elements()],
            'recorded_only': [dict(zip(IDENTITY, k)) for k in (cb - ca).elements()]}


def saved_chains(root, modality, membership_hash):
    """Newest pre-validation snapshot first; do not choose a snapshot by its redrawn sample."""
    candidates, rejected = [], []
    boundary = datetime.strptime(CUTOFF, '%Y%m%dT%H%M%S')
    for folder in sorted(Path(root).iterdir()):
        manifest_path = folder / f'repairs-{modality}/manifest.json'
        if not manifest_path.exists():
            continue
        try:
            manifest = json.loads(manifest_path.read_text())
            created = datetime.fromisoformat(manifest['created']).replace(tzinfo=None)
            if created >= boundary or manifest.get('membership_sha256') != membership_hash:
                continue
            membership = folder / f'merged-{modality}/rows.jsonl'
            repairs = folder / f'repairs-{modality}/repairs.jsonl'
            verdicts = folder / f'verdicts-{modality}/row_verdicts.jsonl'
            for path, key in ((membership, 'membership_sha256'), (repairs, 'repairs_sha256'),
                              (verdicts, 'verdicts_sha256')):
                if sha256(path) != manifest[key]:
                    raise ValueError(f'{key} does not match {paths.portable(path)}')
            candidates.append((created, folder, manifest))
        except (OSError, ValueError, KeyError) as error:
            rejected.append({'directory': paths.portable(folder), 'error': str(error)})
    candidates.sort(key=lambda r: (r[0], str(r[1])), reverse=True)
    return candidates, rejected


def weighted_rate(results):
    """Facts loader: sum(N_m/N * k_m/n_m), with its stratified normal interval (no FPC)."""
    total = sum(r['frame_size'] for r in results.values())
    if not total or any(not r['sample_size'] for r in results.values()):
        return None
    weights = {m: r['frame_size'] / total for m, r in results.items()}
    rates = {m: r['decisions'].get('LABEL_INVERTED_CONFIRMED', 0) / r['sample_size'] for m, r in results.items()}
    rate = sum(weights[m] * rates[m] for m in results)
    variance = sum(weights[m] ** 2 * rates[m] * (1 - rates[m]) / (r['sample_size'] - 1)
                   for m, r in results.items() if r['sample_size'] > 1)
    half = 1.96 * math.sqrt(variance)
    return {'rate': rate, 'ci95': [max(0.0, rate - half), rate + half], 'weights': weights,
            'rates': rates, 'computation': 'sum(N_m/sum(N) * inverted_m/sample_m)',
            'variance': variance, 'variance_computation': 'sum((N_m/sum(N))**2 * p_m*(1-p_m)/(n_m-1))',
            'frame_total': total, 'expected_inversions': rate * total}


def run(out_dir, root=paths.VERIFIER_ROOT):
    out_dir, root = Path(out_dir), Path(root)
    out_dir.mkdir(parents=True, exist_ok=False)
    results = {}
    for modality in MODALITIES:
        adjudication = root / f'adjudication-{modality}-001'
        samples = list(adjudication.glob('decisions-supported_sample-20261001T*.jsonl'))
        if len(samples) != 1:
            raise ValueError(f'Expected one October 1 validation sample for {modality}, found {len(samples)}')
        sample_path = samples[0]
        summary_path = sample_path.with_name(sample_path.name.replace('decisions-', 'summary-')).with_suffix('.json')
        summary = json.loads(summary_path.read_text())
        recorded = read_jsonl(sample_path)
        candidates, rejected = saved_chains(root, modality, summary['membership_sha256'])
        if not candidates:
            raise ValueError(f'Cannot reconstruct {modality}: no intact pre-sample chain with recorded membership hash')
        _, chain, manifest = candidates[0]
        folder = out_dir / modality
        folder.mkdir()
        earlier = before_sample(adjudication)
        decision_dir = folder / 'pre-sample-decisions'
        decision_dir.mkdir()
        for path in earlier:
            (decision_dir / path.name).write_bytes(path.read_bytes())
        membership = chain / f'merged-{modality}'
        repairs = chain / f'repairs-{modality}'
        frame, sample, digest = redraw(membership, repairs, decision_dir, summary['rows'])
        comparison = compare_sample(sample, recorded)
        frame.to_json(folder / 'frame.jsonl', orient='records', lines=True)
        sample.to_json(folder / 'redrawn_sample.jsonl', orient='records', lines=True)
        counts = dict(Counter(r['decision'] for r in recorded))
        inputs = [sample_path, summary_path, membership / 'rows.jsonl', repairs / 'repairs.jsonl',
                  repairs / 'manifest.json', chain / f'verdicts-{modality}/row_verdicts.jsonl', *earlier]
        record = {'frame_size': len(frame), 'expected_frame_size': EXPECTED[modality],
                  'frame_size_difference': len(frame) - EXPECTED[modality], 'sample_size': len(recorded),
                  'redrawn_sample_size': len(sample), 'decisions': counts, 'sample_comparison': comparison,
                  'decisions_match_summary': counts == summary['decisions'],
                  'reconstruction_exact': comparison['ordered_equal'],
                  'state_directory': paths.portable(chain), 'state_created': manifest['created'],
                  'candidate_directories': [paths.portable(c[1]) for c in candidates], 'rejected_candidates': rejected,
                  'membership_sha256_matches_summary': True, 'selection_ast_sha256': digest,
                  'heuristic_ast_sha256': HEURISTIC_SHA256, 'heuristic_bundle_sha256': HEURISTIC_BUNDLE_SHA256,
                  'historical_code_revision': 'a19b79204', 'recorded_module_sha256': summary['module_sha256'],
                  'executed_module_sha256': sha256(labels.__file__), 'seed': 20260930,
                  'models': summary['models'], 'decision_cutoff_local': CUTOFF, 'inputs': [provenance(p) for p in inputs],
                  'limitation': 'The sample summary has no repair hash or launch-directory receipt. The selected '
                                'pre-sample snapshot has the recorded membership hash and its sample is compared '
                                'in full identity and order; this does not prove the original directory spelling.'}
        if not comparison['ordered_equal']:
            record['limitation'] += ' Exact frame recovery is unconfirmed; differences are exported without adjustment.'
        (folder / 'summary.json').write_text(json.dumps(record, indent=2) + '\n')
        results[modality] = record
        print(f'{modality}: frame={len(frame)} expected={EXPECTED[modality]} sample_equal={comparison["ordered_equal"]}', flush=True)
    weighted = weighted_rate(results)
    receipts = []
    for path in sorted((root / 'pipeline-runs').glob('*/manifest.json')):
        manifest = json.loads(path.read_text())
        receipts.append({'path': paths.portable(path), 'started': manifest.get('started'), 'sha256': sha256(path)})
    result = {'schema': SCHEMA, 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'modalities': results,
              'weighted_inversion': weighted, 'exact_samples': all(r['reconstruction_exact'] for r in results.values()),
              'directory_provenance': 'pass-2 pipeline plans start after the validation summaries. Recover from '
                                      'intact pre-validation chains rather than treating those later plans as launch receipts.',
              'pipeline_receipts': receipts, 'selection_source': provenance(labels.__file__),
              'heuristic_source': provenance(verdicts.__file__),
              'job_script': provenance(paths.REPO / 'scripts/verifier_v3_openai_jobs_2026_10_01.sh')}
    (out_dir / 'validation_summary.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--out', type=Path, required=True, help='new scratch directory (refuses to overwrite)')
    parser.add_argument('--root', type=Path, default=paths.VERIFIER_ROOT)
    args = parser.parse_args(argv)
    try:
        run(args.out, args.root)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f'Validation frame recovery failed: {error}\n')


if __name__ == '__main__':
    main()
