#!/usr/bin/env python3
"""ClockBase verifier: one command for the whole audit, from the saved agent outputs to the ranked atlas.

  python scripts/clockbase_verifier.py plan     --root ROOT --final NAME
  python scripts/clockbase_verifier.py run      --root ROOT --final NAME [--models cache-only|live] [--force STAGE ...]
  python scripts/clockbase_verifier.py verify   --run ROOT/NAME --reference ROOT/REFERENCE
  python scripts/clockbase_verifier.py selftest --root SANDBOX --reference ROOT/REFERENCE [--per-modality 40]

Stages run in dependency order. Each stage writes one directory and is skipped when its completion file exists,
so a run resumes where it stopped and a finished root is re-verified without recomputation. Model stages send
requests only with --models live (neolink and openai are aliases); with --models cache-only they answer from their ledgers and any
request that is not cached is reported, never sent. Inputs: the three atlas correction tables (CORRECTIONS, or
CLOCKBASE_CORRECTIONS), GEO_data (per-sample predictions, metadata, DNAm matrices) and the saved agent traces.
Locations come from audit_system.paths; each model command uses its saved run's pair and provider.

  upstream    membership recovery, trace replay, replay links (per modality); DNAm input-integrity screen;
              GEO source inventories and stream repairs of empty and detached DNAm studies
  pass 1      merge, verdicts and repairs without model decisions (selects rows for the model stages)
  models 1    design reconstruction of unrecovered rows; orientation adjudication of label-dependent rows
  pass 2      chain with design and adjudication (atlas used to select ranking candidates)
  models 2    adjudication of rows still pending, unconfirmed rows and a validation sample; contrast
              classification; design audit of ranking candidates (and of candidates added by a later rule)
  pass 3      the chain with the design-respecting estimand (individual means, paired arms, age as a covariate,
              overflow values dropped); its atlas selects the candidates of the next two model stages
  models 3    design audit of candidates added by the design estimand; two-reader audit of the individuals behind
              the samples of ranking candidates whose series has no GEO individual field
  final       chain: merge, verdicts, repairs (design estimand, audited individuals), published-output status, error
              census, atlas, ranking, export
  audits      every saved model decision re-derived from the stored raw model responses
  analyses    blood cross-clock structure and GrimAge input check; the sampling frame
              of the orientation-validation sample
  tests       unit tests of the verifier modules

`verify` compares two chain directories row by row (final class and recomputed statistics), the error census,
published-output fates, the ranking and the per-study export, and writes verify.json. Timestamps, module hashes,
run-directory paths and model-call receipts are provenance: differences in them are counted, not failed.
`selftest` plants known defects into copies of verified atlas rows (a disclosed simulation on real records),
runs the deterministic stages on them in a sandbox and reports, per defect type, how many were caught.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from audit_system import paths

PYTHON = os.environ.get('CLOCKBASE_PYTHON') or sys.executable  # the interpreter of this run, for every stage
PY = [PYTHON, '-B', '-m']
MODALITIES = ('dnam', 'hs_rna', 'mm_rna')
EMPTY_LIST = paths.EMPTY_MATRIX_LIST
MODEL_PAIRS = {'tokenswitch': 'gp-6-as,cl-op-5', 'openai': 'openai:gpt-6-astra,openai:gpt-6-sol'}
VERIFIER_TESTS = ['test_row_verdicts_2026_09_27', 'test_row_repairs_2026_09_28', 'test_error_census_2026_09_29',
                  'test_comparison_ranking_2026_09_30', 'test_design_audit_2026_09_30', 'test_input_integrity_2026_09_30',
                  'test_blood_clock_structure_2026_09_30', 'test_membership_recovery_2026_09_27',
                  'test_design_reconstruction_2026_09_28', 'test_intervention_classification_2026_09_29',
                  'test_dnam_stream_repair_2026_09_26', 'test_clockbase_verifier_2026_10_04',
                  'test_v3_public_export_2026_09_30', 'test_model_decision_audit_2026_10_04',
                  'test_prompts_release_2026_10_04', 'test_validation_frame_2026_10_05',
                  'test_design_estimand_2026_10_10', 'test_individuals_audit_2026_10_10', 'test_final_verification_rules_2026_11_13']


class Stage:
    def __init__(self, name, group, done, commands, note=''):
        self.name, self.group, self.done, self.note = name, group, done, note
        self.commands = [(c[0], c[1], c[2] if len(c) > 2 else None) for c in commands]
        self.model = any(provider for _, _, provider in self.commands)


def model_command(command, provider):
    return command, {'AUDIT_MODELS': MODEL_PAIRS[provider]}, provider


def chain(root, name, stream, detached, extra=()):
    env = {'STREAM_ROOT': str(stream), 'STREAM_DETACHED_ROOT': str(detached), 'CLOCKBASE_PYTHON': PYTHON}
    return [(['bash', str(REPO / 'scripts/verifier_v3_final_chain_2026_09_30.sh'), str(root), name, *extra], env)]


def mini_chain(root, out, design=False, adjudication=False):
    """merge, verdicts and repairs (no census, ranking or export): passes 1 and 2."""
    commands = []
    for m in MODALITIES:
        merge = PY + ['audit_system.membership_merge_2026_09_27', '--modality', m, '--recovery', f'{root}/membership-{m}-001',
                      '--link', f'{root}/link-{m}-001', '--replay', f'{root}/replay-{m}-001', '--out', f'{out}/merged-{m}']
        if design:
            merge += ['--design', f'{root}/design-{m}-001']
        verdicts = PY + ['audit_system.row_verdicts_2026_09_27', '--modality', m, '--membership', f'{out}/merged-{m}',
                         '--out', f'{out}/verdicts-{m}']
        if m == 'dnam':
            verdicts += ['--integrity', f'{root}/input-integrity-dnam-001']
        repairs = PY + ['audit_system.row_repairs_2026_09_28', '--modality', m, '--membership', f'{out}/merged-{m}',
                        '--verdicts', f'{out}/verdicts-{m}', '--replay', f'{root}/replay-{m}-001', '--out', f'{out}/repairs-{m}']
        if adjudication:
            repairs += ['--adjudication', f'{root}/adjudication-{m}-001']
        commands += [(merge, {}), (verdicts, {}), (repairs, {})]
    return commands


def stages(root, final, stream, detached):
    r = Path(root)
    out = []
    for m in MODALITIES:
        out.append(Stage(f'membership-{m}', 'upstream', r / f'membership-{m}-001/manifest.json',
                         [(PY + ['audit_system.membership_recovery_2026_09_27', 'run', '--modality', m, '--out',
                                 f'{r}/membership-{m}-001', '--workers', '16', '--resume'], {})]))
        out.append(Stage(f'replay-{m}', 'upstream', r / f'replay-{m}-001/analysis_summary.json',
                         [(PY + ['audit_system.trace_replay_2026_09_27', step, '--modality', m, '--out', f'{r}/replay-{m}-001'], {})
                          for step in ('extract', 'split-inputs', 'run', 'analyze')], note='sandboxed with unshare -rn'))
        out.append(Stage(f'link-{m}', 'upstream', r / f'link-{m}-001/manifest.json',
                         [(PY + ['audit_system.replay_link_2026_09_27', '--modality', m, '--replay', f'{r}/replay-{m}-001',
                                 '--out', f'{r}/link-{m}-001'], {})]))
    out.append(Stage('integrity-dnam', 'upstream', r / 'input-integrity-dnam-001/samples.csv',
                     [(PY + ['audit_system.input_integrity_2026_09_30', 'dnam', '--out', f'{r}/input-integrity-dnam-001',
                             '--workers', '12'], {})]))
    inventory = [PYTHON, '-B']
    out.append(Stage('stream-repair-empty', 'upstream', Path(stream) / 'summary.json', [
        (inventory + [str(REPO / 'src/dnam_source_inventory_2026_09_30/geo_inventory.py')],
         {'INVENTORY_OUT': f'{r}/empty-repair-inventory', 'INVENTORY_SOURCE': str(EMPTY_LIST), 'INVENTORY_NAME': 'empty_dnam_geo_repair'}),
        (inventory + [str(REPO / 'src/dnam_source_inventory_2026_09_30/classify.py')],
         {'INVENTORY_OUT': f'{r}/empty-repair-inventory', 'INVENTORY_NAME': 'empty_dnam_geo_repair'}),
        (PY + ['src.dnam_stream_repair_2026_09_26', 'run', '--inventory',
               f'{r}/empty-repair-inventory/empty_dnam_geo_repair_inventory_2026-09-25.csv', '--root', str(stream), '--workers', '2'], {})],
        note='downloads GEO supplementary files'))
    out.append(Stage('stream-repair-detached', 'upstream', Path(detached) / 'summary.json', [
        (PY + ['audit_system.input_integrity_2026_09_30', 'detached-list', '--screen', f'{r}/input-integrity-dnam-001',
               '--out', f'{r}/detached-repair-inventory/detached_studies.jsonl'], {}),
        (inventory + [str(REPO / 'src/dnam_source_inventory_2026_09_30/geo_inventory.py')],
         {'INVENTORY_OUT': f'{r}/detached-repair-inventory',
          'INVENTORY_SOURCE': f'{r}/detached-repair-inventory/detached_studies.jsonl', 'INVENTORY_NAME': 'detached_dnam_geo_repair'}),
        (inventory + [str(REPO / 'src/dnam_source_inventory_2026_09_30/classify.py')],
         {'INVENTORY_OUT': f'{r}/detached-repair-inventory', 'INVENTORY_NAME': 'detached_dnam_geo_repair'}),
        (PY + ['src.dnam_stream_repair_2026_09_26', 'run', '--inventory',
               f'{r}/detached-repair-inventory/detached_dnam_geo_repair_inventory.csv', '--root', str(detached), '--workers', '2'], {})],
        note='downloads GEO supplementary files'))
    out.append(Stage('pass-1', 'pass 1', r / 'pass-1/repairs-mm_rna/manifest.json', mini_chain(r, r / 'pass-1')))
    for m in MODALITIES:
        out.append(Stage(f'design-{m}', 'models 1', r / f'design-{m}-001/summary.json',
                         [model_command(PY + ['audit_system.design_reconstruction_2026_09_28', '--modality', m,
                                              '--repairs', f'{r}/pass-1/repairs-{m}', '--out', f'{r}/design-{m}-001'], 'tokenswitch')]))
        out.append(Stage(f'adjudication-decisive-{m}', 'models 1', r / f'adjudication-{m}-001', [
            model_command(PY + ['audit_system.label_adjudication_2026_09_27', '--modality', m,
                                '--membership', f'{r}/pass-1/merged-{m}', '--out', f'{r}/adjudication-{m}-001',
                                '--select', 'decisive', '--verdicts', f'{r}/pass-1/verdicts-{m}'], 'tokenswitch')]))
    out.append(Stage('pass-2', 'pass 2', r / 'pass-2/repairs-mm_rna/manifest.json',
                     mini_chain(r, r / 'pass-2', design=True, adjudication=True)))
    for m in MODALITIES:
        out.append(Stage(f'adjudication-followup-{m}', 'models 2', r / f'adjudication-{m}-001/.followup-done', [
            model_command(PY + ['audit_system.label_adjudication_2026_09_27', '--modality', m,
                                '--membership', f'{r}/pass-2/merged-{m}', '--out', f'{r}/adjudication-{m}-001',
                                '--select', select, '--verdicts', f'{r}/pass-2/repairs-{m}']
                          + (['--limit', '300'] if select == 'supported_sample' else []),
                          'openai' if select == 'supported_sample' else 'tokenswitch')
            for select in ('pending', 'unconfirmed', 'supported_sample')]
            + [(['touch', f'{r}/adjudication-{m}-001/.followup-done'], {})]))
    rep = [x for m in MODALITIES for x in ('--repairs', f'{m}={r}/pass-2/repairs-{m}')]
    out.append(Stage('pass-2-atlas', 'pass 2', r / 'pass-2/atlas-v3/atlas_v3_combined.csv',
                     [(PY + ['audit_system.atlas_rebuild_2026_09_29', *rep, '--out', f'{r}/pass-2/atlas-v3'], {})]))
    out.append(Stage('classification', 'models 2', r / 'classification-001/summary.json',
                     [model_command(PY + ['audit_system.intervention_classification_2026_09_29', '--atlas',
                                          f'{r}/pass-2/atlas-v3/atlas_v3_combined.csv', '--out', f'{r}/classification-001',
                                          *rep, '--batch', '8', '--max-tokens', '10000'], 'tokenswitch')]))
    out.append(Stage('design-audit', 'models 2', r / 'design-audit-001/summary.json',
                     [model_command(PY + ['audit_system.design_audit_2026_09_30', '--atlas', f'{r}/pass-2/atlas-v3/atlas_v3_combined.csv',
                                          '--classification', f'{r}/classification-001', *rep,
                                          '--out', f'{r}/design-audit-001'], 'openai')]))
    out.append(Stage('design-audit-added', 'models 2', r / 'design-audit-002/summary.json',
                     [model_command(PY + ['audit_system.design_audit_2026_09_30', '--atlas', f'{r}/pass-2/atlas-v3/atlas_v3_combined.csv',
                                          '--classification', f'{r}/classification-001', *rep, '--out', f'{r}/design-audit-002',
                                          '--exclude', f'{r}/design-audit-001'], 'openai')],
                     note='candidates under the contrast-type rule that design-audit-001 did not cover'))
    # Pass 3: the chain under the design-respecting estimand, before the individuals of ranking candidates are read.
    pre = f'pre{final}'
    rep3 = [x for m in MODALITIES for x in ('--repairs', f'{m}={r}/{pre}/repairs-{m}')]
    atlas3 = f'{r}/{pre}/atlas-v3/atlas_v3_combined.csv'
    out.append(Stage('pass-3', 'pass 3', r / pre / 'public-export/manifest.json', chain(r, pre, stream, detached)))
    out.append(Stage('design-audit-design-estimand', 'models 3', r / 'design-audit-003/summary.json',
                     [model_command(PY + ['audit_system.design_audit_2026_09_30', '--atlas', atlas3,
                                          '--classification', f'{r}/classification-001', *rep3, '--out', f'{r}/design-audit-003',
                                          '--exclude', f'{r}/design-audit-001', '--exclude', f'{r}/design-audit-002'], 'openai')],
                     note='candidates that reach the ranking only under the design estimand'))
    out.append(Stage('individuals-audit', 'models 3', r / 'individuals-audit-001/summary.json',
                     [model_command(PY + ['audit_system.individuals_audit_2026_10_10', '--atlas', atlas3,
                                          '--classification', f'{r}/classification-001',
                                          *[x for n in ('001', '002', '003') for x in ('--design-audit', f'{r}/design-audit-{n}')],
                                          *rep3, '--out', f'{r}/individuals-audit-001'], 'openai')]))
    out.append(Stage('final-chain', 'final', r / final / 'public-export/manifest.json', chain(r, final, stream, detached)))
    # Audit of the model layer: every saved model decision re-derived from the stored raw responses (no requests).
    out.append(Stage('decision-audit', 'audits', r / final / 'decision-audit/audit.json',
                     [(PY + ['audit_system.model_decision_audit_2026_10_04', '--verifier', str(r),
                             '--out', f'{r}/{final}/decision-audit'], {})]))
    blood = r / 'r6-6-blood'
    out.append(Stage('blood-structure-cohort', 'analyses', blood / 'cohort-002/summary.json',
                     [(PY + ['audit_system.blood_clock_structure_2026_09_30', '--integrity', f'{r}/input-integrity-dnam-001',
                             '--out', f'{blood}/cohort-002'], {})]))
    out.append(Stage('grimage-input', 'analyses', blood / 'grimage-age-input-002/summary.json',
                     [(PY + ['audit_system.grimage_age_input_2026_09_30', '--cohort', f'{blood}/cohort-002/cohort.csv',
                             '--out', f'{blood}/grimage-age-input-002'], {})]))
    out.append(Stage('blood-structure', 'analyses', blood / 'structure-002/summary.json',
                     [(PY + ['audit_system.blood_clock_structure_2026_09_30', '--integrity', f'{r}/input-integrity-dnam-001',
                             '--grimage', f'{blood}/grimage-age-input-002', '--out', f'{blood}/structure-002'], {})]))
    # The frame of the orientation-validation sample, reconstructed with the original selection and redrawn (no models).
    out.append(Stage('validation-frame', 'analyses', r / 'validation-frame-001/validation_summary.json',
                     [(PY + ['audit_system.validation_frame_2026_10_05', '--root', str(r),
                             '--out', f'{r}/validation-frame-001'], {})]))
    out.append(Stage('tests', 'tests', None, [(PY + ['unittest', *[f'audit_system.tests.{t}' for t in VERIFIER_TESTS]],
                                               {'TMPDIR': os.environ.get('TMPDIR', tempfile.gettempdir())})]))
    return out


def model_env(mode, provider):
    if mode == 'cache-only':
        return {'AUDIT_CACHE_ONLY': '1'}, []
    return {'AUDIT_CACHE_ONLY': '0'}, [token.replace('{provider}', provider) for token in paths.MODEL_LAUNCHER]


def run_stages(args):
    root = Path(args.root)
    if not args.plan_only:
        root.mkdir(parents=True, exist_ok=True)
    plan = stages(root, args.final, args.stream, args.stream_detached)
    log_dir = root / 'pipeline-runs' / f"{time.strftime('%Y%m%dT%H%M%S')}-{os.getpid()}"
    if not args.plan_only:
        log_dir.mkdir(parents=True)
    manifest = {'command': sys.argv, 'models': args.models, 'started': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'stages': []}
    # Passes 1 and 2 only select rows for the model stages; skip them when every stage that reads them is done
    # or itself skipped (checked downstream first, so pass-2 is skipped when pass-2-atlas is).
    readers = {'pass-2-atlas': [s for s in plan if s.name in ('classification', 'design-audit', 'design-audit-added')],
               'pass-2': [s for s in plan if s.name.startswith('adjudication-followup-')
                          or s.name in ('pass-2-atlas', 'classification', 'design-audit', 'design-audit-added')],
               'pass-1': [s for s in plan if s.name in [f'design-{m}' for m in MODALITIES]
                          or s.name.startswith('adjudication-decisive-')],
               'pass-3': [s for s in plan if s.name in ('design-audit-design-estimand', 'individuals-audit')]}
    not_needed = set()
    for name, reading in readers.items():
        if not args.force and all(r.name in not_needed or (r.done and Path(r.done).exists()) for r in reading):
            not_needed.add(name)
    for stage in plan:
        entry = {'stage': stage.name, 'group': stage.group, 'done_file': str(stage.done) if stage.done else None}
        if stage.model:
            entry['model_commands'] = [{'argv': command, 'provider': provider, 'models': env['AUDIT_MODELS'].split(',')}
                                       for command, env, provider in stage.commands if provider]
        if stage.name in not_needed:
            entry['status'] = 'NOT_NEEDED'
            manifest['stages'].append(entry)
            print(f'[skip]  {stage.name} (every stage that reads it is done)')
            continue
        forced = stage.name in args.force or (args.force and 'all' in args.force)
        done = stage.done is not None and Path(stage.done).exists()
        if done and not forced:
            entry['status'] = 'REUSED'
            manifest['stages'].append(entry)
            print(f'[reuse] {stage.name}')
            continue
        if args.plan_only:
            entry['status'] = 'WOULD_RUN'
            manifest['stages'].append(entry)
            print(f'[run]   {stage.name}  ({stage.group}{", model" if stage.model else ""}{", " + stage.note if stage.note else ""})')
            continue
        started = time.time()
        entry['commands'] = []
        for command, env, provider in stage.commands:
            menv, wrapper = model_env(args.models, provider) if provider else ({}, [])
            full = wrapper + command
            environment = {**os.environ, **env, **menv}
            with open(log_dir / f'{stage.name}.log', 'a') as handle:
                handle.write(f'$ {" ".join(full)}\n')
                handle.flush()
                result = subprocess.run(full, cwd=REPO, env=environment, stdout=handle, stderr=subprocess.STDOUT)
            entry['commands'].append({'argv': full, 'exit': result.returncode})
            if result.returncode != 0:
                entry.update(status='FAILED', seconds=round(time.time() - started, 1))
                manifest['stages'].append(entry)
                (log_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2))
                raise SystemExit(f'stage {stage.name} failed; see {log_dir / (stage.name + ".log")}')
        entry.update(status='RAN', seconds=round(time.time() - started, 1))
        manifest['stages'].append(entry)
        (log_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        print(f'[done]  {stage.name}  {entry["seconds"]}s')
    if args.plan_only:
        return
    manifest['finished'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    (log_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print('manifest:', log_dir / 'manifest.json')


def load_repairs(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            r = json.loads(line)
            rec = r.get('recomputed') or {}
            fields = ('n_case', 'n_control', 'mean_case', 'mean_control', 'mean_diff',
                      'cohen_d', 't', 'p', 'fdr_family', 'fdr_family_size', 'undefined_reason')
            out[r['row_ordinal']] = (r['final_class'], *(rec.get(k) for k in fields))
    return out


def same(a, b, tolerance=1e-9):
    if isinstance(a, float) or isinstance(b, float):
        if a is None or b is None:
            return a is b
        return abs(a - b) <= tolerance * max(1.0, abs(a), abs(b))
    return a == b


PROVENANCE = ('created', 'module_sha256', 'repairs_sha256', 'dir', 'source', 'path', 'receipts')  # hashes stay compared


def content(obj):
    """Scientific content of a record: drops timestamps, module hashes, run-directory paths and model-call receipts."""
    if isinstance(obj, dict):
        return {k: content(v) for k, v in obj.items() if k not in PROVENANCE}
    if isinstance(obj, list):
        return [content(v) for v in obj]
    return obj


def verify(args):
    new, ref = Path(args.run), Path(args.reference)
    report = {'run': str(new), 'reference': str(ref), 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'checks': {},
              'compared_as_provenance': list(PROVENANCE)}
    for m in MODALITIES:
        a, b = load_repairs(new / f'repairs-{m}/repairs.jsonl'), load_repairs(ref / f'repairs-{m}/repairs.jsonl')
        diffs = [o for o in sorted(set(a) | set(b)) if o not in a or o not in b or not all(same(x, y) for x, y in zip(a[o], b[o]))]
        report['checks'][f'rows-{m}'] = {'rows': len(b), 'differing': len(diffs), 'examples': [
            {'row': o, 'run': a.get(o), 'reference': b.get(o)} for o in diffs[:5]]}
    for name in ('error-census/census.json', 'published-status/summary.json'):
        a, b = json.loads((new / name).read_text()), json.loads((ref / name).read_text())
        report['checks'][name] = {'identical': content(a) == content(b), 'byte_identical_after_timestamps': a == b}
    import csv
    def ranked(path):
        with open(path) as handle:
            return [(r['gse_id'], r['label'], r['clock'], round(float(r['hedges_g']), 9), round(float(r['rank_bound']), 9))
                    for r in csv.DictReader(handle)]
    ra, rb = ranked(new / 'comparison-ranking/ranked_comparisons_all.csv'), ranked(ref / 'comparison-ranking/ranked_comparisons_all.csv')
    report['checks']['ranking'] = {'ranked_run': len(ra), 'ranked_reference': len(rb), 'identical_order_and_values': ra == rb}
    ea, eb = json.loads((new / 'public-export/manifest.json').read_text()), json.loads((ref / 'public-export/manifest.json').read_text())
    study_files = lambda e: {k: v for k, v in e['files'].items() if k.startswith('studies/')}
    fa, fb = study_files(ea), study_files(eb)
    hash_differs = sorted(k for k in set(fa) | set(fb) if fa.get(k) != fb.get(k))
    content_differs, provenance_only = [], []
    for k in hash_differs:
        if k not in fa or k not in fb:
            content_differs.append(k)
            continue
        x, y = json.loads((new / 'public-export' / k).read_text()), json.loads((ref / 'public-export' / k).read_text())
        (provenance_only if content(x) == content(y) else content_differs).append(k)
    report['checks']['export'] = {'study_files_run': len(fa), 'study_files_reference': len(fb),
                                  'byte_identical': len(fb) - len(hash_differs), 'provenance_only': len(provenance_only),
                                  'provenance_only_examples': provenance_only[:20], 'differing_files': content_differs[:50],
                                  'differing_count': len(content_differs)}
    ok = (all(c['differing'] == 0 for k, c in report['checks'].items() if k.startswith('rows-'))
          and all(report['checks'][n]['identical'] for n in ('error-census/census.json', 'published-status/summary.json'))
          and report['checks']['ranking']['identical_order_and_values'] and not content_differs)
    report['reproduced'] = ok
    (new / 'verify.json').write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(report, indent=1, default=str)[:4000])
    raise SystemExit(0 if ok else 1)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('plan', 'run'):
        p = sub.add_parser(name)
        p.add_argument('--root', type=Path, default=paths.VERIFIER_ROOT)
        p.add_argument('--final', required=True, help='name of the final chain directory under ROOT')
        p.add_argument('--models', choices=['cache-only', 'live', 'neolink', 'openai'], default='cache-only')
        p.add_argument('--stream', type=Path, default=paths.STREAM_EMPTY)
        p.add_argument('--stream-detached', type=Path, default=paths.STREAM_DETACHED)
        p.add_argument('--force', nargs='*', default=[], help="stage names to rerun ('all' for every stage)")
    v = sub.add_parser('verify')
    v.add_argument('--run', type=Path, required=True)
    v.add_argument('--reference', type=Path, required=True)
    s = sub.add_parser('selftest')
    s.add_argument('--root', type=Path, required=True)
    s.add_argument('--reference', type=Path, required=True)
    s.add_argument('--per-modality', type=int, default=40)
    args = parser.parse_args(argv)
    if args.command in ('plan', 'run'):
        args.plan_only = args.command == 'plan'
        run_stages(args)
    elif args.command == 'verify':
        verify(args)
    else:
        sys.path.insert(0, str(REPO))
        from audit_system.verifier_selftest_2026_10_01 import run as selftest
        selftest(args.root, args.reference, args.per_modality)


if __name__ == '__main__':
    main()
