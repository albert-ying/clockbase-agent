"""Dual-model design audit of ranking candidates: specimen context, tumour material and arm confounds.

The comparison ranking needs three properties that the reproduction checks cannot establish:
  1. the case/control difference arose in living organisms and the samples were profiled without
     culture. GPT-4 per-sample in vivo labels call cultured primary cells in vivo (GSE128909: bone
     marrow neutrophils cultured for 4 h in GM-CSF were labelled in_vivo), and replicate cultures have
     so little variance that they dominate any standardized ranking;
  2. no arm is tumour material; the clocks were trained on non-malignant tissue;
  3. the arms differ only by the labelled condition (GSE74738 compares miscarriage placentas at 6-10
     weeks with control placentas up to 38 weeks, so the clock difference tracks gestational age).
GPT-6 Astra and Claude Opus 5 independently read the series title and summary and the literal GEO
records of the samples in each arm (title, source, characteristics, treatment, growth and extraction
protocols, description, platform), never P values or effect sizes. Each returns per comparison:
  orientation  CASE_IS_CONDITION | REVERSED | UNCLEAR  (is the case arm the labelled condition?)
  specimen     ORGANISM_DIRECT | EX_VIVO_CULTURE | CELL_LINE | STEM_CELL_DERIVED | UNCLEAR
  tumor        YES | NO | UNCLEAR
  confounds    subset of CONFOUNDS, each only when the literal fields show it
  evidence     up to 3 short verbatim field values
Only agreement admits: a comparison is RANKABLE when both models return CASE_IS_CONDITION, ORGANISM_DIRECT
and tumor NO and neither lists a confound. Every other comparison is kept with both answers and its
reasons. The orientation question also checks the half of usable rows whose orientation rested on the
label-token heuristic alone.

Candidates are the comparisons the ranking could use: final class VERIFIED or PROVISIONAL, family FDR
below 0.05, at least 3 samples per arm, a rankable category and an agreed matching contrast type.

  agent-api run --provider tokenswitch -- python -m audit_system.design_audit_2026_09_30 \
      --atlas atlas_v3_combined.csv --classification DIR --repairs dnam=DIR --repairs hs_rna=DIR \
      --repairs mm_rna=DIR --out DIR [--dry-run]
"""
import argparse
import json
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from audit_system.comparison_ranking_2026_09_30 import CATEGORY_CONTRAST, FDR, MIN_ARM, UNRANKED_CLOCKS, contrast_agreed
from audit_system.intervention_classification_2026_09_29 import MODALITY_KEYS, OVERVIEW
from audit_system.label_adjudication_2026_09_27 import Ledger, call as ledger_call
from audit_system.membership_recovery_2026_09_27 import GEO, MODALITIES, sha256

SCHEMA = 'clockbase-design-audit-v1'
MODELS = tuple(os.environ.get('AUDIT_MODELS', 'gp-6-as,cl-op-5').split(','))  # see label_adjudication_2026_09_27
ORIENTATIONS = ['CASE_IS_CONDITION', 'REVERSED', 'UNCLEAR']
SPECIMENS = ['ORGANISM_DIRECT', 'EX_VIVO_CULTURE', 'CELL_LINE', 'STEM_CELL_DERIVED', 'UNCLEAR']
TUMOR = ['YES', 'NO', 'UNCLEAR']
CONFOUNDS = ['TISSUE_OR_CELL_TYPE', 'AGE_OR_DEVELOPMENTAL_STAGE', 'TIME_POINT', 'SEX', 'BATCH_OR_PLATFORM',
             'ADDITIONAL_TREATMENT', 'OTHER']
KEY = ['modality', 'gse_id', 'comparison_name', 'condition_name']
FIELD_PREFIXES = ('title', 'source_name', 'organism', 'characteristics', 'treatment_protocol', 'growth_protocol',
                  'extract_protocol', 'description', 'molecule', 'platform_id', 'instrument_model', 'library_strategy')
LONG_FIELDS = ('treatment_protocol', 'growth_protocol', 'extract_protocol', 'description')
SAMPLES_SHOWN, VALUE_CHARS, LONG_CHARS = 8, 160, 400
BATCH_ITEMS, BATCH_CHARS = 6, 24000

PROMPT = """You audit the design of comparisons from an aging-clock atlas. For each comparison you see the study
title and summary and the literal GEO records of the samples in the case arm and the control arm. No results are
shown; do not guess them. Judge only from the records.

orientation: CASE_IS_CONDITION if the case arm is the labelled condition and the control arm its comparator,
  REVERSED if the control arm is the labelled condition, UNCLEAR if the records do not say.
specimen (the samples that were profiled):
  ORGANISM_DIRECT: the case/control difference arose in living organisms (genotype, or a drug, diet, exposure,
    infection or disease present in the animals or people) and the tissue, blood or cells were profiled directly
    after collection; sorting or dissociation is allowed, culture is not.
  EX_VIVO_CULTURE: primary cells or tissue were cultured, stimulated or treated after isolation, including cells
    from treated animals that were cultured before profiling, organ explants, and differences introduced in culture.
  CELL_LINE: immortalized, transformed or cancer cell lines.
  STEM_CELL_DERIVED: pluripotent or embryonic stem cells, cells derived from them, or organoids.
  UNCLEAR: the records do not say.
tumor: YES if either arm consists of tumour or cancer tissue or cancer cells, NO if neither does, UNCLEAR.
confounds: each factor that differs systematically between the arms besides the labelled condition, only when the
  records show it: TISSUE_OR_CELL_TYPE (different tissues or cell types), AGE_OR_DEVELOPMENTAL_STAGE (different age
  groups or stages), TIME_POINT (different collection times or treatment durations), SEX (arms mostly of different
  sex), BATCH_OR_PLATFORM (different platforms, batches or processing), ADDITIONAL_TREATMENT (a second treatment in one
  arm only), OTHER. Return [] when the arms differ only by the labelled condition or the records show no difference.
evidence: up to 3 short verbatim field values supporting your answers.

Comparisons:
{items}

Return JSON only: {{"comparisons": [{{"id": <id>, "orientation": "...", "specimen": "...", "tumor": "...",
"confounds": [...], "evidence": ["...", "..."]}}]}}"""


def meta_candidates(atlas_path, classification_dir, exclude_dirs=(), context='in_vivo'):
    """Comparisons of the intervention meta-analysis, minus those already audited in exclude_dirs.

    Added 2026-10-16. The ranking excludes tumour material and cultured samples, the meta-analysis did not, so
    its Methods claim of samples profiled directly from organisms rested on the per-sample context label alone.
    The meta-analysis admits comparisons the ranking never considered (no significance filter, two samples per
    arm), so auditing its own comparisons is what lets the same exclusion be applied to every intervention.
    """
    from audit_system.intervention_meta_2026_09_29 import usable_effects
    rows = usable_effects(atlas_path, classification_dir, context)
    rows = rows.sort_values(['v3_n_case', 'v3_n_control'], ascending=False).drop_duplicates(KEY)
    for folder in exclude_dirs:
        done = pd.read_json(Path(folder) / 'comparisons.jsonl', lines=True)[KEY]
        rows = rows.merge(done.assign(_done=True), on=KEY, how='left')
        rows = rows[rows._done.isna()].drop(columns='_done')
    return rows[KEY + ['condition_category', 'atlas_row_ordinal']].reset_index(drop=True)


def candidates(atlas_path, classification_dir, exclude_dirs=()):
    """Ranking candidates (the ranking's eligibility before the design audit), minus comparisons audited in exclude_dirs."""
    atlas = pd.read_csv(atlas_path, low_memory=False)
    classes = pd.read_json(Path(classification_dir) / 'comparisons.jsonl', lines=True)
    atlas = atlas.merge(contrast_agreed(classes)[KEY + ['contrast_type']], on=KEY, how='left')
    base = (atlas.v3_final_class.isin(['VERIFIED', 'PROVISIONAL']) & (pd.to_numeric(atlas.v3_fdr, errors='coerce') < FDR)
            & ~atlas.clock.isin(UNRANKED_CLOCKS)
            & (atlas.v3_n_case >= MIN_ARM) & (atlas.v3_n_control >= MIN_ARM) & atlas.condition_category.isin(CATEGORY_CONTRAST)
            & (atlas.contrast_type == atlas.condition_category.map(CATEGORY_CONTRAST)))
    rows = atlas[base].sort_values(['v3_n_case', 'v3_n_control'], ascending=False).drop_duplicates(KEY)
    for folder in exclude_dirs:
        done = pd.read_json(Path(folder) / 'comparisons.jsonl', lines=True)[KEY]
        rows = rows.merge(done.assign(_done=True), on=KEY, how='left')
        rows = rows[rows._done.isna()].drop(columns='_done')
    return rows[KEY + ['condition_category', 'atlas_row_ordinal']].reset_index(drop=True)


def sample_records(modality_key, gse, cache):
    if gse not in cache:
        path = GEO / MODALITIES[modality_key]['folder'] / 'metadata' / f'{gse}.csv'
        frame = None
        if path.exists():
            frame = pd.read_csv(path, dtype=str)
            index = 'geo_accession' if 'geo_accession' in frame.columns else frame.columns[0]
            frame = frame.drop_duplicates(index).set_index(index)
            frame = frame[[c for c in frame.columns if c.startswith(FIELD_PREFIXES)]]
        if len(cache) > 500:
            cache.clear()
        cache[gse] = frame
    return cache[gse]


def clip(column, value):
    limit = LONG_CHARS if column.startswith(LONG_FIELDS) else VALUE_CHARS
    text = re.sub(r'\s+', ' ', str(value)).strip()
    return text if len(text) <= limit else text[:limit] + '...'


def describe_arms(records, case, control):
    """Fields identical across every sample of both arms once, then the varying fields per sample."""
    if records is None:
        return 'GEO sample records unavailable.'
    both = records.reindex(list(case) + list(control))
    shared, varying = [], []
    for column in both.columns:
        values = both[column].dropna().astype(str)
        if values.empty:
            continue
        (shared if values.nunique() == 1 and len(values) == len(both) else varying).append(column)
    lines = ['Fields identical in every sample of both arms:']
    lines += [f'  {c}: {clip(c, both[c].dropna().iloc[0])}' for c in shared] or ['  none']
    for name, gsms in (('Case arm', case), ('Control arm', control)):
        lines.append(f'{name}, {len(gsms)} samples:')
        for gsm in list(gsms)[:SAMPLES_SHOWN]:
            row = records.reindex([gsm]).iloc[0] if gsm in records.index else None
            values = ('; '.join(f'{c}={clip(c, row[c])}' for c in varying if pd.notna(row[c]))
                      if row is not None else 'record missing')
            lines.append(f'  - {gsm}: {values}')
        if len(gsms) > SAMPLES_SHOWN:
            lines.append(f'  - ... {len(gsms) - SAMPLES_SHOWN} more samples')
    return '\n'.join(lines)


def parse(text):
    """Whole JSON first; otherwise salvage complete comparison objects (answers cut at the token limit)."""
    match = re.search(r'\{.*\}', text or '', re.S)
    if match:
        try:
            return {int(c['id']): c for c in json.loads(match.group(0)).get('comparisons', []) if 'id' in c}
        except (json.JSONDecodeError, TypeError, ValueError, AttributeError):
            pass
    salvaged = {}
    for chunk in re.split(r'(?=\{\s*"id"\s*:)', text or ''):
        found = re.match(r'\{\s*"id"\s*:\s*(\d+)', chunk)
        if not found:
            continue
        item = {'id': int(found.group(1))}
        for field, options in (('orientation', ORIENTATIONS), ('specimen', SPECIMENS), ('tumor', TUMOR)):
            value = re.search(rf'"{field}"\s*:\s*"({"|".join(options)})"', chunk)
            if value:
                item[field] = value.group(1)
        listed = re.search(r'"confounds"\s*:\s*\[([^\]]*)\]', chunk)
        if 'orientation' in item and 'specimen' in item and 'tumor' in item and listed:
            item['confounds'] = [c for c in CONFOUNDS if f'"{c}"' in listed.group(1)]
            item['parse'] = 'field_extraction'
            salvaged[item['id']] = item
    return salvaged or None


def normalize(answer):
    orientation = str(answer.get('orientation', '')).upper()
    specimen = str(answer.get('specimen', '')).upper()
    tumor = str(answer.get('tumor', '')).upper()
    confounds = answer.get('confounds') if isinstance(answer.get('confounds'), list) else None
    return {'orientation': orientation if orientation in ORIENTATIONS else 'UNCLEAR',
            'specimen': specimen if specimen in SPECIMENS else 'UNCLEAR', 'tumor': tumor if tumor in TUMOR else 'UNCLEAR',
            'confounds': sorted({str(c).upper() for c in confounds} & set(CONFOUNDS)) if confounds is not None else None,
            'evidence': answer.get('evidence')}


def reconcile(answers):
    a, b = (normalize(x) for x in answers)
    reasons = []
    if a['orientation'] != b['orientation']:
        reasons.append('ORIENTATION_DISAGREE')
    elif a['orientation'] != 'CASE_IS_CONDITION':
        reasons.append('ORIENTATION_' + a['orientation'])
    if a['specimen'] != b['specimen']:
        reasons.append('SPECIMEN_DISAGREE')
    elif a['specimen'] != 'ORGANISM_DIRECT':
        reasons.append('SPECIMEN_' + a['specimen'])
    if a['tumor'] != b['tumor']:
        reasons.append('TUMOR_DISAGREE')
    elif a['tumor'] != 'NO':
        reasons.append('TUMOR_' + a['tumor'])
    if a['confounds'] is None or b['confounds'] is None:
        reasons.append('CONFOUNDS_MISSING')
    else:
        for c in sorted(set(a['confounds']) | set(b['confounds'])):
            reasons.append(('CONFOUND_AGREED_' if c in a['confounds'] and c in b['confounds'] else 'CONFOUND_ONE_MODEL_') + c)
    return {'status': 'NOT_RANKABLE' if reasons else 'RANKABLE', 'reasons': reasons,
            'orientation': a['orientation'] if a['orientation'] == b['orientation'] else None,
            'specimen': a['specimen'] if a['specimen'] == b['specimen'] else None,
            'tumor': a['tumor'] if a['tumor'] == b['tumor'] else None, 'answers': [a, b]}


def arms_from_repairs(repairs):
    arms = {}
    for mkey, run_dir in repairs.items():
        for line in open(Path(run_dir) / 'repairs.jsonl'):
            r = json.loads(line)
            inputs = r.get('inputs') or {}
            arms[(mkey, r['row_ordinal'])] = (inputs.get('case_gsms') or [], inputs.get('control_gsms') or [])
    return arms


def build_items(frame, repairs):
    arms = arms_from_repairs(repairs)
    titles = {}
    for mkey, name in OVERVIEW.items():
        overview = pd.read_csv(GEO / name, usecols=['accession', 'title', 'description'], dtype=str).drop_duplicates('accession')
        titles[mkey] = overview.set_index('accession')[['title', 'description']].to_dict('index')
    cache, items = {}, []
    for i, c in frame.iterrows():
        mkey = MODALITY_KEYS[c.modality]
        case, control = arms.get((mkey, int(c.atlas_row_ordinal)), ([], []))
        series = titles[mkey].get(c.gse_id, {})
        text = (f'Comparison {i}. [{c.modality}] study {c.gse_id}: "{clip("title", series.get("title", ""))}". '
                f'Summary: {clip("description", series.get("description", ""))}\n'
                f'Atlas label: comparison "{c.comparison_name}", condition "{c.condition_name}", category {c.condition_category}.\n'
                + describe_arms(sample_records(mkey, c.gse_id, cache), case, control))
        items.append({'id': int(i), **{k: c[k] for k in KEY + ['condition_category']}, 'text': text,
                      'n_case': len(case), 'n_control': len(control)})
    return items


def pack(items):
    batches, current, size = [], [], 0
    for item in items:  # candidates are ordered by study, so a study's comparisons usually share a request
        if current and (len(current) >= BATCH_ITEMS or size + len(item['text']) > BATCH_CHARS):
            batches.append(current)
            current, size = [], 0
        current.append(item)
        size += len(item['text'])
    return batches + ([current] if current else [])


def run(atlas_path, classification_dir, repairs, out_dir, workers=12, max_requests=6000, dry_run=False, max_tokens=4000,
        exclude_dirs=(), meta=False):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    select = meta_candidates if meta else candidates
    frame = select(atlas_path, classification_dir, exclude_dirs).sort_values(['modality', 'gse_id']).reset_index(drop=True)
    frame.to_csv(out_dir / 'candidates.csv', index=False)
    items = build_items(frame, repairs)
    batches = pack(items)
    prompts = [PROMPT.format(items='\n\n'.join(i['text'] for i in b)) for b in batches]
    if dry_run:
        report = {'candidates': len(items), 'requests_per_model': len(batches),
                  'prompt_chars_total': sum(map(len, prompts)), 'prompt_chars_max': max(map(len, prompts), default=0)}
        (out_dir / 'dry_run.json').write_text(json.dumps(report, indent=2))
        (out_dir / 'dry_run_example_prompt.txt').write_text(prompts[0] if prompts else '')
        print(json.dumps(report, indent=1))
        return
    ledger = Ledger(out_dir / 'ledger.jsonl')
    budget = {'max': max_requests, 'sent': 0, 'lock': threading.Lock()}

    def audit(job):
        batch, prompt = job
        answers, receipts = {}, []
        for model in MODELS:
            record, cached = ledger_call(model, prompt, ledger, budget, parse=parse, max_tokens=max_tokens)
            receipts.append({'model': model, 'status': record.get('status'), 'request_sha256': record.get('request_sha256'),
                             'cached': cached})
            parsed = record.get('parsed') if record.get('status') == 'OK' else None
            answers[model] = {int(k): v for k, v in parsed.items()} if parsed else None
        out = []
        for item in batch:
            if not all(answers.values()):
                result = {'status': 'MODEL_CALL_FAILED'}
            elif not all(item['id'] in answers[m] for m in MODELS):
                result = {'status': 'MODEL_ANSWER_MISSING'}
            else:
                result = reconcile([answers[m][item['id']] for m in MODELS])
            out.append({**{k: item[k] for k in KEY + ['condition_category', 'n_case', 'n_control']}, **result,
                        'receipts': receipts})
        return out

    results = []
    with ThreadPoolExecutor(workers) as pool:
        for chunk in pool.map(audit, zip(batches, prompts)):
            results.extend(chunk)
    with open(out_dir / 'comparisons.jsonl', 'w') as sink:
        for r in results:
            sink.write(json.dumps(r, default=str) + '\n')
    table = pd.DataFrame(results)
    reasons = pd.Series([x for r in results for x in r.get('reasons', [])]).value_counts().to_dict()
    empty = table.empty  # every candidate already audited
    summary = {'schema': SCHEMA, 'candidates': len(table), 'status': {} if empty else table.status.value_counts().to_dict(),
               'rankable_by_category': {} if empty else
               table[table.status == 'RANKABLE'].condition_category.value_counts().to_dict(),
               'reasons': reasons, 'requests_sent': budget['sent'], 'models': MODELS,
               'atlas_sha256': sha256(atlas_path), 'classification_sha256': sha256(Path(classification_dir) / 'comparisons.jsonl'),
               'excluded_already_audited': [str(d) for d in exclude_dirs],
               'module_sha256': sha256(__file__), 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=1, default=str))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--atlas', type=Path, required=True)
    parser.add_argument('--classification', type=Path, required=True)
    parser.add_argument('--repairs', action='append', required=True, help='modality=DIR (repairs.jsonl gives the arms)')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=12)
    parser.add_argument('--max-requests', type=int, default=6000)
    parser.add_argument('--max-tokens', type=int, default=4000)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--exclude', type=Path, action='append', default=[],
                        help='earlier design-audit run whose comparisons are skipped')
    parser.add_argument('--meta', action='store_true',
                        help='audit the comparisons of the intervention meta-analysis instead of the ranking '
                             'candidates, so tumour material can be excluded from every intervention')
    args = parser.parse_args(argv)
    run(args.atlas, args.classification, dict(item.split('=', 1) for item in args.repairs), args.out, args.workers,
        args.max_requests, args.dry_run, args.max_tokens, args.exclude, args.meta)


if __name__ == '__main__':
    main()
