"""Dual-model classification of atlas comparisons, canonical intervention identity and a relevance score.

For every usable comparison (verifier-v3 VERIFIED or PROVISIONAL), GPT-6 Astra and Claude Opus 5
independently read the study title, the comparison labels and the literal metadata values that
distinguish the case arm from the control arm. P values and effect sizes are never shown, so the
relevance score is independent of the statistical evidence. Each model returns:

  contrast_type      INTERVENTION | DISEASE_VS_HEALTHY | TISSUE_OR_CELL_TYPE | DEVELOPMENT_OR_TIME |
                     CELL_STATE | AGE | OTHER
  intervention_class GENE_LOSS | GENE_GAIN | GENE_MODULATION | DRUG | DIET | ENVIRONMENT | PATHOGEN |
                     CELL_OR_TISSUE_TRANSFER | OTHER | NA
  canonical_target   gene symbol (upper case) or drug/diet/exposure name (lower case); NA when not an intervention
  manipulation       DECREASE | INCREASE | EXPOSURE | NA
  relevance_score    0-10 by the published rubric below
  relevance_codes    subset of RELEVANCE_CODES

A comparison is AGREED when both models give the same contrast type, class, manipulation and
normalized target. Its relevance score is the mean of the two scores; their absolute difference is kept
as the reliability measure. Disagreements are kept with both answers.

  agent-api run --provider tokenswitch -- python -m audit_system.intervention_classification_2026_09_29 \
      --atlas atlas_v3_combined.csv --out DIR [--workers 12] [--max-requests N]
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

from audit_system.label_adjudication_2026_09_27 import Ledger, call as ledger_call
from audit_system.membership_recovery_2026_09_27 import GEO, sha256

SCHEMA = 'clockbase-intervention-classification-v1'
MODELS = tuple(os.environ.get('AUDIT_MODELS', 'gp-6-as,cl-op-5').split(','))  # see label_adjudication_2026_09_27
BATCH = 20
MODALITY_KEYS = {'human DNAm': 'dnam', 'human RNA-seq': 'hs_rna', 'mouse RNA-seq': 'mm_rna'}
FOLDERS = {'dnam': 'GEO_hs_dnam', 'hs_rna': 'GEO_hs_rna', 'mm_rna': 'GEO_mm_rna'}
OVERVIEW = {'dnam': 'gse_me_hs.csv', 'hs_rna': 'gse_rna_hs.csv', 'mm_rna': 'gse_rna_mm.csv'}
CONTRAST_TYPES = ['INTERVENTION', 'DISEASE_VS_HEALTHY', 'TISSUE_OR_CELL_TYPE', 'DEVELOPMENT_OR_TIME', 'CELL_STATE',
                  'AGE', 'OTHER']
CLASSES = ['GENE_LOSS', 'GENE_GAIN', 'GENE_MODULATION', 'DRUG', 'DIET', 'ENVIRONMENT', 'PATHOGEN',
           'CELL_OR_TISSUE_TRANSFER', 'OTHER', 'NA']
MANIPULATIONS = ['DECREASE', 'INCREASE', 'EXPOSURE', 'NA']
RELEVANCE_CODES = ['IN_VIVO', 'PRIMARY_TISSUE_OR_CELLS', 'CELL_LINE', 'PLURIPOTENT_OR_REPROGRAMMING', 'TUMOR',
                   'AGING_RELEVANT_TISSUE', 'ESTABLISHED_AGING_PATHWAY', 'TRANSLATABLE_INTERVENTION',
                   'CONTROL_MATCHED', 'CONFOUNDED_DESIGN']

RUBRIC = """Relevance score (0-10) for using this comparison as evidence about an intervention that changes the pace
of organismal aging. Judge only study design and biological context, never the outcome:
  +3 in vivo organism (or +1 primary human/animal tissue or cells ex vivo); 0 for immortalized or tumor cell lines
  +2 the case arm differs from a matched control only by the intervention (same tissue, age, sex, time point)
  +2 the intervention acts on an established aging pathway (e.g. mTOR, insulin/IGF-1, sirtuins, AMPK, senescence,
     inflammation, mitochondrial function, NAD+, proteostasis, epigenetic regulators)
  +2 translatable intervention (approved or clinical-stage drug, diet, exercise, or a druggable gene target)
  +1 aging-relevant tissue or model (e.g. aged animals, progeroid or long-lived models, tissues that decline with age)
  Score 0-2 for comparisons that are not interventions (disease vs healthy, tissue vs tissue, developmental stages,
  reprogramming or cell states)."""

PROMPT = """You classify comparisons from an aging-clock atlas. For each comparison you see the study title, the
comparison labels, and the literal sample-metadata values that distinguish the case arm from the control arm.
No results are shown; do not guess them.

Definitions:
contrast_type: INTERVENTION (a deliberate genetic, drug, diet, environmental or pathogen manipulation of the case arm),
DISEASE_VS_HEALTHY, TISSUE_OR_CELL_TYPE, DEVELOPMENT_OR_TIME, CELL_STATE (e.g. reprogramming, differentiation,
activation), AGE, OTHER.
intervention_class: {classes} (NA unless contrast_type is INTERVENTION).
canonical_target: for genes the official symbol in upper case without species prefix (e.g. SIRT6); for drugs the
generic name in lower case (e.g. rapamycin); for diets/exposures a short lower-case name (e.g. caloric restriction,
cigarette smoke); NA unless contrast_type is INTERVENTION.
manipulation: DECREASE (knockout, knockdown, inhibition, deletion), INCREASE (overexpression, activation, agonist),
EXPOSURE (drug, diet or environment given), NA.
{rubric}
relevance_codes: any of {codes}.

Comparisons:
{items}

Return JSON only: {{"comparisons": [{{"id": <id>, "contrast_type": "...", "intervention_class": "...",
"canonical_target": "...", "manipulation": "...", "relevance_score": <0-10>, "relevance_codes": [...]}}]}}"""


def normalize_target(value):
    text = str(value or 'NA').strip()
    if text.upper() in ('NA', 'N/A', 'NONE', ''):
        return 'NA'
    text = re.sub(r'^(human|mouse|murine|h|m)[\s_-]+', '', text, flags=re.I)
    return re.sub(r'[\s_\-]+', ' ', text).strip().casefold()


def arm_contrast(meta, case, control, limit=6):
    """Metadata fields whose values differ between arms, as literal value sets."""
    if meta is None:
        return ''
    a, b = meta.reindex(case), meta.reindex(control)
    lines = []
    for column in meta.columns:
        va = sorted(set(a[column].dropna().astype(str)))[:4]
        vb = sorted(set(b[column].dropna().astype(str)))[:4]
        if va and vb and set(va) != set(vb):
            lines.append(f'{column}: case={va} control={vb}')
        if len(lines) >= limit:
            break
    return '; '.join(lines)


def load_meta(key, gse, cache):
    if (key, gse) not in cache:
        frame = None
        for sub in ('processed_metadata_preharmonize', 'processed_metadata'):
            path = GEO / FOLDERS[key] / sub / f'{gse}.csv'
            if path.exists():
                frame = pd.read_csv(path, dtype=str)
                frame = frame.drop_duplicates(frame.columns[0]).set_index(frame.columns[0])
                break
        if len(cache) > 2000:
            cache.clear()
        cache[(key, gse)] = frame
    return cache[(key, gse)]


def parse(text):
    """Whole JSON first; otherwise salvage every complete comparison object (answers cut at the token limit)."""
    match = re.search(r'\{.*\}', text or '', re.S)
    if match:
        try:
            return {int(c['id']): c for c in json.loads(match.group(0)).get('comparisons', []) if 'id' in c}
        except (json.JSONDecodeError, TypeError, ValueError, AttributeError):
            pass
    salvaged = {}
    for chunk in re.findall(r'\{[^{}]*"id"\s*:[^{}]*\}', text or ''):
        try:
            item = json.loads(chunk)
            salvaged[int(item['id'])] = item
        except (json.JSONDecodeError, TypeError, ValueError, KeyError):
            continue
    return salvaged or None


GENE_CLASSES = {'GENE_LOSS', 'GENE_GAIN', 'GENE_MODULATION'}


def intervention_key(answer):
    """Meta-analysis grouping key. Tolerates schema granularity, not substance:
    genetic interventions need the same targets and direction (loss vs gain); drug, diet, environmental,
    pathogen and other exposures need the same target set; combinations compare as sets."""
    contrast = str(answer.get('contrast_type', '')).upper()
    if contrast != 'INTERVENTION':
        return (contrast, 'NA', 'NA', 'NA')
    targets = sorted({normalize_target(t) for t in re.split(r'\s*(?:[;,+/&]|\band\b|\bplus\b)\s*', str(answer.get('canonical_target') or ''))
                      if normalize_target(t) != 'NA'})
    cls, manipulation = str(answer.get('intervention_class', 'NA')).upper(), str(answer.get('manipulation', 'NA')).upper()
    if cls in GENE_CLASSES:
        if (cls, manipulation) in (('GENE_LOSS', 'INCREASE'), ('GENE_GAIN', 'DECREASE')):
            return (contrast, 'GENE', 'CONFLICT', str(id(answer)))  # self-contradictory answer never agrees
        direction = ('DECREASE' if cls == 'GENE_LOSS' or manipulation == 'DECREASE' else
                     'INCREASE' if cls == 'GENE_GAIN' or manipulation == 'INCREASE' else 'NA')
        return (contrast, 'GENE', direction, ';'.join(targets) or 'NA')
    return (contrast, 'EXPOSURE', 'EXPOSURE', ';'.join(targets) or 'NA')


def reconcile(answers):
    parts = []
    for answer in answers:
        if not answer:
            return {'status': 'MODEL_ANSWER_MISSING'}
        parts.append(intervention_key(answer))
    scores = [a.get('relevance_score') for a in answers]
    scores = [float(s) for s in scores if isinstance(s, (int, float))]
    result = {'answers': answers, 'relevance_score': sum(scores) / len(scores) if len(scores) == 2 else None,
              'relevance_score_abs_diff': abs(scores[0] - scores[1]) if len(scores) == 2 else None,
              'relevance_codes_agreed': sorted(set(answers[0].get('relevance_codes') or []) & set(answers[1].get('relevance_codes') or []))}
    if parts[0] == parts[1]:
        result.update(status='AGREED', contrast_type=parts[0][0], intervention_class=parts[0][1],
                      manipulation=parts[0][2], canonical_target=parts[0][3])
    else:
        result.update(status='MODELS_DISAGREE', contrast_type=parts[0][0] if parts[0][0] == parts[1][0] else None)
    return result


def run(atlas_path, out_dir, workers=12, max_requests=20000, limit=None, repairs=None, batch=BATCH, max_tokens=6000):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    atlas = pd.read_csv(atlas_path, low_memory=False)
    usable = atlas[atlas.v3_final_class.isin(['VERIFIED', 'PROVISIONAL'])]
    key = ['modality', 'gse_id', 'comparison_name', 'condition_name']
    comparisons = usable.drop_duplicates(key)[key + ['condition_category', 'atlas_row_ordinal']].reset_index(drop=True)
    if limit:
        comparisons = comparisons.head(limit)
    titles = {}
    for mkey, name in OVERVIEW.items():
        overview = pd.read_csv(GEO / name, usecols=['accession', 'title'], dtype=str)
        titles[mkey] = overview.drop_duplicates('accession').set_index('accession').title.to_dict()
    members = {}
    for mkey, run_dir in (repairs or {}).items():
        path = Path(run_dir) / 'repairs.jsonl'
        if path.exists():
            for line in open(path):
                r = json.loads(line)
                inputs = r.get('inputs') or {}
                members[(mkey, r['row_ordinal'])] = (inputs.get('case_gsms') or [], inputs.get('control_gsms') or [])
    previous = {}
    if (out_dir / 'comparisons.jsonl').exists():
        # Incremental: comparisons already answered by both models are kept (re-reconciled under the current
        # rule, no new request); only the rest are sent.
        for line in open(out_dir / 'comparisons.jsonl'):
            r = json.loads(line)
            if r.get('status') in ('AGREED', 'MODELS_DISAGREE') and r.get('answers'):
                r.update({k: None for k in ('contrast_type', 'intervention_class', 'manipulation', 'canonical_target')})
                r.update(reconcile(r['answers']))
                previous[(r['modality'], r['gse_id'], r['comparison_name'], str(r['condition_name']))] = r
    meta_cache, items = {}, []
    for i, c in comparisons.iterrows():
        if (c.modality, c.gse_id, c.comparison_name, str(c.condition_name)) in previous:
            continue
        mkey = MODALITY_KEYS[c.modality]
        case, control = members.get((mkey, int(c.atlas_row_ordinal)), ([], []))
        items.append({'id': int(i), 'modality': c.modality, 'gse_id': c.gse_id, 'comparison_name': c.comparison_name,
                      'condition_name': c.condition_name, 'category': c.condition_category,
                      'title': titles[mkey].get(c.gse_id, ''),
                      'contrast': arm_contrast(load_meta(mkey, c.gse_id, meta_cache), case, control)})
    ledger = Ledger(out_dir / 'ledger.jsonl')
    budget = {'max': max_requests, 'sent': 0, 'lock': threading.Lock()}

    def batch_prompt(batch):
        lines = [f'{b["id"]}. [{b["modality"]}; {b["gse_id"]}: "{b["title"][:160]}"] comparison "{b["comparison_name"]}", '
                 f'condition "{b["condition_name"]}", category {b["category"]}. Arms differ in: {b["contrast"] or "not recorded"}'
                 for b in batch]
        return PROMPT.format(classes=CLASSES, rubric=RUBRIC, codes=RELEVANCE_CODES, items='\n'.join(lines))

    def classify(batch):
        prompt = batch_prompt(batch)
        answers = {}
        receipts = []
        for model in MODELS:
            record, cached = ledger_call(model, prompt, ledger, budget, parse=parse, max_tokens=max_tokens)
            receipts.append({'model': model, 'status': record.get('status'), 'request_sha256': record.get('request_sha256'), 'cached': cached})
            parsed = record.get('parsed') if record.get('status') == 'OK' else None
            answers[model] = {int(k): v for k, v in parsed.items()} if parsed else None
        out = []
        for b in batch:
            if not all(answers.values()):
                result = {'status': 'MODEL_CALL_FAILED'}
            elif not all(b['id'] in answers[m] for m in MODELS):
                result = {'status': 'MODEL_ANSWER_TRUNCATED'}  # re-requested in a smaller batch on the next run
            else:
                result = reconcile([answers[m].get(b['id']) for m in MODELS])
            out.append({**{k: b[k] for k in ('modality', 'gse_id', 'comparison_name', 'condition_name', 'category')},
                        **result, 'receipts': receipts})
        return out

    batches = [items[i:i + batch] for i in range(0, len(items), batch)]
    results = list(previous.values())
    with ThreadPoolExecutor(workers) as pool:
        for chunk in pool.map(classify, batches):
            results.extend(chunk)
    with open(out_dir / 'comparisons.jsonl', 'w') as sink:
        for r in results:
            sink.write(json.dumps(r, default=str) + '\n')
    frame = pd.DataFrame(results)
    agreed = frame[frame.status == 'AGREED']
    summary = {'schema': SCHEMA, 'comparisons': len(frame), 'reused_answers': len(previous),
               'status': frame.status.value_counts().to_dict(),
               'agreed_contrast_types': agreed.contrast_type.value_counts().to_dict() if len(agreed) else {},
               'relevance_abs_diff_mean': float(frame.relevance_score_abs_diff.dropna().mean()) if 'relevance_score_abs_diff' in frame else None,
               'models': MODELS, 'atlas_sha256': sha256(atlas_path), 'module_sha256': sha256(__file__),
               'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=1, default=str))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--atlas', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=12)
    parser.add_argument('--max-requests', type=int, default=20000)
    parser.add_argument('--limit', type=int)
    parser.add_argument('--repairs', action='append', required=True, help='modality=DIR (repairs.jsonl gives the arms)')
    parser.add_argument('--batch', type=int, default=BATCH)
    parser.add_argument('--max-tokens', type=int, default=6000)
    args = parser.parse_args(argv)
    run(args.atlas, args.out, args.workers, args.max_requests, args.limit,
        dict(item.split('=', 1) for item in args.repairs), args.batch, args.max_tokens)


if __name__ == '__main__':
    main()
