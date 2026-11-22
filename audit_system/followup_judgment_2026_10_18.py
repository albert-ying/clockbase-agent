"""Model judgement of the ranked comparisons: the questions a scientist asks before choosing an intervention to test.

The comparison ranking (comparison_ranking_2026_09_30) orders verified comparisons by the conservative bound of their
effect; models enter it only to establish validity (contrast type, design audit). This module adds a
judgement version beside the statistical one. Two models independently read the series title and summary and the
literal GEO records of both arms (oriented so that the case arm is the labelled condition), never P values or effect
sizes, and answer three questions per ranked comparison:

  age_at_start  OLD | ADULT | YOUNG | GERMLINE_OR_PRENATAL | NOT_STATED
                age of the subjects when the condition of the case arm began (treatment, diet, exposure, induction
                of a disease or of a conditional genetic change). OLD: mice or rats aged 18 months or more, humans
                aged 60 years or more, other species in the last third of their lifespan. GERMLINE_OR_PRENATAL:
                present from conception or before birth (constitutive mutations, prenatal exposures).
  background    HEALTHY | DISEASE_MODEL | NOT_STATED
                apart from the labelled condition, were the subjects of both arms otherwise healthy? DISEASE_MODEL
                when both arms carry an induced or genetic disease, injury, infection, tumour, diet-induced obesity
                or another pathological perturbation (for example treated and untreated mice with induced colitis).
  control       UNTREATED | ACTIVE_COMPARATOR | NOT_STATED
                is the control arm untreated, vehicle-, sham- or placebo-treated, wild type or healthy (UNTREATED),
                or does it receive a different active treatment, vaccine or adjuvant (ACTIVE_COMPARATOR)?

Only agreement counts. The judgement version of each table keeps the ranked comparisons for which both models give the
admitting answers and orders them by the unchanged statistical bound:
  clock outputs lower (Table 2): age_at_start OLD, background HEALTHY, control UNTREATED: interventions that lower
    clock outputs when started in old, otherwise healthy organisms, against untreated controls.
  clock outputs higher (Table 1): background HEALTHY, control UNTREATED.
Comparisons that a model leaves out of its answer to a batch are re-sent once in a new request. The questions, the
two models and the admission rules were committed before any model call (2026-10-18). If the
provider refuses the cross-provider pair for quota, the run uses GPT-6 Astra and GPT-6 Sol through the OpenAI API, as
the design audit did (AUDIT_MODELS=openai:gpt-6-astra,openai:gpt-6-sol). The statistical ranking is not changed.

  agent-api run --provider tokenswitch -- python -m audit_system.followup_judgment_2026_10_18 \
      --final final-2026-10-10 --out DIR [--dry-run]

Use --cached-only SOURCE to rebuild the tables from SOURCE/comparisons.jsonl (or
the JSONL file itself), matching comparison identities without making model calls.
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

from audit_system import paths
from audit_system.design_audit_2026_09_30 import KEY, clip, describe_arms, pack, sample_records
from audit_system.intervention_classification_2026_09_29 import MODALITY_KEYS, OVERVIEW
from audit_system.label_adjudication_2026_09_27 import Ledger, call as ledger_call
from audit_system.membership_recovery_2026_09_27 import GEO, sha256

SCHEMA = 'clockbase-followup-judgment-v1'
MODELS = tuple(os.environ.get('AUDIT_MODELS', 'gp-6-as,cl-op-5').split(','))
AGES = ['OLD', 'ADULT', 'YOUNG', 'GERMLINE_OR_PRENATAL', 'NOT_STATED']
BACKGROUNDS = ['HEALTHY', 'DISEASE_MODEL', 'NOT_STATED']
CONTROLS = ['UNTREATED', 'ACTIVE_COMPARATOR', 'NOT_STATED']
FIELDS = {'age_at_start': AGES, 'background': BACKGROUNDS, 'control': CONTROLS}
ADMIT = {'age-lowering': {'age_at_start': 'OLD', 'background': 'HEALTHY', 'control': 'UNTREATED'},
         'age-raising': {'background': 'HEALTHY', 'control': 'UNTREATED'}}

PROMPT = """You help choose interventions for aging experiments. For each comparison from an aging-clock atlas you see
the study title and summary and the literal GEO records of the samples in the case arm (the labelled condition) and
the control arm. No results are shown; do not guess them. Judge only from the records and the study description.

age_at_start: age of the subjects when the condition of the case arm began (treatment, diet, exposure, induction of a
  disease or of a conditional genetic change).
  OLD: mice or rats aged 18 months or more, humans aged 60 years or more, other species in the last third of life.
  ADULT: mature but not old. YOUNG: before maturity (mice younger than 2 months).
  GERMLINE_OR_PRENATAL: present from conception or before birth (constitutive mutations, prenatal exposures).
  NOT_STATED: the records and description do not give it.
background: apart from the labelled condition, were the subjects of both arms otherwise healthy?
  HEALTHY: yes.
  DISEASE_MODEL: both arms carry an induced or genetic disease, injury, infection, tumour, diet-induced obesity or
    another pathological perturbation (for example treated and untreated mice with induced colitis).
  NOT_STATED: the records do not say.
control: UNTREATED if the control arm is untreated, vehicle-, sham- or placebo-treated, wild type or healthy;
  ACTIVE_COMPARATOR if it receives a different active treatment, vaccine or adjuvant; NOT_STATED.
evidence: up to 3 short verbatim field values supporting your answers.

Comparisons:
{items}

Return JSON only: {{"comparisons": [{{"id": <id>, "age_at_start": "...", "background": "...", "control": "...",
"evidence": ["...", "..."]}}]}}"""


def ranked_rows(final):
    """The ranked comparisons (one row per comparison) with their direction and statistical rank."""
    table = pd.read_csv(Path(final) / 'comparison-ranking' / 'ranked_comparisons_all.csv', low_memory=False)
    table['statistical_rank'] = table.groupby('direction_v3').cumcount() + 1
    return table.drop_duplicates(KEY).reset_index(drop=True)


def oriented_arms(final, modality_key, gse, comparison, condition, cache):
    """Case and control samples of a comparison as the verifier tested it (arms swapped where the labels were inverted)."""
    if (modality_key, gse) not in cache:
        cache[(modality_key, gse)] = json.loads((Path(final) / 'public-export' / 'studies' / modality_key / f'{gse}.json').read_text())
    for row in cache[(modality_key, gse)]['rows']:
        if row['comparison_name'] == comparison and row['condition_name'] == condition and row.get('membership'):
            m = row['membership']
            case, control = list(m.get('case_gsms') or []), list(m.get('control_gsms') or [])
            if any(a.get('action') == 'SWAP_ARMS_CONFIRMED_INVERSION' for a in row.get('actions') or []):
                case, control = control, case
            return case, control
    return [], []


def build_items(final, frame):
    titles = {}
    for mkey, name in OVERVIEW.items():
        overview = pd.read_csv(GEO / name, usecols=['accession', 'title', 'description'], dtype=str).drop_duplicates('accession')
        titles[mkey] = overview.set_index('accession')[['title', 'description']].to_dict('index')
    records, studies, items = {}, {}, []
    for i, c in frame.iterrows():
        mkey = MODALITY_KEYS[c.modality]
        case, control = oriented_arms(final, mkey, c.gse_id, c.comparison_name, c.condition_name, studies)
        series = titles[mkey].get(c.gse_id, {})
        text = (f'Comparison {i}. [{c.modality}] study {c.gse_id}: "{clip("title", series.get("title", ""))}". '
                f'Summary: {clip("description", series.get("description", ""))}\n'
                f'Atlas label: comparison "{c.comparison_name}", condition "{c.condition_name}", category {c.condition_category}.\n'
                + describe_arms(sample_records(mkey, c.gse_id, records), case, control))
        items.append({'id': int(i), **{k: c[k] for k in KEY}, 'text': text, 'n_case': len(case), 'n_control': len(control)})
    return items


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
        item = {'id': int(found.group(1)), 'parse': 'field_extraction'}
        for field, options in FIELDS.items():
            value = re.search(rf'"{field}"\s*:\s*"({"|".join(options)})"', chunk)
            if value:
                item[field] = value.group(1)
        if all(field in item for field in FIELDS):
            salvaged[item['id']] = item
    return salvaged or None


def normalize(answer):
    out = {field: str(answer.get(field, '')).upper() for field in FIELDS}
    out = {field: value if value in FIELDS[field] else 'NOT_STATED' for field, value in out.items()}
    return {**out, 'evidence': answer.get('evidence')}


def reconcile(answers):
    """Agreed answer per question (None where the models differ)."""
    a, b = (normalize(x) for x in answers)
    agreed = {field: a[field] if a[field] == b[field] else None for field in FIELDS}
    return {**agreed, 'answers': [a, b]}


def admitted(direction, judgement):
    """True when both models gave every admitting answer for the comparison's direction."""
    rule = ADMIT.get(direction)
    return bool(rule) and all(judgement.get(field) == value for field, value in rule.items())


def judged_tables(ranked, results):
    """The ranked table with the agreed answers, and the judgement version of each direction in statistical order."""
    answers = pd.DataFrame([{**{k: r[k] for k in KEY}, **{f: r.get(f) for f in FIELDS}, 'status': r['status']}
                            for r in results])
    table = ranked.merge(answers, on=KEY, how='left')
    table['judgement_admitted'] = [admitted(d, row) for d, row in zip(table.direction_v3, table[list(FIELDS)].to_dict('records'))]
    versions = {}
    for direction in ADMIT:
        kept = table[(table.direction_v3 == direction) & table.judgement_admitted].copy()
        kept['judgement_rank'] = range(1, len(kept) + 1)
        versions[direction] = kept
    return table, versions


def write_tables(final, out_dir, ranked, results, models, requests_sent, provenance=None):
    """Write the ranking views with the same admission rules for fresh and saved answers."""
    with open(out_dir / 'comparisons.jsonl', 'w') as sink:
        for r in results:
            sink.write(json.dumps(r, default=str) + '\n')
    table, versions = judged_tables(ranked, results)
    table.to_csv(out_dir / 'ranked_with_judgement.csv', index=False)
    for direction, kept in versions.items():
        kept.to_csv(out_dir / f'judgement_{direction}.csv', index=False)
    judged = table[table.status == 'JUDGED']
    summary = {'schema': SCHEMA, 'final': final.name, 'comparisons': len(table), 'models': models,
               'status': table.status.value_counts(dropna=False).to_dict(),
               'agreed': {f: judged[f].value_counts(dropna=False).to_dict() for f in FIELDS},
               'admission_rules': ADMIT,
               'admitted': {d: {'comparisons': len(v), 'by_category': v.condition_category.value_counts().to_dict()}
                            for d, v in versions.items()},
               'requests_sent': requests_sent,
               'ranking_sha256': sha256(final / 'comparison-ranking' / 'ranked_comparisons_all.csv'),
               'module_sha256': sha256(__file__), 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    if provenance:
        summary.update(provenance)
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=1, default=str))


def run_cached(final, out_dir, source):
    """Reuse saved answers by comparison identity; never enter the model-call path."""
    final, out_dir, source = Path(final), Path(out_dir), Path(source)
    if source.is_dir():
        source = source / 'comparisons.jsonl'
    if not source.is_file():
        raise SystemExit(f'Cannot proceed - required data file not found: {source}')
    if source.resolve() == (out_dir / 'comparisons.jsonl').resolve():
        raise SystemExit(f'refusing to overwrite judgement source: {source}')
    identity = lambda record: tuple(None if pd.isna(record[k]) else record[k] for k in KEY)
    try:
        saved = [json.loads(line) for line in source.read_text().splitlines() if line.strip()]
        by_key = {identity(record): record for record in saved}
        if len(by_key) != len(saved):
            raise ValueError(f'duplicate comparison identities in {source}')
        ranked = ranked_rows(final)
        missing = [identity(row) for row in ranked.to_dict('records') if identity(row) not in by_key]
        if missing:
            raise ValueError(f'cached judgement missing for comparison identities: {missing}')
        results = [by_key[identity(row)] for row in ranked.to_dict('records')]
        models = sorted({receipt['model'] for record in results for receipt in record.get('receipts', [])})
        out_dir.mkdir(parents=True, exist_ok=True)
        write_tables(final, out_dir, ranked, results, models, 0,
                     {'cached_only': True, 'source': str(source.resolve()), 'source_sha256': sha256(source),
                      'source_comparisons': len(saved), 'source_comparisons_not_ranked': len(saved) - len(results)})
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(f'Cannot rebuild judgement tables from {source}: {error}') from error


def run(final, out_dir, workers=12, max_requests=2000, dry_run=False, max_tokens=4000):
    final, out_dir = Path(final), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ranked = ranked_rows(final)
    frame = ranked.sort_values(['modality', 'gse_id']).reset_index(drop=True)
    items = build_items(final, frame)
    batches = pack(items)
    prompts = [PROMPT.format(items='\n\n'.join(i['text'] for i in b)) for b in batches]
    if dry_run:
        report = {'comparisons': len(items), 'requests_per_model': len(batches), 'models': MODELS,
                  'prompt_chars_total': sum(map(len, prompts)), 'prompt_chars_max': max(map(len, prompts), default=0),
                  'arms_missing': sum(1 for i in items if not i['n_case'] or not i['n_control'])}
        (out_dir / 'dry_run.json').write_text(json.dumps(report, indent=2, default=str))
        (out_dir / 'dry_run_example_prompt.txt').write_text(prompts[0] if prompts else '')
        print(json.dumps(report, indent=1, default=str))
        return
    ledger = Ledger(out_dir / 'ledger.jsonl')
    budget = {'max': max_requests, 'sent': 0, 'lock': threading.Lock()}

    def judge(job):
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
                result = {'status': 'JUDGED', **reconcile([answers[m][item['id']] for m in MODELS])}
            out.append({**{k: item[k] for k in KEY}, 'n_case': item['n_case'], 'n_control': item['n_control'],
                        **result, 'receipts': receipts})
        return out

    results = []
    with ThreadPoolExecutor(workers) as pool:
        for chunk in pool.map(judge, zip(batches, prompts)):
            results.extend(chunk)
    # One retry, in a new request, for comparisons that a model left out of its answer to a batch (added after the
    # first pass of 2026-10-18 left 5 of 946 unanswered; questions and admission rules unchanged).
    key = lambda record: tuple(record[k] for k in KEY)
    by_key = {key(item): item for item in items}
    missing = [by_key[key(r)] for r in results if r['status'] == 'MODEL_ANSWER_MISSING']
    if missing:
        retry_batches = pack(missing)
        retry_prompts = [PROMPT.format(items='\n\n'.join(i['text'] for i in b)) for b in retry_batches]
        retried = {}
        with ThreadPoolExecutor(workers) as pool:
            for chunk in pool.map(judge, zip(retry_batches, retry_prompts)):
                retried.update({key(r): {**r, 'retried': True} for r in chunk})
        results = [retried.get(key(r), r) for r in results]
    write_tables(final, out_dir, ranked, results, MODELS, budget['sent'])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--final', required=True, help='verifier final run (name under the verifier root, or a path)')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=12)
    parser.add_argument('--max-requests', type=int, default=2000)
    parser.add_argument('--max-tokens', type=int, default=4000)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--cached-only', type=Path,
                        help='saved comparisons JSONL or directory; reuse answers by identity without model calls')
    args = parser.parse_args(argv)
    final = Path(args.final) if Path(args.final).is_dir() else paths.VERIFIER_ROOT / args.final
    if args.cached_only:
        run_cached(final, args.out, args.cached_only)
    else:
        run(final, args.out, args.workers, args.max_requests, args.dry_run, args.max_tokens)


if __name__ == '__main__':
    main()
