"""Two-reader audit of the individuals behind the samples of ranking candidates.

The design-respecting estimand needs to know which samples come from the same individual. A GEO sample field says so
for about one series in fourteen; elsewhere the individual is written into sample titles or the series text (patient
codes at two visits, "biological_rep_1_technical_rep_3", one library on four sequencing lanes). For the comparisons
the ranking could use and whose series has no GEO individual field, two models read the series text, the comparison
names and the GEO fields of the samples in the two arms, and no result. Each states which samples come from the same
individual. Their answers are turned into findings with the rule of design_estimand_2026_10_10 (repeated individual,
paired arms):

  SHARED_AGREED     both readers report shared individuals and give the same partition of the arms' samples; the
                    partition is saved and row repairs recompute the comparison on individual means
  INDEPENDENT       both readers find no finding (no shared individual, or established cell lines)
  SHARED_DISAGREE   both report shared individuals with different partitions or findings; nothing is changed
  UNDETERMINED      a reader cannot tell; nothing is changed
  NOT_READ_ARMS_ABOVE_LIMIT, METADATA_UNAVAILABLE, MODEL_CALL_FAILED

Candidates are the ranking candidates of the design audit (final class VERIFIED or PROVISIONAL, family FDR below 0.05,
at least 3 per arm, a rankable category, an agreed matching contrast type) that the design audit admitted as RANKABLE.
The prompt is the one validated in benchmark/published_clock_2026_10_05/groupbench (285 comparisons outside the
benchmark: the two readers gave the same finding in 99%; 89 of 94 flags from GEO fields confirmed).

  python -m audit_system.individuals_audit_2026_10_10 --atlas atlas_v3_combined.csv --classification DIR \
      --design-audit DIR [--design-audit DIR ...] --repairs dnam=DIR --repairs hs_rna=DIR --repairs mm_rna=DIR \
      --out DIR [--exclude DIR ...] [--cache LEDGER ...] [--dry-run]
"""
import argparse
import gzip
import hashlib
import json
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from audit_system import design_audit_2026_09_30 as design_audit
from audit_system import design_estimand_2026_10_10 as design_estimand
from audit_system.comparison_ranking_2026_09_30 import design_audit_table
from audit_system.intervention_classification_2026_09_29 import MODALITY_KEYS
from audit_system.label_adjudication_2026_09_27 import Ledger, call as ledger_call
from audit_system.membership_recovery_2026_09_27 import GEO, MODALITIES, sha256

SCHEMA = 'clockbase-individuals-audit-v1'
MODELS = tuple(os.environ.get('AUDIT_MODELS', 'openai:gpt-6-astra,openai:gpt-6-sol').split(','))
KEY = ['modality', 'gse_id', 'comparison_name', 'condition_name']
MAX_SAMPLES = 120
MAX_TOKENS = 2000
PROMPT = """You are reading the public GEO record of one series to describe two groups of samples that were compared.

Series {gse}: {title}
Summary: {summary}
Overall design: {design}

Comparison: {comparison}
Condition: {condition}

Fields with one value for all listed samples: {constant}

Group A samples (GEO sample, then its fields):
{group_a}

Group B samples:
{group_b}

Answer three questions from this record only.

1. Individuals. An individual is the person or animal the material was taken from; for primary cells it is the donor.
An established cell line is not an individual. Do any two of the listed samples come from the same individual (for
example the same patient at two visits, two tissues or cell types of one donor, replicate cultures of one donor's
cells, technical replicates of one specimen)? Use a field, a pattern in the titles, or the series text. Give each
sample an individual identifier when the record lets you (the same identifier means the same individual); use null for
a sample whose individual cannot be determined.

2. Ages. Give the age at sampling of each sample's individual when recorded, with one unit for all samples; null when
it is not recorded. Do not use gestational age, age at diagnosis or culture time.

3. Whether age or time is what the comparison contrasts (for example old versus young, or one time point versus
another of a time course), as opposed to a contrast of condition, treatment, genotype, tissue or cell type.

Return one JSON object and nothing else:
{{"material": "PRIMARY_OR_IN_VIVO" | "CELL_LINE" | "UNCLEAR",
 "shared_individuals": "YES" | "NO" | "CANNOT_TELL",
 "individual_evidence": "the field, title pattern or sentence used, in a few words",
 "individuals": {{"GSM...": "identifier or null", ...}},
 "age_unit": "years" | "months" | "weeks" | "days" | null,
 "ages": {{"GSM...": number or null, ...}},
 "age_or_time_is_the_contrast": true | false}}
The "individuals" and "ages" objects must list every sample of both groups."""


def series_text(folder, gse):
    """Title, summary and overall design from the GEO series matrix."""
    out = {'title': '', 'summary': '', 'design': ''}
    keys = {'!Series_title': 'title', '!Series_summary': 'summary', '!Series_overall_design': 'design'}
    files = sorted((folder / 'series_matrix').glob(f'{gse}[_-]*series_matrix.txt.gz'))
    if not files:
        return out
    try:
        with gzip.open(files[0], 'rt', errors='replace') as handle:
            for line in handle:
                if not line.startswith('!Series'):
                    if line.startswith('!Sample') or line.startswith('"ID_REF"'):
                        break
                    continue
                name, _, value = line.rstrip('\n').partition('\t')
                if name in keys:
                    out[keys[name]] = (out[keys[name]] + ' ' + value.strip('"')).strip()
    except OSError as error:
        raise ValueError(f'Cannot read series matrix {files[0]}: {error}') from error
    return {k: v[:1800] for k, v in out.items()}


def sample_lines(raw, gsms):
    """One line per sample with its varying GEO fields, and the fields constant across the listed samples."""
    frame = raw.loc[[g for g in gsms if g in raw.index]]
    constant, varying = [], []
    for column in frame.columns:
        values = {design_estimand.clean(v) for v in frame[column]}
        if len(values) == 1:
            value = next(iter(values))
            if value is not None:
                constant.append(f'{column}={value[:80]}')
        else:
            varying.append(column)
    return constant, varying


def render(raw, gsms, varying):
    lines = []
    for gsm in gsms:
        if gsm not in raw.index:
            lines.append(f'{gsm}: (no GEO fields)')
            continue
        fields = '; '.join(f'{c}={str(raw.at[gsm, c])[:80]}' for c in varying if design_estimand.clean(raw.at[gsm, c]) is not None)
        lines.append(f'{gsm}: {fields}')
    return '\n'.join(lines)


def package(unit, raw, text):
    gsms = unit['case'] + unit['control']
    constant, varying = sample_lines(raw, gsms)
    return PROMPT.format(gse=unit['gse'], title=text['title'] or '(not available)', summary=text['summary'] or '(not available)',
                         design=text['design'] or '(not available)', comparison=unit['comparison_name'],
                         condition=unit['condition_name'], constant='; '.join(constant) or '(none)',
                         group_a=render(raw, unit['case'], varying), group_b=render(raw, unit['control'], varying))


def parse(text):
    """The reader's JSON object, or None."""
    try:
        value = json.loads(text)
    except (TypeError, ValueError):
        start, end = str(text).find('{'), str(text).rfind('}')
        if start < 0 or end <= start:
            return None
        try:
            value = json.loads(str(text)[start:end + 1])
        except ValueError:
            return None
    return value if isinstance(value, dict) else None


def reader_view(case, control, answer):
    """One reader's answer as {'determined', 'findings', 'partition'} over the samples of the two arms."""
    if not isinstance(answer, dict) or not isinstance(answer.get('individuals'), dict):
        return None
    shared = answer.get('shared_individuals')
    cell_line = answer.get('material') == 'CELL_LINE'  # an established cell line is not an individual
    subjects = {str(g): str(v) for g, v in answer['individuals'].items()
                if v is not None and design_estimand.clean(v) is not None and str(g) in set(case) | set(control)}
    if shared == 'NO' or cell_line:
        return {'determined': True, 'findings': [], 'partition': []}
    covered = len(subjects) >= design_estimand.COVERAGE * (len(case) + len(control))
    if shared != 'YES' or not covered:
        return {'determined': False, 'findings': [], 'partition': []}
    groups = {}
    for gsm, individual in subjects.items():
        groups.setdefault(individual, []).append(gsm)
    partition = sorted(sorted(members) for members in groups.values())
    return {'determined': True, 'findings': design_estimand.subject_findings(case, control, subjects), 'partition': partition}


def reconcile(case, control, answers):
    """The decision for one comparison from the two readers' answers."""
    views = [reader_view(case, control, answer) for answer in answers]
    if len(views) != 2 or any(v is None for v in views):
        return {'status': 'MODEL_CALL_FAILED'}
    if not all(v['determined'] for v in views):
        return {'status': 'UNDETERMINED'}
    first, second = views
    if not first['findings'] and not second['findings']:
        return {'status': 'INDEPENDENT'}
    if first['findings'] == second['findings'] and first['partition'] == second['partition']:
        individuals = {gsm: f'I{index:03d}' for index, members in enumerate(first['partition'], 1) for gsm in members}
        return {'status': 'SHARED_AGREED', 'findings': first['findings'], 'individuals': individuals}
    return {'status': 'SHARED_DISAGREE', 'findings_by_reader': [first['findings'], second['findings']]}


def candidates(atlas_path, classification_dir, design_audit_dirs, exclude_dirs=()):
    """Ranking candidates the design audit admitted, minus comparisons read in exclude_dirs."""
    frame = design_audit.candidates(atlas_path, classification_dir, exclude_dirs)
    audit = design_audit_table(design_audit_dirs)
    if audit is None or frame.empty:
        return frame.iloc[0:0] if audit is None else frame
    admitted = audit[audit.status == 'RANKABLE'][KEY].assign(_admitted=True)
    frame = frame.merge(admitted, on=KEY, how='left')
    return frame[frame._admitted.notna()].drop(columns='_admitted').reset_index(drop=True)


def build_items(frame, repairs):
    arms = design_audit.arms_from_repairs(repairs)
    series_cache, items = {}, []
    for _, c in frame.iterrows():
        mkey = MODALITY_KEYS[c.modality]
        case, control = arms.get((mkey, int(c.atlas_row_ordinal)), ([], []))
        item = {**{k: c[k] for k in KEY}, 'case_gsms': [str(g) for g in case], 'control_gsms': [str(g) for g in control]}
        folder = GEO / MODALITIES[mkey]['folder']
        if (mkey, c.gse_id) not in series_cache:
            raw = design_estimand.read_frame(folder / 'processed_metadata_preharmonize' / f'{c.gse_id}.csv')
            series_cache[(mkey, c.gse_id)] = (raw, design_estimand.subject_table(raw)[0] if raw is not None else None)
        raw, field = series_cache[(mkey, c.gse_id)]
        if raw is None:
            item['skip'] = 'METADATA_UNAVAILABLE'
        elif field:
            item['skip'] = 'GEO_INDIVIDUAL_FIELD'
            item['subject_field'] = field
        elif len(case) + len(control) > MAX_SAMPLES:
            item['skip'] = 'NOT_READ_ARMS_ABOVE_LIMIT'
        else:
            unit = {'gse': c.gse_id, 'comparison_name': c.comparison_name, 'condition_name': c.condition_name,
                    'case': item['case_gsms'], 'control': item['control_gsms']}
            item['prompt'] = package(unit, raw, series_text(folder, c.gse_id))
            item['prompt_sha256'] = hashlib.sha256(item['prompt'].encode()).hexdigest()
        items.append(item)
    return items


def run(atlas_path, classification_dir, design_audit_dirs, repairs, out_dir, exclude_dirs=(), cache_ledgers=(), workers=8,
        max_requests=4000, dry_run=False):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = candidates(atlas_path, classification_dir, design_audit_dirs, exclude_dirs)
    frame = frame.sort_values(['modality', 'gse_id', 'comparison_name', 'condition_name']).reset_index(drop=True)
    frame.to_csv(out_dir / 'candidates.csv', index=False)
    items = build_items(frame, repairs)
    readable = [i for i in items if 'prompt' in i]
    if dry_run:
        report = {'candidates': len(items), 'to_read': len(readable),
                  'skipped': pd.Series([i['skip'] for i in items if 'skip' in i]).value_counts().to_dict(),
                  'prompt_chars_total': sum(len(i['prompt']) for i in readable)}
        (out_dir / 'dry_run.json').write_text(json.dumps(report, indent=2))
        print(json.dumps(report, indent=1))
        return
    ledger = Ledger(out_dir / 'ledger.jsonl')
    for path in cache_ledgers:  # answers saved by earlier reads of the same prompt (request hash) are copied in
        for digest, record in Ledger(path).cache.items():
            if digest not in ledger.cache:
                ledger.append(record)
    budget = {'max': max_requests, 'sent': 0, 'lock': threading.Lock()}

    def read(item):
        base = {k: item[k] for k in KEY + ['case_gsms', 'control_gsms']}
        if 'skip' in item:
            return {**base, 'status': item['skip'], **({'subject_field': item['subject_field']} if 'subject_field' in item else {}),
                    'receipts': []}
        answers, receipts = [], []
        for model in MODELS:
            record, cached = ledger_call(model, item['prompt'], ledger, budget, parse=parse, max_tokens=MAX_TOKENS)
            receipts.append({'model': model, 'status': record.get('status'), 'request_sha256': record.get('request_sha256'),
                             'cached': cached})
            answers.append(record.get('parsed') if record.get('status') == 'OK' else None)
        kept = [None if a is None else {k: a.get(k) for k in ('material', 'shared_individuals', 'individual_evidence', 'individuals')}
                for a in answers]
        return {**base, **reconcile(item['case_gsms'], item['control_gsms'], answers), 'answers': kept,
                'prompt_sha256': item['prompt_sha256'], 'receipts': receipts}

    with ThreadPoolExecutor(workers) as pool:
        results = list(pool.map(read, items))
    with open(out_dir / 'comparisons.jsonl', 'w') as sink:
        for r in results:
            sink.write(json.dumps(r, default=str) + '\n')
    statuses = pd.Series([r['status'] for r in results]).value_counts().to_dict() if results else {}
    summary = {'schema': SCHEMA, 'candidates': len(results), 'status': statuses, 'requests_sent': budget['sent'], 'models': MODELS,
               'max_samples': MAX_SAMPLES, 'atlas_sha256': sha256(atlas_path),
               'classification_sha256': sha256(Path(classification_dir) / 'comparisons.jsonl'),
               'design_audit': [str(d) for d in design_audit_dirs], 'excluded_already_read': [str(d) for d in exclude_dirs],
               'cache_ledgers': [str(p) for p in cache_ledgers], 'module_sha256': sha256(__file__),
               'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=1, default=str))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--atlas', type=Path, required=True)
    parser.add_argument('--classification', type=Path, required=True)
    parser.add_argument('--design-audit', type=Path, action='append', required=True)
    parser.add_argument('--repairs', action='append', required=True, help='modality=DIR (repairs.jsonl gives the arms)')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--exclude', type=Path, action='append', default=[], help='earlier run whose comparisons are skipped')
    parser.add_argument('--cache', type=Path, action='append', default=[], help='ledger of earlier reads of the same prompts')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--max-requests', type=int, default=4000)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args(argv)
    run(args.atlas, args.classification, args.design_audit, dict(item.split('=', 1) for item in args.repairs), args.out,
        args.exclude, args.cache, args.workers, args.max_requests, args.dry_run)


if __name__ == '__main__':
    main()
