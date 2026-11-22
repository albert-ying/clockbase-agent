"""Compare two independent extractions of the POSSIBLE papers, then assemble the expanded ground truth.

compare: every record is checked against the cached text (each quote fragment, split at '...', occurs in one text block;
every copied statistic occurs in the quote). Records that pass in both extractions with the same series, clock,
measure, direction, P value, group sizes and effect value are agreed. All other records of a paper go to
review/<paper>.json for a reviewer, who reads the text and returns the records it confirms.

assemble: the earlier ground truth, the agreed records and the reviewed records, with ids assigned per paper.

  python -m benchmark.published_clock_2026_10_05.extraction.reconcile compare --possible DIR --first codex --second codex2
  python -m benchmark.published_clock_2026_10_05.extraction.reconcile collect --possible DIR --reviewer codex-review
  python -m benchmark.published_clock_2026_10_05.extraction.reconcile assemble --possible DIR --base GT.jsonl --out NEW.jsonl
"""
import argparse
import json
import re
from collections import Counter
from decimal import Decimal
from pathlib import Path

from .check_quotes import cached_texts, normalized, numbers
from .possible_inputs import safe_name

RESULT_FIELDS = ('tier', 'gse', 'modality', 'tissue', 'clock_reported', 'clock_clockbase', 'measure', 'case_label',
                 'control_label', 'case_definition', 'control_definition', 'n_case', 'n_control', 'direction',
                 'significant', 'threshold', 'p_value', 'effect', 'test', 'covariates')
MEASURES = {'age', 'acceleration_residual', 'acceleration_difference', 'rate', 'other'}
CLOCKS = {'Horvathv1', 'Hannum', 'PhenoAge', 'GrimAgeV1', 'GrimAgeV2', 'DunedinPACE', 'Horvathv2', 'Lin', 'Zhang_10',
          'YingDamAge', 'YingCausAge', 'YingAdaptAge', 'Hazard_Mouse_All_All_scaleddiff_EN_tAge',
          'Hazard_Rodents_All_All_scaleddiff_EN_tAge', 'Hazard2_All_All_All_scaleddiff_EN_tAge'}


def clockbase_clock(value):
    """The ClockBase clock name of an extracted mapping ('Horvathv2 (skin and blood)' -> 'Horvathv2'), else None."""
    if not value:
        return None
    name = value.split(' (')[0].strip()
    return name if name in CLOCKS else None


def read_output(path):
    """Return the records of one extraction output, or None when it is missing or not one JSON object."""
    try:
        text = Path(path).read_text().strip()
        text = re.sub(r'^```(?:json)?\s*|\s*```$', '', text)
        value = json.loads(text)
        return value if isinstance(value, dict) and isinstance(value.get('records'), list) else None
    except (OSError, ValueError):
        return None


TYPOGRAPHY = str.maketrans({'\u2018': "'", '\u2019': "'", '\u201c': '"', '\u201d': '"', '\u2032': "'",
                            '\u2013': '-', '\u2014': '-', '\u2212': '-', '\u2010': '-', '\u2011': '-',
                            '\u00a0': ' ', '\u2009': ' ', '\u202f': ' ', '\u00ad': None, '\u200b': None})


def canonical(text):
    """Whitespace, quote and dash variants compared as equal."""
    return normalized(text.translate(TYPOGRAPHY))


MIN_PIECE = 8


def verbatim_pieces(fragment, document):
    """Cover the fragment greedily with the longest pieces that occur verbatim in the document; None if a piece shorter
    than MIN_PIECE characters (other than trailing punctuation or space) is needed."""
    pieces, rest = [], fragment
    while rest.strip(' .,;:'):
        low, high = 0, len(rest)
        while low < high:
            middle = (low + high + 1) // 2
            if rest[:middle] in document:
                low = middle
            else:
                high = middle - 1
        if low < MIN_PIECE and low < len(rest):
            return None
        pieces.append(rest[:low])
        rest = rest[low:]
    return pieces


def problems(record, texts):
    """List the reasons a record fails the schema, quote or statistic checks (empty when it passes).

    Extractors join non-adjacent verbatim text (sentences, table cells and their footnotes, caption panels), so a quote
    passes when it is covered by verbatim pieces of the document (the blocks joined in document order) of at least
    MIN_PIECE characters each."""
    document = ' '.join(canonical(text) for text in texts)
    found = []
    if record.get('direction') not in ('higher', 'lower', 'none'):
        found.append('direction')
    if record.get('measure') not in MEASURES:
        found.append('measure')
    quote = (record.get('source') or {}).get('quote') or ''
    fragments = [canonical(f) for f in re.split(r'\.\.\.|…|(?<=[.;!?])\s+(?=[A-Z(\[])', quote) if canonical(f)]
    if not fragments:
        found.append('no quote')
    for fragment in fragments:
        if fragment not in document and verbatim_pieces(fragment, document) is None:
            found.append(f'quote fragment absent: {fragment[:80]}')
    values = numbers(quote)
    effect = record.get('effect') or {}
    copied = [record.get('n_case'), record.get('n_control'), record.get('p_value'), effect.get('value')]
    for value in copied:
        if value is None:
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            found.append(f'non-numeric statistic {value!r}')
        elif Decimal(str(value)) not in values and Decimal(str(abs(value))) not in values:
            found.append(f'statistic {value} absent from quote')
    return found


def key(record):
    clock = clockbase_clock(record.get('clock_clockbase')) or re.sub(r'[^a-z0-9]', '', (record.get('clock_reported') or '').lower())
    effect = (record.get('effect') or {}).get('value')
    return (record.get('gse'), clock, record.get('measure'), record.get('direction'), record.get('p_value'),
            record.get('n_case'), record.get('n_control'), effect)


def compare(possible, first, second):
    possible = Path(possible)
    agreed, summary = [], Counter()
    (possible / 'review').mkdir(exist_ok=True)
    for stale in (possible / 'review').glob('*.json'):  # this step owns review/; a rerun replaces it
        stale.unlink()
    for input_path in sorted((possible / 'inputs').glob('*.json')):
        name, paper = input_path.stem, json.loads(input_path.read_text())
        outputs = {label: read_output(possible / label / f'{name}.json') for label in (first, second)}
        records = {label: (out or {}).get('records', []) for label, out in outputs.items()}
        if not any(records.values()):
            summary['papers_without_records'] += 1
            continue
        texts = cached_texts(Path(paper['paper']['text_path']))
        checked = {label: [(r, problems(r, texts)) for r in rows] for label, rows in records.items()}
        remaining = {label: list(rows) for label, rows in checked.items()}
        paper_agreed = []
        for record, issues in checked[first]:
            if issues:
                continue
            match = next((pair for pair in remaining[second] if not pair[1] and key(pair[0]) == key(record)), None)
            if match is None:
                continue
            remaining[second].remove(match)
            remaining[first].remove((record, issues))
            paper_agreed.append(record)
        for record in paper_agreed:
            agreed.append({'paper': paper['paper'], 'record': record, 'extractors': [first, second]})
        summary['agreed_records'] += len(paper_agreed)
        if any(remaining.values()):
            review = {'paper': paper['paper'], 'series': paper['series'], 'agreed': paper_agreed,
                      'unmatched': {label: [{'record': r, 'problems': p} for r, p in rows]
                                    for label, rows in remaining.items()}}
            (possible / 'review' / f'{name}.json').write_text(json.dumps(review, indent=1) + '\n')
            summary['papers_for_review'] += 1
            summary['unmatched_records'] += sum(map(len, remaining.values()))
    with open(possible / 'agreed.jsonl', 'w') as handle:
        for row in agreed:
            handle.write(json.dumps(row) + '\n')
    summary = dict(summary, first=first, second=second)
    (possible / 'compare_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


def collect(possible, reviewer):
    """Check each reviewer record against the cached text; write reviewed/<paper>.json with the records that pass."""
    possible = Path(possible)
    (possible / 'reviewed').mkdir(exist_ok=True)
    summary, dropped = Counter(), []
    for review_path in sorted((possible / 'review').glob('*.json')):
        paper = json.loads(review_path.read_text())['paper']
        output = read_output(possible / reviewer / review_path.name)
        if output is None:
            summary['papers_without_review'] += 1
            dropped.append({'paper': review_path.stem, 'reason': 'no readable reviewer output'})
            continue
        texts, kept = cached_texts(Path(paper['text_path'])), []
        for record in output['records']:
            issues = problems(record, texts)
            if record.get('agreement') not in ('resolved', 'single') or not isinstance(record.get('extractors'), list):
                issues.append('agreement or extractors missing')
            if issues:
                dropped.append({'paper': review_path.stem, 'gse': record.get('gse'), 'problems': issues})
                summary['records_failing_checks'] += 1
            else:
                kept.append(record)
                summary[f"records_{record['agreement']}"] += 1
        summary['records_rejected_by_reviewer'] += len(output.get('rejected', []))
        (possible / 'reviewed' / review_path.name).write_text(json.dumps({'records': kept}, indent=1) + '\n')
    with open(possible / 'review_dropped.jsonl', 'w') as handle:
        for row in dropped:
            handle.write(json.dumps(row) + '\n')
    summary = dict(summary, reviewer=reviewer)
    (possible / 'collect_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


def ground_truth_record(paper, record, extraction):
    ids = paper['ids']
    source = {'paper_id': paper['paper_id'], **{k: ids.get(k) for k in ('id', 'source', 'pmid', 'pmcid', 'doi')},
              'citation': paper['title'], 'year': paper['year'],
              'location': (record.get('source') or {}).get('location'), 'quote': record['source']['quote'],
              'text_path': paper['text_path']}
    row = {field: record.get(field) for field in RESULT_FIELDS}
    row['tier'] = 'T1'
    row['clock_clockbase'] = clockbase_clock(row['clock_clockbase'])
    row['covariates'] = row['covariates'] or []
    return {**row, 'source': source, 'extraction': extraction}


def assemble(possible, base, out):
    """Write base + agreed + reviewed records; reviewed/<paper>.json holds {"records": [...]} with extraction notes."""
    possible = Path(possible)
    rows = [json.loads(line) for line in Path(base).read_text().splitlines()]
    ids = {row['id'] for row in rows}
    in_base = {((row.get('source') or {}).get('paper_id'), row.get('gse')) for row in rows}
    added, skipped = Counter(), Counter()
    per_paper = {}

    def add(paper, record, extraction):
        listed = {s['gse']: s['modality'] for s in json.loads(
            (possible / 'inputs' / f"{safe_name(paper['paper_id'])}.json").read_text())['series']}
        if record.get('gse') not in listed:
            skipped['series_not_listed_for_paper'] += 1
        elif record.get('modality') != listed[record['gse']]:
            skipped['modality_differs_from_atlas_series'] += 1
        elif (paper['paper_id'], record['gse']) in in_base:
            skipped['paper_series_already_in_base'] += 1
        else:
            per_paper.setdefault(paper['paper_id'], []).append(ground_truth_record(paper, record, extraction))
            added[extraction['agreement']] += 1

    for line in (possible / 'agreed.jsonl').read_text().splitlines():
        item = json.loads(line)
        add(item['paper'], item['record'], {'extractors': item['extractors'], 'agreement': 'agreed', 'notes': ''})
    for path in sorted((possible / 'reviewed').glob('*.json')):
        paper = json.loads((possible / 'inputs' / path.name).read_text())['paper']
        for record in json.loads(path.read_text())['records']:
            add(paper, record, {'extractors': record.pop('extractors'), 'agreement': record.pop('agreement'),
                                'notes': record.pop('notes', '')})
    for paper_id, records in sorted(per_paper.items()):
        prefix = paper_id.replace(':', '-')
        for n, row in enumerate(records):
            row['id'] = f"{prefix}-{row['gse']}-{row['clock_clockbase'] or 'other'}-{n}"
            if row['id'] in ids:
                raise ValueError(f"Duplicate ground-truth id: {row['id']}")
            ids.add(row['id'])
            rows.append(row)
    with open(out, 'w') as handle:
        for row in rows:
            handle.write(json.dumps(row) + '\n')
    return {'records': len(rows), 'added': dict(added), 'skipped': dict(skipped), 'papers_added': len(per_paper)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    commands = parser.add_subparsers(dest='command', required=True)
    command = commands.add_parser('compare')
    command.add_argument('--possible', type=Path, required=True)
    command.add_argument('--first', required=True)
    command.add_argument('--second', required=True)
    command = commands.add_parser('collect')
    command.add_argument('--possible', type=Path, required=True)
    command.add_argument('--reviewer', required=True)
    command = commands.add_parser('assemble')
    command.add_argument('--possible', type=Path, required=True)
    command.add_argument('--base', type=Path, required=True)
    command.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'compare':
        result = compare(args.possible, args.first, args.second)
    elif args.command == 'collect':
        result = collect(args.possible, args.reviewer)
    else:
        if args.out.exists():
            raise SystemExit(f'refusing to overwrite {args.out}')
        result = assemble(args.possible, args.base, args.out)
    print(json.dumps(result, indent=1))


if __name__ == '__main__':
    main()
