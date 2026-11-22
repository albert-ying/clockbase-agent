"""Write one extraction input per paper for the POSSIBLE paper-series pairs of the retrieval.

Each input lists the paper, its cached open-access text and the POSSIBLE series with their atlas clocks. Pairs without
cached full text are not extracted; they are listed in skipped.jsonl with the retrieval's full-text status.

  python -m benchmark.published_clock_2026_10_05.extraction.possible_inputs --root DATA --out DATA/extraction/possible
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path


def safe_name(paper_id):
    return paper_id.replace(':', '_').replace('/', '_')


def build(candidates, out):
    out, papers, skipped = Path(out), defaultdict(list), []
    with open(candidates) as handle:
        for line in handle:
            pair = json.loads(line)
            if pair['triage_class'] != 'POSSIBLE':
                continue
            if pair['fulltext_status'] != 'ok' or not pair.get('text_path'):
                skipped.append({'paper_id': pair['paper_id'], 'gse': pair['gse_id'],
                                'reason': f"no cached full text ({pair['fulltext_status']})"})
                continue
            papers[pair['paper_id']].append(pair)
    (out / 'inputs').mkdir(parents=True, exist_ok=True)
    for paper_id, pairs in sorted(papers.items()):
        first = pairs[0]
        record = {'paper': {'paper_id': paper_id, 'ids': first['ids'], 'title': first['title'], 'year': first['year'],
                            'text_path': first['text_path']},
                  'series': [{'gse': p['gse_id'], 'modality': p['modality'], 'atlas_clocks': p['atlas_clocks']}
                             for p in sorted(pairs, key=lambda p: p['gse_id'])]}
        (out / 'inputs' / f'{safe_name(paper_id)}.json').write_text(json.dumps(record, indent=1) + '\n')
    with open(out / 'skipped.jsonl', 'w') as handle:
        for row in skipped:
            handle.write(json.dumps(row) + '\n')
    summary = {'papers': len(papers), 'pairs': sum(map(len, papers.values())), 'skipped_pairs': len(skipped)}
    (out / 'inputs_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True, help='benchmark data directory')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.root / 'retrieval/candidates.jsonl', args.out), indent=1))


if __name__ == '__main__':
    main()
