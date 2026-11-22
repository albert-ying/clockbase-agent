"""Expand each T2 expectation (no clock named) into one record per ClockBase clock of its series.

ComputAgeBench states an expected direction for a condition against healthy controls without naming a clock. Each
expanded record keeps that expectation and the group definitions, names one clock that the verified atlas holds for
the series, and compares age acceleration (residual on chronological age), the quantity ComputAgeBench evaluates.

  python -m benchmark.published_clock_2026_10_05.curated.expand_t2 --ground-truth T2.jsonl --atlas ATLAS.csv --out OUT.jsonl
"""
import argparse
import json
from pathlib import Path

import pandas as pd


def expand(records, atlas):
    clocks = atlas.groupby(atlas.gse_id.str.upper().str.strip())['clock'].agg(lambda s: sorted(set(s)))
    expanded, missing = [], []
    for record in records:
        found = clocks.get(record['gse'].upper(), [])
        if not found:
            missing.append(record['id'])
        for clock in found:
            row = dict(record, id=f"{record['id']}-{clock}", clock_clockbase=clock, measure='acceleration_residual')
            notes = record['extraction'].get('notes') or []
            notes = (notes if isinstance(notes, list) else [notes]) + [
                f'Expanded to ClockBase clock {clock}; measure set to acceleration_residual for the comparison.']
            row['extraction'] = dict(record['extraction'], notes=notes)
            expanded.append(row)
    return expanded, missing


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--ground-truth', type=Path, required=True)
    parser.add_argument('--atlas', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f'refusing to overwrite {args.out}')
    records = [json.loads(line) for line in args.ground_truth.read_text().splitlines()]
    if any(r['tier'] != 'T2' or r['clock_clockbase'] for r in records):
        raise SystemExit('expected T2 records without a ClockBase clock')
    atlas = pd.read_csv(args.atlas, usecols=['gse_id', 'clock'])
    expanded, missing = expand(records, atlas)
    args.out.write_text(''.join(json.dumps(r) + '\n' for r in expanded))
    print(json.dumps({'records': len(records), 'expanded': len(expanded), 'series_without_atlas_clock': missing}, indent=1))


if __name__ == '__main__':
    main()
