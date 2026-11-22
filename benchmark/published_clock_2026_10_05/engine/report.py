"""Write per-contrast receipts, a flat comparison table and an aggregate report."""
import csv
import json
from pathlib import Path

from .io import clean, write_json


def write_report(out, records, summary, manifest):
    out = Path(out)
    write_json(out / 'summary.json', summary)
    write_json(out / 'manifest.json', manifest)
    try:
        with open(out / 'records.jsonl', 'w') as handle:
            for record in records:
                handle.write(json.dumps(clean(record), allow_nan=False) + '\n')
        systems = sorted({'S1', 'S2', 'S3'} | {s for r in records for s in r.get('systems', {})})
        rows = []
        for record in records:
            row = {k: record.get(k) for k in ('id', 'tier', 'gse', 'modality', 'clock_clockbase', 'measure',
                                             'case_label', 'control_label')}
            row.update({f'published_{k}': record.get(k) for k in ('direction', 'significant', 'p_value', 'threshold')})
            row.update(published_effect=json.dumps(record.get('effect')), source_url=record.get('source', {}).get('url'),
                       groups_status=record['groups']['status'], match_status=record['match']['status'])
            for system in systems:
                value = record.get('systems', {}).get(system, {})
                for key in ('status', 'covered', 'direction', 'significant', 'p', 'cohen_d', 'mean_diff', 'n_case',
                            'n_control', 'row_ordinal', 'swapped', 'measure', 'estimand_mismatch', 'confidence'):
                    row[f'{system}_{key}'] = value.get(key)
            rows.append(row)
        with open(out / 'per-record.csv', 'w', newline='') as handle:
            fields = list(rows[0]) if rows else ['id']
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        lines = ['# Published-clock comparison', '',
                 f"{len(records)} contrasts; bootstrap unit: paper; seed: {summary['bootstrap']['seed']}.", '',
                 'S1 is the matched historical row. S2 uses admitted verified or provisional rows and descriptive ',
                 'results. S3 tests the published groups with the saved ClockBase predictions. S4 entries, when ',
                 'provided, are scored separately for each model.', '',
                 'Direction agreement uses published-significant contrasts with an available direction. ',
                 'Unknown significance thresholds and descriptive contrasts have no significance decision. ',
                 'Coverage counts matched historical rows for S1 and available results for S2/S3. ',
                 'T1 and T2 are reported separately. Intervals are unavailable with fewer than two papers.', '']
        for tier, sections in summary['tiers'].items():
            lines.extend([f'## {tier}', '', '| System | Covered / records | Direction agreement (n) | Significance agreement (n) |',
                          '| --- | --- | --- | --- |'])
            for system, values in sections['overall']['all'].items():
                def fmt(key, denominator):
                    point = values[key]
                    ci = values['intervals'][key]
                    estimate = f'{point:.3f}' if point is not None else 'undefined'
                    if ci['low'] is not None:
                        estimate += f" [{ci['low']:.3f}, {ci['high']:.3f}]"
                    return f'{estimate} ({values[denominator]})'
                lines.append(f"| {system} | {values['covered']} / {values['records']} | "
                             f"{fmt('direction_agreement', 'direction_denominator')} | "
                             f"{fmt('significance_agreement', 'significance_denominator')} |")
            lines.append('')
        lines.extend(['## Limits', '',
                      'Matching a contrast does not make its historical or verified estimand identical to the ',
                      'paper. Atlas age comparisons against published acceleration contrasts carry an ',
                      '`estimand_mismatch` flag. S3 residuals are fit on all series samples with recorded age, ',
                      'which can differ from control-only residuals in a publication. Covariates, pairing, ',
                      'family effects and repeated measurements are not modeled by the Welch test. ',
                      'Published statistics remain quoted source values; the engine does not fill missing values.', '',
                      'All matches, model selections, source hashes, exclusions and residual-fit samples are in ',
                      '`records.jsonl`. Breakdowns by modality, clock family and source, effect correlations, ',
                      'non-significant agreement and optional confidence scores are in `summary.json`.'])
        (out / 'report.md').write_text('\n'.join(lines) + '\n')
    except OSError as error:
        raise OSError(f'Cannot write benchmark report in {out}: {error}') from error
