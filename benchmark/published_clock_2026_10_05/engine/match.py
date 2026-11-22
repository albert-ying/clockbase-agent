"""Match atlas contrasts by sample overlap or independent textual judgments."""
import json

from .stats import direction, threshold


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 0.0


def rank_rows(record, rows, groups, adjudicator=None):
    candidates = []
    resolved = groups.get('status') in ('AGREED', 'MANUAL')
    for row in rows:
        if row.get('clock') != record.get('clock_clockbase') or row.get('gse_id', record['gse']) != record['gse']:
            continue
        item = {'row_ordinal': row['row_ordinal'], 'row': row, 'matched': False, 'swapped': None,
                'jaccard_case': None, 'jaccard_control': None, 'orientations': []}
        if resolved:
            for swapped in (False, True):
                # Atlas rows without recovered samples store None for both arms.
                case, control = (row.get('control_gsms') or [], row.get('case_gsms') or []) if swapped else (
                    row.get('case_gsms') or [], row.get('control_gsms') or [])
                a, b = jaccard(groups['case_gsms'], case), jaccard(groups['control_gsms'], control)
                item['orientations'].append({'swapped': swapped, 'jaccard_case': a, 'jaccard_control': b,
                                             'matched': a >= 0.5 and b >= 0.5})
            best = max(item['orientations'], key=lambda x: (min(x['jaccard_case'], x['jaccard_control']),
                       x['jaccard_case'] + x['jaccard_control'], not x['swapped']))
            item.update(best, method='sample_overlap')
        candidates.append(item)
    if not resolved and candidates and adjudicator is not None:
        descriptions = [{'row_ordinal': c['row_ordinal'], 'comparison': c['row'].get('comparison_name'),
                         'case': c['row'].get('condition_name'), 'clock': c['row']['clock']} for c in candidates]
        prompt = ('Determine which atlas descriptions express the supplied published contrast. '
                  'Judge only the descriptions, without inferring results. Return JSON with candidates: '
                  'a list of objects with row_ordinal and orientation (SAME, SWAPPED, DIFFERENT, UNCLEAR). '
                  'Omitted rows are unresolved. SAME means atlas case is published case; SWAPPED means '
                  'atlas case is published control.\n' + json.dumps(
                      {'gse': record['gse'], 'clock': record['clock_clockbase'],
                       'case_definition': record['case_definition'], 'control_definition': record['control_definition'],
                       'atlas': descriptions}, sort_keys=True))
        answers = adjudicator.ask(prompt)
        judgments = []
        for answer in answers:
            mapping = {}
            entries = (answer.get('parsed') or {}).get('candidates', [])
            if answer.get('status') == 'OK' and isinstance(entries, list):
                for entry in entries:
                    if isinstance(entry, dict) and isinstance(entry.get('row_ordinal'), int):
                        ordinal = entry['row_ordinal']
                        mapping[ordinal] = 'UNCLEAR' if ordinal in mapping else entry.get('orientation')
            judgments.append(mapping)
        for item in candidates:
            orientations = [j.get(item['row_ordinal'], 'UNCLEAR') for j in judgments]
            agreed = len(orientations) == 2 and orientations[0] == orientations[1] and orientations[0] in ('SAME', 'SWAPPED')
            item.update(matched=agreed, swapped=(orientations[0] == 'SWAPPED') if agreed else None,
                        method='text_agreement', judgments=orientations, answers=answers)
    candidates.sort(key=lambda c: (not c['matched'],
                    -min(c['jaccard_case'] or 0, c['jaccard_control'] or 0),
                    -((c['jaccard_case'] or 0) + (c['jaccard_control'] or 0)),
                    c['row'].get('final_class') == 'DUPLICATE', c['row_ordinal']))
    selected = next((c for c in candidates if c['matched']), None)
    return {'status': 'MATCHED' if selected else 'UNRESOLVED',
            'selected_row_ordinal': selected['row_ordinal'] if selected else None, 'candidates': candidates}


def atlas_result(record, candidate, system):
    if candidate is None:
        return {'status': 'UNMATCHED', 'covered': False}
    row, alpha = candidate['row'], threshold(record)
    swapped = candidate['swapped']
    accepted = row.get('final_class') in ('VERIFIED', 'PROVISIONAL', 'DESCRIPTIVE')
    data = (row.get('historical') or {}) if system == 'S1' else (row.get('recomputed') or {})
    p = data.get('pvalue') if system == 'S1' else data.get('p')
    d = data.get('effect_size') if system == 'S1' else data.get('cohen_d')
    diff = data.get('mean_diff')
    # Exported arms follow repairs; historical statistics precede any orientation swap.
    historical_swaps = sum(a.get('action') == 'SWAP_ARMS_CONFIRMED_INVERSION'
                           for a in row.get('actions', []))
    if system == 'S1' and historical_swaps % 2:
        swapped = not swapped
    orient = -1 if swapped else 1
    d, diff = (orient * d if d is not None else None), (orient * diff if diff is not None else None)
    reported_direction = data.get('direction')
    sign = direction(diff if diff is not None else d)
    if sign is None and reported_direction in ('higher', 'lower', 'none'):
        sign = {'higher': 'lower', 'lower': 'higher', 'none': 'none'}[reported_direction] if swapped else reported_direction
    covered = system == 'S1' or accepted
    descriptive = row.get('final_class') == 'DESCRIPTIVE' and system == 'S2'
    return {'status': 'HISTORICAL' if system == 'S1' else row.get('final_class'), 'covered': covered,
            'row_ordinal': row['row_ordinal'], 'swapped': swapped, 'p': p if covered and not descriptive else None,
            'cohen_d': d if covered else None, 'mean_diff': diff if covered else None,
            'direction': sign if covered else None,
            'significant': bool(p < alpha) if covered and not descriptive and p is not None and alpha is not None else None,
            'threshold': alpha, 'n_case': data.get('n_control' if swapped else 'n_case'),
            'n_control': data.get('n_case' if swapped else 'n_control'), 'measure': row.get('measure', 'age'),
            'estimand_mismatch': row.get('measure', 'age') != record['measure'],
            'effect_type': 'historical unspecified' if system == 'S1' else "Cohen's d",
            'effect_definitions_matched': row.get('evidence', {}).get('effect_definitions_matched', []),
            'checks': row.get('checks', {}), 'soft_flags': row.get('soft_flags', []),
            'duplicate_of': row.get('duplicate_of')}
