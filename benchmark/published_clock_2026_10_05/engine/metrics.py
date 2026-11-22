"""Coverage and agreement, with seeded bootstrap intervals resampling whole papers."""
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr

from .io import clean

METRICS = ('coverage', 'direction_agreement', 'significance_agreement', 'nonsignificant_agreement',
           'effect_spearman', 'direction_accuracy', 'confidence_brier', 'confidence_accuracy_gap')


def paper_key(record):
    source = record.get('source') or {}
    return str(next((source[k] for k in ('doi', 'pmid', 'pmcid', 'citation', 'url') if source.get(k)), record['id']))


def source_key(record):
    source = record.get('source') or {}
    return str(record.get('source_family') or source.get('collection') or source.get('citation') or 'unspecified')


def clock_family(record):
    if record.get('clock_family'):
        return record['clock_family']
    name = (record.get('clock_clockbase') or record.get('clock_reported') or '').lower()
    for family in ('horvath', 'hannum', 'grimage', 'phenoage', 'dunedin', 'tage', 'ying'):
        if family in name:
            return family
    return record.get('clock_clockbase') or 'unmapped'


def effect_pair(record, result, system):
    """Only compare the same estimand and effect definition; do not infer missing SDs."""
    effect = record.get('effect') or {}
    if effect.get('value') is None or result.get('measure') != record.get('measure'):
        return None
    unit, kind = str(effect.get('unit', '')).lower(), str(effect.get('type', '')).lower()
    if unit in ('sd', 'standard deviations') and kind in ("cohen's d", 'cohen d', 'cohens d', 'cohen_d'):
        if system == 'S1' and 'cohen_d_pooled' not in result.get('effect_definitions_matched', []):
            return None
        value = result.get('cohen_d')
        key = ('cohen_d', record['measure'], 'SD')
    elif kind in ('mean difference', 'mean_diff') and unit in ('years', 'year', 'rate', 'years/year'):
        value = result.get('mean_diff')
        key = ('mean_diff', record['measure'], unit.replace('year', 'years') if unit == 'year' else unit)
    else:
        return None
    if value is None or not np.isfinite(value) or not np.isfinite(effect['value']):
        return None
    return key, float(effect['value']), float(value)


def point_metrics(items, system):
    available = [(r, r.get('systems', {}).get(system, {})) for r in items]
    covered = [(r, s) for r, s in available if s.get('covered') is True]
    direction_pairs = [(r['direction'], s['direction']) for r, s in covered
                       if r.get('significant') is True and s.get('direction') in ('higher', 'lower', 'none')]
    significance_pairs = [(r['significant'], s['significant']) for r, s in covered
                          if isinstance(r.get('significant'), bool) and isinstance(s.get('significant'), bool)]
    nonsignificant = [(a, b) for a, b in significance_pairs if a is False]
    all_directions = [(r.get('direction'), s['direction']) for r, s in covered
                      if r.get('direction') in ('higher', 'lower', 'none') and s.get('direction') in ('higher', 'lower', 'none')]
    buckets = defaultdict(list)
    for record, result in covered:
        pair = effect_pair(record, result, system)
        if pair is not None:
            key, published, predicted = pair
            buckets[key].append((published, predicted))
    correlations = []
    for key, pairs in sorted(buckets.items()):
        x, y = np.asarray(pairs).T
        rho = float(spearmanr(x, y).statistic) if len(x) >= 3 and len(set(x)) > 1 and len(set(y)) > 1 else None
        correlations.append({'unit_group': list(key), 'n': len(pairs), 'rho': rho})
    # Mixing mean differences and SDs would change ranks; overall rho exists only for one unit group.
    rho = correlations[0]['rho'] if len(correlations) == 1 else None
    confidence = [(float(s['confidence']), float(r.get('direction') == s.get('direction'))) for r, s in covered
                  if isinstance(s.get('confidence'), (int, float)) and not isinstance(s.get('confidence'), bool)
                  and 0 <= s['confidence'] <= 1 and r.get('direction') in ('higher', 'lower', 'none')
                  and s.get('direction') in ('higher', 'lower', 'none')]
    agreement = lambda pairs: float(np.mean([a == b for a, b in pairs])) if pairs else None
    return {'records': len(items), 'papers': len({paper_key(r) for r in items}), 'covered': len(covered),
            'coverage': len(covered) / len(items) if items else None,
            'direction_agreement': agreement(direction_pairs), 'direction_denominator': len(direction_pairs),
            'significance_agreement': agreement(significance_pairs), 'significance_denominator': len(significance_pairs),
            'nonsignificant_agreement': agreement(nonsignificant), 'nonsignificant_denominator': len(nonsignificant),
            'direction_accuracy': agreement(all_directions), 'direction_accuracy_denominator': len(all_directions),
            'effect_spearman': rho, 'effect_groups': correlations, 'effect_denominator': sum(c['n'] for c in correlations),
            'confidence_brier': float(np.mean([(c - correct) ** 2 for c, correct in confidence])) if confidence else None,
            'confidence_accuracy_gap': float(np.mean([c - correct for c, correct in confidence])) if confidence else None,
            'confidence_denominator': len(confidence)}


def summarize(items, system, bootstrap=1000, seed=20261005):
    point = point_metrics(items, system)
    papers = defaultdict(list)
    for record in items:
        papers[paper_key(record)].append(record)
    keys, estimates = sorted(papers), {m: [] for m in METRICS}
    effect_estimates = {tuple(group['unit_group']): [] for group in point['effect_groups']}
    rng = np.random.default_rng(seed)
    # A single paper cannot estimate between-paper uncertainty.
    if len(keys) >= 2:
        for _ in range(bootstrap):
            sampled = [r for i in rng.integers(0, len(keys), len(keys)) for r in papers[keys[i]]]
            value = point_metrics(sampled, system)
            for metric in METRICS:
                if value[metric] is not None and np.isfinite(value[metric]):
                    estimates[metric].append(value[metric])
            for group in value['effect_groups']:
                key = tuple(group['unit_group'])
                if group['rho'] is not None:
                    effect_estimates[key].append(group['rho'])
    point['intervals'] = {m: {'low': float(np.quantile(v, 0.025)) if v else None,
                               'high': float(np.quantile(v, 0.975)) if v else None,
                               'valid_resamples': len(v)} for m, v in estimates.items()}
    for group in point['effect_groups']:
        values = effect_estimates[tuple(group['unit_group'])]
        group['interval'] = {'low': float(np.quantile(values, 0.025)) if values else None,
                             'high': float(np.quantile(values, 0.975)) if values else None,
                             'valid_resamples': len(values)}
    return clean(point)


def metrics(records, bootstrap=1000, seed=20261005):
    systems = sorted({'S1', 'S2', 'S3'} | {s for r in records for s in r.get('systems', {})})
    output = {'bootstrap': {'unit': 'paper', 'replicates': bootstrap, 'seed': seed, 'level': 0.95}, 'tiers': {}}
    for tier in sorted({r.get('tier', 'T1') for r in records}):
        items = [r for r in records if r.get('tier', 'T1') == tier]
        subsets = {'overall': {'all': items}}
        for dimension, key in [('modality', lambda r: r['modality']), ('clock_family', clock_family), ('source', source_key)]:
            values = sorted({key(r) for r in items})
            subsets[dimension] = {v: [r for r in items if key(r) == v] for v in values}
        post_cutoff = [r for r in items if r.get('post_cutoff') is True]
        if post_cutoff:
            subsets['post_cutoff'] = {'known': post_cutoff}
        output['tiers'][tier] = {dim: {value: {s: summarize(subset, s, bootstrap, seed) for s in systems}
                                                    for value, subset in values.items()} for dim, values in subsets.items()}
    return output
