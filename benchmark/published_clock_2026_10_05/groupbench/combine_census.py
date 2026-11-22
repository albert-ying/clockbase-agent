"""Combine rule outcomes with the three census reviews (PROTOCOL.md v2, Adjudication).

Per unit: a rule error is overturned only when all three reviewers find no error and give the same reason category
(the unit then becomes correct, or not assessable when the shared reason is an assessability reason); it is confirmed
when at least two reviewers find an error; otherwise unresolved. A rule-correct unit becomes a confirmed error when at
least two reviewers find one; it becomes not assessable when all three give the same assessability reason. Fleiss'
kappa is computed on the reviewers' unit verdicts (error, correct, not assessable).
"""
from collections import Counter

ASSESSABILITY = {'unencoded_factor', 'annotation_defect', 'ambiguous_label', 'unassignable'}


def unit_verdict(review):
    """A reviewer's unit verdict: error if any variant is an error; not assessable if every variant is; else correct."""
    verdicts = [v['verdict'] for v in review['variants']]
    if 'error' in verdicts:
        return 'error'
    if verdicts and all(v == 'not_assessable' for v in verdicts):
        return 'not_assessable'
    return 'correct'


def reason(review):
    reasons = {v.get('not_assessable_reason') or v.get('overturn_reason_category') or ''
               for v in review['variants'] if v['verdict'] != 'error'}
    reasons.discard('')
    return tuple(sorted(reasons))


def combine(rule_class, reviews):
    verdicts = [unit_verdict(r) for r in reviews]
    errors = verdicts.count('error')
    types = Counter(t for r in reviews for v in r['variants'] if v['verdict'] == 'error' for t in v.get('error_types', []))
    shared = {reason(r) for r in reviews}
    if rule_class == 'ERROR':
        if len(reviews) == 3 and errors == 0 and len(shared) == 1 and next(iter(shared)):
            shared_reason = next(iter(shared))
            outcome = 'NOT_ASSESSABLE' if set(shared_reason) & ASSESSABILITY and all(v == 'not_assessable' for v in verdicts) else 'OVERTURNED'
        elif errors >= 2:
            outcome = 'CONFIRMED'
        else:
            outcome = 'UNRESOLVED'
    else:
        if errors >= 2:
            outcome = 'CONFIRMED_MISSED_BY_RULES'
        elif len(reviews) == 3 and all(v == 'not_assessable' for v in verdicts) and len(shared) == 1:
            outcome = 'NOT_ASSESSABLE'
        else:
            outcome = 'CORRECT'
    return {'outcome': outcome, 'reviewer_verdicts': verdicts, 'reviewer_error_types': dict(types),
            'shared_reason': sorted(next(iter(shared))) if len(shared) == 1 else None}


def fleiss_kappa(table):
    """table: list of per-unit category counts (each row sums to the number of raters)."""
    if not table:
        return None
    raters = sum(table[0])
    n = len(table)
    categories = len(table[0])
    p_j = [sum(row[j] for row in table) / (n * raters) for j in range(categories)]
    p_i = [(sum(c * c for c in row) - raters) / (raters * (raters - 1)) for row in table]
    p_bar, p_e = sum(p_i) / n, sum(p * p for p in p_j)
    return (p_bar - p_e) / (1 - p_e) if p_e < 1 else None


def kappa(all_reviews):
    order = ('error', 'correct', 'not_assessable')
    table = []
    for reviews in all_reviews:
        if len(reviews) == 3:
            verdicts = Counter(unit_verdict(r) for r in reviews)
            table.append([verdicts[c] for c in order])
    return fleiss_kappa(table)
