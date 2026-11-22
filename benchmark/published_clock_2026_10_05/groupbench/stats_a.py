"""Error-rate statistics of PROTOCOL.md v2 for clustered binary outcomes (units within series).

- Random-intercept logistic model fitted by maximum likelihood with Gauss-Hermite quadrature; the marginal error rate
  E[expit(b0 + u)] with a profile-likelihood 95% interval (b0 profiled with sigma re-optimized, each bound mapped to
  the marginal rate at its sigma).
- Fallback when there are fewer than 5 errors or the fit does not converge: Clopper-Pearson on the effective sample
  n / DEFF, DEFF = 1 + (m - 1) ICC, m the size-weighted mean cluster size, ICC the ANOVA estimate floored at 0.
- Co-primary: the proportion of series with at least one error (Clopper-Pearson).
- Secondary: unadjusted Clopper-Pearson, series bootstrap (BCa, 10,000 replicates, seed 20261006), leave-one-series-out
  range, exclusion bounds.
"""
import math

import numpy as np
from scipy import optimize, stats

SEED = 20261006
NODES, WEIGHTS = np.polynomial.hermite.hermgauss(30)


def clopper_pearson(k, n, alpha=0.05):
    if n <= 0:
        return (None, None)
    lower = stats.beta.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    upper = stats.beta.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return (float(lower), float(upper))


def icc_anova(clusters):
    """ANOVA intraclass correlation for binary outcomes; clusters is a list of 0/1 lists."""
    clusters = [c for c in clusters if c]
    k, n = len(clusters), sum(len(c) for c in clusters)
    if k < 2 or n <= k:
        return 0.0
    grand = sum(sum(c) for c in clusters) / n
    ssb = sum(len(c) * (np.mean(c) - grand) ** 2 for c in clusters)
    ssw = sum(sum((x - np.mean(c)) ** 2 for x in c) for c in clusters)
    msb, msw = ssb / (k - 1), ssw / (n - k)
    m0 = (n - sum(len(c) ** 2 for c in clusters) / n) / (k - 1)
    denominator = msb + (m0 - 1) * msw
    return float(max(0.0, (msb - msw) / denominator)) if denominator > 0 else 0.0


def effective_interval(clusters, alpha=0.05):
    n = sum(len(c) for c in clusters)
    k = sum(sum(c) for c in clusters)
    icc = icc_anova(clusters)
    m = sum(len(c) ** 2 for c in clusters) / n if n else 0
    deff = 1 + (m - 1) * icc
    n_eff, k_eff = n / deff, k / deff
    lower = stats.beta.ppf(alpha / 2, k_eff, n_eff - k_eff + 1) if k > 0 else 0.0
    upper = stats.beta.ppf(1 - alpha / 2, k_eff + 1, n_eff - k_eff) if k < n else 1.0
    return {'rate': k / n if n else None, 'interval': [float(lower), float(upper)], 'icc': icc, 'deff': deff,
            'n_effective': n_eff, 'size_weighted_mean_cluster': m}


def _loglik(b0, log_sigma, clusters):
    sigma = math.exp(log_sigma)
    total = 0.0
    for c in clusters:
        y, n = sum(c), len(c)
        u = math.sqrt(2) * sigma * NODES
        p = 1 / (1 + np.exp(-(b0 + u)))
        p = np.clip(p, 1e-300, 1 - 1e-16)
        total += math.log(max(np.sum(WEIGHTS * p ** y * (1 - p) ** (n - y)) / math.sqrt(math.pi), 1e-300))
    return total


def marginal_rate(b0, log_sigma):
    sigma = math.exp(log_sigma)
    u = math.sqrt(2) * sigma * NODES
    return float(np.sum(WEIGHTS / (1 + np.exp(-(b0 + u)))) / math.sqrt(math.pi))


def glmm_rate(clusters):
    """Maximum likelihood random-intercept logistic fit; None when it fails or is degenerate."""
    clusters = [c for c in clusters if c]
    k = sum(sum(c) for c in clusters)
    n = sum(len(c) for c in clusters)
    if k == 0 or k == n:
        return None
    objective = lambda x: -_loglik(x[0], x[1], clusters)
    start = [math.log((k + 0.5) / (n - k + 0.5)), math.log(0.5)]
    fit = optimize.minimize(objective, start, method='Nelder-Mead', options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 4000})
    if not fit.success:
        return None
    b0, log_sigma = fit.x

    def profile(b):
        res = optimize.minimize_scalar(lambda s: -_loglik(b, s, clusters), bounds=(-8, 4), method='bounded')
        return -res.fun, res.x

    # The bounded profile can sit marginally below the unconstrained optimum; use the larger as the reference.
    best = max(-fit.fun, profile(b0)[0])
    threshold = best - stats.chi2.ppf(0.95, 1) / 2

    def bound(direction):
        step, b = 0.25, b0
        for _ in range(200):
            b_next = b + direction * step
            value, _ = profile(b_next)
            if value < threshold:
                f_low, f_high = profile(b)[0] - threshold, value - threshold
                if f_low <= 0:  # no sign change: the profile is flat near the optimum
                    return None, None
                root = optimize.brentq(lambda x: profile(x)[0] - threshold, min(b, b_next), max(b, b_next))
                return root, profile(root)[1]
            b = b_next
        return None, None

    try:
        lower_b, lower_s = bound(-1)
        upper_b, upper_s = bound(1)
    except ValueError:
        return None
    interval = [marginal_rate(lower_b, lower_s) if lower_b is not None else None,
                marginal_rate(upper_b, upper_s) if upper_b is not None else None]
    return {'rate': marginal_rate(b0, log_sigma), 'interval': interval, 'b0': float(b0), 'sigma': float(math.exp(log_sigma))}


def bootstrap(clusters, replicates=10000, seed=SEED, alpha=0.05):
    """Series bootstrap of the pooled rate: percentile and BCa intervals."""
    clusters = [c for c in clusters if c]
    if len(clusters) < 2:
        return None
    rng = np.random.default_rng(seed)
    sums = np.array([sum(c) for c in clusters], dtype=float)
    sizes = np.array([len(c) for c in clusters], dtype=float)
    estimate = sums.sum() / sizes.sum()
    draws = rng.integers(0, len(clusters), size=(replicates, len(clusters)))
    rates = sums[draws].sum(axis=1) / sizes[draws].sum(axis=1)
    percentile = [float(np.quantile(rates, alpha / 2)), float(np.quantile(rates, 1 - alpha / 2))]
    proportion = np.mean(rates < estimate)
    if proportion in (0.0, 1.0):
        return {'percentile': percentile, 'bca': None}
    z0 = stats.norm.ppf(proportion)
    jack = np.array([(sums.sum() - sums[i]) / (sizes.sum() - sizes[i]) for i in range(len(clusters))])
    diff = jack.mean() - jack
    acceleration = (diff ** 3).sum() / (6 * ((diff ** 2).sum()) ** 1.5) if (diff ** 2).sum() > 0 else 0.0
    def adjusted(q):
        z = stats.norm.ppf(q)
        return float(np.quantile(rates, stats.norm.cdf(z0 + (z0 + z) / (1 - acceleration * (z0 + z)))))
    return {'percentile': percentile, 'bca': [adjusted(alpha / 2), adjusted(1 - alpha / 2)]}


def estimate(units, key='gse', outcome='error'):
    """All primary and secondary statistics for units: dicts with a cluster key and a 0/1 outcome."""
    by_cluster = {}
    for unit in units:
        by_cluster.setdefault(unit[key], []).append(int(unit[outcome]))
    clusters = list(by_cluster.values())
    n, k = sum(len(c) for c in clusters), sum(sum(c) for c in clusters)
    result = {'units': n, 'errors': k, 'series': len(clusters), 'rate': k / n if n else None,
              'unadjusted_interval': clopper_pearson(k, n), 'effective': effective_interval(clusters) if n else None}
    model = glmm_rate(clusters) if k >= 5 else None
    result['model'] = model
    valid = (model is not None and None not in model['interval']
             and model['interval'][0] < model['rate'] < model['interval'][1]
             and model['interval'][1] - model['interval'][0] > 1e-6)
    result['model_valid'] = valid
    result['primary'] = ({'method': 'random-intercept logistic, profile likelihood', **model} if valid
                         else {'method': 'Clopper-Pearson on effective n (DEFF)', 'rate': result['rate'],
                               'interval': result['effective']['interval'] if result['effective'] else None})
    series_with_error = sum(1 for c in clusters if any(c))
    result['series_with_error'] = {'count': series_with_error, 'of': len(clusters),
                                   'proportion': series_with_error / len(clusters) if clusters else None,
                                   'interval': clopper_pearson(series_with_error, len(clusters))}
    result['bootstrap'] = bootstrap(clusters)
    loso = []
    for name in by_cluster:
        rest = [u for u in units if u[key] != name]
        if rest:
            loso.append(sum(int(u[outcome]) for u in rest) / len(rest))
    result['leave_one_series_out'] = [min(loso), max(loso)] if loso else None
    return result


def exclusion_bounds(errors, assessable, excluded):
    """The rate if every excluded comparison were an error, and if none were."""
    total = assessable + excluded
    return {'all_excluded_errors': (errors + excluded) / total if total else None,
            'no_excluded_errors': errors / total if total else None}
