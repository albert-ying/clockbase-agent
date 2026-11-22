"""Intervention-level random-effects meta-analysis over verifier-v3 comparisons.

Unit of evidence: one comparison x clock with a recomputed Cohen's d (pooled SD, case minus control)
and arm sizes. Each is converted to Hedges' g with its large-sample variance. Within a study, all
effects for the same intervention (several clocks, tissues or time points) are averaged into one
study effect whose variance is the mean of their variances; this treats the effects as perfectly
correlated and so never counts one study as several. Studies are then pooled by DerSimonian-Laird
random effects, separately per species.

Interventions come from the dual-model classification (intervention_classification_2026_09_29.py):
only comparisons both models called INTERVENTION with the same canonical target, class and direction
of manipulation enter. All usable comparisons enter regardless of significance.

Independent projects (added 2026-10-14). Several GEO series can hold one experiment: SubSeries of one
SuperSeries, fractions or tissues of the same animals, or series of one paper. Counted as separate studies,
they make one experiment look replicated. With --projects, the unit of replication is the project: series
linked, directly or through other series, by a shared SuperSeries (SubSeries of the same series, or one the
SuperSeries of the other), a shared PubMed ID, or the same contact name and institute in their public GEO
records (`geo` subcommand). Within a project, effects are averaged as within a study above. Series whose GEO
record could not be read stay their own project and are counted.

Correction of 2026-10-17: exclude tumour material and cultured samples, as the ranking does. The ranking admits
only comparisons that a dual-model design audit reads as profiled directly from organisms with no arm of tumour
material; this analysis filtered on the per-sample context label alone, so its Methods claim of samples profiled
directly from organisms did not hold. The audit was therefore run on this analysis's own comparisons
(design_audit_2026_09_30.py --meta; 8,673 of them had not been read for the ranking), and a comparison is dropped
when both models call it tumour material or call its samples a cell line, an ex vivo culture or stem-cell derived.
Only agreement drops a comparison; where the two models disagree it stays and is flagged. It was noticed on
BLZ945, the largest pooled decrease, whose two projects both profiled tumour-bearing mice. --keep-tumour
reproduces the earlier runs.

Correction of 2026-10-15: pool within one atlas category. An intervention's comparisons were pooled across the
four atlas categories (genetic, drug, environmental, disease), so a comparison the atlas does not call a drug
contrast could be averaged into a drug. It was noticed on PLX5622, whose five projects include two that contrast
microglial repopulation after the drug was withdrawn rather than drug exposure: GSE146526 (`Saline_MG-REP vs
Saline_CTRL`, `MIA_MG-REP vs MIA_CTRL`, condition `Saline`/`MIA`) and GSE94042 (`Aged microglia repopulated vs
Aged control`, condition `Microglia repopulation`). The classifier mapped both to plx5622 because PLX5622 produced
the depletion, and the atlas category of all three comparisons is Environment or Other, not Drug. Grouping within
the category is mechanical, applies to every intervention, and is blind to the outcome: here it keeps the ten
drug-exposure comparisons, including one whose effect is positive, and separates the three others. Library-wide it
splits interventions that span categories (3,719 groups instead of 3,547), lowers the replicated count (805 from
828) because some interventions no longer have two projects in one category, and leaves every intervention that
was already robust unchanged. `--pool-across-categories` reproduces the 2026-10-14 grouping.

Specification fixed on 2026-10-14 before this module was run on the atlas (the motivating question, PLX5622,
had been explored on one clock and the drug category; see the Methods): rows of the final run supporting
inference (verified or provisional), excluding GrimAge2 (wrong age input; not ranked); samples profiled directly
from organisms (context in_vivo, as for the ranking); projects as above; tiers and the robust set as below
(replicated: at least two projects; robust: BH-FDR < 0.05 over interventions and the 95% interval excludes
zero); ranked by the bound nearest zero. Sensitivities: series as the unit, all contexts, projects without the
contact link, and the Hartung-Knapp interval (t quantile on k - 1 degrees of freedom), reported beside the
normal interval.

  python -m audit_system.intervention_meta_2026_09_29 geo --atlas atlas_v3_combined.csv \
      --classification DIR --cache GEO_DIR --out projects.json [--offline]
  python -m audit_system.intervention_meta_2026_09_29 run --atlas atlas_v3_combined.csv \
      --classification DIR --out DIR [--projects projects.json] [--context in_vivo]

Lifespan databases (added 2026-10-14, rules fixed before the run; `concordance` subcommand). Mouse interventions of
the primary analysis are matched to GenAge model-organism genes (Mus musculus, longevity influence Pro- or
Anti-Longevity; gene interventions by symbol, case-insensitive) and to DrugAge compounds (Mus musculus; exposure
interventions by compound name, case-insensitive). A compound is lifespan-extending when DrugAge lists a significant
positive average lifespan change and no significant negative one, lifespan-shortening in the reverse case, and is
otherwise not used. Expected direction of the pooled effect: lower clock outputs for loss of an anti-longevity gene,
gain of a pro-longevity gene or a lifespan-extending compound; higher outputs for the reverse. Reported: the share of
interventions whose pooled effect has the expected sign, for all matched interventions, those with a 95% interval
excluding zero, and those with at least two projects, each with a two-sided exact binomial test against one half.

  python -m audit_system.intervention_meta_2026_09_29 concordance --meta DIR/intervention_meta.csv \
      --genage genage_models.csv --drugage drugage.csv --out DIR
"""
import argparse
import concurrent.futures
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from audit_system.membership_recovery_2026_09_27 import sha256

SPECIES = {'human DNAm': 'human', 'human RNA-seq': 'human', 'mouse RNA-seq': 'mouse'}
Z = stats.norm.ppf(0.975)
RELEVANCE_THRESHOLD = 5
EXCLUDED_CLOCKS = ('GrimAgeV2',)  # the historical pipeline supplied the Horvath prediction as its age input
GROUP_KEYS = ['species', 'canonical_target', 'intervention_class', 'manipulation']
GEO_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}&targ=self&form=text&view=brief'
ACCESSION = re.compile(r'^GSE\d+$')


def hedges_g(d, n1, n2):
    """Small-sample corrected standardized mean difference and its variance (Hedges & Olkin 1985)."""
    df = n1 + n2 - 2
    j = 1 - 3 / (4 * df - 1)
    g = j * d
    var = j ** 2 * ((n1 + n2) / (n1 * n2) + d ** 2 / (2 * (n1 + n2)))
    return g, var


def dersimonian_laird(effects, variances):
    """Random-effects pooled estimate; returns dict with estimate, se, CI, tau2, I2, Q, p."""
    y, v = np.asarray(effects, float), np.asarray(variances, float)
    k = len(y)
    w = 1 / v
    fixed = np.sum(w * y) / np.sum(w)
    q = float(np.sum(w * (y - fixed) ** 2))
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    w_star = 1 / (v + tau2)
    estimate = float(np.sum(w_star * y) / np.sum(w_star))
    se = float(np.sqrt(1 / np.sum(w_star)))
    i2 = max(0.0, (q - (k - 1)) / q) if k > 1 and q > 0 else 0.0
    p = float(2 * stats.norm.sf(abs(estimate / se)))
    # Hartung-Knapp: weighted residual variance with a t quantile; undefined for a single study.
    if k > 1:
        se_hk = float(np.sqrt(np.sum(w_star * (y - estimate) ** 2) / ((k - 1) * np.sum(w_star))))
        half_hk = float(stats.t.ppf(0.975, k - 1) * se_hk)
        hk = (estimate - half_hk, estimate + half_hk)
    else:
        hk = (np.nan, np.nan)
    return {'k': k, 'estimate': estimate, 'se': se, 'ci_low': estimate - Z * se, 'ci_high': estimate + Z * se,
            'ci_low_hk': hk[0], 'ci_high_hk': hk[1],
            'tau2': float(tau2), 'i2': float(i2), 'q': q, 'p': p,
            'same_sign_fraction': float(np.mean(np.sign(y) == np.sign(estimate)))}


def study_effects(rows, unit='gse_id'):
    """One effect per study (or project): mean g and mean variance across its comparisons and clocks."""
    out = []
    for key, group in rows.groupby(unit):
        out.append({unit: key, 'g': float(group.g.mean()), 'var': float(group.var_g.mean()),
                    'series': ';'.join(sorted(group.gse_id.unique())),
                    'comparisons': int(group[['gse_id', 'comparison_name', 'condition_name']].drop_duplicates().shape[0]),
                    'effects': int(len(group)), 'n_total': int((group.v3_n_case + group.v3_n_control).max())})
    return pd.DataFrame(out)


def parse_geo_record(text):
    """Links of one GEO series record (brief SOFT text): SuperSeries relations, PubMed IDs, contact and institute."""
    def values(field):
        return [v.strip() for v in re.findall(rf'^!Series_{field} = (.*)$', text, flags=re.M)]
    related = []
    for relation in values('relation'):
        match = re.match(r'(SubSeries of|SuperSeries of): (GSE\d+)', relation)
        if match:
            related.append(match.group(2))
    contact = (values('contact_name') or [''])[0]
    parts = [part.strip().lower() for part in contact.split(',') if part.strip()]
    institute = (values('contact_institute') or [''])[0].lower()
    return {'accession': (values('geo_accession') or [''])[0], 'title': (values('title') or [''])[0],
            'related': sorted(set(related)), 'pubmed': sorted(set(values('pubmed_id'))),
            'contact': f'{parts[0]} {parts[-1]} | {institute}' if parts and institute else ''}


def projects_from_records(records, link_contact=True):
    """Group series into projects: linked by a shared SuperSeries, a shared PubMed ID or (optionally) the same contact
    name and institute, transitively. `records` maps an accession to parse_geo_record() output or None (not read).
    Returns accession -> project ID (the smallest member accession) and, per project, the kinds of link that joined it."""
    parent = {s: s for s in records}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def number(accession):
        return int(accession[3:])
    owner, links = {}, []
    for series in sorted(records, key=number):
        record = records[series]
        keys = {('series', series)}
        if record:
            keys |= {('series', s) for s in record['related']} | {('pubmed', p) for p in record['pubmed']}
            if link_contact and record['contact']:
                keys.add(('contact', record['contact']))
        for key in sorted(keys):
            if key in owner and find(owner[key]) != find(series):
                links.append((series, key[0]))
                a, b = find(owner[key]), find(series)
                parent[max(a, b, key=number)] = min(a, b, key=number)
            owner.setdefault(key, series)
    project = {s: find(s) for s in records}
    kinds = {}
    for series, kind in links:
        kinds.setdefault(project[series], set()).add(kind)
    return project, {p: sorted(k) for p, k in kinds.items()}


def fetch_geo_record(series, cache, offline=False):
    """The public GEO record of one series, cached as text; None when it cannot be read."""
    path = Path(cache) / f'{series}.txt'
    if path.exists() and '!Series_title' in path.read_text(errors='replace'):
        return path.read_text(errors='replace')
    if offline:
        return None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(GEO_URL.format(series), timeout=40) as response:
                text = response.read().decode('utf-8', 'replace')
            if '!Series_title' in text:
                path.write_text(text)
                return text
        except Exception:  # network errors are retried, then the series is reported as not read
            pass
        time.sleep(2 + 2 * attempt)
    return None


def usable_effects(atlas_path, classification_dir, context=None):
    """Comparison x clock effects that enter the meta-analysis, with their intervention keys and Hedges' g."""
    atlas = pd.read_csv(atlas_path, low_memory=False)
    usable = atlas[atlas.v3_final_class.isin(['VERIFIED', 'PROVISIONAL'])].copy()
    usable = usable[(usable.v3_n_case >= 2) & (usable.v3_n_control >= 2) & usable.v3_cohen_d.notna()]
    usable = usable[~usable.clock.isin(EXCLUDED_CLOCKS)]
    if context:
        usable = usable[usable.v3_context == context]  # e.g. in vivo effects only, as for the ranking
    classes = pd.read_json(Path(classification_dir) / 'comparisons.jsonl', lines=True)
    agreed = classes[(classes.status == 'AGREED') & (classes.contrast_type == 'INTERVENTION')]
    key = ['modality', 'gse_id', 'comparison_name', 'condition_name']
    usable = usable.merge(agreed[key + ['canonical_target', 'intervention_class', 'manipulation', 'relevance_score']],
                          on=key, how='inner')
    usable['species'] = usable.modality.map(SPECIES)
    g, var = hedges_g(usable.v3_cohen_d.to_numpy(float), usable.v3_n_case.to_numpy(float), usable.v3_n_control.to_numpy(float))
    usable['g'], usable['var_g'] = g, var
    return usable


def geo(atlas_path, classification_dir, cache, out_path, offline=False, workers=3):
    """Read the GEO records of every series in an intervention with two or more series and write the projects."""
    cache, out_path = Path(cache), Path(out_path)
    cache.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        raise SystemExit(f'{out_path} exists')
    usable = usable_effects(atlas_path, classification_dir)
    counts = usable.groupby(GROUP_KEYS).gse_id.nunique()
    multi = usable.set_index(GROUP_KEYS).loc[counts[counts >= 2].index]
    series = sorted({s for s in multi.gse_id.unique() if ACCESSION.match(str(s))}, key=lambda s: int(s[3:]))
    with concurrent.futures.ThreadPoolExecutor(1 if offline else workers) as pool:
        texts = dict(zip(series, pool.map(lambda s: fetch_geo_record(s, cache, offline), series)))
    records = {s: (parse_geo_record(t) if t else None) for s, t in texts.items()}
    with_contact, contact_kinds = projects_from_records(records, link_contact=True)
    without_contact, _ = projects_from_records(records, link_contact=False)
    out = {'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'source': GEO_URL.format('GSE...'),
           'rule': 'shared SuperSeries, shared PubMed ID, or same contact name and institute (transitive)',
           'series': len(series), 'not_read': sorted(s for s, r in records.items() if r is None),
           'projects': len(set(with_contact.values())), 'projects_without_contact_link': len(set(without_contact.values())),
           'project': with_contact, 'project_without_contact_link': without_contact, 'project_links': contact_kinds,
           'records': {s: r for s, r in records.items() if r},
           'cache_sha256': {f'{s}.txt': sha256(cache / f'{s}.txt') for s, r in records.items() if r}}
    out_path.write_text(json.dumps(out, indent=1, sort_keys=True))
    print(json.dumps({k: out[k] for k in ('series', 'projects', 'projects_without_contact_link')}
                     | {'not_read': len(out['not_read'])}))
    if out['not_read']:
        print(f"GEO record not read for {len(out['not_read'])} series; each stays its own project", file=sys.stderr)


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(p)
    out[order] = np.minimum(adjusted, 1)
    return out


def design_flags(usable, design_audit_dirs):
    """Per effect: what the dual-model design audit said about the specimen, where it read that comparison.

    The audit was run for the comparison ranking and covers the comparisons eligible there, so most effects of
    this analysis are unaudited; the counts below therefore carry their own denominator and are a disclosure, not
    a filter. `tumor` marks comparisons both models called tumour material, as for BLZ945, whose pooled effect
    comes from tumour-bearing mice.
    """
    from audit_system.comparison_ranking_2026_09_30 import KEY, design_audit_table
    audit = design_audit_table(design_audit_dirs)
    if audit is None:
        return usable.assign(design_audited=False, design_tumor=False)
    audit = audit.drop_duplicates(KEY)[KEY + ['tumor', 'specimen']]
    merged = usable.merge(audit, on=KEY, how='left', suffixes=('', '_audit'))
    if len(merged) != len(usable):
        raise SystemExit('the design audit matched an effect more than once')
    merged['design_audited'] = merged.tumor.notna()
    merged['design_tumor'] = merged.tumor.eq('YES')
    merged['design_not_from_organism'] = merged.specimen.isin(
        ['CELL_LINE', 'EX_VIVO_CULTURE', 'STEM_CELL_DERIVED'])
    return merged.drop(columns=['tumor', 'specimen'])


def run(atlas_path, classification_dir, out_dir, context=None, projects_path=None, link_contact=True,
        pool_across_categories=False, design_audit_dirs=(), keep_tumour=False):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    usable = usable_effects(atlas_path, classification_dir, context)
    usable = design_flags(usable, design_audit_dirs)
    # Correction of 2026-10-17 (docstring): drop the comparisons the design audit excludes from the ranking, so
    # that the same standard applies here. Only agreement excludes; where the two models disagree the comparison
    # stays and is flagged. --keep-tumour reproduces the earlier runs.
    excluded = {'tumour': int(usable.design_tumor.sum()),
                'not_from_organism': int(usable.design_not_from_organism.sum()),
                'either': int((usable.design_tumor | usable.design_not_from_organism).sum()),
                'audited': int(usable.design_audited.sum()), 'effects_before': int(len(usable))}
    if not keep_tumour:
        usable = usable[~(usable.design_tumor | usable.design_not_from_organism)].copy()
    # Correction of 2026-10-15 (docstring): one group per intervention and atlas category. A missing category is
    # its own group, so that no effect is dropped by the grouping.
    usable['condition_category'] = usable.condition_category.fillna('Unclassified')
    unit = 'gse_id'
    if projects_path:
        projects = json.loads(Path(projects_path).read_text())
        mapping = projects['project' if link_contact else 'project_without_contact_link']
        usable['project'] = usable.gse_id.map(mapping).fillna(usable.gse_id)  # series outside the file: own project
        unit = 'project'
    usable.to_csv(out_dir / 'meta_input_effects.csv', index=False)
    results, units = [], []
    group_keys = GROUP_KEYS if pool_across_categories else GROUP_KEYS + ['condition_category']
    grouped = usable.groupby(group_keys, dropna=False)
    if sum(len(rows) for _, rows in grouped) != len(usable):
        raise SystemExit('grouping dropped effects')
    for keys, rows in grouped:
        studies = study_effects(rows, unit)
        pooled = dersimonian_laird(studies.g, studies['var'])
        units.append(studies.rename(columns={unit: 'unit'}).assign(**dict(zip(group_keys, keys))))
        results.append({**dict(zip(group_keys, keys)), **pooled, 'studies': ';'.join(studies[unit]),
                        'series': ';'.join(sorted(rows.gse_id.unique())), 'n_series': int(rows.gse_id.nunique()),
                        'modalities': ';'.join(sorted(rows.modality.unique())), 'effects': int(len(rows)),
                        'clocks': ';'.join(sorted(rows.clock.unique())),
                        'design_audited_comparisons': int(rows[rows.design_audited]
                                                          .drop_duplicates(['gse_id', 'comparison_name']).shape[0]),
                        'design_tumor_comparisons': int(rows[rows.design_tumor]
                                                        .drop_duplicates(['gse_id', 'comparison_name']).shape[0]),
                        'comparisons': int(rows.drop_duplicates(['gse_id', 'comparison_name']).shape[0]),
                        'provisional_share': float((rows.v3_final_class == 'PROVISIONAL').mean()),
                        'relevance_mean': float(rows.drop_duplicates(['gse_id', 'comparison_name']).relevance_score.mean())})
    table = pd.DataFrame(results)
    if len(table):
        table['fdr'] = bh(table.p)
        table['direction'] = np.where(table.estimate < 0, 'age-lowering', 'age-raising')
        table['rank_bound'] = np.where(table.estimate < 0, -table.ci_high, table.ci_low)  # magnitude of bound nearest zero
        table['tier'] = np.where(table.k >= 2, 'replicated', 'single-study')
        table = table.sort_values('rank_bound', ascending=False)
    table.to_csv(out_dir / 'intervention_meta.csv', index=False)
    pd.concat(units).to_csv(out_dir / 'unit_effects.csv', index=False)  # one row per intervention and study or project
    robust = table[(table.fdr < 0.05) & (table.rank_bound > 0)] if len(table) else table
    for tier in ('replicated', 'single-study'):
        for direction in ('age-lowering', 'age-raising'):
            part = robust[(robust.tier == tier) & (robust.direction == direction)]
            part.to_csv(out_dir / f'ranked_{tier}_{direction}.csv', index=False)
            # Sensitivity: the same ranking restricted to interventions the models rated relevant (blind to results).
            part[part.relevance_mean >= RELEVANCE_THRESHOLD].to_csv(
                out_dir / f'ranked_{tier}_{direction}_relevance{RELEVANCE_THRESHOLD}.csv', index=False)
    summary = {'effects_in': int(len(usable)), 'interventions': int(len(table)),
               'by_tier': table.tier.value_counts().to_dict() if len(table) else {},
               'robust_by_tier_direction': {' '.join(k): int(v) for k, v in robust.groupby(['tier', 'direction']).size().items()} if len(robust) else {},
               'robust_relevant_by_tier_direction': {' '.join(k): int(v) for k, v in robust[robust.relevance_mean >= RELEVANCE_THRESHOLD]
                                                     .groupby(['tier', 'direction']).size().items()} if len(robust) else {},
               'method': ('Hedges g per comparison x clock; within-study mean (variance = mean variance); '
                          'DerSimonian-Laird across studies per species, target, class and manipulation; BH over interventions'),
               'group_keys': list(group_keys), 'pool_across_categories': bool(pool_across_categories),
               'design_audit': [str(d) for d in design_audit_dirs], 'design_excluded': excluded,
               'keep_tumour': bool(keep_tumour),
               'design_audited_effects': int(usable.design_audited.sum()),
               'design_tumor_effects': int(usable.design_tumor.sum()),
               'unit': unit, 'context': context or 'all', 'excluded_clocks': list(EXCLUDED_CLOCKS),
               'link_contact': bool(projects_path) and link_contact,
               'projects_sha256': sha256(projects_path) if projects_path else None,
               'atlas_sha256': sha256(atlas_path), 'classification_sha256': sha256(Path(classification_dir) / 'comparisons.jsonl'),
               'module_sha256': sha256(__file__), 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=1, default=str))


def lifespan_expectations(genage_path, drugage_path):
    """{(class, target, manipulation): expected sign} from GenAge models and DrugAge (mouse entries only)."""
    genage = pd.read_csv(genage_path)
    genage = genage[(genage.organism == 'Mus musculus')
                    & genage['longevity influence'].isin(['Pro-Longevity', 'Anti-Longevity'])]
    expected = {}
    for symbol, influence in zip(genage.symbol.astype(str).str.strip().str.lower(), genage['longevity influence']):
        pro = influence == 'Pro-Longevity'
        expected[('GENE', symbol, 'DECREASE')] = 1 if pro else -1
        expected[('GENE', symbol, 'INCREASE')] = -1 if pro else 1
    drugage = pd.read_csv(drugage_path)
    drugage = drugage[(drugage.species == 'Mus musculus') & (drugage.avg_lifespan_significance == 'S')]
    for compound, part in drugage.groupby(drugage.compound_name.astype(str).str.strip().str.lower()):
        up, down = (part.avg_lifespan_change_percent > 0).any(), (part.avg_lifespan_change_percent < 0).any()
        if up != down:
            expected[('EXPOSURE', compound, 'EXPOSURE')] = -1 if up else 1
    return expected


def concordance(meta_path, genage_path, drugage_path, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    table = pd.read_csv(meta_path)
    table = table[table.species == 'mouse'].copy()
    expected = lifespan_expectations(genage_path, drugage_path)
    keys = zip(table.intervention_class, table.canonical_target.astype(str).str.strip().str.lower(), table.manipulation)
    table['expected'] = [expected.get(key) for key in keys]
    matched = table[table.expected.notna()].copy()
    matched['source'] = np.where(matched.intervention_class == 'GENE', 'GenAge', 'DrugAge')
    matched['concordant'] = np.sign(matched.estimate) == matched.expected
    matched['interval_excludes_zero'] = (matched.ci_low > 0) | (matched.ci_high < 0)
    sets = {'all matched': matched, 'interval excludes zero': matched[matched.interval_excludes_zero],
            'two or more projects': matched[matched.k >= 2]}
    summary = {'genage_sha256': sha256(genage_path), 'drugage_sha256': sha256(drugage_path), 'meta_sha256': sha256(meta_path),
               'expectations': len(expected), 'sets': {}}
    for name, part in sets.items():
        for source in ('GenAge', 'DrugAge', 'both'):
            rows = part if source == 'both' else part[part.source == source]
            k, n = int(rows.concordant.sum()), int(len(rows))
            summary['sets'].setdefault(name, {})[source] = {
                'concordant': k, 'n': n, 'share': k / n if n else None,
                'p_binomial': float(stats.binomtest(k, n, 0.5).pvalue) if n else None}
    # Descriptive, added after the first run (2026-10-14): most interventions raise clock outputs, so the share of
    # matched interventions expected to raise them is what a rule calling every effect 'higher' would reach, and the
    # shares are also given by expected direction.
    summary['post_hoc'] = {
        name: {'expected_higher': int((part.expected > 0).sum()),
               'constant_higher_rule_share': float((part.expected > 0).mean()) if len(part) else None,
               'concordant_when_expected_higher': int(part[part.expected > 0].concordant.sum()),
               'concordant_when_expected_lower': int(part[part.expected < 0].concordant.sum()),
               'expected_lower': int((part.expected < 0).sum())}
        for name, part in sets.items()}
    matched.to_csv(out_dir / 'lifespan_concordance.csv', index=False)
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary['sets'], indent=1))


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == 'concordance':
        parser = argparse.ArgumentParser(description='lifespan-database concordance of the pooled effects')
        parser.add_argument('command')
        parser.add_argument('--meta', type=Path, required=True)
        parser.add_argument('--genage', type=Path, required=True)
        parser.add_argument('--drugage', type=Path, required=True)
        parser.add_argument('--out', type=Path, required=True)
        args = parser.parse_args(argv)
        return concordance(args.meta, args.genage, args.drugage, args.out)
    if not argv or argv[0] not in ('geo', 'run'):
        argv = ['run'] + argv  # the earlier flat command line
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('geo', 'run'):
        p = sub.add_parser(name)
        p.add_argument('--atlas', type=Path, required=True)
        p.add_argument('--classification', type=Path, required=True)
        p.add_argument('--out', type=Path, required=True)
    sub.choices['geo'].add_argument('--cache', type=Path, required=True)
    sub.choices['geo'].add_argument('--offline', action='store_true', help='use cached records only')
    sub.choices['run'].add_argument('--context', choices=['in_vivo', 'in_vitro'])
    sub.choices['run'].add_argument('--projects', type=Path, help='projects.json from the geo command')
    sub.choices['run'].add_argument('--no-contact-link', action='store_true',
                                    help='sensitivity: projects by SuperSeries and PubMed ID only')
    sub.choices['run'].add_argument('--pool-across-categories', action='store_true',
                                    help='reproduce the 2026-10-14 grouping, which pooled an intervention across '
                                         'the four atlas categories')
    sub.choices['run'].add_argument('--design-audit', type=Path, action='append', default=[],
                                    help='design-audit run whose specimen and tumour verdicts exclude comparisons, '
                                         'as they do for the ranking, and annotate the pooled rows')
    sub.choices['run'].add_argument('--keep-tumour', action='store_true',
                                    help='sensitivity: keep tumour material and cultured samples, as the runs '
                                         'before 2026-10-17 did')
    args = parser.parse_args(argv)
    if args.command == 'geo':
        geo(args.atlas, args.classification, args.cache, args.out, args.offline)
    else:
        run(args.atlas, args.classification, args.out, args.context, args.projects, not args.no_contact_link,
            args.pool_across_categories, args.design_audit, args.keep_tumour)


if __name__ == '__main__':
    main()
