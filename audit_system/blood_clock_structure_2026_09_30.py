"""Cross-clock structure in whole blood: raw predictions versus age-adjusted residuals.

The historical Extended Data Figure 2 correlated raw clock predictions pooled across blood studies. Two
problems make that estimand uninformative about clock biology. Every clock tracks chronological age, so
raw correlations are dominated by the shared age signal; and pooling across studies mixes in study-level
offsets from different preprocessing. The historical matrix also carried constant outputs from studies
whose DNAm matrix was empty. This analysis:
  - takes exact (GSE, GSM) pairs annotated Whole Blood in the GPT-4 metadata (the historical membership
    source) with a recorded age of 18-100 years;
  - uses the stream-repaired predictions (re-derived from the GEO supplementary data, QC-retained samples)
    for studies whose historical matrix was empty, and excludes empty-matrix studies without a repair; a
    re-derived clock whose coverage check withheld it is not used (its values are imputed constants);
  - excludes samples whose historical input profile is detached from its probe identities
    (input_integrity_2026_09_30);
  - keeps studies with at least MIN_STUDY samples with all 17 clocks of the historical PCA;
  - within each study, regresses each clock on age and correlates the residuals (Spearman and Pearson); a
    clock whose outputs are constant within a study contributes no correlation there;
    study matrices are pooled by Fisher z weighted by n - 3. The pooled raw-value correlation (the
    historical estimand) is reported beside it;
  - clusters the clocks by average linkage on 1 - r and reports each clock's nearest neighbour and its
    cluster at k = 2, 3 and 4;
  - with --grimage (grimage_age_input_2026_09_30 output), repeats the within-study residual analysis on
    the stream-repaired studies twice: with the GrimAge outputs of the historical recipe, which received
    the Horvathv1 prediction as age, and with GrimAge recomputed from the recorded chronological age.
Descriptive; no P values.

  python -m audit_system.blood_clock_structure_2026_09_30 --integrity DIR --out DIR [--grimage DIR]
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform

from audit_system import paths
from audit_system.membership_recovery_2026_09_27 import GEO, sha256

PREDICTIONS = GEO / 'result' / 'qced_human_dnam_age.csv'
METADATA = GEO / 'result' / 'gpt4_human_dnam_metadata.csv'
EMPTY_MATRICES = paths.EMPTY_MATRIX_LIST
STREAM_REPAIR = paths.STREAM_EMPTY
CLOCKS = ['Horvathv1', 'Hannum', 'Lin', 'PhenoAge', 'YingCausAge', 'YingDamAge', 'YingAdaptAge', 'Horvathv2', 'PEDBE',
          'Zhang_10', 'DunedinPoAm38', 'DunedinPACE', 'GrimAgeV1', 'GrimAgeV2', 'DNAmTL', 'Knight', 'LeeRobust']
MIN_AGE, MAX_AGE, MIN_STUDY, DETACHED_R = 18, 100, 20, 0.10  # row_verdicts_2026_09_27.PROFILE_DETACHED_R


def repaired_predictions(gse):
    """Stream-repaired clock outputs for QC-retained samples, or None when the study was not repaired."""
    study = STREAM_REPAIR / 'studies' / gse
    status = study / 'status.json'
    if not status.exists() or json.loads(status.read_text()).get('state') != 'REPAIRED':
        return None
    table = pd.read_csv(study / 'clock_predictions.csv', index_col=0)
    table = table.rename(columns={'DNAmGrimAge': 'GrimAgeV1', 'DNAmGrimAge.1': 'GrimAgeV2'})
    withheld = [v['clock'] for v in json.loads(status.read_text()).get('clock_verdicts', []) if v.get('verdict') == 'WITHHELD']
    table[[c for c in withheld if c in table.columns]] = np.nan  # failed coverage: outputs are imputed constants
    retained = pd.read_csv(study / 'sample_status.csv', index_col=0)
    table = table.loc[table.index.isin(retained.index[retained.status == 'RETAINED'])]
    return table.reindex(columns=CLOCKS).assign(gse_id=gse, gsm_id=table.index.astype(str)).reset_index(drop=True)


def cohort(integrity_dir):
    """Eligible samples, with the count left after each step."""
    historical = pd.read_csv(PREDICTIONS, low_memory=False, dtype={'gse_id': str, 'gsm_id': str})
    meta = pd.read_csv(METADATA, low_memory=False, dtype=str, usecols=['gse_id', 'gsm_id', 'tissue', 'age', 'sex'])
    meta = meta.drop_duplicates(['gse_id', 'gsm_id'], keep=False)  # conflicting duplicate annotations are not used
    blood = meta[meta.tissue.eq('Whole Blood')].assign(age=lambda m: pd.to_numeric(m.age, errors='coerce'))
    blood = blood[blood.age.between(MIN_AGE, MAX_AGE)]
    empty = {json.loads(line)['study_id'] for line in open(EMPTY_MATRICES)}
    steps = {'whole_blood_pairs_age_18_100': len(blood)}
    frame = historical[~historical.gse_id.isin(empty)].merge(blood, on=['gse_id', 'gsm_id']).assign(source='historical')
    repaired = [repaired_predictions(g) for g in sorted(set(blood.gse_id) & empty)]
    repaired = pd.concat([r for r in repaired if r is not None] or [pd.DataFrame(columns=CLOCKS + ['gse_id', 'gsm_id'])])
    repaired = repaired.merge(blood, on=['gse_id', 'gsm_id']).assign(source='stream_repaired')
    steps['historical_non_empty'] = len(frame)
    steps['stream_repaired_empty'] = len(repaired)
    steps['empty_matrix_not_repaired'] = int(blood.gse_id.isin(empty).sum() - len(repaired))
    frame = pd.concat([frame, repaired], ignore_index=True)
    screen = pd.read_csv(Path(integrity_dir) / 'samples.csv', usecols=['gse_id', 'gsm', 'r_reference'], dtype={'gsm': str})
    bad = screen[screen.r_reference < DETACHED_R]
    detached = set(zip(bad.gse_id, bad.gsm))
    keep = [(src == 'stream_repaired') or ((g, s) not in detached) for src, g, s in zip(frame.source, frame.gse_id, frame.gsm_id)]
    steps['detached_profiles_removed'] = int(len(frame) - sum(keep))
    frame = frame[keep]
    for clock in CLOCKS:
        frame[clock] = pd.to_numeric(frame[clock], errors='coerce')
    frame = frame.dropna(subset=CLOCKS)
    steps['all_17_clocks_present'] = len(frame)
    sizes = frame.groupby('gse_id').size()
    frame = frame[frame.gse_id.isin(sizes.index[sizes >= MIN_STUDY])]
    steps[f'in_studies_with_at_least_{MIN_STUDY}'] = len(frame)
    return frame.reset_index(drop=True), steps


def age_residuals(frame):
    """Residual of each clock from a linear regression on age within each study."""
    out = frame[CLOCKS].copy()
    for _, index in frame.groupby('gse_id').groups.items():
        age = frame.loc[index, 'age'].to_numpy(float)
        design = np.column_stack([np.ones_like(age), age])
        values = frame.loc[index, CLOCKS].to_numpy(float)
        coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
        residual = values - design @ coefficients
        constant = np.ptp(values, axis=0) == 0  # residuals would be floating-point noise
        residual[:, constant] = np.nan
        out.loc[index, CLOCKS] = residual
    return out


def pooled_within_study(values, studies, method):
    """Fisher-z average of per-study correlation matrices, weighted by n - 3; undefined entries are skipped."""
    total, weight = np.zeros((len(CLOCKS), len(CLOCKS))), np.zeros((len(CLOCKS), len(CLOCKS)))
    for _, index in studies.groupby(studies).groups.items():
        corr = values.loc[index].corr(method=method).to_numpy()
        z = np.arctanh(np.clip(corr, -0.999999, 0.999999))
        ok = np.isfinite(z)
        total[ok] += (len(index) - 3) * z[ok]
        weight[ok] += len(index) - 3
    pooled = np.tanh(np.divide(total, weight, out=np.full_like(total, np.nan), where=weight > 0))
    np.fill_diagonal(pooled, 1.0)
    return pd.DataFrame(pooled, index=CLOCKS, columns=CLOCKS)


def structure(corr):
    distance = squareform(1 - corr.to_numpy(), checks=False)
    tree = linkage(distance, 'average')
    info = {'order': [CLOCKS[i] for i in leaves_list(tree)],
            'nearest_neighbour': {c: corr[c].drop(c).idxmax() for c in CLOCKS}, 'linkage': tree.tolist()}
    for k in (2, 3, 4):
        info[f'clusters_k{k}'] = dict(zip(CLOCKS, map(int, fcluster(tree, k, 'maxclust'))))
    return info


def grimage_contrast(frame, grimage_dir, out_dir):
    """Within-study residual structure on repaired studies with historical versus chronological-age GrimAge."""
    recomputed = pd.read_csv(Path(grimage_dir) / 'grimage_age_input.csv', dtype={'gsm_id': str})
    subset = frame[frame.source == 'stream_repaired'].merge(recomputed.drop(columns=['age', 'Horvathv1', 'predicted_sex']),
                                                             on=['gse_id', 'gsm_id'], how='inner')
    sizes = subset.groupby('gse_id').size()
    subset = subset[subset.gse_id.isin(sizes.index[sizes >= MIN_STUDY])].reset_index(drop=True)
    variants = {'historical_horvath_age_input': subset[CLOCKS],
                'sex_corrected_horvath_age_input': subset[CLOCKS].assign(GrimAgeV1=subset.GrimAgeV1_sex_corrected,
                                                                          GrimAgeV2=subset.GrimAgeV2_sex_corrected),
                'chronological_age_input': subset[CLOCKS].assign(GrimAgeV1=subset.GrimAgeV1_chronological,
                                                                  GrimAgeV2=subset.GrimAgeV2_chronological)}
    out = {'samples': len(subset), 'studies': int(subset.gse_id.nunique()),
           'stored_vs_recomputed_historical_max_abs': {
               v: float(np.nanmax(np.abs(subset[v] - subset[f'{v}_historical']))) for v in ('GrimAgeV1', 'GrimAgeV2')}}
    for name, values in variants.items():
        residuals = age_residuals(values.assign(gse_id=subset.gse_id, age=subset.age))
        corr = pooled_within_study(residuals, subset.gse_id, 'spearman')
        corr.to_csv(out_dir / f'correlation_repaired_residual_spearman_{name}.csv')
        info = structure(corr)
        info['horvath_grimage_residual_r'] = {v: float(corr.loc['Horvathv1', v]) for v in ('GrimAgeV1', 'GrimAgeV2')}
        out[name] = info
    return out


def run(integrity_dir, out_dir, grimage_dir=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    frame, steps = cohort(integrity_dir)
    frame[['gse_id', 'gsm_id', 'source', 'age', 'sex'] + CLOCKS].to_csv(out_dir / 'cohort.csv', index=False)
    residuals = age_residuals(frame)
    residuals.assign(gse_id=frame.gse_id, gsm_id=frame.gsm_id).to_csv(out_dir / 'age_residuals.csv', index=False)
    matrices = {'raw_pooled_pearson': frame[CLOCKS].corr(),  # the historical estimand, on the cleaned cohort
                'residual_within_study_spearman': pooled_within_study(residuals, frame.gse_id, 'spearman'),
                'residual_within_study_pearson': pooled_within_study(residuals, frame.gse_id, 'pearson'),
                'raw_within_study_spearman': pooled_within_study(frame[CLOCKS], frame.gse_id, 'spearman')}
    results = {}
    for name, corr in matrices.items():
        corr.to_csv(out_dir / f'correlation_{name}.csv')
        results[name] = structure(corr)
    age_r = {}
    for clock in CLOCKS:
        zs, ws = [], []
        for _, x in frame.groupby('gse_id'):
            r = x[clock].corr(x.age, method='spearman')
            if np.isfinite(r):
                zs.append(np.arctanh(np.clip(r, -0.999999, 0.999999)))
                ws.append(len(x) - 3)
        age_r[clock] = float(np.tanh(np.average(zs, weights=ws)))
    if grimage_dir:
        results['grimage_age_input_contrast'] = grimage_contrast(frame, grimage_dir, out_dir)
    summary = {'cohort_steps': steps, 'samples': len(frame), 'studies': int(frame.gse_id.nunique()),
               'samples_by_source': frame.source.value_counts().to_dict(),
               'age_range': [float(frame.age.min()), float(frame.age.max())],
               'within_study_clock_age_spearman': age_r, 'structures': results,
               'inputs': {str(p): sha256(p) for p in (PREDICTIONS, METADATA, EMPTY_MATRICES, Path(integrity_dir) / 'samples.csv')
                   + ((Path(grimage_dir) / 'grimage_age_input.csv',) if grimage_dir else ())},
               'module_sha256': sha256(__file__), 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2))
    for name, info in results.items():
        if name == 'grimage_age_input_contrast':
            for variant in ('historical_horvath_age_input', 'sex_corrected_horvath_age_input', 'chronological_age_input'):
                print(variant, '| Horvathv1 nearest:', info[variant]['nearest_neighbour']['Horvathv1'],
                      '| r(Horvath, GrimAge) residual:', {k: round(v, 3) for k, v in info[variant]['horvath_grimage_residual_r'].items()})
            continue
        print(name, '| Horvathv1 nearest:', info['nearest_neighbour']['Horvathv1'], '| k4:',
              {c: info['clusters_k4'][c] for c in ('Horvathv1', 'Hannum', 'Lin', 'Horvathv2', 'PhenoAge', 'GrimAgeV1',
                                                   'GrimAgeV2', 'DunedinPACE')})
    print(json.dumps(steps, indent=1))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--integrity', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--grimage', type=Path, help='grimage_age_input_2026_09_30 output directory')
    args = parser.parse_args(argv)
    run(args.integrity, args.out, args.grimage)


if __name__ == '__main__':
    main()
