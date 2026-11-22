"""GrimAge with its intended inputs versus the historical inputs.

GrimAge takes chronological age and sex as inputs. The historical recipe (GEO_data/src/4_compute_age.py,
lines 227-235) had no chronological age, so it passed each sample's Horvathv1 prediction as "age" and the
SexEstimation label ('Female'/'Male') as "sex". biolearn's GrimageModel codes Female = 1 only when sex == 0,
so the string labels made every sample male, and the outputs embed the Horvathv1 estimate. This module
recomputes GrimAge V1 and V2 for whole-blood samples of the stream-repaired studies (measured clock CpGs
from the GEO supplementary data, GSE40279-mean imputation as in the historical recipe) three ways:
  historical      Horvathv1 prediction as age, string sex labels: must reproduce the stored repaired output
  sex_corrected   Horvathv1 prediction as age, sex coded 0 = female, 1 = male
  chronological   recorded age, sex coded 0 = female, 1 = male
Only the sites in the GrimAge coefficient files enter: imputation fills each site from its reference mean
independently, so the outputs equal those of the full-matrix recipe (checked by the reproduction column).
Samples and studies follow blood_clock_structure_2026_09_30.cohort.

  python -m audit_system.grimage_age_input_2026_09_30 --cohort DIR/cohort.csv --out DIR
"""
import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system import paths
from audit_system.membership_recovery_2026_09_27 import sha256

STREAM_REPAIR = paths.STREAM_EMPTY
REFERENCE_CACHE = STREAM_REPAIR / 'reference-cache'


def grimage_sites():
    from src.dnam_stream_repair_2026_09_26 import gallery
    sites = set()
    for name in ('GrimAgeV1', 'GrimAgeV2'):
        model = gallery().get(name)
        model = getattr(model, 'model', model)  # imputation decorator wraps the GrimageModel
        sites |= {v for v in model.coefficients['var'] if str(v).startswith(('cg', 'ch'))}
    return sorted(sites)


def predict(beta, age, sex):
    from biolearn.data_library import GeoData
    from biolearn.imputation import impute_from_standard
    from src.dnam_stream_repair_2026_09_26 import gallery, reference_averages
    means = reference_averages(REFERENCE_CACHE)
    index = means.index.union(beta.index).intersection(grimage_sites())
    data = GeoData(metadata=pd.DataFrame({'age': age, 'sex': sex}, index=beta.columns),
                   dnam=impute_from_standard(beta.reindex(index), means.reindex(index)))
    out = {}
    for name in ('GrimAgeV1', 'GrimAgeV2'):
        result = gallery().get(name).predict(GeoData(metadata=data.metadata, dnam=data.dnam.copy()))
        out[name] = result['DNAmGrimAge'] if 'DNAmGrimAge' in result.columns else result.iloc[:, 0]
    return pd.DataFrame(out)


def study(gse, samples):
    folder = STREAM_REPAIR / 'studies' / gse
    beta = pd.read_csv(folder / 'clock_cpg_beta.csv.gz', index_col=0)
    stored = pd.read_csv(folder / 'clock_predictions.csv', index_col=0)
    stored = stored.rename(columns={'DNAmGrimAge': 'GrimAgeV1', 'DNAmGrimAge.1': 'GrimAgeV2'})
    gsms = [g for g in samples.gsm_id if g in beta.columns and g in stored.index]
    beta, stored, samples = beta[gsms], stored.loc[gsms], samples.set_index('gsm_id').loc[gsms]
    labels = stored['predicted_sex'].to_numpy()
    coded = np.where(labels == 'Female', 0, 1)
    variants = {'historical': predict(beta, stored['Horvathv1'].to_numpy(), labels),
                'sex_corrected': predict(beta, stored['Horvathv1'].to_numpy(), coded),
                'chronological': predict(beta, samples['age'].to_numpy(float), coded)}
    out = pd.DataFrame({'gse_id': gse, 'gsm_id': gsms, 'age': samples['age'].to_numpy(float),
                        'predicted_sex': labels, 'Horvathv1': stored['Horvathv1'].to_numpy(),
                        'GrimAgeV1_stored': stored['GrimAgeV1'].to_numpy(), 'GrimAgeV2_stored': stored['GrimAgeV2'].to_numpy()})
    for variant, values in variants.items():
        for version in ('GrimAgeV1', 'GrimAgeV2'):
            out[f'{version}_{variant}'] = values[version].to_numpy()
    return out


def run(cohort_path, out_dir, workers=8):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    cohort = pd.read_csv(cohort_path, dtype={'gsm_id': str})
    repaired = cohort[cohort.source == 'stream_repaired']
    jobs = [(gse, rows[['gsm_id', 'age']]) for gse, rows in repaired.groupby('gse_id')]
    with ProcessPoolExecutor(workers) as pool:
        parts = list(pool.map(study, *zip(*jobs)))
    table = pd.concat(parts, ignore_index=True)
    table.to_csv(out_dir / 'grimage_age_input.csv', index=False)
    check = {v: float(np.nanmax(np.abs(table[f'{v}_historical'] - table[f'{v}_stored']))) for v in ('GrimAgeV1', 'GrimAgeV2')}
    shifts = {f'{v}_{variant}_minus_historical': float((table[f'{v}_{variant}'] - table[f'{v}_historical']).mean())
              for v in ('GrimAgeV1', 'GrimAgeV2') for variant in ('sex_corrected', 'chronological')}
    age_r = {f'{v}_{variant}': float(np.corrcoef(table.age, table[f'{v}_{variant}'])[0, 1])
             for v in ('GrimAgeV1', 'GrimAgeV2') for variant in ('historical', 'sex_corrected', 'chronological')}
    summary = {'samples': len(table), 'studies': int(table.gse_id.nunique()),
               'predicted_sex_counts': table.predicted_sex.value_counts().to_dict(),
               'max_abs_difference_historical_recompute_vs_stored': check,
               'mean_shift_vs_historical': shifts, 'pooled_pearson_with_age': age_r,
               'cohort_sha256': sha256(cohort_path), 'module_sha256': sha256(__file__),
               'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=1))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--cohort', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args(argv)
    run(args.cohort, args.out, args.workers)


if __name__ == '__main__':
    main()
