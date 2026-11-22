"""Execute EN scaled-log RNA clocks on preserved preprocessed molecular inputs.

Checks prediction parity with saved outputs and an independent NumPy evaluation
of persisted parameters. Does not rerun raw-count preprocessing, validate a
biological comparison, or establish historical execution provenance.
"""
import argparse
import csv
from datetime import datetime,timezone
import importlib.util
import json
from pathlib import Path
import re
import time
import warnings

import joblib
import numpy as np
import pandas as pd
import scipy
import sklearn

from .agent_workflow_2026_09_10 import sha,write_new


def independent_prediction(model,frame, *, inverse_transform=True):
    """Read fitted parameters; do not call sklearn transform/predict/support."""
    if [type(step).__name__ for _,step in model.steps]!=['SimpleImputer','StandardScaler','SelectKBest','ElasticNet']:
        raise ValueError('Independent verifier does not support this estimator sequence')
    imputer,scaler,selector,estimator=[step for _,step in model.steps]
    if imputer.strategy not in ('mean','median') or imputer.add_indicator or not np.isnan(imputer.missing_values):
        raise ValueError('Unsupported imputation semantics')
    values=frame[model.feature_names].to_numpy(dtype=float,copy=True)
    means=imputer.statistics_
    keep=~np.isnan(means)
    if getattr(imputer,'keep_empty_features',False) and not keep.all():
        raise ValueError('Unsupported retained empty training features')
    values=values[:,keep];means=means[keep]
    values=np.where(np.isnan(values),means,values)
    if scaler.with_mean:values=values-scaler.mean_
    if scaler.with_std:values=values/scaler.scale_
    scores=np.asarray(selector.scores_).copy()
    scores[np.isnan(scores)]=np.finfo(scores.dtype).min
    mask=np.ones(len(scores),dtype=bool) if selector.k=='all' else np.zeros(len(scores),dtype=bool)
    if selector.k!='all' and selector.k>0:mask[np.argsort(scores,kind='mergesort')[-selector.k:]]=True
    linear=values[:,mask]@estimator.coef_+estimator.intercept_
    return 10**(-(10**(-linear))) if inverse_transform else linear


def run(geo_root,study,clock_names,output):
    if not re.fullmatch(r'GSE\d+',study):raise ValueError('Invalid study identity')
    if sklearn.__version__!='1.3.2':raise ValueError('Use the pinned sklearn 1.3.2 replay environment')
    output.mkdir(parents=True,exist_ok=False)
    source=geo_root/'src_rna/2_run_rna_clock.py'
    spec=importlib.util.spec_from_file_location('clockbase_preserved_producer',source)
    producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)
    matrix=geo_root/f'GEO_mm_rna/processed_data/{study}/{study}_scaled_data.csv'
    saved=geo_root/f'GEO_mm_rna/clock_result/{study}.csv'
    metadata=geo_root/f'GEO_mm_rna/metadata/{study}.csv'
    paths=[source,matrix,saved,metadata]
    inputs={str(p):{'sha256':sha(p),'size_bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in paths}
    with matrix.open() as handle:header=next(csv.reader(handle))
    if len(set(header[1:]))!=len(header[1:]):raise ValueError('Duplicate molecular feature columns')
    frame=pd.read_csv(matrix,index_col=0)
    original=pd.read_csv(saved).set_index('sample_id',verify_integrity=True)
    meta=pd.read_csv(metadata)
    if frame.index.has_duplicates or set(frame.index)!=set(original.index) or set(frame.index)!=set(meta['GSM']):
        raise ValueError('Molecular/source/prediction sample identity differs')
    predictions=pd.DataFrame(index=frame.index);results=[];started=time.monotonic()
    for name in clock_names:
        if not name.endswith('_scaledlog_EN'):raise ValueError('Only scaled-log EN clock replay is supported')
        matches=[p for folder in producer.clock_folders for p in Path(folder).glob('*.pkl') if producer.generate_clock_name(str(p))==name]
        if len(matches)!=1:raise ValueError('Clock model identity absent or ambiguous: '+name)
        model_path=matches[0];inputs[str(model_path)]={'sha256':sha(model_path),'size_bytes':model_path.stat().st_size,'mtime_ns':model_path.stat().st_mtime_ns}
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always');model=producer.load_clock(model_path)
            predicted=np.asarray(producer.predict_tAge(frame,model,name)['tAge'])
        if any('version' in str(w.message).lower() for w in caught):raise ValueError('Persisted estimator version differs from runtime')
        verified=independent_prediction(model,frame)
        column=name+'_tAge'
        historical=original.loc[frame.index,column].to_numpy(dtype=float)
        finite=np.isfinite(frame[model.feature_names].to_numpy(dtype=float))
        checks={'independent_parameter_calculation':bool(np.allclose(predicted,verified,rtol=1e-10,atol=1e-14)),
                'saved_prediction_parity':bool(np.allclose(predicted,historical,rtol=1e-10,atol=1e-14))}
        results.append({'clock':name,'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
            'imputation_strategy':model.steps[0][1].strategy,
            'sample_count':len(frame),'input_clock_features':len(model.feature_names),
            'observed_clock_features_per_sample_min':int(finite.sum(axis=1).min()),
            'observed_clock_features_per_sample_max':int(finite.sum(axis=1).max()),
            'max_abs_difference_independent':float(np.max(np.abs(predicted-verified))),
            'max_abs_difference_saved':float(np.max(np.abs(predicted-historical))),
            'warnings':[str(w.message) for w in caught]})
        predictions[column]=predicted;predictions[column+'_independent']=verified
    for path,record in inputs.items():
        if sha(path)!=record['sha256']:raise ValueError('Replay input changed: '+path)
    predictions.rename_axis('sample_id').to_csv(output/'predictions.csv',mode='x')
    write_new(output/'result.json',{'study_id':study,'recorded_at_utc':datetime.now(timezone.utc).isoformat(),
        'elapsed_seconds':time.monotonic()-started,'results':results,'inputs':inputs,
        'implementation_sha256':sha(__file__),'prediction_output_sha256':sha(output/'predictions.csv'),
        'versions':{'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__,'joblib':joblib.__version__},
        'clock_predictions_recomputed':True,'raw_count_preprocessing_recomputed':False,
        'scientific_status':'CANNOT VERIFY','development_only':True,'api_calls':0,'api_cost_usd':0,
        'scope':'Two numerical parity checks per clock from existing preprocessed inputs; no full-pipeline or biological validity claim'})
    print(json.dumps({'study_id':study,'results':results,'scientific_status':'CANNOT VERIFY'}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--geo-root',type=Path,required=True);p.add_argument('--study',required=True)
    p.add_argument('--clock',dest='clock_names',action='append',required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    try:run(args.geo_root,args.study,args.clock_names,args.output)
    except Exception as error:
        if args.output.is_dir() and not (args.output/'failure.json').exists() and not (args.output/'result.json').exists():
            write_new(args.output/'failure.json',{'status':'CANNOT VERIFY',
                'error':f'{type(error).__name__}: {error}','implementation_sha256':sha(__file__),
                'recorded_at_utc':datetime.now(timezone.utc).isoformat()})
        raise SystemExit(f'Clock replay stopped: {type(error).__name__}: {error}')
