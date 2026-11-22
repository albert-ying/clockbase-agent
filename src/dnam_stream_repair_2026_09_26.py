"""Streaming repair of human DNAm studies whose historical input matrix was empty.

For each study: download one GEO source (author beta/M-value table first, raw
IDATs as fallback), build a GSM-keyed beta matrix, run the historical biolearn
clock recipe, apply the verifier-v3 input contract, then delete the raw
download. Studies run concurrently under a disk reservation budget.

  python -m src.dnam_stream_repair_2026_09_26 run --inventory INV.csv --root ROOT
  python -m src.dnam_stream_repair_2026_09_26 study --root ROOT --gse GSE179325
  python -m src.dnam_stream_repair_2026_09_26 summary --root ROOT

The historical recipe (GEO_data/src/4_compute_age.py) is reproduced exactly:
GSE40279 population-mean imputation via impute_from_standard, age from
Horvathv1 and sex from SexEstimation(averaging) as GrimAge metadata, and the
25 clocks recorded in the historical missing-percentage tables. No model calls.
"""
import argparse
import csv
import fcntl
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from audit_system import paths  # noqa: E402
from audit_system.agent_workflow_2026_09_10 import sha  # noqa: E402

SCHEMA = 'clockbase-dnam-stream-repair-v1'
GEO_DNAM = paths.GEO / 'GEO_hs_dnam'
HISTORICAL_CLOCKS = [
    'Horvathv1', 'Hannum', 'Lin', 'PhenoAge', 'YingCausAge', 'YingDamAge', 'YingAdaptAge',
    'Horvathv2', 'PEDBE', 'Zhang_10', 'DunedinPoAm38', 'DunedinPACE', 'GrimAgeV1', 'GrimAgeV2',
    'AlcoholMcCartney', 'BMI_McCartney', 'DNAmTL', 'HRSInCHPhenoAge', 'Knight', 'LeeControl',
    'SexEstimation', 'LeeRefinedRobust', 'LeeRobust', 'SmokingMcCartney', 'DownSyndrome']
# Clocks whose verdict is not assessed: SexEstimation only supplies GrimAge metadata.
UTILITY_CLOCKS = {'SexEstimation'}
# Verifier-v3 stage-1 thresholds (minimum per-sample fraction of required CpGs measured).
COVERAGE_VERIFIED = 0.80
COVERAGE_PROVISIONAL = 0.60
DETECTION_P = 0.01
# Samples with more than this fraction of array features missing (detection-P or NA) fail sample QC.
SAMPLE_MAX_MISSING = 0.10
# Sample-level assay screen: enrichment or oxidative chemistries are not bisulfite beta values.
# Underscore counts as a separator (titles such as tissue_1_oxBS, GSE73895).
NON_STANDARD_ASSAY = re.compile(
    r'5hmc|hmedip|medip|hmc[-_ ]?ip|immuno[- ]?precipitat|(?<![a-z0-9])ip[-_ ]?(enriched|fraction)(?![a-z0-9])|'
    r'(?<![a-z0-9])oxbs(?![a-z0-9])|oxidative bisulfite|(?<![a-z0-9])tab[-_ ]?(seq|array)(?![a-z0-9])|'
    r'(?<![a-z0-9])mbd[-_ ]?(seq|cap|enrich)|(?<![a-z0-9])chip[-_ ]?(seq|on[-_ ]chip|enrich[a-z]*|input|dna|ip)(?![a-z0-9])', re.I)
# ChIP only with assay context; bare 'ChIP: 123'/'chip' fields name the array slide (e.g. GSE56046).
CHIP_ASSAY = re.compile(r'(?<![A-Za-z0-9])ChIP[-_ ]?(seq|on[-_ ]chip|enrich\w*|input|DNA|IP)(?![A-Za-z0-9])', re.I)
PREDICTION_CHUNK = 256
# Clocks whose preprocessing reads beyond their own coefficient sites; they get the full matrix.
# Every other clock gives identical output on its own sites (max |diff| = 0, checked
# on all 25 clocks; evidence in the run root, subset-equivalence-2026-09-26.json).
FULL_MATRIX_CLOCKS = {'DunedinPACE', 'SexEstimation'}
IDAT_CHUNK = 96
NA_TOKENS = ['', ' ', '  ', 'NA', 'NaN', 'nan', 'null', 'NULL', 'N/A', '.', '-']
DETECTION_COLUMN = re.compile(r'(detection|det)[ ._-]*p|p[ ._-]?val|pval', re.I)
ANNOTATION_COLUMNS = {
    'id_ref', 'id', 'probe', 'probeid', 'probe_id', 'targetid', 'target_id', 'ilmnid', 'name',
    'chr', 'chromosome', 'mapinfo', 'position', 'pos', 'start', 'end', 'strand', 'index', 'row',
    'gene', 'genes', 'ucsc_refgene_name', 'infinium_design_type', 'probe_type', 'type', 'color',
    'color_channel', 'relation_to_ucsc_cpg_island', 'cpg_island', 'annotation'}
ALIAS_KEYS = {'sample id', 'sample identifier', 'sample code', 'sample name', 'sample', 'id',
              'sentrix id', 'sentrix_id', 'sentrix barcode', 'barcode', 'array id', 'basename',
              'sample_id', 'sample_name', 'subject id', 'patient id', 'individual'}
BETA_SUFFIX = re.compile(r'([ ._-]*(avg[ ._-]?beta|beta[ ._-]?values?|beta|m[ ._-]?values?|methylation|meth))$', re.I)
AUTHOR_DNAM_AGE_KEY = re.compile(r'^\s*(dnam ?age|dna methylation age|horvath[ _-]?(dnam)?[ _-]?age)\s*$', re.I)
EPICV2_SUFFIX = re.compile(r'^(cg\d{8}|ch\.[^_]+|rs\d+)_[A-Z]{2}\d+$')
AGE_KEY = re.compile(r'^\s*(age|age \(y(ea)?rs?\)|age_years|age years|age at [a-z ]+)\s*$', re.I)


class RouteFailure(Exception):
    """A source route cannot produce a valid beta matrix; the reason is recorded."""

    def __init__(self, code, detail):
        super().__init__(f'{code}: {detail}')
        self.code, self.detail = code, detail


def now():
    return time.strftime('%Y-%m-%dT%H:%M:%S%z')


def json_safe(value):
    """NaN and infinities become null; numpy scalars become Python numbers."""
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    return value


def write_json(path, value):
    tmp = Path(str(path) + '.tmp')
    tmp.write_text(json.dumps(json_safe(value), indent=2, allow_nan=False, default=str) + '\n')
    tmp.replace(path)


def append_event(root, event):
    event = {'time': now(), **event}
    with open(Path(root) / 'events.jsonl', 'a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        handle.write(json.dumps(event, default=str) + '\n')
        fcntl.flock(handle, fcntl.LOCK_UN)


# ---------------------------------------------------------------- source metadata

def series_samples(gse):
    """Per-GSM header fields from every local GEO Series Matrix file of the study."""
    paths = sorted(GEO_DNAM.joinpath('series_matrix').glob(f'{gse}_series_matrix*.txt.gz')) + \
        sorted(GEO_DNAM.joinpath('series_matrix').glob(f'{gse}-GPL*_series_matrix*.txt.gz'))
    if not paths:
        raise FileNotFoundError(f'Cannot proceed - required data file not found: {GEO_DNAM}/series_matrix/{gse}*')
    samples, files = {}, []
    for path in paths:
        fields = {}
        with gzip.open(path, 'rt', newline='') as handle:
            for record in csv.reader(handle, delimiter='\t'):
                if record and record[0] == '!series_matrix_table_begin':
                    break
                if record and record[0].startswith('!Sample_'):
                    fields.setdefault(record[0][8:], []).append(record[1:])
        ids = fields.get('geo_accession', [[]])[0]
        for i, gsm in enumerate(ids):
            entry = samples.setdefault(gsm, {})
            for key, rows in fields.items():
                entry.setdefault(key, []).extend(row[i] for row in rows if i < len(row) and row[i] != '')
        files.append({'path': str(path), 'sha256': sha(path)})
    if not samples:
        raise RouteFailure('SERIES_WITHOUT_SAMPLES', 'no GSM records in local series matrix')
    return samples, files


def normalize_token(value):
    token = re.sub(r'__dup\d+$', '', str(value).strip().strip('"').strip())
    token = BETA_SUFFIX.sub('', token)
    return re.sub(r'[\s_.\-]+', '_', token).strip('_').casefold()


def sample_aliases(samples):
    """Candidate column tokens per GSM. Tokens owned by several GSMs are dropped."""
    owners = {}

    def add(token, gsm, rule):
        token = normalize_token(token)
        if token:
            owners.setdefault(token, {}).setdefault(gsm, rule)

    for gsm, fields in samples.items():
        add(gsm, gsm, 'gsm')
        for title in fields.get('title', []):
            add(title, gsm, 'title')
            words = re.split(r'[\s_]+', title.strip())
            if len(words) > 1:
                add(words[0], gsm, 'title_leading_token')
                add(words[-1], gsm, 'title_trailing_token')
            if ':' in title:
                for part in title.split(':'):
                    add(part, gsm, 'title_part')
            bracket = re.search(r'\[([^\[\]]+)\]\s*$', title)
            if bracket:
                add(bracket[1], gsm, 'title_bracket')
        for value in fields.get('description', []) + fields.get('source_name_ch1', []):
            add(value, gsm, 'description_or_source')
        sentrix = {}
        for value in fields.get('characteristics_ch1', []):
            key, sep, token = value.partition(':')
            key = key.strip().casefold()
            if not sep:
                continue
            if key in ALIAS_KEYS:
                add(token, gsm, f'characteristic:{key}')
            if 'sentrix' in key and ('position' in key or 'pos' in key):
                sentrix['pos'] = token.strip()
            elif 'sentrix' in key or key in {'slide', 'chip', 'array'}:
                sentrix['id'] = token.strip()
        if {'id', 'pos'} <= set(sentrix):
            add(f"{sentrix['id']}_{sentrix['pos']}", gsm, 'sentrix_id_position')
        for url in fields.get('supplementary_file', []):
            match = re.search(r'/(GSM\d+)_([^/]+)_(?:Grn|Red)\.idat(?:\.gz)?$', url)
            if match and match[1] == gsm:
                add(match[2], gsm, 'idat_basename')
    return {token: next(iter(gsms.items())) for token, gsms in owners.items() if len(gsms) == 1}


# ---------------------------------------------------------------- beta-table route

def table_layout(path):
    """Header, delimiter and data offset for plain or GEO-formatted tables."""
    if str(path).lower().endswith(('.xlsx', '.xls')):
        raise RouteFailure('UNSUPPORTED_SPREADSHEET', Path(path).name)
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='', errors='replace') as handle:
        skip, line = 0, handle.readline()
        while line and (line.startswith(('!', '#', '^')) or not line.strip()):
            skip += 1
            if line.strip() == '!series_matrix_table_begin':
                line = handle.readline()
                break
            line = handle.readline()
        if not line:
            raise RouteFailure('EMPTY_TABLE', Path(path).name)
        delimiter = '\t' if '\t' in line else ',' if ',' in line else None
        header = next(csv.reader([line], delimiter=delimiter)) if delimiter else line.split()
        first = handle.readline()
    width = len(next(csv.reader([first], delimiter=delimiter)) if delimiter else first.split())
    feature_column = width == len(header)
    if width == len(header) + 1:
        header = ['ID_REF'] + header
    elif width != len(header):
        raise RouteFailure('HEADER_WIDTH_MISMATCH', f'header {len(header)} vs row {width}')
    names, seen = [], {}
    for h in (h.strip().strip('"') for h in header):
        seen[h] = seen.get(h, 0) + 1
        names.append(h if seen[h] == 1 else f'{h}__dup{seen[h] - 1}')  # pandas needs unique names
    return {'header': names, 'delimiter': delimiter,
            'skiprows': skip, 'feature_column_labelled': feature_column}


def resolve_unlabelled_detection_pairs(path, layout, mapping):
    """Same-named column pairs: keep the beta column when exactly one of the pair looks like detection P.

    Detection-P-like: at least 95% of the first 5,000 non-missing tokens are numbers below 0.05
    or textual thresholds such as 'p<0.01' (GSE112596).
    """
    owners = {}
    for column, (gsm, _) in mapping.items():
        owners.setdefault(gsm, []).append(column)
    pairs = {g: c for g, c in owners.items() if len(c) > 1}
    columns = [layout['header'][0]] + [c for cs in pairs.values() for c in cs]
    head = pd.read_csv(path, sep=layout['delimiter'] or r'\s+', skiprows=layout['skiprows'] + 1,
                       names=layout['header'], header=None, usecols=columns, index_col=0, nrows=5000,
                       na_values=NA_TOKENS, keep_default_na=False, dtype=str,
                       engine='c' if layout['delimiter'] else 'python')

    def detection_like(tokens):
        tokens = tokens.dropna().str.strip()
        if tokens.empty:
            return False
        numeric = pd.to_numeric(tokens, errors='coerce')
        textual = tokens.str.fullmatch(r'p\s*[<≤]=?\s*0?\.0\d*', case=False)
        return bool(((numeric < 0.05) | textual).mean() >= 0.95)

    dropped = []
    for gsm, cols in pairs.items():
        plike = [c for c in cols if detection_like(head[c])]
        if not plike:
            # e.g. 'X_Meth' + 'X_Beta' (GSE52270): drop the intensity column, keep the one bounded by [0, 1].
            numeric = {c: pd.to_numeric(head[c], errors='coerce').dropna() for c in cols}
            plike = [c for c in cols if len(numeric[c]) and numeric[c].max() > 1.5]
            if len(plike) == 1 and all(numeric[c].between(0, 1).all() for c in cols if c not in plike):
                dropped.append(plike[0])
                continue
        if len(cols) != 2 or len(plike) != 1:
            raise RouteFailure('DUPLICATE_SAMPLE_COLUMNS', f'{gsm}: cannot separate beta from detection P in {cols}')
        dropped.append(plike[0])
    return {c: v for c, v in mapping.items() if c not in dropped}, dropped


def classify_columns(header, aliases, gsms, allow_duplicates=False):
    """Map each value column to a GSM, or to detection-P / annotation. Unmapped numeric columns fail."""
    mapping, detection, annotation, unmapped = {}, [], [], []
    for column in header[1:]:
        token = normalize_token(column)
        if DETECTION_COLUMN.search(column):
            detection.append(column)
        elif token in aliases:
            gsm, rule = aliases[token]
            mapping[column] = (gsm, rule)
        elif token in ANNOTATION_COLUMNS or token.replace(' ', '_') in ANNOTATION_COLUMNS:
            annotation.append(column)
        else:
            unmapped.append(column)
    owners = {}
    for column, (gsm, _) in mapping.items():
        owners.setdefault(gsm, []).append(column)
    duplicated = {gsm: cols for gsm, cols in owners.items() if len(cols) > 1}
    if duplicated and not allow_duplicates:
        raise RouteFailure('DUPLICATE_SAMPLE_COLUMNS', json.dumps(dict(list(duplicated.items())[:5])))
    if not mapping:
        raise RouteFailure('NO_SAMPLE_COLUMNS_MAPPED', f'{len(unmapped)} unmapped columns, e.g. {unmapped[:5]}')
    if set(owners) - set(gsms):
        raise RouteFailure('FOREIGN_SAMPLES', str(sorted(set(owners) - set(gsms))[:5]))
    return mapping, detection, annotation, unmapped


def collapse_feature_ids(frame):
    """EPICv2 replicate probes (cgXXXXXXXX_TC21) are averaged to their cg identifier."""
    ids = frame.index.to_series().str.replace(EPICV2_SUFFIX, r'\1', regex=True)
    if ids.equals(frame.index.to_series()):
        return frame, 0
    collapsed = frame.groupby(ids.values).mean()
    return collapsed, int(len(frame) - len(collapsed))


def read_beta_table(path, samples, *, kind):
    """GSM-keyed beta matrix (features x GSM) from an author table, with its evidence."""
    layout = table_layout(path)
    aliases = sample_aliases(samples)
    mapping, detection, annotation, unmapped = classify_columns(layout['header'], aliases, list(samples),
                                                               allow_duplicates=True)
    inferred_detection = []
    if len({g for g, _ in mapping.values()}) < len(mapping):
        mapping, inferred_detection = resolve_unlabelled_detection_pairs(path, layout, mapping)
    columns = [layout['header'][0], *mapping]
    frames = []
    try:
        reader = pd.read_csv(path, sep=layout['delimiter'] or r'\s+', skiprows=layout['skiprows'] + 1,
                         names=layout['header'], header=None, usecols=columns, index_col=0,
                         na_values=NA_TOKENS, keep_default_na=False, chunksize=100_000,
                         dtype={c: 'float64' for c in mapping} | {layout['header'][0]: str},
                         comment=None, engine='c' if layout['delimiter'] else 'python',
                         on_bad_lines='error')
        for chunk in reader:
            chunk = chunk[~chunk.index.isin(['!series_matrix_table_end'])]
            frames.append(chunk.astype('float32'))
    except (ValueError, pd.errors.ParserError) as error:
        raise RouteFailure('NON_NUMERIC_OR_MALFORMED_VALUES', str(error)[:300]) from error
    frame = pd.concat(frames)
    frame.index = frame.index.astype(str).str.strip().str.strip('"')
    frame.columns = [mapping[c][0] for c in frame.columns]
    if frame.index.has_duplicates:
        frame, _ = collapse_feature_ids(frame)
        if frame.index.has_duplicates:
            raise RouteFailure('DUPLICATE_FEATURES', str(frame.index[frame.index.duplicated()][:5].tolist()))
    frame, collapsed = collapse_feature_ids(frame)
    values = frame.to_numpy()
    finite = np.isfinite(values)
    if not finite.any():
        raise RouteFailure('NO_FINITE_VALUES', Path(path).name)
    lo, hi = float(np.nanmin(values)), float(np.nanmax(values))
    scale = 'beta'
    if lo < 0 or kind == 'mvalue':
        if lo < -20 or hi > 20:
            raise RouteFailure('UNRECOGNISED_VALUE_SCALE', f'range [{lo:.3g}, {hi:.3g}]')
        frame = (2.0 ** frame) / (1.0 + 2.0 ** frame)
        scale = 'mvalue_converted'
    elif hi > 1.0 + 1e-6:
        raise RouteFailure('UNRECOGNISED_VALUE_SCALE', f'range [{lo:.3g}, {hi:.3g}] (not beta or M)')
    frame = frame.clip(0.0, 1.0) if scale == 'beta' else frame
    empty = frame.columns[~np.isfinite(frame.to_numpy()).any(axis=0)]
    if len(empty):
        raise RouteFailure('SAMPLE_WITHOUT_VALUES', str(list(empty)[:5]))
    evidence = {
        'route': 'author_table', 'file': Path(path).name, 'layout': {k: v for k, v in layout.items() if k != 'header'},
        'value_scale': scale, 'features': int(len(frame)), 'epicv2_replicates_collapsed': collapsed,
        'mapped_columns': len(mapping), 'detection_p_columns_ignored': len(detection),
        'unlabelled_detection_p_columns_inferred': len(inferred_detection),
        'annotation_columns': annotation, 'unmapped_columns': unmapped,
        'mapping_rules': pd.Series([r for _, r in mapping.values()]).value_counts().to_dict(),
        'series_samples': len(samples), 'samples_in_table': frame.shape[1],
        'series_samples_absent_from_table': sorted(set(samples) - set(frame.columns))}
    if unmapped:
        # Numeric columns that match no source sample make the membership ambiguous.
        raise RouteFailure('UNMAPPED_COLUMNS', f'{len(unmapped)} columns, e.g. {unmapped[:5]}')
    sample_map = pd.DataFrame([{'column': c, 'gsm': g, 'rule': r} for c, (g, r) in mapping.items()])
    return frame, evidence, sample_map


# ---------------------------------------------------------------- IDAT route

def idat_pairs(tar_path, work, samples):
    """Extract GSM-prefixed Grn/Red IDAT members; every kept GSM needs exactly one pair."""
    pairs = {}
    work.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as archive:
        for member in archive:
            match = re.fullmatch(r'(GSM\d+)_?(.*)_(Grn|Red|green|red)\.idat(\.gz)?', Path(member.name).name, re.I)
            if not member.isfile() or not match:
                continue
            # Some submitters name channels green/red; minfi requires _Grn/_Red.
            gsm, stem, channel = match[1], match[2], 'Grn' if match[3].lower() in {'grn', 'green'} else 'Red'
            if gsm not in samples:
                continue
            base = work / f'{gsm}_{stem}' if stem else work / gsm
            pairs.setdefault(gsm, {}).setdefault(str(base), set()).add(channel)
            target = Path(f'{base}_{channel}.idat' + ('.gz' if match[4] else ''))
            with archive.extractfile(member) as source, open(target, 'wb') as sink:
                shutil.copyfileobj(source, sink, 16 * 1024 * 1024)
    complete, problems = {}, {}
    for gsm, bases in pairs.items():
        good = [b for b, ch in bases.items() if ch == {'Grn', 'Red'}]
        if len(good) == 1:
            complete[gsm] = good[0]
        else:
            problems[gsm] = f'{len(good)} complete pairs'
    if not complete:
        raise RouteFailure('NO_COMPLETE_IDAT_PAIRS', f'{len(pairs)} GSMs with IDAT members')
    return complete, problems


def process_idats(pairs, samples, work, log):
    """Chunked minfi ssNoob per platform group; returns features x GSM beta and per-sample QC."""
    groups = {}
    for gsm, base in pairs.items():
        platform = (samples[gsm].get('platform_id') or ['unknown'])[0]
        groups.setdefault(platform, []).append((gsm, base))
    worker = REPO / 'src' / 'dnam_idat_chunk_2026_09_26.R'
    betas, qcs, receipts, unreadable = [], [], [], {}

    def run_chunk(chunk, tag):
        config = {'basenames': [b for _, b in chunk], 'sample_ids': [g for g, _ in chunk],
                  'detection_p_threshold': DETECTION_P,
                  'beta_output': str(work / f'beta-{tag}.tsv.gz'),
                  'qc_output': str(work / f'qc-{tag}.tsv'),
                  'receipt_output': str(work / f'receipt-{tag}.json')}
        path = work / f'config-{tag}.json'
        path.write_text(json.dumps(config))
        result = subprocess.run(['Rscript', '--vanilla', str(worker), str(path)],
                                capture_output=True, text=True, timeout=6 * 3600)
        log.write(f'--- R chunk {tag}\n{result.stdout[-4000:]}\n{result.stderr[-8000:]}\n')
        if result.returncode:
            raise RouteFailure('IDAT_PROCESSING_FAILED', f'{tag}: {result.stderr.strip()[-400:]}')
        betas.append(pd.read_csv(config['beta_output'], sep='\t', index_col=0, na_values=['NA'],
                                 dtype={g: 'float32' for g, _ in chunk}))
        qcs.append(pd.read_csv(config['qc_output'], sep='\t'))
        receipts.append(json.loads(Path(config['receipt_output']).read_text()))

    for platform, members in groups.items():
        for start in range(0, len(members), IDAT_CHUNK):
            chunk = members[start:start + IDAT_CHUNK]
            tag = f'{platform}-{start // IDAT_CHUNK:04d}'
            try:
                run_chunk(chunk, tag)
            except RouteFailure as failure:
                # One unreadable IDAT fails the whole chunk; isolate it by processing samples singly.
                log.write(f'chunk {tag} failed; isolating samples: {failure.detail[-200:]}\n')
                for i, member in enumerate(chunk):
                    try:
                        run_chunk([member], f'{tag}-s{i:03d}')
                    except RouteFailure as single:
                        unreadable[member[0]] = single.detail[-200:]
    if not betas:
        raise RouteFailure('IDAT_PROCESSING_FAILED', f'no readable IDAT pairs; {len(unreadable)} unreadable')
    frame = pd.concat(betas, axis=1)
    frame, collapsed = collapse_feature_ids(frame)
    evidence = {'route': 'idat', 'platform_groups': {p: len(m) for p, m in groups.items()},
                'chunks': receipts, 'features': int(len(frame)), 'epicv2_replicates_collapsed': collapsed,
                'series_samples': len(samples), 'samples_processed': frame.shape[1],
                'series_samples_absent': sorted(set(samples) - set(frame.columns)),
                'unreadable_idat_gsms': unreadable}
    return frame, evidence, pd.concat(qcs, ignore_index=True)


# ---------------------------------------------------------------- clocks

_GALLERY = None
_MEANS = None


def gallery():
    global _GALLERY
    if _GALLERY is None:
        from biolearn.model_gallery import ModelGallery
        _GALLERY = ModelGallery()
    return _GALLERY


def reference_averages(cache):
    """GSE40279 population means, reproduced by the verified pipeline reader."""
    global _MEANS
    if _MEANS is None:
        from src.pipeline_inputs_2026_09_10 import reference_means
        _MEANS, record = reference_means(None, Path(cache))
        _MEANS.attrs['record'] = record
    return _MEANS


def clock_sites():
    sites = {}
    for name in HISTORICAL_CLOCKS:
        sites[name] = [s for s in gallery().get(name).methylation_sites() if s != 'intercept']
    return sites


def predict_historical(beta, means):
    """The historical 4_compute_age.py recipe, run in sample chunks."""
    from biolearn.data_library import GeoData
    from biolearn.imputation import impute_from_standard
    outputs, errors = [], {}
    index = means.index.union(beta.index)
    for start in range(0, beta.shape[1], PREDICTION_CHUNK):
        chunk = beta.iloc[:, start:start + PREDICTION_CHUNK].astype('float64').reindex(index)
        data = GeoData(metadata=pd.DataFrame(index=chunk.columns), dnam=impute_from_standard(chunk, means))
        horvath = gallery().get('Horvathv1').predict(data)
        sex = gallery().get('SexEstimation', imputation_method='averaging').predict(data)
        data.metadata = pd.concat([horvath, sex[['predicted_sex']]], axis=1)
        data.metadata.columns = ['age', 'sex']
        results, primary = [], {}
        for name in HISTORICAL_CLOCKS:
            model = gallery().get(name, imputation_method='averaging') if name == 'SexEstimation' else gallery().get(name)
            view = data
            if name not in FULL_MATRIX_CLOCKS:
                sites = data.dnam.index.intersection(model.methylation_sites())
                view = GeoData(metadata=data.metadata, dnam=data.dnam.loc[sites])
            try:
                result = model.predict(view)
            except ValueError as error:
                # biolearn refuses clocks whose sites are neither measured nor in the reference.
                errors.setdefault(name, str(error)[:300])
                result = pd.DataFrame({name: np.nan}, index=data.dnam.columns)
            if len(result.columns) == 1:
                result.columns = [name]
            column = next((c for c in (name, 'DNAmGrimAge', 'Predicted') if c in result.columns), result.columns[0])
            primary[name] = result[column]
            results.append(result)
        outputs.append((pd.concat(results, axis=1), pd.DataFrame(primary), data.dnam))
    combined = pd.concat([o for o, _, _ in outputs])
    primary = pd.concat([p for _, p, _ in outputs])
    return combined, primary, outputs, errors


def independent_linear_check(outputs):
    """Coefficient arithmetic for the three linear clocks, on the same imputed values."""
    worst = {}
    for name in ['Horvathv1', 'Hannum', 'PhenoAge']:
        coefficients = gallery().get(name, imputation_method='none').coefficients.iloc[:, 0]
        sites = [s for s in coefficients.index if s != 'intercept']
        diffs = []
        for _, predictions, dnam in outputs:
            values = dnam.reindex(sites).to_numpy(float)
            if not np.isfinite(values).all():
                diffs.append(np.inf)
                continue
            linear = coefficients.loc[sites].to_numpy(float) @ values + float(coefficients.get('intercept', 0))
            if name == 'Horvathv1':
                t = linear + .696
                linear = np.where(t < 0, 21 * np.exp(t) - 1, 21 * t + 20)
            diffs.append(float(np.max(np.abs(predictions[name].to_numpy(float) - linear))))
        worst[name] = max(diffs)
    return worst


def recorded_ages(samples):
    ages = {}
    for gsm, fields in samples.items():
        for value in fields.get('characteristics_ch1', []):
            key, sep, token = value.partition(':')
            if sep and AGE_KEY.match(key):
                try:
                    ages[gsm] = float(token.strip().split()[0])
                except (ValueError, IndexError):
                    pass
    return pd.Series(ages, dtype=float)


def author_dnam_ages(samples):
    """Submitter-reported DNAm age (e.g. 'dnamage: 47.1'), an external check on Horvathv1."""
    values = {}
    for gsm, fields in samples.items():
        for value in fields.get('characteristics_ch1', []):
            key, sep, token = value.partition(':')
            if sep and AUTHOR_DNAM_AGE_KEY.match(key):
                try:
                    values[gsm] = float(token.strip())
                except ValueError:
                    pass
    return pd.Series(values, dtype=float)


def author_age_agreement(samples, primary, retained):
    reported = author_dnam_ages(samples)
    shared = reported.index.intersection(pd.Index(retained))
    if len(shared) < 5 or 'Horvathv1' not in primary.columns:
        return None
    ours = primary.loc[shared, 'Horvathv1'].astype(float)
    return {'n': int(len(shared)), 'pearson_r': float(np.corrcoef(reported.loc[shared], ours)[0, 1]),
            'mean_abs_diff_years': float((reported.loc[shared] - ours).abs().mean()),
            'mean_diff_ours_minus_reported': float((ours - reported.loc[shared]).mean())}


def assay_flags(samples):
    """GSMs whose own title/source/characteristics/description name an enrichment or oxidative assay."""
    flags = {}
    for gsm, fields in samples.items():
        text = ' | '.join(fields.get('title', []) + fields.get('source_name_ch1', []) +
                          fields.get('characteristics_ch1', []) + fields.get('description', []))
        match = NON_STANDARD_ASSAY.search(text) or CHIP_ASSAY.search(text)
        if match:
            flags[gsm] = match.group(0)
    return flags


def sample_qc(beta, samples=None):
    missing = beta.isna().mean(axis=0)
    flags = assay_flags(samples or {})
    reasons = []
    for gsm, fraction in missing.items():
        reason = []
        if gsm in flags:
            reason.append(f'NON_STANDARD_ASSAY:{flags[gsm]}')
        if fraction > SAMPLE_MAX_MISSING:
            reason.append(f'MISSING_FRACTION_ABOVE_{SAMPLE_MAX_MISSING:.2f}')
        reasons.append(';'.join(reason))
    status = pd.DataFrame({'missing_fraction': missing, 'reason': reasons}, index=missing.index)
    status.insert(1, 'status', np.where(status.reason != '', 'WITHHELD', 'RETAINED'))
    return status.rename_axis('sample_id')


def verdicts(beta, predictions, sites, ages, retained=None, errors=None):
    """Verifier-v3 stage 1 per clock over QC-retained samples: coverage, finiteness, non-degeneracy."""
    rows, coverage, errors = [], {}, errors or {}
    present = beta.notna()
    retained = list(beta.columns) if retained is None else list(retained)
    for name, required in sites.items():
        observed = present.reindex(required, fill_value=False).sum(axis=0) / max(len(required), 1)
        coverage[name] = observed
        if name in UTILITY_CLOCKS:
            continue
        values = predictions[name].reindex(retained) if name in predictions.columns else None
        reasons = []
        minimum = float(observed.reindex(retained).min())
        if name in errors:
            reasons.append('MISSING_REQUIRED_SITES')
        elif values is None or not np.isfinite(values.to_numpy(float)).all():
            reasons.append('NONFINITE_PREDICTION')
        elif len(values) >= 3 and float(np.nanstd(values.to_numpy(float))) < 1e-8:
            reasons.append('DEGENERATE_OUTPUT')
        if minimum < COVERAGE_PROVISIONAL:
            reasons.append(f'COVERAGE_BELOW_{COVERAGE_PROVISIONAL:.2f}')
        verdict = 'WITHHELD' if reasons else 'VERIFIED' if minimum >= COVERAGE_VERIFIED else 'PROVISIONAL'
        if verdict == 'PROVISIONAL':
            reasons.append(f'COVERAGE_{COVERAGE_PROVISIONAL:.2f}_TO_{COVERAGE_VERIFIED:.2f}')
        row = {'clock': name, 'verdict': verdict, 'reasons': reasons, 'required_sites': len(required),
               'coverage_min': minimum, 'coverage_median': float(observed.reindex(retained).median()),
               'samples_assessed': len(retained)}
        if name in errors:
            row['model_error'] = errors[name]
        shared = ages.index.intersection(values.index) if values is not None else []
        if values is not None and len(shared) >= 10 and ages.loc[shared].std() > 0:
            row['pearson_r_recorded_age'] = float(np.corrcoef(ages.loc[shared], values.loc[shared])[0, 1])
            row['age_n'] = int(len(shared))
        rows.append(row)
    return rows, pd.DataFrame(coverage)


def historical_comparison(gse, predictions):
    path = GEO_DNAM / 'clock_result' / f'{gse}.csv'
    if not path.exists():
        return {'historical_file': None}
    old = pd.read_csv(path, index_col=0)
    old = old.loc[:, ~old.columns.duplicated()]
    out = {'historical_file': str(path), 'historical_sha256': sha(path), 'clocks': {}}
    for name in ['Horvathv1', 'Hannum', 'PhenoAge', 'DunedinPACE', 'GrimAgeV2']:
        if name in old.columns and name in predictions.columns:
            shared = old.index.intersection(predictions.index)
            a, b = old.loc[shared, name].astype(float), predictions.loc[shared, name].astype(float)
            out['clocks'][name] = {'n': int(len(shared)), 'historical_sd': float(a.std()), 'repaired_sd': float(b.std()),
                                   'historical_constant': bool(a.std() < 1e-8)}
    return out


# ---------------------------------------------------------------- one study

def download(url, target, expected, log):
    part = Path(str(target) + '.partial')
    for attempt in range(1, 6):
        result = subprocess.run(['curl', '-fsSL', '--retry', '5', '--retry-delay', '20', '--connect-timeout', '60',
                                 '-C', '-', '-o', str(part), url], capture_output=True, text=True)
        log.write(f'curl attempt {attempt} rc={result.returncode} {result.stderr[-500:]}\n')
        if result.returncode == 0:
            break
        time.sleep(30 * attempt)
    else:
        raise RouteFailure('DOWNLOAD_FAILED', f'{url}: {result.stderr[-300:]}')
    size = part.stat().st_size
    if expected and size != expected:
        raise RouteFailure('DOWNLOAD_SIZE_MISMATCH', f'{size} bytes, listing declared {expected}')
    part.rename(target)
    return {'url': url, 'bytes': size, 'sha256': sha(target), 'downloaded': now()}


def routes_for(row):
    """Author tables first (largest beta, then M-value), then the IDAT archive."""
    tables, idat = [], None
    for item in str(row['files']).split(' | '):
        match = re.fullmatch(r'(.+?)\[(\d+)(?:,([^\]]+))?\]', item.strip())
        if not match:
            continue
        name, size, kind = match[1], int(match[2]), (match[3] or '').strip()
        if kind in {'beta', 'mvalue'}:
            tables.append((0 if kind == 'beta' else 1, -size, name, size, kind))
        elif name.endswith('_RAW.tar'):
            idat = (name, size)
    routes = [{'type': 'table', 'file': n, 'bytes': s, 'kind': k} for _, _, n, s, k in sorted(tables)]
    if idat and row['class'] == 'IDAT_AVAILABLE':
        routes.append({'type': 'idat', 'file': idat[0], 'bytes': idat[1], 'kind': 'idat'})
    return routes


def suppl_url(gse, name):
    return f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/suppl/{name}'


def run_study(root, gse, keep_raw=False):
    root = Path(root)
    row = pd.read_csv(root / 'inventory.csv').set_index('gse').loc[gse]
    out = root / 'studies' / gse
    if (out / 'status.json').exists():
        return json.loads((out / 'status.json').read_text())
    if out.exists():
        attempts = root / 'attempts' / gse
        attempts.mkdir(parents=True, exist_ok=True)
        out.rename(attempts / f'attempt-{len(list(attempts.iterdir())) + 1:03d}')
    out.mkdir(parents=True)
    work = root / 'work' / gse
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    started = time.time()
    status = {'schema': SCHEMA, 'gse': gse, 'class': row['class'], 'platform': row['platform'],
              'n_samples_geo': int(row['n_samples']), 'table_linked_32': bool(row['table_linked_32']),
              'historical_table_rows': int(row['historical_table_rows']), 'attempted_routes': []}
    with open(out / 'log.txt', 'w') as log:
        try:
            samples, series_files = series_samples(gse)
            status['series_matrix'] = series_files
            beta = None
            for route in routes_for(row):
                attempt = {**route, 'url': suppl_url(gse, route['file'])}
                status['attempted_routes'].append(attempt)
                target = work / route['file']
                try:
                    t0 = time.time()
                    attempt['download'] = download(attempt['url'], target, route['bytes'], log)
                    attempt['download']['seconds'] = round(time.time() - t0, 1)
                    t0 = time.time()
                    if route['type'] == 'table':
                        beta, evidence, sample_map = read_beta_table(target, samples, kind=route['kind'])
                        qc = None
                    else:
                        pairs, problems = idat_pairs(target, work / 'idat', samples)
                        if not keep_raw:
                            target.unlink()
                        beta, evidence, qc = process_idats(pairs, samples, work / 'idat', log)
                        evidence['incomplete_idat_gsms'] = problems
                        sample_map = pd.DataFrame([{'column': Path(b).name, 'gsm': g, 'rule': 'idat_gsm_prefix'}
                                                   for g, b in pairs.items()])
                    attempt['result'] = 'PASS'
                    attempt['parse_seconds'] = round(time.time() - t0, 1)
                    break
                except RouteFailure as failure:
                    attempt['result'] = 'FAIL'
                    attempt['failure'] = {'code': failure.code, 'detail': failure.detail}
                    log.write(f'route failed {failure}\n')
                    beta = None
                finally:
                    if not keep_raw:
                        for leftover in [target, Path(str(target) + '.partial')]:
                            leftover.unlink(missing_ok=True)
                    shutil.rmtree(work / 'idat', ignore_errors=True)
            if beta is None:
                codes = [a.get('failure', {}).get('code') for a in status['attempted_routes']]
                status.update(state='WITHHELD', reason='NO_ROUTE_PRODUCED_BETA' if codes else 'NO_USABLE_SOURCE',
                              route_failures=codes)
                return status
            status['input'] = evidence
            sample_map.to_csv(out / 'sample_map.csv', index=False)
            if qc is not None:
                qc.to_csv(out / 'idat_qc.csv', index=False)
            missing = beta.isna().mean(axis=0)
            status['input']['missing_fraction'] = {'min': float(missing.min()), 'median': float(missing.median()),
                                                   'max': float(missing.max())}
            means = reference_averages(root / 'reference-cache')
            status['reference'] = {k: v for k, v in means.attrs['record'].items() if k != 'created_unix'}
            qc_status = sample_qc(beta, samples)
            qc_status.to_csv(out / 'sample_status.csv')
            retained = qc_status.index[qc_status.status == 'RETAINED']
            status['sample_qc'] = {'threshold_missing_fraction': SAMPLE_MAX_MISSING,
                                   'retained': int(len(retained)), 'withheld': int((qc_status.status == 'WITHHELD').sum()),
                                   'non_standard_assay': int(qc_status.reason.str.contains('NON_STANDARD_ASSAY').sum())}
            if not len(retained):
                raise RouteFailure('ALL_SAMPLES_FAILED_QC', f'every sample exceeds {SAMPLE_MAX_MISSING} missing')
            t0 = time.time()
            predictions, primary, outputs, clock_errors = predict_historical(beta, means)
            status['predict_seconds'] = round(time.time() - t0, 1)
            status['clock_errors'] = clock_errors
            predictions.rename_axis('sample_id').to_csv(out / 'clock_predictions.csv')
            status['independent_linear_max_abs_diff'] = independent_linear_check(outputs)
            sites = clock_sites()
            rows, coverage = verdicts(beta, primary, sites, recorded_ages(samples), retained, clock_errors)
            coverage.rename_axis('sample_id').to_csv(out / 'clock_coverage.csv')
            status['clock_verdicts'] = rows
            status['historical_comparison'] = historical_comparison(gse, primary)
            status['author_reported_dnam_age_vs_horvathv1'] = author_age_agreement(samples, primary, retained)
            keep = sorted(set().union(*[set(v) for k, v in sites.items() if k not in UTILITY_CLOCKS]))
            beta.reindex(keep).dropna(how='all').to_csv(out / 'clock_cpg_beta.csv.gz', float_format='%.6g',
                                                         compression={'method': 'gzip', 'compresslevel': 6})
            linear_ok = all(v <= 1e-6 for v in status['independent_linear_max_abs_diff'].values())
            counts = pd.Series([r['verdict'] for r in rows]).value_counts().to_dict()
            status.update(state='REPAIRED' if linear_ok else 'WITHHELD',
                          reason=None if linear_ok else 'INDEPENDENT_LINEAR_CHECK_FAILED',
                          verdict_counts=counts, samples_predicted=int(predictions.shape[0]))
            return status
        except RouteFailure as failure:
            status.update(state='WITHHELD', reason=failure.code, detail=failure.detail)
            return status
        except Exception as error:
            status.update(state='ERROR', reason=type(error).__name__, detail=str(error)[:1000],
                          traceback=traceback.format_exc()[-4000:])
            return status
        finally:
            status['elapsed_seconds'] = round(time.time() - started, 1)
            status['finished'] = now()
            status['software'] = {'python': sys.version.split()[0], 'biolearn': _biolearn_version(),
                                  'module_sha256': sha(__file__),
                                  'r_worker_sha256': sha(REPO / 'src' / 'dnam_idat_chunk_2026_09_26.R')}
            for path in sorted(out.iterdir()):
                if path.is_file() and path.name not in {'status.json', 'log.txt'}:
                    status.setdefault('outputs', {})[path.name] = sha(path)
            if status.get('state') != 'ERROR' or not keep_raw:
                shutil.rmtree(work, ignore_errors=True)
            try:
                write_json(out / 'status.json', status)
            except Exception as error:
                write_json(out / 'status.json', {'schema': SCHEMA, 'gse': gse, 'state': 'ERROR',
                                                 'reason': 'STATUS_WRITE_FAILED', 'detail': repr(error)[:500],
                                                 'partial_state': status.get('state')})


def _biolearn_version():
    try:
        from importlib.metadata import version
        return version('biolearn')
    except Exception:
        return None


# ---------------------------------------------------------------- orchestrator

class DiskBudget:
    """Reserve bytes before a study starts; refuse when the filesystem is near full."""

    def __init__(self, root, max_inflight, min_free):
        self.root, self.max_inflight, self.min_free = Path(root), max_inflight, min_free
        self.reserved, self.lock = 0, threading.Condition()

    def need(self, row):
        routes = routes_for(row)
        largest = max((r['bytes'] * (2 if r['type'] == 'idat' else 1) for r in routes), default=0)
        return largest + 2 * 1024 ** 3

    def acquire(self, amount, stop):
        with self.lock:
            while not stop.is_set():
                free = shutil.disk_usage(self.root).free
                fits = self.reserved == 0 or self.reserved + amount <= self.max_inflight
                if fits and free - amount >= self.min_free:
                    self.reserved += amount
                    return True
                self.lock.wait(timeout=60)
            return False

    def release(self, amount):
        with self.lock:
            self.reserved -= amount
            self.lock.notify_all()


def orchestrate(args):
    root = Path(args.root)
    root.mkdir(parents=True, exist_ok=True)
    inventory = root / 'inventory.csv'
    if not inventory.exists():
        shutil.copyfile(args.inventory, inventory)
        write_json(root / 'run.json', {'schema': SCHEMA, 'created': now(), 'inventory_source': str(args.inventory),
                                       'inventory_sha256': sha(inventory), 'module_sha256': sha(__file__),
                                       'thresholds': {'coverage_verified': COVERAGE_VERIFIED,
                                                      'coverage_provisional': COVERAGE_PROVISIONAL,
                                                      'detection_p': DETECTION_P}})
    elif sha(inventory) != json.loads((root / 'run.json').read_text())['inventory_sha256']:
        raise SystemExit('Inventory in the run root changed; use a new root')
    rows = pd.read_csv(inventory)
    rows = rows[rows['class'].isin(['IDAT_AVAILABLE', 'BETA_TABLE_AVAILABLE'])]
    if args.only:
        rows = rows[rows.gse.isin(args.only.split(','))]
    rows = rows.assign(priority=~rows.table_linked_32.astype(bool)).sort_values(['priority', 'best_route_bytes'])
    if args.limit:
        rows = rows.head(args.limit)
    for gse in (args.retry.split(',') if args.retry else []):
        previous = root / 'studies' / gse
        if previous.exists():
            attempts = root / 'attempts' / gse
            attempts.mkdir(parents=True, exist_ok=True)
            previous.rename(attempts / f'attempt-{len(list(attempts.iterdir())) + 1:03d}')
            append_event(root, {'event': 'study_archived_for_retry', 'gse': gse})
    done = {p.parent.name for p in (root / 'studies').glob('*/status.json')}
    queue = [r for _, r in rows.iterrows() if r.gse not in done]
    means = reference_averages(root / 'reference-cache')
    equivalence = root / 'subset-equivalence-2026-09-26.json'
    if not equivalence.exists() or json.loads(equivalence.read_text()).get('module_sha256') != sha(__file__):
        if check_subset_equivalence(root)['status'] != 'PASS':
            raise SystemExit('Per-clock subset equivalence failed; see subset-equivalence-2026-09-26.json')
    append_event(root, {'event': 'run_start', 'pending': len(queue), 'already_done': len(done),
                        'reference_means_sha256': means.attrs['record']['means_sha256']})
    budget = DiskBudget(root, args.max_inflight_gb * 1024 ** 3, args.min_free_gb * 1024 ** 3)
    stop = threading.Event()

    def one(row):
        amount = budget.need(row)
        if not budget.acquire(amount, stop):
            return
        append_event(root, {'event': 'study_start', 'gse': row.gse, 'reserved_bytes': amount})
        try:
            command = [sys.executable, '-m', 'src.dnam_stream_repair_2026_09_26', 'study', '--root', str(root),
                       '--gse', row.gse] + (['--keep-raw'] if args.keep_raw else [])
            result = subprocess.run(command, cwd=REPO, capture_output=True, text=True, timeout=args.study_timeout_hours * 3600)
            status = root / 'studies' / row.gse / 'status.json'
            state = json.loads(status.read_text()).get('state') if status.exists() else 'CRASHED'
            if not status.exists():
                (root / 'crashes').mkdir(exist_ok=True)
                (root / 'crashes' / f'{row.gse}-{int(time.time())}.txt').write_text(result.stdout[-5000:] + result.stderr[-20000:])
            append_event(root, {'event': 'study_end', 'gse': row.gse, 'state': state, 'rc': result.returncode})
        except subprocess.TimeoutExpired:
            append_event(root, {'event': 'study_timeout', 'gse': row.gse})
        finally:
            budget.release(amount)

    try:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(one, queue))
    except KeyboardInterrupt:
        stop.set()
        raise
    append_event(root, {'event': 'run_end'})
    print(json.dumps(summarize(root), indent=2))


def summarize(root):
    root = Path(root)
    records = [json.loads(p.read_text()) for p in sorted((root / 'studies').glob('*/status.json'))]
    table = []
    for s in records:
        counts = s.get('verdict_counts', {})
        table.append({'gse': s['gse'], 'state': s.get('state'), 'reason': s.get('reason'),
                      'route': (s.get('input') or {}).get('route'), 'table_linked_32': s.get('table_linked_32'),
                      'historical_table_rows': s.get('historical_table_rows'),
                      'samples_predicted': s.get('samples_predicted'), 'VERIFIED': counts.get('VERIFIED', 0),
                      'PROVISIONAL': counts.get('PROVISIONAL', 0), 'WITHHELD': counts.get('WITHHELD', 0),
                      'elapsed_seconds': s.get('elapsed_seconds')})
    frame = pd.DataFrame(table)
    if not frame.empty:
        frame.to_csv(root / 'summary.csv', index=False)
    summary = {'studies_finished': len(records),
               'states': frame['state'].value_counts().to_dict() if not frame.empty else {},
               'withheld_reasons': frame.loc[frame.state != 'REPAIRED', 'reason'].value_counts().to_dict() if not frame.empty else {},
               'clock_verdicts': {k: int(frame[k].sum()) for k in ('VERIFIED', 'PROVISIONAL', 'WITHHELD')} if not frame.empty else {},
               'linked_32_finished': int(frame['table_linked_32'].sum()) if not frame.empty else 0,
               'samples_predicted': int(frame['samples_predicted'].fillna(0).sum()) if not frame.empty else 0}
    write_json(root / 'summary.json', summary)
    return summary


def check_subset_equivalence(root, samples=32, seed=1):
    """Software check on SYNTHETIC beta values: per-clock site subsets must equal full-matrix output."""
    from biolearn.data_library import GeoData
    from biolearn.imputation import impute_from_standard
    means = reference_averages(Path(root) / 'reference-cache')
    rng = np.random.default_rng(seed)
    index = means.index[:460000]
    beta = pd.DataFrame(rng.beta(2, 2, (len(index), samples)), index=index, columns=[f'S{i}' for i in range(samples)])
    beta[rng.random(beta.shape) < 0.03] = np.nan
    chunk = beta.reindex(means.index.union(beta.index))
    data = GeoData(metadata=pd.DataFrame(index=chunk.columns), dnam=impute_from_standard(chunk, means))
    age = gallery().get('Horvathv1').predict(data)
    sex = gallery().get('SexEstimation', imputation_method='averaging').predict(data)
    data.metadata = pd.concat([age, sex[['predicted_sex']]], axis=1)
    data.metadata.columns = ['age', 'sex']
    rows = []
    for name in HISTORICAL_CLOCKS:
        model = gallery().get(name, imputation_method='averaging') if name == 'SexEstimation' else gallery().get(name)
        view = GeoData(metadata=data.metadata, dnam=data.dnam.loc[data.dnam.index.intersection(model.methylation_sites())])
        try:
            full, part = model.predict(data), model.predict(view)
        except ValueError as error:
            rows.append({'clock': name, 'status': 'MODEL_ERROR', 'error': str(error)[:200]})
            continue
        diff = float(np.nanmax(np.abs(full.select_dtypes('number').to_numpy(float) - part.select_dtypes('number').to_numpy(float))))
        rows.append({'clock': name, 'max_abs_diff': diff, 'identical': diff == 0.0,
                     'uses_full_matrix_in_runner': name in FULL_MATRIX_CLOCKS})
    unsafe = [r['clock'] for r in rows if r.get('identical') is False and r['clock'] not in FULL_MATRIX_CLOCKS]
    value = {'check': 'per-clock subset equivalence', 'data': f'SYNTHETIC beta(2,2), {samples} samples, 3% missing, seed {seed}',
             'biolearn': _biolearn_version(), 'module_sha256': sha(__file__), 'clocks': rows,
             'status': 'PASS' if not unsafe else 'FAIL', 'unsafe_subset_clocks': unsafe}
    write_json(Path(root) / 'subset-equivalence-2026-09-26.json', value)
    return value


def skip_studies(root, gses, reason, detail):
    """Record studies deliberately not processed, with the reason, instead of leaving them silently absent."""
    root = Path(root)
    inventory = pd.read_csv(root / 'inventory.csv').set_index('gse')
    for gse in gses:
        out = root / 'studies' / gse
        if (out / 'status.json').exists():
            raise SystemExit(f'{gse} already has a final status; archive it with --retry before skipping')
        if out.exists():
            attempts = root / 'attempts' / gse
            attempts.mkdir(parents=True, exist_ok=True)
            out.rename(attempts / f'attempt-{len(list(attempts.iterdir())) + 1:03d}')
        out.mkdir(parents=True)
        row = inventory.loc[gse]
        status = {'schema': SCHEMA, 'gse': gse, 'class': row['class'], 'platform': row['platform'],
                  'n_samples_geo': int(row['n_samples']), 'table_linked_32': bool(row['table_linked_32']),
                  'historical_table_rows': int(row['historical_table_rows']), 'best_route': row['best_route'],
                  'best_route_bytes': int(row['best_route_bytes']), 'state': 'WITHHELD', 'reason': reason,
                  'detail': detail, 'attempted_routes': [], 'finished': now()}
        write_json(out / 'status.json', status)
        append_event(root, {'event': 'study_skipped', 'gse': gse, 'reason': reason})
        print(gse, reason)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    commands = parser.add_subparsers(dest='command', required=True)
    run = commands.add_parser('run')
    run.add_argument('--inventory', type=Path, required=True)
    run.add_argument('--root', type=Path, required=True)
    run.add_argument('--only')
    run.add_argument('--limit', type=int)
    run.add_argument('--workers', type=int, default=4)
    run.add_argument('--max-inflight-gb', type=int, default=200)
    run.add_argument('--min-free-gb', type=int, default=150)
    run.add_argument('--study-timeout-hours', type=float, default=12)
    run.add_argument('--keep-raw', action='store_true')
    run.add_argument('--retry', help='comma-separated GSEs whose finished result is archived and rerun')
    study = commands.add_parser('study')
    study.add_argument('--root', type=Path, required=True)
    study.add_argument('--gse', required=True)
    study.add_argument('--keep-raw', action='store_true')
    check = commands.add_parser('check-subset-equivalence')
    check.add_argument('--root', type=Path, required=True)
    skip = commands.add_parser('skip', help='record studies deliberately not processed')
    skip.add_argument('--root', type=Path, required=True)
    skip.add_argument('--gse', required=True, help='comma-separated GSEs')
    skip.add_argument('--reason', required=True)
    skip.add_argument('--detail', required=True)
    summary = commands.add_parser('summary')
    summary.add_argument('--root', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == 'run':
        orchestrate(args)
    elif args.command == 'study':
        status = run_study(args.root, args.gse, keep_raw=args.keep_raw)
        print(json.dumps({k: status.get(k) for k in ('gse', 'state', 'reason', 'verdict_counts', 'elapsed_seconds')}))
    elif args.command == 'skip':
        skip_studies(args.root, args.gse.split(','), args.reason, args.detail)
    elif args.command == 'check-subset-equivalence':
        value = check_subset_equivalence(args.root)
        print(value['status'], value['unsafe_subset_clocks'])
        if value['status'] != 'PASS':
            raise SystemExit(1)
    else:
        print(json.dumps(summarize(args.root), indent=2))


if __name__ == '__main__':
    main()
