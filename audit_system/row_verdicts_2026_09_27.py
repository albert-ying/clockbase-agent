"""Verifier v3 row verdicts: VERIFIED / PROVISIONAL / WITHHELD for every atlas row.

Inputs: a modality correction table and the membership-recovery rows.jsonl for it. Every check is
deterministic and every triggered check leaves a reason code. Missing evidence is itself a
reason code, so there is no undecided state. See ROW_VERDICT_SPEC_2026-09-27.md.

  python -m audit_system.row_verdicts_2026_09_27 --modality dnam --membership DIR --out DIR
"""
from audit_system.adjudication_files import decision_files
import argparse
import csv
import gzip
import hashlib
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system import paths
from audit_system.membership_recovery_2026_09_27 import CORRECTIONS, GEO, MODALITIES, sha256

SCHEMA = 'clockbase-row-verdict-v1'
MIN_ARM = 3
COVERAGE_HARD, COVERAGE_SOFT = 0.60, 0.80
# Correlation of a DNAm sample's beta profile with the GSE40279 population-mean methylome
# (input_integrity_2026_09_30.py; 218,969 screened samples). Below 0.10 (757 samples in 69 studies; a random
# permutation of 450K values gives |r| < 0.01) the values are detached from their probe identities or inverted
# (index case GSE73950): hard. From 0.10 to 0.50 (407 samples) the profiles are real but atypical for clock
# inputs: TAB-array and oxBS assays, GoldenGate arrays, and demethylated or DNMT-knockout cell lines: soft.
PROFILE_DETACHED_R, PROFILE_ATYPICAL_R = 0.10, 0.50
# DNAm clocks in the atlas: training tissue domain and output target.
DNAM_CLOCKS = {
    'Horvathv1': {'domain': 'pan-tissue', 'target': 'chronological age (years)'},
    'PhenoAge': {'domain': 'blood', 'target': 'phenotypic age (years)'},
    # GrimAge takes chronological age and sex as inputs; the historical recipe (GEO_data/src/4_compute_age.py lines
    # 227-235) passed each sample's Horvathv1 prediction as age, so its outputs embed the Horvathv1 estimate, and
    # passed sex as 'Female'/'Male' strings, which biolearn's GrimageModel reads as male for every sample.
    'GrimAgeV2': {'domain': 'blood', 'target': 'mortality-risk surrogate (years)', 'age_input': 'Horvathv1 prediction'},
    'DunedinPACE': {'domain': 'blood', 'target': 'pace of aging (years per year)'},
    'YingDamAge': {'domain': 'blood', 'target': 'damage-associated age (years)'},
}
# Output of the RNA mortality clocks. The mouse model's unit is established; the human multi-species model's published
# identity and unit are not (evidence/rna-application-profile-2026-09-17), so its records say so.
RNA_TARGETS = {'mm_rna': 'mouse multi-tissue mortality clock: cohort-centered log10 mortality-hazard score (not years)',
               'hs_rna': 'multi-species multi-tissue mortality clock: score unit not established'}
BLOOD = re.compile(r'blood|pbmc|leukocyte|lymphocyte|monocyte|neutrophil|granulocyte|t[- ]?cell|b[- ]?cell|'
                   r'cd4|cd8|cd14|cd19|nk cell|buffy|whole blood|plasma cell|erythro|platelet|macrophage', re.I)
# Underscore counts as a separator (titles such as tissue_1_oxBS, GSE73895).
NON_STANDARD_ASSAY = re.compile(
    r'5hmc|hmedip|medip|hmc[-_ ]?ip|immuno[- ]?precipitat|(?<![a-z0-9])ip[-_ ]?(enriched|fraction)(?![a-z0-9])|'
    r'(?<![a-z0-9])oxbs(?![a-z0-9])|oxidative bisulfite|(?<![a-z0-9])tab[-_ ]?(seq|array)(?![a-z0-9])|'
    r'(?<![a-z0-9])mbd[-_ ]?(seq|cap|enrich)|(?<![a-z0-9])chip[-_ ]?(seq|on[-_ ]chip|enrich[a-z]*|input|dna|ip)(?![a-z0-9])', re.I)
STOP = {'vs', 'versus', 'v', 'and', 'the', 'of', 'in', 'with', 'without', 'group', 'control', 'controls', 'ctrl',
        'sample', 'samples', 'cells', 'cell', 'treated', 'treatment', 'untreated', 'a', 'an', 'to', 'for', 'at',
        'on', 'by', 'from', 'nan', 'none', 'na'}


def tokens(text):
    return {t for t in re.findall(r'[a-z0-9]+', str(text).lower()) if t not in STOP and len(t) > 1}


def label_orientation(comparison, condition, case_level, control_level, stratum):
    """SUPPORTED / INVERTED / UNCONFIRMED from literal token overlap of recovered levels and labels."""
    parts = re.split(r'\s+(?:vs\.?|versus|v\.?)\s+|_vs_|_VS_', str(comparison), maxsplit=1, flags=re.I)
    left = tokens(parts[0]) | tokens(condition)
    right = tokens(parts[1]) if len(parts) == 2 else set()
    case, control = tokens(case_level), tokens(control_level)
    if control_level == '<rest>':
        control = set()
    forward = len(case & left) + len(control & right)
    backward = len(case & right) + len(control & left)
    if forward > backward:
        return 'SUPPORTED'
    if backward > forward:
        return 'INVERTED'
    return 'UNCONFIRMED'


def dnam_coverage(gse, gsms, clock):
    path = GEO / 'GEO_hs_dnam/missing_percentage' / f'{gse}.csv'
    if not path.exists():
        return None
    missing = pd.read_csv(path, index_col=0)
    if clock not in missing.columns:
        return None
    values = 1 - pd.to_numeric(missing[clock].reindex(gsms), errors='coerce')
    return float(values.min()) if values.notna().all() else None


def series_sample_text(gse):
    """GSM -> concatenated title/source/characteristics/description from local series matrices."""
    text = {}
    for path in sorted((GEO / 'GEO_hs_dnam/series_matrix').glob(f'{gse}*_series_matrix*.txt.gz')):
        fields = {}
        with gzip.open(path, 'rt', newline='', errors='replace') as handle:
            for record in csv.reader(handle, delimiter='\t'):
                if record and record[0] == '!series_matrix_table_begin':
                    break
                if record and record[0] in ('!Sample_geo_accession', '!Sample_title', '!Sample_source_name_ch1',
                                            '!Sample_characteristics_ch1', '!Sample_description'):
                    fields.setdefault(record[0], []).append(record[1:])
        ids = fields.get('!Sample_geo_accession', [[]])[0]
        for i, gsm in enumerate(ids):
            text[gsm] = ' | '.join(row[i] for key, rows in fields.items() if key != '!Sample_geo_accession'
                                   for row in rows if i < len(row))
    return text


def sample_tissue(modality, gse, gsms):
    folder = MODALITIES[modality]['folder']
    path = GEO / folder / 'processed_metadata' / f'{gse}.csv'
    if not path.exists():
        return ''
    meta = pd.read_csv(path, dtype=str).set_index('GSM').reindex(gsms)
    columns = [c for c in ('source_name', 'harmonized.specimen_part', 'harmonized.cell_line') if c in meta.columns]
    return ' | '.join(meta[columns].fillna('').astype(str).agg(' '.join, axis=1))


def verdict_for(row, rec, context):
    """Return (verdict, reasons, evidence) for one atlas row."""
    hard, soft, info, evidence = [], [], [], {}
    modality = context['modality']
    gse = row['gse_id']
    if gse in context['input_artifact']:
        hard.append(context['input_artifact'][gse])
    if row.get('fdr_below_raw_p'):
        hard.append('FDR_BELOW_RAW_P')
    status = (rec or {}).get('status')
    evidence['membership_status'] = status
    evidence['membership_source'] = (rec or {}).get('source')
    if status != 'UNIQUE':
        hard.append({'NOT_TWO_SAMPLE_T': 'TEST_NOT_RECONSTRUCTIBLE', 'STORED_P_UNUSABLE': 'STORED_P_UNUSABLE',
                     'NO_PREDICTIONS': 'PREDICTIONS_UNAVAILABLE', 'AMBIGUOUS': 'MEMBERSHIP_AMBIGUOUS',
                     'ONE_SAMPLE_TEST': 'ONE_SAMPLE_TEST_NOT_CONTRAST',
                     'CONTINUOUS_TERM': 'CONTINUOUS_TERM_NOT_GROUP_CONTRAST',
                     'ARMS_UNRESOLVED': 'ARMS_UNRESOLVED'}
                    .get(status, 'MEMBERSHIP_UNRECOVERED'))
    else:
        case, control = rec['case_gsms'], rec['control_gsms']
        members = case + control
        evidence.update(n_case=rec['n_case'], n_control=rec['n_control'], test=rec.get('test'),
                        split=[rec.get('field'), rec.get('case_level'), rec.get('control_level'), rec.get('stratum')],
                        prediction_source=rec.get('prediction_source'), link_traces=rec.get('link_traces'))
        info.append('MEMBERSHIP_RECOVERED')
        if rec.get('kind') == 'model':
            info.append('COVARIATE_MODEL')
            evidence['covariates'] = rec.get('covariates')
        if rec.get('misfiled'):
            info.append('OUTPUT_MISFILED')
        if min(rec['n_case'], rec['n_control']) < MIN_ARM:
            hard.append('ARM_TOO_SMALL')
        designed = rec.get('source') == 'design_reconstruction'
        effect_matches = set(rec.get('effect_definitions_matched') or [])
        evidence['effect_definitions_matched'] = sorted(effect_matches)
        if designed:
            # New arms from dual-model design reconstruction: there is no historical value to reproduce.
            info.append('DESIGN_RECONSTRUCTED')
        elif not effect_matches:
            hard.append('EFFECT_NOT_REPRODUCED')
        elif effect_matches == {'model_direction_sign'}:
            soft.append('EFFECT_SIGN_ONLY')
        adjudicated = context.get('adjudication', {}).get(row.get('row_ordinal'))
        if designed:
            orientation = 'SUPPORTED'  # both models assigned the arms to the labels
        elif adjudicated in ('LABEL_SUPPORTED', 'LABEL_INVERTED_CONFIRMED'):
            orientation = 'SUPPORTED' if adjudicated == 'LABEL_SUPPORTED' else 'INVERTED'
            info.append('LABEL_ADJUDICATED_DUAL_MODEL')
        elif rec.get('case_level') is None:
            orientation = 'UNCONFIRMED'
        else:
            orientation = label_orientation(row['comparison_name'], row['condition_name'], rec['case_level'],
                                            rec['control_level'], rec.get('stratum'))
        evidence['label_orientation'] = orientation
        if orientation == 'INVERTED':
            hard.append('LABEL_ORIENTATION_INVERTED')
        elif orientation == 'UNCONFIRMED':
            soft.append('LABEL_UNCONFIRMED')
        if modality == 'dnam':
            texts = context['sample_text'](gse)
            flagged = [g for g in members if NON_STANDARD_ASSAY.search(texts.get(g, ''))]
            if flagged:
                hard.append('NON_STANDARD_ASSAY')
                evidence['non_standard_assay_gsms'] = flagged[:20]
            integrity = context.get('integrity')
            if integrity is not None:
                screened = {g: integrity[g] for g in members if g in integrity}
                evidence['profile_r_reference_min'] = min(screened.values()) if screened else None
                detached = [g for g, r in screened.items() if r < PROFILE_DETACHED_R]
                if detached:
                    hard.append('INPUT_PROFILE_DETACHED')
                    evidence['detached_gsms'] = detached
                atypical = [g for g, r in screened.items() if PROFILE_DETACHED_R <= r < PROFILE_ATYPICAL_R]
                if atypical:
                    soft.append('PROFILE_ATYPICAL')
                    evidence['atypical_profile_gsms'] = atypical
                if len(screened) < len(members):
                    info.append('PROFILE_UNSCREENED')
            coverage = dnam_coverage(gse, members, row['clock'])
            evidence['clock_coverage_min'] = coverage
            if coverage is None:
                hard.append('COVERAGE_UNAVAILABLE')
            elif coverage < COVERAGE_HARD:
                hard.append('CLOCK_COVERAGE_LOW')
            elif coverage < COVERAGE_SOFT:
                soft.append('CLOCK_COVERAGE_MID')
            clock = DNAM_CLOCKS.get(row['clock'])
            if clock is None:
                hard.append('CLOCK_UNREGISTERED')
            else:
                evidence['clock_target'] = clock['target']
                if clock['domain'] == 'blood' and not BLOOD.search(sample_tissue(modality, gse, members)):
                    soft.append('TISSUE_OUT_OF_DOMAIN')
                if clock.get('age_input'):
                    soft.append('CLOCK_AGE_INPUT_FROM_HORVATH')
        else:
            evidence['clock_target'] = RNA_TARGETS[modality]
    reasons = hard + soft + info
    verdict = 'WITHHELD' if hard else 'PROVISIONAL' if soft else 'VERIFIED'
    return verdict, reasons, evidence


def run(modality, membership_dir, out_dir, adjudication_dirs=(), integrity_dir=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    table = CORRECTIONS / f"2026-09-08-{MODALITIES[modality]['table']}-atlas_WORKING_NOT_VALIDATED.csv"
    atlas = pd.read_csv(table, low_memory=False)
    p = pd.to_numeric(atlas.pvalue, errors='coerce')
    fdr = pd.to_numeric(atlas.get('fdr'), errors='coerce')
    atlas['fdr_below_raw_p'] = (fdr < p * (1 - 1e-10)) & fdr.notna() & p.notna()
    membership_path = Path(membership_dir) / 'rows.jsonl'
    recovered = {}
    with open(membership_path) as handle:
        for line in handle:
            record = json.loads(line)
            recovered[record['row_ordinal']] = record
    input_artifact = {}
    if modality == 'dnam':
        census = paths.EMPTY_MATRIX_LIST
        for line in open(census):
            input_artifact[json.loads(line)['study_id']] = 'INPUT_EMPTY'
    else:
        qc = pd.read_csv(GEO / 'result/geo_rna_imputation_qc_all_2026-06-14.csv', low_memory=False)
        species = 'human' if modality == 'hs_rna' else 'mouse'
        flagged = qc[(qc.species == species) & qc.artifact_flag.astype(str).eq('True')]
        input_artifact = {g: 'RNA_INPUT_ARTIFACT' for g in flagged.GSE_ID}
    cache = {}

    def sample_text(gse):
        if gse not in cache:
            cache.clear() if len(cache) > 256 else None
            cache[gse] = series_sample_text(gse)
        return cache[gse]

    adjudication = {}
    for folder in adjudication_dirs:
        for path in decision_files(folder):
            for line in open(path):
                decision = json.loads(line)
                if decision.get('decision') in ('LABEL_SUPPORTED', 'LABEL_INVERTED_CONFIRMED', 'LABEL_UNCONFIRMED'):
                    adjudication[decision['row_ordinal']] = decision['decision']
    integrity = None
    if modality == 'dnam' and integrity_dir:
        screen = pd.read_csv(Path(integrity_dir) / 'samples.csv', usecols=['gsm', 'r_reference']).dropna()
        integrity = dict(zip(screen.gsm.astype(str), screen.r_reference))  # samples without values: coverage checks
    context = {'modality': modality, 'input_artifact': input_artifact, 'sample_text': sample_text,
               'adjudication': adjudication, 'integrity': integrity}
    records = []
    for ordinal, row in enumerate(atlas.to_dict('records')):
        verdict, reasons, evidence = verdict_for({**row, 'row_ordinal': ordinal}, recovered.get(ordinal), context)
        records.append({'row_ordinal': ordinal, 'gse_id': row['gse_id'], 'comparison_name': row['comparison_name'],
                        'condition_name': row['condition_name'], 'clock': row['clock'], 'verdict': verdict,
                        'reasons': reasons, 'evidence': evidence})
    frame = pd.DataFrame(records)
    frame.to_json(out_dir / 'row_verdicts.jsonl', orient='records', lines=True)
    reason_counts = frame.explode('reasons').reasons.value_counts().to_dict()
    manifest = {'schema': SCHEMA, 'modality': modality, 'table': str(table), 'table_sha256': sha256(table),
                'membership_rows': str(membership_path), 'membership_sha256': sha256(membership_path),
                'module_sha256': sha256(__file__), 'rows': len(frame),
                'verdicts': frame.verdict.value_counts().to_dict(), 'reason_counts': reason_counts,
                'thresholds': {'min_arm': MIN_ARM, 'coverage_hard': COVERAGE_HARD, 'coverage_soft': COVERAGE_SOFT,
                               'profile_detached_r': PROFILE_DETACHED_R, 'profile_atypical_r': PROFILE_ATYPICAL_R},
                'integrity_dir': str(integrity_dir) if integrity_dir else None,
                'integrity_sha256': sha256(Path(integrity_dir) / 'samples.csv') if integrity_dir else None,
                'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                'verdicts_sha256': sha256(out_dir / 'row_verdicts.jsonl')}
    (out_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps({'verdicts': manifest['verdicts'], 'reasons': reason_counts}, indent=1))
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--modality', choices=sorted(MODALITIES), required=True)
    parser.add_argument('--membership', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--adjudication', type=Path, action='append', default=[])
    parser.add_argument('--integrity', type=Path, help='input_integrity_2026_09_30 output (DNAm)')
    args = parser.parse_args(argv)
    run(args.modality, args.membership, args.out, args.adjudication, args.integrity)


if __name__ == '__main__':
    main()
