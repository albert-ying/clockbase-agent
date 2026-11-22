"""Public, per-study export of every verifier-v3 decision and repair path.

For each study the export lists every atlas row with: the historical values, the verdict and its
reason codes, the repair actions in order, the case and control GSMs actually used, the recomputed
statistics, the membership source (replayed agent test, numerical recovery or dual-model design
reconstruction), the replay provenance (trace file and logged test), the design-reconstruction rules
and the adjudicating models' verbatim judgments where they exist. A manifest hashes every release
file except the manifest itself. A reader can follow any published number back to its GSMs.

  python -m audit_system.v3_public_export_2026_09_30 --final DIR --verifier DIR --out DIR
"""
from audit_system.adjudication_files import decision_files
import argparse
import csv
import json
import math
import shutil
import time
from collections import defaultdict
from pathlib import Path

from audit_system.membership_recovery_2026_09_27 import sha256
from audit_system.paths import ROOTS, portable

MODALITIES = ('dnam', 'hs_rna', 'mm_rna')
MODALITY_LABELS = {'dnam': 'human DNAm', 'hs_rna': 'human RNA-seq', 'mm_rna': 'mouse RNA-seq'}
CLASSIFICATION_FIELDS = ('status', 'contrast_type', 'intervention_class', 'manipulation', 'canonical_target',
                         'relevance_score', 'relevance_score_abs_diff', 'answers', 'receipts', 'relevance_codes_agreed')
DESIGN_AUDIT_FIELDS = ('status', 'reasons', 'orientation', 'specimen', 'tumor', 'answers', 'receipts')
INDIVIDUALS_AUDIT_FIELDS = ('status', 'findings', 'findings_by_reader', 'subject_field', 'receipts')
INDIVIDUALS_STATUS = {
    'SHARED_AGREED': 'both readers report shared individuals with the same partition of the arms; the comparison is recomputed on individual means',
    'INDEPENDENT': 'both readers find no shared individual (or established cell lines)',
    'SHARED_DISAGREE': 'both readers report shared individuals with different partitions or findings; nothing is changed',
    'UNDETERMINED': 'a reader cannot tell from the records; nothing is changed',
    'GEO_INDIVIDUAL_FIELD': 'the series has a GEO individual field; row repairs use it and the comparison is not read',
    'NOT_READ_ARMS_ABOVE_LIMIT': 'the arms hold more than 120 samples; not read',
    'METADATA_UNAVAILABLE': 'no GEO sample table for the series; not read',
    'MODEL_CALL_FAILED': 'a reader returned no usable answer',
}
RANKING_FIELDS = ('rank', 'rank_within_direction', 'hedges_g', 'g_ci_low', 'g_ci_high', 'ci_level', 'selection_k',
                  'direction_v3', 'clock')
MODEL_ALIASES = {
    'gp-6-as': 'GPT-6 Astra (OpenAI), called through the NeoLink gateway',
    'cl-op-5': 'Claude Opus 5 (Anthropic), called through the NeoLink gateway',
    'openai:gpt-6-astra': 'GPT-6 Astra (OpenAI API)',
    'openai:gpt-6-sol': 'GPT-6 Sol (OpenAI API)',
}
REASON_CODES = {
    'CLOCK_OUTPUTS_OUTSIDE_PLAUSIBLE_RANGE': 'clock outputs outside the plausible range',
    'clock outputs outside the plausible range': 'an age-clock arm mean is outside the screening bounds or saturates a defined lower transform bound',
    'SINGLE_CELLS_AS_INDEPENDENT_SAMPLES': 'single cells require tests on recorded individuals',
    'samples are single cells; individuals not recorded': 'single-cell samples cannot establish independent biological replication without recorded individuals',
    'ARMS_DO_NOT_MATCH_GEO_SAMPLE_RECORDS': 'arms do not match the GEO sample records',
    'arms do not match the GEO sample records': 'arms do not match the GEO sample records',
    'ARMS_DEFINED_BY_CLOCK_OUTPUT': 'arms defined by a clock output',
    'arms defined by a clock output': 'arms defined by a clock output',
    'MEMBERSHIP_RECOVERED': 'case/control GSMs recovered (replayed agent test or exact numerical reproduction)',
    'DESIGN_RECONSTRUCTED': 'arms reconstructed by two models from literal metadata rules; identical GSM sets required',
    'MEMBERSHIP_UNRECOVERED': 'historical membership could not be recovered',
    'TEST_NOT_RECONSTRUCTIBLE': 'historical test family cannot be reproduced from saved predictions',
    'STORED_P_UNUSABLE': 'stored P value missing, 0 or 1',
    'MEMBERSHIP_AMBIGUOUS': 'several distinct memberships reproduce the stored statistics',
    'ONE_SAMPLE_TEST_NOT_CONTRAST': 'historical one-sample test reported as a group comparison',
    'CONTINUOUS_TERM_NOT_GROUP_CONTRAST': 'historical continuous regression term reported as a group comparison',
    'ARMS_UNRESOLVED': 'arms of a linked model fit could not be resolved',
    'ARM_TOO_SMALL': 'an arm has fewer than 3 samples',
    'OUTPUTS_CONSTANT_WITHIN_BOTH_ARMS': 'outputs constant within both arms',
    'outputs constant within both arms': 'test undefined because each arm is constant or pooled SD is at the round-off floor',
    'EFFECT_NOT_REPRODUCED': 'stored effect not reproduced by any standard definition',
    'EFFECT_SIGN_ONLY': 'only the direction of a model coefficient could be checked (soft)',
    'LABEL_ORIENTATION_INVERTED': 'arms swapped relative to the comparison label',
    'LABEL_UNCONFIRMED': 'label orientation not confirmed (soft)',
    'LABEL_ADJUDICATED_DUAL_MODEL': 'orientation decided by agreement of GPT-6 Astra and Claude Opus 5',
    'NON_STANDARD_ASSAY': 'enrichment or oxidative-assay samples among the arms',
    'INPUT_EMPTY': 'historical DNAm matrix had no probe rows',
    'INPUT_PROFILE_DETACHED': 'a member DNAm profile correlates with the reference methylome below 0.10 (values detached from probes or inverted)',
    'PROFILE_ATYPICAL': 'a member DNAm profile correlates with the reference methylome at 0.10-0.50 (non-bisulfite assay, GoldenGate or demethylated cells) (soft)',
    'PROFILE_UNSCREENED': 'some members absent from the input-integrity screen',
    'RNA_INPUT_ARTIFACT': 'historical RNA clock input nearly all missing or outputs constant',
    'CLOCK_COVERAGE_LOW': 'fewer than 60% of required CpGs measured in some member',
    'CLOCK_COVERAGE_MID': '60-80% of required CpGs measured (soft)',
    'COVERAGE_UNAVAILABLE': 'per-sample coverage file unavailable',
    'TISSUE_OUT_OF_DOMAIN': 'blood-trained clock applied to non-blood samples (soft)',
    'CLOCK_AGE_INPUT_FROM_HORVATH': 'GrimAge computed with the Horvathv1 prediction as its chronological-age input and every sample coded male (soft)',
    'FDR_BELOW_RAW_P': 'stored FDR smaller than the raw P',
    'COVARIATE_MODEL': 'historical test was a covariate-adjusted model',
    'OUTPUT_MISFILED': "producing trace stored under another study's key",
    'PREDICTIONS_UNAVAILABLE': 'per-sample predictions unavailable',
    'CLOCK_UNREGISTERED': 'DNAm clock absent from the registered clock targets',
    'ARM_ENTIRELY_NON_STANDARD_ASSAY': 'dropping non-standard assay samples empties an arm; the comparison contrasts assays, not biological groups',
    'REPAIRED_CLOCK_WITHHELD': 'the stream-repaired clock verdict is WITHHELD; its reasons are carried into pending_reason',
    'REPAIRED_CLOCK_PROVISIONAL': 'the stream-repaired clock verdict is PROVISIONAL (soft)',
    'MISSING_REQUIRED_SITES': 'stream-repaired clock model reported an error for required sites',
    'NONFINITE_PREDICTION': 'stream-repaired clock predictions absent or not all finite on QC-retained samples',
    'DEGENERATE_OUTPUT': 'stream-repaired clock has at least 3 retained predictions with standard deviation below 1e-8',
    'CLOCK_VALUES_OVERFLOW': 'every clock value of an arm is an overflow artefact (absolute value above 1,000)',
    'AGE_SEPARATES_ARMS': 'recorded ages of the two arms do not overlap, so age cannot be fitted; the statistic is not age-adjusted (soft)',
    'COVERAGE_BELOW_0.60': 'stream-repaired clock minimum retained-sample coverage below 0.60',
    'COVERAGE_0.60_TO_0.80': 'stream-repaired clock minimum retained-sample coverage at least 0.60 and below 0.80 (soft)',
}
ADJUDICATION_ISSUES = {
    'CASE_NOT_CONDITION': 'arm A does not match the labelled case condition',
    'CONTROL_INAPPROPRIATE': 'arm B is not an appropriate control',
    'TISSUE_OR_CELLTYPE_CONFOUNDED': 'tissue or cell type confounds the comparison',
    'TIMEPOINT_CONFOUNDED': 'time point confounds the comparison',
    'TECHNICAL_REPLICATES': 'the arms include technical replicates',
    'NON_STANDARD_ASSAY': REASON_CODES['NON_STANDARD_ASSAY'],
    'CELL_LINE_ONLY': 'the comparison uses cell lines only',
    'MIXED_PLATFORM_OR_BATCH': 'the comparison mixes platforms or batches',
    'OTHER': 'another concern reported in the model judgment',
}
REASON_CODES.update({'ADJUDICATED_' + code: 'both label-adjudicating models reported: ' + meaning + ' (soft)'
                     for code, meaning in ADJUDICATION_ISSUES.items()
                     if code not in ('CASE_NOT_CONDITION', 'OTHER', 'NON_STANDARD_ASSAY')})
ACTIONS = {
    'WITHHOLD_IMPLAUSIBLE_CLOCK_OUTPUTS': 'withhold implausible age-clock outputs from inference and the BH family',
    'SINGLE_CELL_INDIVIDUAL_UNITS': 'aggregate single cells per recorded individual; report descriptively when individuals are unrecorded',
    'PAIR_SUBJECT_TOKENS_IN_TITLES': 'paired test on exact one-to-one subject tokens with opposite temporal labels in GEO titles',
    'WITHHOLD_GEO_ARM_MISMATCH': 'withhold arms that do not match the GEO sample records from inference and the BH family',
    'WITHHOLD_CLOCK_OUTPUT_ARMS': 'withhold arms defined by a clock output from inference and the BH family',
    'DESIGN_RECONSTRUCTED_DUAL_MODEL': 'use arms assigned to the labels by both models from the saved design rules',
    'USE_REPAIRED_INPUT': 'replace historical empty or detached DNAm input with QC-retained stream-repaired predictions; the saved clock verdict governs coverage',
    'DROP_DETACHED_INPUT_SAMPLES': 'drop the recorded detached GSMs from both arms',
    'DROP_NON_STANDARD_ASSAY_SAMPLES': 'drop flagged non-standard assay GSMs from both arms',
    'SWAP_ARMS_CONFIRMED_INVERSION': 'swap arms when both models place the labelled condition in the stored control arm',
    'KEEP_ORIENTATION_ADJUDICATED_SUPPORTED': 'keep the arms when adjudication supports the label despite an inverted heuristic',
    'RECOMPUTE_DECLARED_ESTIMAND': "recompute Welch two-sample t-test; Cohen d with pooled SD; case minus control",
    'DROP_OVERFLOW_CLOCK_VALUES': 'drop GSMs whose clock value is an overflow artefact (absolute value above 1,000) from both arms',
    'RECOMPUTE_DESIGN_ESTIMAND': ('replace the Welch result by the design-respecting statistic named in estimand (individual means, '
                                  'paired t-test or least squares on arm and age); the Welch result stays in recomputed_declared'),
    'DESCRIPTIVE_ONLY_UNITS_BELOW_3': 'fewer than 3 independent individuals per arm (or fewer than 3 pairs): effect reported without a P value; excluded from inference',
    'REPLACE_UNREPRODUCIBLE_EFFECT': 'replace the unreproduced stored effect by the recomputed effect under the declared estimand',
    'DESCRIPTIVE_ONLY_ARM_BELOW_3': 'report effect without a P value; exclude from inference (also emitted when recomputation yields no finite P)',
    'DESCRIPTIVE_ONLY_CONSTANT_OUTPUTS': 'outputs constant within both arms: report mean difference without d, t, P or FDR; exclude from inference',
    'MARK_DUPLICATE_TEST': 'point to the representative of the same study, clock and repaired arms; enter the BH family once',
}
DESIGN_CODES = {
    'ORIENTATION_DISAGREE': 'design-audit models disagree on arm orientation',
    'ORIENTATION_REVERSED': 'both design-audit models place the labelled condition in the control arm',
    'ORIENTATION_UNCLEAR': 'both design-audit models say the records do not establish arm orientation',
    'SPECIMEN_DISAGREE': 'design-audit models disagree on specimen context',
    'SPECIMEN_EX_VIVO_CULTURE': 'both models report primary cells or tissue cultured, stimulated or treated after isolation',
    'SPECIMEN_CELL_LINE': 'both models report immortalized, transformed or cancer cell lines',
    'SPECIMEN_STEM_CELL_DERIVED': 'both models report pluripotent or embryonic stem cells, their derivatives or organoids',
    'SPECIMEN_UNCLEAR': 'both models say the records do not establish specimen context',
    'TUMOR_DISAGREE': 'design-audit models disagree on tumour material',
    'TUMOR_YES': 'both models report tumour or cancer tissue or cells in either arm',
    'TUMOR_UNCLEAR': 'both models say the records do not establish tumour status',
    'CONFOUNDS_MISSING': 'at least one design-audit answer lacks a confounds list',
}
CONFOUND_CODES = {
    'TISSUE_OR_CELL_TYPE': 'different tissues or cell types',
    'AGE_OR_DEVELOPMENTAL_STAGE': 'different age groups or stages',
    'TIME_POINT': 'different collection times or treatment durations',
    'SEX': 'arms mostly of different sex',
    'BATCH_OR_PLATFORM': 'different platforms, batches or processing',
    'ADDITIONAL_TREATMENT': 'a second treatment in one arm only',
    'OTHER': 'another factor differs systematically between the arms besides the labelled condition',
}
for prefix, meaning in (('CONFOUND_AGREED_', 'both design-audit models report '),
                        ('CONFOUND_ONE_MODEL_', 'one design-audit model reports ')):
    DESIGN_CODES.update({prefix + code: meaning + text for code, text in CONFOUND_CODES.items()})
RELEVANCE_CODES = {
    'IN_VIVO': 'in vivo material', 'PRIMARY_TISSUE_OR_CELLS': 'primary tissue or cells',
    'CELL_LINE': 'cell-line material', 'PLURIPOTENT_OR_REPROGRAMMING': 'pluripotent or reprogramming material',
    'TUMOR': 'tumour material', 'AGING_RELEVANT_TISSUE': 'aging-relevant tissue or model',
    'ESTABLISHED_AGING_PATHWAY': 'intervention acts on an established aging pathway',
    'TRANSLATABLE_INTERVENTION': 'approved or clinical-stage drug, diet, exercise, or a druggable gene target',
    'CONTROL_MATCHED': 'case differs from a matched control only by the intervention (same tissue, age, sex, time point)',
    'CONFOUNDED_DESIGN': 'confounded design',
}
# The emitters preserve verbatim model lists, including tokens outside their declared schemas.
RAW_MODEL_CODES = {code: 'unregistered token in a saved model code list; the emitter declares no meaning for this token; retained verbatim, not interpreted as a verifier reason'
                   for code in ('CONFORED_DESIGN', 'IN_VVO', 'CONTROL_INAPPROPRIATE_CONTROL_PLACEHOLDER',
                                'CONFOUNDED_DESIXN', 'CONTROL_INAPPROUNIATE', 'CONFOUNDED_DESICESIGN', 'BELOW',
                                'TRANAPOLATABLE_INTERVENTION', 'AGING_REVENT_TISSUE', 'EPIGENETIC_REGULATORS',
                                'CONFOUNDS_DESIGN', 'AGING_RElevant_TISSUE', 'OT_VIVO', 'EPIGENETIC_REGULATOR',
                                'AGING_RELEVANTANT_TISSUE', 'CONFOULDED_DESIGN', 'ENVIRONMENTAL_EXPOSURE',
                                'CASE_NOT_CONVENTION', 'CASE_NOT_CONDTION', 'CONTROL_INAPPROPRIATE_CONTROL', 'EPIGENETIC')}
RAW_MODEL_CODES.update({
    'CELL_STATE': 'classification contrast type (reprogramming, differentiation or activation), also emitted outside its field in saved model code lists',
    'DEVELOPMENT_OR_TIME': 'classification contrast type for development or time, also emitted outside its field in saved model code lists',
    'CELL_OR_TISSUE_TRANSFER': 'classification intervention class for cell or tissue transfer, also emitted outside its field in saved model code lists',
    'PATHOGEN': 'classification intervention class for pathogen manipulation, also emitted outside its field in saved model code lists',
})
FINAL_CLASSES = {
    'VERIFIED': 'all checks pass after repair; inferential',
    'PROVISIONAL': 'all hard checks pass; soft flags shown; inferential',
    'DESCRIPTIVE': 'test undefined or fewer than 3 samples, independent individuals or pairs per arm: mean difference reported, no P or FDR, excluded from inference; constant outputs have no d or t',
    'PENDING_DESIGN': 'arms not recovered and not reconstructed by two agreeing models',
    'PENDING_ADJUDICATION': 'orientation contradicted by the label heuristic and not yet adjudicated',
    'PENDING_INPUT': 'input or clock coverage insufficient even after repair',
    'NOT_REPAIRABLE': 'the comparison does not contrast biological groups (assays, or arms defined by a clock output)',
    'DUPLICATE': 'the same study, clock and arms as its representative row (duplicate_of); counted once',
}
CHECK_OUTCOMES = {
    'PASS': 'the saved check evaluated this row and recorded no failure or soft concern',
    'FAIL': 'a saved hard reason or duplicate designation failed this check',
    'FLAG': 'a saved soft concern applies; it does not withhold the row by itself',
    'NOT_APPLICABLE': 'the check does not apply to this modality, clock or reconstructed historical estimand',
    'NOT_EVALUATED': 'the saved evidence does not establish an outcome, or an earlier check stopped evaluation',
}
# Codes and evidence refer to the historical verdict; repair outcomes are recorded separately.
CHECKS = {
    'biological_group_contrast': ('ARMS_DEFINED_BY_CLOCK_OUTPUT', 'clock_output_arms',
                                  'Output-defined arms fail the biological group comparison requirement.'),
    'geo_arm_agreement': ('ARMS_DO_NOT_MATCH_GEO_SAMPLE_RECORDS', 'comparison_reclassification',
                          'Arms contradicted by the checked GEO sample records fail.'),
    'membership_recovery': ('MEMBERSHIP_RECOVERED DESIGN_RECONSTRUCTED MEMBERSHIP_UNRECOVERED MEMBERSHIP_AMBIGUOUS ARMS_UNRESOLVED PREDICTIONS_UNAVAILABLE',
                            'membership_status membership_source', 'UNIQUE passes; a recorded recovery failure fails; other stopped recovery is not evaluated.'),
    'test_type': ('TEST_NOT_RECONSTRUCTIBLE ONE_SAMPLE_TEST_NOT_CONTRAST CONTINUOUS_TERM_NOT_GROUP_CONTRAST COVARIATE_MODEL',
                  'test covariates', 'Recovered historical membership passes; recorded unsupported test fails. Historical status is retained for reconstructed arms.'),
    'stored_p_value': ('STORED_P_UNUSABLE', '', 'Numerical recovery past the stored-P gate passes; STORED_P_UNUSABLE fails. Replay or design alone does not validate stored P.'),
    'effect_reproduction': ('EFFECT_NOT_REPRODUCED EFFECT_SIGN_ONLY DESIGN_RECONSTRUCTED', 'effect_definitions_matched',
                            'Recorded matched definitions pass; sign-only is soft; no match fails. Reconstructed designs have no historical effect to reproduce.'),
    'arm_orientation': ('LABEL_ORIENTATION_INVERTED LABEL_UNCONFIRMED LABEL_ADJUDICATED_DUAL_MODEL DESIGN_RECONSTRUCTED', 'label_orientation split',
                        'Saved supported orientation passes, inverted fails, unconfirmed is soft; saved dual-model decisions and swaps describe repair.'),
    'arm_size': ('ARM_TOO_SMALL', 'n_case n_control', 'Saved arm-size evaluation passes unless ARM_TOO_SMALL; repaired DESCRIPTIVE rows fail the inferential arm check.'),
    'assay': ('NON_STANDARD_ASSAY ARM_ENTIRELY_NON_STANDARD_ASSAY', 'non_standard_assay_gsms',
              'DNAm UNIQUE membership reaches the assay screen; absence of its failure code passes. RNA is not applicable.'),
    'dnam_input_integrity': ('INPUT_EMPTY INPUT_PROFILE_DETACHED PROFILE_ATYPICAL PROFILE_UNSCREENED',
                             'profile_r_reference_min detached_gsms atypical_profile_gsms',
                             'Empty or detached fails; atypical is soft; incomplete or absent profile screening is not evaluated. A recorded complete screen without defects passes.'),
    'dnam_clock_coverage': ('CLOCK_COVERAGE_LOW CLOCK_COVERAGE_MID COVERAGE_UNAVAILABLE REPAIRED_CLOCK_WITHHELD REPAIRED_CLOCK_PROVISIONAL COVERAGE_BELOW_0.60 COVERAGE_0.60_TO_0.80',
                           'clock_coverage_min', 'Saved coverage without a coverage reason passes; low or unavailable fails, mid is soft. Stream clock verdict supplies the repair outcome.'),
    'clock_target_age_input': ('CLOCK_UNREGISTERED CLOCK_AGE_INPUT_FROM_HORVATH', 'clock_target',
                               'Recorded target passes; unregistered fails; historical Horvath age input is soft. Target text retains the RNA unit limitation.'),
    'tissue_domain': ('TISSUE_OUT_OF_DOMAIN', 'clock_target',
                      'For a registered blood-trained DNAm clock with UNIQUE membership, no domain flag passes; out of domain is soft. Pan-tissue and RNA clocks are not applicable.'),
    'rna_input': ('RNA_INPUT_ARTIFACT', '', 'The modality-wide RNA input census applied to every verdict row; absence of RNA_INPUT_ARTIFACT passes. DNAm is not applicable.'),
    'fdr_vs_p': ('FDR_BELOW_RAW_P', '', 'Saved verdict comparison passes when historical FDR and P are present without the failure code; missing values are not evaluated. Repaired family FDR is separate.'),
    'duplicate': ('', '', 'Saved duplicate_of fails uniqueness; completed VERIFIED, PROVISIONAL or DESCRIPTIVE rows without it pass. Pending rows were not evaluated.'),
    'design_concerns': (' '.join(k for k in REASON_CODES if k.startswith('ADJUDICATED_')), '',
                        'Adjudicated soft flags or design-audit reasons are soft. Completed label adjudication with no concerns or RANKABLE design audit passes; absent/incomplete judgments are not evaluated.'),
}
HARD_CHECK_CODES = set(REASON_CODES) - {
    'MEMBERSHIP_RECOVERED', 'DESIGN_RECONSTRUCTED', 'LABEL_ADJUDICATED_DUAL_MODEL', 'COVARIATE_MODEL',
    'OUTPUT_MISFILED', 'EFFECT_SIGN_ONLY', 'LABEL_UNCONFIRMED', 'PROFILE_ATYPICAL', 'PROFILE_UNSCREENED',
    'CLOCK_COVERAGE_MID', 'TISSUE_OUT_OF_DOMAIN', 'CLOCK_AGE_INPUT_FROM_HORVATH',
    'REPAIRED_CLOCK_PROVISIONAL', 'COVERAGE_0.60_TO_0.80',
} - {k for k in REASON_CODES if k.startswith('ADJUDICATED_')}


def row_checks(entry, member, modality):
    """Describe saved evaluations and repairs without running a scientific check again."""
    evidence = entry.get('evidence') or {}
    reasons = set(entry['reasons'])
    codes = reasons | set(entry.get('pending_reason') or []) | set(entry.get('soft_flags') or [])
    unique = evidence.get('membership_status') == 'UNIQUE'
    designed = member.get('source') == 'design_reconstruction'
    evaluated = {
        'membership_recovery': unique,
        'test_type': unique and not designed,
        'stored_p_value': member.get('source') in ('recovery', 'replay+recovery') or
                          member.get('status') in ('NONE', 'AMBIGUOUS', 'NO_PREDICTIONS'),
        'effect_reproduction': unique and not designed and bool(evidence.get('effect_definitions_matched')),
        'arm_orientation': evidence.get('label_orientation') == 'SUPPORTED',
        'arm_size': unique and 'n_case' in evidence and 'n_control' in evidence,
        'assay': unique,
        'dnam_input_integrity': evidence.get('profile_r_reference_min') is not None and 'PROFILE_UNSCREENED' not in codes,
        'dnam_clock_coverage': evidence.get('clock_coverage_min') is not None,
        'clock_target_age_input': bool(evidence.get('clock_target')),
        'tissue_domain': unique and bool(evidence.get('clock_target')),
        'rna_input': entry.get('evidence') is not None,
        'fdr_vs_p': entry.get('evidence') is not None and all(release_value(entry['historical'].get(k)) is not None
                                                          for k in ('fdr', 'pvalue')),
        'duplicate': entry['final_class'] in ('VERIFIED', 'PROVISIONAL', 'DESCRIPTIVE', 'DUPLICATE'),
    }
    adjudication = entry.get('adjudication') or {}
    audit = entry['ranking_evidence'].get('design_audit') or {}
    evaluated['design_concerns'] = (adjudication.get('decision') in
                                   ('LABEL_SUPPORTED', 'LABEL_INVERTED_CONFIRMED', 'LABEL_UNCONFIRMED') or
                                   audit.get('status') == 'RANKABLE')
    not_applicable = set()
    if modality != 'dnam':
        not_applicable.update(('assay', 'dnam_input_integrity', 'dnam_clock_coverage', 'tissue_domain'))
    else:
        not_applicable.add('rna_input')
        if entry['clock'] == 'Horvathv1':
            not_applicable.add('tissue_domain')
    if designed:
        not_applicable.add('effect_reproduction')
    checks = {}
    for name, (listed, fields, _) in CHECKS.items():
        triggered = sorted(codes & set(listed.split()))
        support = {f'evidence.{k}': evidence[k] for k in fields.split() if k in evidence}
        support.update({'reasons': sorted(reasons & set(listed.split()))})
        if name in ('membership_recovery', 'test_type', 'stored_p_value'):
            support['membership'] = {k: member.get(k) for k in ('status', 'source', 'historical_status', 'kind', 'stored_test', 'stored_p')}
        if name == 'test_type' and member.get('historical_status') in ('NOT_TWO_SAMPLE_T', 'ONE_SAMPLE_TEST', 'CONTINUOUS_TERM'):
            triggered.append({'NOT_TWO_SAMPLE_T': 'TEST_NOT_RECONSTRUCTIBLE', 'ONE_SAMPLE_TEST': 'ONE_SAMPLE_TEST_NOT_CONTRAST',
                              'CONTINUOUS_TERM': 'CONTINUOUS_TERM_NOT_GROUP_CONTRAST'}[member['historical_status']])
        if name == 'stored_p_value' and member.get('historical_status') == 'STORED_P_UNUSABLE':
            triggered.append('STORED_P_UNUSABLE')
        if name == 'arm_orientation':
            support['adjudication_fields'] = ['adjudication.decision', 'adjudication.judgments', 'adjudication.receipts']
        if name == 'fdr_vs_p':
            support['historical'] = {k: entry['historical'].get(k) for k in ('fdr', 'pvalue')}
        if name == 'duplicate':
            support.update(duplicate_of=entry.get('duplicate_of'), final_class=entry['final_class'])
        if name == 'design_concerns':
            triggered += audit.get('reasons') or []
            support.update(adjudication_fields=['adjudication.issues_agreed', 'adjudication.judgments'],
                           design_audit_fields=['ranking_evidence.design_audit.status', 'ranking_evidence.design_audit.reasons',
                                                'ranking_evidence.design_audit.answers'], soft_flags=entry.get('soft_flags'))
        soft = set(triggered) - HARD_CHECK_CODES - {'MEMBERSHIP_RECOVERED', 'DESIGN_RECONSTRUCTED',
                                                   'LABEL_ADJUDICATED_DUAL_MODEL', 'COVARIATE_MODEL', 'PROFILE_UNSCREENED'}
        outcome = ('NOT_APPLICABLE' if name in not_applicable else
                   'FAIL' if set(triggered) & HARD_CHECK_CODES or (name == 'duplicate' and entry.get('duplicate_of') is not None) else
                   'FLAG' if soft else 'PASS' if evaluated.get(name) else 'NOT_EVALUATED')
        checks[name] = {'outcome': outcome, 'reason_codes': sorted(set(triggered)), 'evidence': support}
    actions = {a['action']: a for a in entry['actions']}
    repairs = {
        'test_type': ('RECOMPUTE_DECLARED_ESTIMAND',),
        'arm_orientation': ('SWAP_ARMS_CONFIRMED_INVERSION', 'KEEP_ORIENTATION_ADJUDICATED_SUPPORTED'),
        'assay': ('DROP_NON_STANDARD_ASSAY_SAMPLES',),
        'dnam_input_integrity': ('USE_REPAIRED_INPUT', 'DROP_DETACHED_INPUT_SAMPLES'),
        'dnam_clock_coverage': ('USE_REPAIRED_INPUT',),
    }
    for name, selected in repairs.items():
        saved = [actions[a] for a in selected if a in actions]
        if not saved:
            continue
        result = 'NOT_EVALUATED'
        if entry.get('recomputed'):
            result = 'PASS'
        if name == 'dnam_clock_coverage':
            verdict = saved[0].get('clock_verdict') or {}
            result = {'VERIFIED': 'PASS', 'PROVISIONAL': 'FLAG', 'WITHHELD': 'FAIL'}.get(verdict.get('verdict'), 'NOT_EVALUATED')
        if name == 'assay' and 'ARM_ENTIRELY_NON_STANDARD_ASSAY' in codes:
            result = 'FAIL'
        if name == 'dnam_input_integrity' and 'INPUT_PROFILE_DETACHED' in (entry.get('pending_reason') or []):
            result = 'FAIL'
        checks[name]['after_repair'] = {'outcome': result, 'actions': saved}
    if entry.get('recomputed'):
        checks['arm_size']['after_repair'] = {'outcome': 'FAIL' if 'DESCRIPTIVE_ONLY_ARM_BELOW_3' in actions else 'PASS',
                                             'evidence': {k: entry['recomputed'].get(k) for k in ('n_case', 'n_control', 'p')}}
        checks['effect_reproduction']['replacement'] = {'actions': [actions[a] for a in
                                                       ('RECOMPUTE_DECLARED_ESTIMAND', 'REPLACE_UNREPRODUCIBLE_EFFECT') if a in actions],
                                                       'recomputed': entry['recomputed']}
        if entry['recomputed'].get('fdr_family') is not None:
            checks['fdr_vs_p']['after_repair'] = {'outcome': 'PASS', 'evidence': {k: entry['recomputed'].get(k)
                                                   for k in ('p', 'fdr_family', 'fdr_family_size')}}
    return checks


def code_definitions():
    return {'reason_codes': REASON_CODES, 'actions': ACTIONS, 'adjudication_issues': ADJUDICATION_ISSUES,
            'design_audit_codes': DESIGN_CODES, 'design_confounds': CONFOUND_CODES,
            'relevance_codes': RELEVANCE_CODES, 'raw_model_codes': RAW_MODEL_CODES}


def undefined_codes(value, definitions=None):
    """Find undefined emitted codes, including codes inside model answers and repair evidence."""
    definitions = definitions or code_definitions()
    defined = {code for mapping in definitions.values() for code in mapping}
    missing = set()
    code_fields = {'reasons', 'reason_codes', 'pending_reason', 'soft_flags', 'issues', 'issues_agreed',
                   'confounds', 'relevance_codes', 'relevance_codes_agreed'}

    def visit(record):
        if isinstance(record, dict):
            for field, item in record.items():
                if field in code_fields and isinstance(item, list):
                    missing.update(code for code in item if isinstance(code, str) and code not in defined)
                elif field == 'action' and isinstance(item, str) and item not in definitions['actions']:
                    missing.add(item)
                visit(item)
        elif isinstance(record, list):
            for item in record:
                visit(item)

    visit(value)
    return missing


IDENTIFIER_COLUMNS = {
    'gse_id', 'comparison_name', 'condition_name', 'condition_category', 'clock', 'modality', 'atlas_row_ordinal',
    'intervention', 'dataset_id', 'case_id', 'case_id_x', 'case_id_y', 'tissue', 'cell_type', 'control_group',
}
V3_COLUMNS = {
    'v3_final_class': 'Final repair class (codes.json).',
    'v3_original_verdict': 'Verdict before repair.', 'v3_reasons': 'Historical verdict reason codes, separated by semicolons.',
    'v3_actions': 'Ordered repair action codes, separated by semicolons.',
    'v3_soft_flags': 'Soft concerns after repair, separated by semicolons.', 'v3_pending': 'Withholding reason codes, separated by semicolons.',
    'v3_n_case': 'Number of case units used after repair (samples; individuals or pairs under a design estimand).',
    'v3_n_control': 'Number of control units used after repair (samples; individuals or pairs under a design estimand).',
    'v3_mean_diff': 'Repaired case-minus-control estimate in the clock target unit (mean difference; arm coefficient under age adjustment).',
    'v3_cohen_d': 'Repaired case-minus-control Cohen d (pooled SD; residual SD under age adjustment).',
    'v3_p': 'Repaired two-sided P of the row estimand (Welch two-sample unless v3_design_method is set); empty for DESCRIPTIVE and unresolved rows.',
    'v3_design_method': 'Design-respecting estimand used instead of the Welch test: UNIT_MEANS, PAIRED, AGE_ADJUSTED or AGE_ADJUSTED_UNIT_MEANS; empty otherwise.',
    'v3_fdr': 'BH FDR over unique inferential repaired rows within the modality; duplicates inherit the representative FDR.',
    'v3_membership_source': 'Replayed test, numerical recovery or design reconstruction membership source.',
    'v3_duplicate_of': 'Representative row ordinal within the same modality.', 'v3_representative_class': 'Final class of a duplicate representative.',
    'v3_in_vivo_share': 'Share of annotated repaired members labelled in vivo by historical GPT-4 annotations.',
    'v3_d_se': 'Hedges-Olkin standard error of repaired Cohen d.', 'v3_d_ci_low': 'Lower normal 95% confidence bound of repaired Cohen d.',
    'v3_d_ci_high': 'Upper normal 95% confidence bound of repaired Cohen d.', 'v3_direction': 'Sign of repaired Cohen d.',
    'v3_robust': 'VERIFIED/PROVISIONAL, FDR below .05 and normal 95% d interval excludes zero.',
    'v3_rank_bound': 'Absolute normal 95% d interval bound nearest zero; comparison-ranking tables use adjusted exact Hedges g intervals instead.',
    'v3_label': 'Intervention or condition label truncated to 80 characters.', 'v3_context': 'in_vivo if annotated share is at least .5, in_vitro otherwise; unknown if absent.',
    'v3_rankable_category': 'Category is Genetic, Drug, Disease or Environment; further ranking checks still apply.',
}
HISTORICAL_COLUMNS = {
    'statistical_test': 'Stored test name or family.', 'pvalue': 'Stored raw P value.',
    'fdr': 'Stored FDR.', 'effect_size': 'Stored effect; its definition depends on the historical test.',
    'direction': 'Stored effect direction.', 'log2FoldChange': 'Stored log2FoldChange field (not necessarily a fold change for clock outputs).',
    'abs_effect': 'Stored absolute effect.', 'abs_log2FoldChange': 'Stored absolute log2FoldChange field.',
    'neg_log10_pvalue': 'Stored negative log10 P.', 'trust_score': 'Historical model trust score.',
    'original_trust_score': 'Trust score before historical score processing.', 'model_version': 'Historical model version.',
    'final_score': 'Historical composite ranking score.', 'final_score_pre_weighted_audit': 'Historical score before the weighted audit.',
    'covariate': 'Stored covariate annotation.', 'covariates': 'Stored model covariates.',
    'outlier_status': 'Historical effect/P outlier label.', 'note': 'Historical row note.', 'notes': 'Historical row notes.',
    'trust_score_justification': 'Historical rationale for the trust score.', 'audit_decision': 'Historical audit decision.',
    'control_is_healthy': 'Historical judgment of control health status.', 'comparison_inverted': 'Historical judgment of comparison inversion.',
    'audit_rationale': 'Historical audit rationale.', 'sample_size': 'Stored sample-size annotation.',
    'n_samples': 'Stored number-of-samples annotation.',
}
for arm, names in {
    'case/condition': ('condition', 'case', 'test', 'treatment', 'intervention', 'group1'),
    'control': ('control', 'group2'),
    'first stored group (historical test definition)': ('1',),
    'second stored group (historical test definition)': ('2',),
}.items():
    for name in names:
        for column in ('n_' + name, name + '_n', 'sample_size_' + name, 'n' + name):
            HISTORICAL_COLUMNS[column] = f'Stored {arm} sample count.'
        for column in ('mean_' + name, name + '_mean', 'mean' + name):
            HISTORICAL_COLUMNS[column] = f'Stored {arm} mean.'
        for column in ('std_' + name, name + '_std', 'std' + name):
            HISTORICAL_COLUMNS[column] = f'Stored {arm} standard deviation.'


def export_atlas(source, target):
    """Rename historical columns only in the exported atlas, preserving its values and identifiers."""
    columns = {}
    try:
        with source.open(newline='') as reader, target.open('w', newline='') as writer:
            rows = csv.DictReader(reader)
            names = {c: c if c.startswith('v3_') or c in IDENTIFIER_COLUMNS else 'historical_' + c for c in rows.fieldnames}
            if len(set(names.values())) != len(names):
                raise ValueError('Historical column renaming would produce duplicate columns')
            for original, name in names.items():
                description = (V3_COLUMNS.get(original) if original.startswith('v3_') else
                               'Identifier or label carried from the source atlas: ' + original + '.' if original in IDENTIFIER_COLUMNS else
                               HISTORICAL_COLUMNS.get(original, 'Historical source field ' + original + '.') +
                               ' Retained for provenance; not a verifier-v3 result.')
                if description is None:
                    raise ValueError(f'Undocumented v3 column: {original}')
                columns[name] = description
            output = csv.DictWriter(writer, fieldnames=list(names.values()))
            output.writeheader()
            for row in rows:
                output.writerow({names[k]: v for k, v in row.items()})
    except (OSError, ValueError) as exc:
        raise RuntimeError(f'Cannot export atlas {source} to {target}: {exc}') from exc
    return columns


def load_jsonl(path):
    if not Path(path).exists():
        return []
    with open(path) as source:
        return [json.loads(line) for line in source if line.strip()]


def release_value(value):
    """Replace non-finite numbers with null and local paths with bundle-relative references."""
    if isinstance(value, dict):
        return {k: release_value(value[k]) for k in sorted(value)}
    if isinstance(value, (list, tuple)):
        return [release_value(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return portable(value.absolute())
    if isinstance(value, str) and (value.startswith(('/', './', '../', 'clockbase_data/'))):
        return portable(Path(value).absolute())
    return value


def write_json(path, value, indent=None):
    try:
        path.write_text(json.dumps(release_value(value), allow_nan=False, indent=indent, default=str))
    except (OSError, ValueError, TypeError) as exc:
        raise RuntimeError(f'Cannot write release JSON {path}: {exc}') from exc


def export_file(source, target):
    """Copy a release file, normalizing JSON and JSONL records before writing."""
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        if source.suffix == '.json':
            write_json(target, json.loads(source.read_text()), indent=2)
        elif source.suffix == '.jsonl':
            with source.open() as reader, target.open('w') as writer:
                for line in reader:
                    if line.strip():
                        writer.write(json.dumps(release_value(json.loads(line)), allow_nan=False) + '\n')
        else:
            shutil.copyfile(source, target)
    except (OSError, ValueError) as exc:
        raise RuntimeError(f'Cannot export {source} to {target}: {exc}') from exc


def study_id(value):
    value = release_value(value)
    return None if value is None or str(value).strip().lower() in ('', 'nan', 'none', 'null') else str(value)


def comparison_key(row, modality=None):
    values = [modality or row.get('modality'), study_id(row.get('gse_id')),
              row.get('comparison_name'), row.get('condition_name')]
    values[0] = MODALITY_LABELS.get(values[0], values[0])
    return tuple(None if v == '' else release_value(v) for v in values)


def ranking_evidence(final_dir, verifier_dir):
    classes = {comparison_key(r): {k: r.get(k) for k in CLASSIFICATION_FIELDS}
               for r in load_jsonl(verifier_dir / 'classification-001' / 'comparisons.jsonl')}
    audits = {}
    for path in sorted(verifier_dir.glob('design-audit-0*/comparisons.jsonl')):
        for row in load_jsonl(path):
            audits.setdefault(comparison_key(row), {'run': path.parent.name,
                                                    **{k: row.get(k) for k in DESIGN_AUDIT_FIELDS}})
    ranked = {}
    ranked_sources = {}
    ranking_path = final_dir / 'comparison-ranking' / 'ranked_comparisons_all.csv'
    if ranking_path.exists():
        directions = defaultdict(int)
        with ranking_path.open(newline='') as source:
            for rank, row in enumerate(csv.DictReader(source), 1):
                direction = row.get('direction_v3') or None
                directions[direction] += 1
                result = {'ranked': True, 'rank': rank, 'rank_within_direction': directions[direction]}
                for field in ('hedges_g', 'g_ci_low', 'g_ci_high', 'ci_level', 'selection_k'):
                    result[field] = float(row[field]) if row.get(field) else None
                result.update(direction_v3=direction, clock=row.get('clock') or None)
                key = (comparison_key(row), result['clock'])
                ranked[key] = result
                ranked_sources[key] = row
    return classes, audits, ranked, ranked_sources, ranking_path.exists()


def individuals_evidence(verifier_dir):
    """{comparison: the two-reader individuals audit record (first run that read it)}."""
    out = {}
    for path in sorted(verifier_dir.glob('individuals-audit-0*/comparisons.jsonl')):
        for row in load_jsonl(path):
            out.setdefault(comparison_key(row), {'run': path.parent.name,
                                                 **{k: row.get(k) for k in INDIVIDUALS_AUDIT_FIELDS if row.get(k) is not None}})
    return out


def ranking_row_matches(row, source):
    """Resolve repeated comparison labels to the selected repaired row, allowing CSV float rounding."""
    if row['final_class'] not in ('VERIFIED', 'PROVISIONAL') or row['final_class'] != source.get('v3_final_class'):
        return False
    recomputed = row.get('recomputed') or {}
    for field in ('n_case', 'n_control', 'mean_diff', 'p'):
        value, stored = recomputed.get(field), source.get('v3_' + field)
        if stored:
            if value is None or not math.isclose(float(value), float(stored), rel_tol=1e-12, abs_tol=1e-15):
                return False
    return True


def readme(columns=None):
    overview = '''# ClockBase verifier-v3 public export

`studies/<modality>/<study>.json` contains the verification records for one study.
Modalities are dnam (human DNAm), hs_rna (human RNA-seq) and mm_rna (mouse RNA-seq).
Joint study IDs use underscores in filenames. Missing IDs use NO_STUDY_ID.json and
have gse_id null. The top-level fields are schema, modality, gse_id, rows and final_classes
(row counts by final class).

Each row contains row_ordinal (the source row index), comparison_name, condition_name,
clock, historical (stored statistics), verdict_before_repair, reasons, evidence,
actions (ordered repairs), final_class, duplicate_of, soft_flags, pending_reason,
estimand, recomputed (statistics after repair), membership_source, case_gsms,
control_gsms and predictions (prediction path and SHA-256). When the sample records
show a repeated individual, paired arms or recorded ages, recomputed holds the
design-respecting statistic named in estimand, design describes it (findings, method,
individual and age fields, units per arm, sample tables read) and recomputed_declared
keeps the Welch result on the same arms.
When present, stored_case_gsms and stored_control_gsms preserve the arms before repair;
replay contains traces, the logged test and output_misfiled; design_rules and
historical_status describe design reconstruction; adjudication contains model decisions,
judgments, agreed issues and receipts.

ranking_evidence joins on modality, gse_id, comparison_name and condition_name.
classification contains status, contrast_type, intervention_class, manipulation,
canonical_target, relevance_score, relevance_score_abs_diff, both model answers,
receipts and relevance_codes_agreed. design_audit contains the first matching run
(runs sorted by name), status, reasons, orientation, specimen, tumor, both answers
and receipts. ranking contains ranked, rank, rank_within_direction, hedges_g,
g_ci_low, g_ci_high, ci_level, selection_k, direction_v3 and clock. Ranks are 1-based
positions in ranked_comparisons_all.csv; direction ranks count positions within each
direction. The selected row is matched by clock, final class, arm sizes, mean difference
and P value (allowing CSV rounding); if these tie, the first matching repaired row wins.
Duplicate and pending rows are not marked ranked. If the ranking table exists,
absent rows have ranked false and null ranking fields; if it is missing, ranking is null.
Missing classification or design_audit is null. Non-finite numbers are null throughout JSON.
individuals_audit is present for ranking candidates read by the two-reader individuals audit
(run, status, findings, receipts); individuals-audit/<run>.jsonl holds the complete records.

The root tables contain the combined atlas, error census, published status summary
and per-publication status CSVs. ranking/ contains every comparison-ranking file,
including the all-comparison table, category lists, direction lists, relevance-filtered
versions and summary. classification/comparisons.jsonl contains the classification
records; design-audit/<run>.jsonl preserves every design-audit run.

codes.json defines verification reason codes, pending reasons, soft flags, actions,
adjudication issues, design-audit reasons/confounds, classification relevance codes,
final classes, model aliases and the checks derivation table.
Design-audit reasons and classification judgments are recorded beside both answers
in ranking_evidence and the accompanying JSONL files.

Filesystem paths use root:relative/path, for example verifier:replay-dnam-001/inputs/GSE1.csv.
The root labels are listed in codes.json under portable_paths; audit_system.paths.resolve()
maps these references to a reader's configured bundle roots. Relative paths use forward slashes.

manifest.json records the source run creation time when available (otherwise the export time),
study and row counts, the exporter module SHA-256,
and SHA-256 hashes keyed by export-relative filename for every exported file except
manifest.json itself (a manifest cannot include its own hash). Hashes describe the
normalized release bytes; source datasets and verifier decisions are not modified.
'''
    checks = '\n## Checks derived from saved records\n\n'
    checks += ('Every row has the same checks object. No check is recomputed by this exporter. '
               'outcome describes the saved historical verdict, with reason_codes and evidence values or field references. '
               'after_repair, where present, describes a recorded repair and its outcome; replacement preserves the '
               'recomputed effect without claiming that a failed historical reproduction subsequently passed. '
               'Missing evidence never becomes a pass. Membership failures stop the per-arm checks. '
               'Stored-P validity is only established by the numerical recovery gate, not by replay or reconstructed arms. '
               'Historical input/coverage failures remain visible even when USE_REPAIRED_INPUT resolves them. '
               'The recorded stream clock verdict covers QC-retained study samples; no per-arm rescreen is inferred. '
               'Design-audit concerns affect ranking and are FLAG here; they do not change the repair class.\n\n')
    checks += ('codes.json also catalogs raw_model_codes found in the saved model lists. These tokens '
               'fall outside the emitter schemas (including spelling errors); they remain verbatim, '
               'with no inferred correction or verification meaning. The reason/action definitions come '
               'from row_verdicts_2026_09_27, row_repairs_2026_09_28 and dnam_stream_repair_2026_09_26; '
               'model vocabularies come from label_adjudication_2026_09_27, design_audit_2026_09_30 '
               'and intervention_classification_2026_09_29.\n\n')
    checks += '\n'.join(f'- {name}: {meaning}' for name, meaning in CHECK_OUTCOMES.items()) + '\n\n'
    checks += '| Check | Reason codes | Saved evidence | Derivation |\n| --- | --- | --- | --- |\n'
    for name, (codes, fields, derivation) in CHECKS.items():
        checks += f'| {name} | {codes or "duplicate_of (no emitted reason code)"} | {fields or "membership, historical, actions or model judgments as recorded"} | {derivation} |\n'
    checks += ('\ncase_gsms/control_gsms are the arms recorded in repair inputs. When repair stopped before a test, '
               'they are not evidence that a test ran; recomputed is null. membership preserves recovered arms and status '
               'even when repair stopped before recording inputs. DESCRIPTIVE rows have no repaired P or FDR; '
               'historical P values remain explicitly labelled historical and are excluded from inference. '
               'Every withheld row retains pending_reason; nothing is dropped.\n')
    checks += ('\n## Exported atlas columns\n\nOnly identifiers, labels and v3_* columns retain their names. '
               'Every other column has a historical_ prefix. The internal atlas is unchanged. '
               'Ranking tables retain their original schemas; ranking_evidence is joined before atlas export. '
               'Consumers of historical atlas fields should use the historical_ names, or rename them back '
               'when passing this exported copy to an internal verifier module.\n\n')
    checks += '| Column | Meaning |\n| --- | --- |\n'
    checks += ''.join(f'| {name} | {meaning} |\n' for name, meaning in (columns or {}).items())
    return overview + checks


def run(final_dir, verifier_dir, out_dir):
    final_dir, verifier_dir, out_dir = Path(final_dir), Path(verifier_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    counts = {}
    classes, audits, ranked, ranked_sources, has_ranking = ranking_evidence(final_dir, verifier_dir)
    individuals = individuals_evidence(verifier_dir)
    ranked_used = set()
    for modality in MODALITIES:
        repairs = load_jsonl(final_dir / f'repairs-{modality}' / 'repairs.jsonl')
        members = {r['row_ordinal']: r for r in load_jsonl(final_dir / f'merged-{modality}' / 'rows.jsonl')}
        verdicts = {r['row_ordinal']: r for r in load_jsonl(final_dir / f'verdicts-{modality}' / 'row_verdicts.jsonl')}
        links = {r['row_ordinal']: r for r in load_jsonl(verifier_dir / f'link-{modality}-001' / 'row_links.jsonl')}
        adjudications = {}
        for path in decision_files(verifier_dir / f'adjudication-{modality}-001'):
            for d in load_jsonl(path):
                if d.get('decision') in ('LABEL_SUPPORTED', 'LABEL_INVERTED_CONFIRMED', 'LABEL_UNCONFIRMED') \
                        or d['row_ordinal'] not in adjudications:  # an incomplete later run never erases a decision
                    adjudications[d['row_ordinal']] = {k: d.get(k) for k in ('decision', 'judgments', 'issues_agreed', 'receipts')}
        by_study = defaultdict(list)
        for r in repairs:
            ordinal = r['row_ordinal']
            member, link, verdict = members.get(ordinal, {}), links.get(ordinal, {}), verdicts.get(ordinal, {})
            entry = {'row_ordinal': ordinal, 'comparison_name': r['comparison_name'], 'condition_name': r['condition_name'],
                     'clock': r['clock'], 'historical': r['original'], 'verdict_before_repair': r['original_verdict'],
                     'reasons': r['reasons'], 'evidence': verdict.get('evidence'), 'actions': r['actions'],
                     'final_class': r['final_class'], 'duplicate_of': r.get('duplicate_of'), 'soft_flags': r.get('soft_flags', []),
                     'pending_reason': r.get('pending_reason'), 'estimand': r.get('estimand'),
                     'recomputed': r.get('recomputed'), 'membership_source': member.get('source'),
                     **({'design': r['design']} if r.get('design') else {}),
                     **({'verification_rules': r['verification_rules']} if r.get('verification_rules') else {}),
                     **({'recomputed_before_verification_rules': r['recomputed_before_verification_rules']}
                        if r.get('recomputed_before_verification_rules') else {}),
                     **({'recomputed_declared': r['recomputed_declared']} if r.get('recomputed_declared') else {}),
                     'case_gsms': (r.get('inputs') or {}).get('case_gsms'), 'control_gsms': (r.get('inputs') or {}).get('control_gsms'),
                     'predictions': (r.get('inputs') or {}).get('predictions')}
            entry['membership'] = {k: member.get(k) for k in ('status', 'source', 'historical_status', 'kind', 'test',
                                   'stored_test', 'stored_p', 'case_gsms', 'control_gsms', 'n_case', 'n_control')}
            entry['representative_class'] = r.get('representative_class')
            inputs = r.get('inputs') or {}
            if (inputs.get('stored_case_gsms'), inputs.get('stored_control_gsms')) != (entry['case_gsms'], entry['control_gsms']):
                entry['stored_case_gsms'] = inputs.get('stored_case_gsms')  # arms before a swap or sample drop
                entry['stored_control_gsms'] = inputs.get('stored_control_gsms')
            if link.get('link_status') == 'LINKED':
                entry['replay'] = {'traces': link.get('link_traces'), 'test': {k: link['link'].get(k) for k in
                                   ('fn', 'kind', 'term', 'formula', 'trace', 'test_line', 'trace_prompt_gse')},
                                   'output_misfiled': link.get('misfiled')}
            if member.get('source') == 'design_reconstruction':
                entry['design_rules'] = member.get('design_rules')
                entry['historical_status'] = member.get('historical_status')
            if ordinal in adjudications:
                entry['adjudication'] = adjudications[ordinal]
            key = comparison_key(r, modality)
            ranking_key = (key, r['clock'])
            ranking = ranked.get(ranking_key)
            if ranking is not None:
                if ranking_key in ranked_used or not ranking_row_matches(r, ranked_sources[ranking_key]):
                    ranking = None
                else:
                    ranked_used.add(ranking_key)
            if ranking is None and has_ranking:
                ranking = {'ranked': False, **{k: None for k in RANKING_FIELDS}}
            entry['ranking_evidence'] = {'classification': classes.get(key), 'design_audit': audits.get(key),
                                         'ranking': ranking}
            if key in individuals:
                entry['ranking_evidence']['individuals_audit'] = individuals[key]
            entry['checks'] = row_checks(entry, member, modality)
            missing = undefined_codes(entry)
            if missing:
                raise ValueError(f'Undefined release codes in {modality} row {ordinal}: {sorted(missing)}')
            by_study[study_id(r.get('gse_id'))].append(entry)
        target = out_dir / 'studies' / modality
        target.mkdir(parents=True)
        for gse, rows in by_study.items():
            final_classes = defaultdict(int)
            for row in rows:
                final_classes[row['final_class']] += 1
            name = 'NO_STUDY_ID' if gse is None else str(gse).replace('/', '_')
            path = target / (name + '.json')  # joint records such as GSE85385/GSE85388
            write_json(path, {'schema': 'clockbase-v3-study-trace-v1', 'modality': modality, 'gse_id': gse,
                              'rows': rows, 'final_classes': final_classes})
        counts[modality] = {'studies': len(by_study), 'rows': len(repairs)}
    atlas_source = final_dir / 'atlas-v3' / 'atlas_v3_combined.csv'
    columns = export_atlas(atlas_source, out_dir / 'atlas_v3_combined.csv') if atlas_source.exists() else {}
    for name, source in (('error_census.json', final_dir / 'error-census' / 'census.json'),
                         ('published_status_summary.json', final_dir / 'published-status' / 'summary.json')):
        if Path(source).exists():
            export_file(source, out_dir / name)
    for path in sorted((final_dir / 'published-status').glob('*_v3_status.csv')):
        export_file(path, out_dir / path.name)
    for path in sorted((final_dir / 'comparison-ranking').rglob('*')):
        if path.is_file():
            export_file(path, out_dir / 'ranking' / path.relative_to(final_dir / 'comparison-ranking'))
    classification_path = verifier_dir / 'classification-001' / 'comparisons.jsonl'
    if classification_path.exists():
        export_file(classification_path, out_dir / 'classification' / 'comparisons.jsonl')
    for path in sorted(verifier_dir.glob('design-audit-0*/comparisons.jsonl')):
        export_file(path, out_dir / 'design-audit' / (path.parent.name + '.jsonl'))
    for path in sorted(verifier_dir.glob('individuals-audit-0*/comparisons.jsonl')):
        export_file(path, out_dir / 'individuals-audit' / (path.parent.name + '.jsonl'))
    write_json(out_dir / 'codes.json', {**code_definitions(), 'final_classes': FINAL_CLASSES,
                                      'individuals_audit_status': INDIVIDUALS_STATUS,
                                      'checks': {name: {'reason_codes': codes.split() + (list(DESIGN_CODES) if name == 'design_concerns' else []),
                                                       'evidence_fields': fields.split(), 'derivation': derivation,
                                                       'outcomes': CHECK_OUTCOMES}
                                                 for name, (codes, fields, derivation) in CHECKS.items()},
                                      'model_aliases': MODEL_ALIASES,
                                      'portable_paths': {'format': 'root:relative/path (forward slashes); resolve the root '
                                                                    'using audit_system.paths.ROOTS',
                                                         'roots': list(ROOTS)}}, indent=2)
    (out_dir / 'README.md').write_text(readme(columns))
    files = {path.relative_to(out_dir).as_posix(): sha256(path) for path in sorted(out_dir.rglob('*')) if path.is_file()}
    provenance = final_dir / 'reclassification_provenance.json'
    created = json.loads(provenance.read_text()).get('created') if provenance.is_file() else None
    manifest = {'schema': 'clockbase-v3-public-export-v1', 'created': created or time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                'counts': counts, 'module_sha256': sha256(__file__), 'files': files}
    write_json(out_dir / 'manifest.json', manifest, indent=1)
    print(json.dumps({'counts': counts, 'files': len(files)}, allow_nan=False, indent=1))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--final', type=Path, required=True)
    parser.add_argument('--verifier', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    run(args.final, args.verifier, args.out)


if __name__ == '__main__':
    main()
