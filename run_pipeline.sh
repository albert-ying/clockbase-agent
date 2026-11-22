#!/usr/bin/env bash
#
# ClockBase Agent: versioned source ingestion, original producer and full audit.
#
# Canonical usage (OpenAI key supplied through the environment):
#   bash run_pipeline.sh --refresh REFRESH_CONFIG.json --output REFRESH_ROOT
#   bash run_pipeline.sh --plan UPDATE_PLAN.json --output REGISTRY --until report
#   bash run_pipeline.sh --discover DISCOVERY_PLAN.json --output DISCOVERY
#   bash run_pipeline.sh --audit-corpus run --output CORPUS_RUN
#
# --update is an alias for the default versioned update controller.
# --producer-only runs the existing eight-stage single-study producer and its
# basic numerical audit; it does not run the comprehensive sixteen-check audit.
# --legacy retains the historical agent/score/audit/plot stage dispatcher.
#
# Plans explicitly declare sources, models, tariffs, cost/call limits and input
# limits. Changed sources or policies use new versions; incomplete results and
# failures are retained. Execution completion does not establish scientific PASS.
# --refresh generates intermediate plans for pinned-local-ARCHS4 human/mouse RNA
# and supported human-DNAm inputs; unsupported source formats stay in its ledger.
# PYTHON overrides the interpreter; the installed .venv-pipeline is preferred.
# Refresh uses the released src/verify_ingestion_2026_09_20.py source snapshot.
# The recorded DNAm R environment is in
# src/dnam_idat_preprocess_2026_09_24.sessionInfo.txt; see its provenance caveat.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# Complete future updates are the default. Narrower producer and historical
# modes require explicit selection.
if [[ "${1:-}" != "--legacy" ]]; then
    MODULE="src.update_pipeline_2026_09_20"
    if [[ "${1:-}" == "--audit-corpus" ]]; then
        MODULE="audit_system.corpus_audit_2026_09_20"
        shift
    elif [[ "${1:-}" == "--update" ]]; then
        shift
    elif [[ "${1:-}" == "--producer-only" ]]; then
        MODULE="src.pipeline_2026_09_10"
        shift
    fi
    if [[ -z "${PYTHON:-}" && -x "$HERE/.venv-pipeline/bin/python" ]]; then
        exec "$HERE/.venv-pipeline/bin/python" -m "$MODULE" "$@"
    fi
    exec "${PYTHON:-python3}" -m "$MODULE" "$@"
fi
shift

PYTHON="${PYTHON:-python3}"
MODALITY="mouse_rna"
STAGES="agent,score,audit,plot"

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------

usage() {
    cat <<'USAGE'
Historical pipeline stage dispatcher:
  bash run_pipeline.sh --legacy --modality mouse_rna|human_rna|human_dnam
  bash run_pipeline.sh --legacy --stages agent,score,audit,plot
Historical agent/audit stages require ANTHROPIC_API_KEY; score uses OPENAI_API_KEY.
Use the default --plan mode for versioned ingestion and comprehensive auditing.
USAGE
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --modality)
            MODALITY="$2"; shift 2 ;;
        --stages)
            STAGES="$2"; shift 2 ;;
        -h|--help)
            usage 0 ;;
        *)
            echo "Unknown argument: $1" >&2
            usage 2 ;;
    esac
done

case "$MODALITY" in
    mouse_rna|human_rna|human_dnam) ;;
    *)
        echo "ERROR: --modality must be one of mouse_rna|human_rna|human_dnam (got: $MODALITY)" >&2
        exit 2 ;;
esac

run_stage() {
    local needle="$1"
    [[ ",${STAGES}," == *",${needle},"* ]]
}

log() {
    printf '\n=== [%s] %s ===\n' "$(date +'%Y-%m-%d %H:%M:%S')" "$*"
}

# -----------------------------------------------------------------------------
# Stage 1 — Agent: primary CodeAct loop (Sonnet 3.7)
# -----------------------------------------------------------------------------

if run_stage agent; then
    log "Stage 1 / agent — modality=${MODALITY}"
    case "$MODALITY" in
        mouse_rna)
            (cd src && "$PYTHON" parallel_mm_rna.py)
            ;;
        human_rna)
            (cd src && "$PYTHON" parallel_hs_rna.py)
            ;;
        human_dnam)
            # No standalone src/parallel_hs_dnam.py is in the working tree;
            # the DNAm modality is run via the minipaper_fork (which combines
            # agent and scoring). The agent stage is therefore a no-op for
            # human_dnam — proceed straight to the score stage.
            echo "INFO: human_dnam has no separate agent stage; see minipaper_fork."
            ;;
    esac
fi

# -----------------------------------------------------------------------------
# Stage 2 — Score: minipaper_fork (gpt-4o-mini) + cheap-RAG
# -----------------------------------------------------------------------------

if run_stage score; then
    log "Stage 2 / score — modality=${MODALITY}"
    # minipaper_fork reads the modality from the CSV_MODALITY env var.
    case "$MODALITY" in
        mouse_rna)
            (cd minipaper_fork && CSV_MODALITY=mm_rna NUM_WORKERS="${NUM_WORKERS:-4}" \
                "$PYTHON" src/run.py)
            ;;
        human_rna)
            (cd minipaper_fork && CSV_MODALITY=hs_rna NUM_WORKERS="${NUM_WORKERS:-12}" \
                "$PYTHON" src/run.py)
            ;;
        human_dnam)
            (cd minipaper_fork && CSV_MODALITY=hs_dnam NUM_WORKERS="${NUM_WORKERS:-4}" \
                "$PYTHON" src/run.py)
            ;;
    esac

    # Score-merging / scoring-input preparation steps (these are the
    # `s0[23]_*` scripts at repo root). Mouse runs use the bare-numbered
    # files; human runs use the *_human.py variants.
    case "$MODALITY" in
        mouse_rna)
            "$PYTHON" summarize_scores.py
            ;;
        human_rna|human_dnam)
            "$PYTHON" s02_prepare_scoring_input_human.py
            "$PYTHON" s03_merge_scores_human.py
            "$PYTHON" summarize_scores_human.py
            ;;
    esac
fi

# -----------------------------------------------------------------------------
# Stage 3 — Audit: autonomous Opus 4.7 plausibility re-audit
# -----------------------------------------------------------------------------

if run_stage audit; then
    log "Stage 3 / audit — modality=${MODALITY}"
    case "$MODALITY" in
        mouse_rna)
            INPUT_CSV="result/main_results_result_mm_rna_with_scores.csv"
            ;;
        human_rna)
            INPUT_CSV="result/main_results_hs_rna_with_scores.csv"
            ;;
        human_dnam)
            INPUT_CSV="result/main_results_hs_dnam_with_scores.csv"
            ;;
    esac

    if [[ ! -f "$INPUT_CSV" ]]; then
        echo "ERROR: expected scored CSV not found at $INPUT_CSV. Run the score stage first." >&2
        exit 3
    fi

    "$PYTHON" src/autonomous_post_audit.py \
        --mode plausibility \
        --input "$INPUT_CSV" \
        --output-dir result/case_audit_results \
        --aggregate "result/dir_audit_reaudit_${MODALITY}.json"
fi

# -----------------------------------------------------------------------------
# Stage 4 — Plot: paper figures from the audited result tables
# -----------------------------------------------------------------------------

if run_stage plot; then
    log "Stage 4 / plot — modality=${MODALITY}"
    case "$MODALITY" in
        mouse_rna)
            "$PYTHON" analysis_1_reproducibility.py
            "$PYTHON" analysis_2_genetic_intervention_types.py
            "$PYTHON" analysis_3_pathway_enrichment.py
            "$PYTHON" analysis_4_fda_approved_drugs.py
            ;;
        human_rna|human_dnam)
            "$PYTHON" analysis_1_reproducibility_human.py
            "$PYTHON" analysis_2_genetic_intervention_types_human.py
            "$PYTHON" analysis_3_pathway_enrichment_human.py
            "$PYTHON" analysis_4_fda_approved_drugs_human.py
            ;;
    esac
    "$PYTHON" create_combined_modality_volcano.py || true
fi

log "Pipeline complete for modality=${MODALITY} stages=${STAGES}"
