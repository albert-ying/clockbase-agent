# ClockBase Agent

ClockBase Agent applies aging clocks to the experimental and clinical contrasts deposited in GEO. A CodingAgent designs
and runs case-control comparisons for each study (mouse RNA-seq, human RNA-seq and human DNA methylation); a verification
system recovered and checked samples for 97,233 rows and withheld the 7,139 rows without
recoverable samples; it repairs what can be repaired and records each check, repair and final class. This repository contains the agent, its prompts, the preprocessing and clock code, the verifier with its
tests, and the scripts that rebuild the verified atlas and its ranking.

## Contents

| Component | Location |
|---|---|
| CodingAgent (LangGraph CodeAct loop; a self-contained version of the BioAgentOS CodeAct agent that ran the historical analysis) | `src/coding_agent.py` |
| Study drivers (the human DNAm driver reconstructed from the saved traces) | `src/parallel_mm_rna.py`, `src/parallel_hs_rna.py`, `src/parallel_hs_dnam.py`, `src/parallel_processing.py` |
| Prompts: the task prompts of the saved runs (recovered from the traces and tested against every trace) and the other prompts, as plain-text copies regenerated from the code | `prompts/` (see `prompts/README.md`) |
| Biological interpreter, scorer and manuscript agent | `minipaper_fork/`, `human_prompts_for_minipaper.py` |
| Scripts of the historical GEO pipeline: Git versions last committed before the outputs' modification times, plus later compatible candidates where those versions cannot have written the observed outputs (`as_run/`); the 2026-10-04 snapshot (`current/`); and recovered atlas analysis scripts (`analysis/`, see its README) | `src/historical_geo_pipeline/` (with `MANIFEST.tsv` and the limits of the reconstruction in its README) |
| Revision preprocessing, clock computation and DNAm re-derivation from GEO | `src/pipeline_inputs_2026_09_10.py`, `pipeline_assets/`, `src/dnam_stream_repair_2026_09_26.py` |
| Update controller for new GEO studies (bounded batches) | `src/update_pipeline_2026_09_20.py`, `PIPELINE_2026_09_10.md` |
| Verifier (membership recovery, trace replay, verdicts, repairs with the design-respecting estimand, individuals audit, census, ranking, export) | `audit_system/` |
| One command for the verification | `scripts/clockbase_verifier.py` |
| Benchmark of the verified atlas against published clock results and expert sample annotation (AgingAgentBench, https://github.com/albert-ying/AgingAgentBench; files written before its renaming call it ClockEffectBench), with the agent reruns of the benchmark series and of the ranked comparisons: code, protocol with amendments, and result tables (its data root is described in the benchmark README; the extracted published results, per-row records and every rerun are in the records deposit) | `benchmark/published_clock_2026_10_05/` |
| Replication of intervention effects across independent projects (GEO series grouped from their public records) and agreement with GenAge and DrugAge | `audit_system/intervention_meta_2026_09_29.py` |
| Development notebooks of the historical analysis (01, 03, 04 and 05 read hard-coded `/dfs/user/` server paths; not runnable from the release) | `notebooks/` |
| Historical plausibility filter and scoring rubric, preserved as used (heading errata below) | `questionaire.md` |
| Data locations | `audit_system/paths.py` |
| Release builder and sandboxed release check | `scripts/build_release_2026_10_04.py`, `scripts/release_check_2026_10_04.sh` |
| Builders of the revision's Figure 3 (with Tables 1, 2 and S17) and Extended Data Figures 2 and 21, run from the bundle, and of Figure 4 and Extended Data Figures 22 and 23, run from the records deposit and the benchmark results | `scripts/figures_2026_10_03/` |
| Figure and analysis scripts of the original manuscript (27 of the 34 read historical result tables that are not released; six need only `figures/combined_3modality_plot_data.csv`, which the records deposit holds under `supporting/historical-atlas/figures/`; `create_gse106161_workflow_visualization.py` is an empty file; see the [dependency inventory](PIPELINE_2026_09_10.md#release-runtime-dependency-inventory-2026-10-05)) | `analysis_*.py`, `create_*.py` at the repository root |

## Install

Use Linux x86-64 and Python 3.12; the sandboxed release check also needs curl and bubblewrap. The release check has
also been run on macOS 27 (Apple M5) with Python 3.12. No non-standard hardware is needed.
Install `uv` for the commands below. For new-study preprocessing, R and the full pipeline dependencies, see
[PIPELINE Install](PIPELINE_2026_09_10.md#install). The exact versions of the environment that produced the verifier runs are pinned in
`requirements-verifier.txt`:

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements-verifier.txt \
    --extra-index-url https://download.pytorch.org/whl/cpu --index-strategy unsafe-best-match
```

The verifier environment has 87 packages (about 0.8 GB); a fresh install takes about one minute on a fast connection.
Activate it before the `python` commands below:

```bash
source .venv/bin/activate
```

The pipeline setup creates `.venv-pipeline` from `pipeline_requirements_2026_09_10.txt`; its interpreter can also
run the verifier. No environment record of the Biolearn version that produced the historical DNA methylation clock predictions
survives; the DNA methylation re-derivation and the verifier use Biolearn 0.9.1 (pinned here), and each re-derived
study records the versions it used in its `status.json`.

`questionaire.md` is preserved as used. Two headings understate their options: Q2.2 reads '(7 points)' and Q3.1 reads '(5 points)', but both offer options up to 10, matching the section maxima (25 and 10), the saved scores and Methods/Extended Data Figure 8 (0–10).

## Data

The verifier reads its inputs from `clockbase_data/` (or the locations given by the environment variables in
`audit_system/paths.py`). The rebuild bundle unpacks to this layout:

```
clockbase_data/
  GEO_data/                     sample metadata, DNAm missing-probe tables, RNA imputation QC and series-matrix
                                headers that the rebuild reads (the per-sample predictions it uses are in
                                verifier/replay-*-001/inputs), and the self-test's inputs: per-sample DNAm prediction tables
                                result/qced_human_dnam_age.csv and result/hs_dnam_age.csv.gz, RNA clock outputs of
                                134 human and 137 mouse studies (GEO_hs_rna/clock_result, GEO_mm_rna/clock_result),
                                and three DNAm matrices (GEO_hs_dnam/processed_data)
  verifier/                     the saved stage outputs (membership, replay, links, model decisions with their ledgers,
                                classification, design audits, individuals audit, input-integrity screen) and the reference final run
  dnam-stream-repair-empty/     DNAm re-derivation of the 453 of 513 empty-matrix studies with a beta-table or IDAT source
                                (status, sample QC, clock outputs and coverage; 402 re-derived, 51 withheld);
                                inventory.csv classifies all 513, and the other 60 (non-Infinium or sequencing
                                platforms, signal-only tables or no usable file) were not attempted
  dnam-stream-repair-detached/  DNAm re-derivation of studies whose historical profiles were detached from probes
  corrections/                  the three atlas correction tables
  inputs/                       the empty-matrix census and the published items checked against the atlas
  agent-traces/                 the saved CodingAgent traces and per-study outputs
```

The bundle holds what the rebuild below reads: from each saved stage, the files the rebuild reads and the stage's
completion file; the model-call ledgers; the inputs of the final chain, of the self-test and of the figure builders
(the builders themselves are in the code release); and the reference final run `final-2026-11-13`. The intermediate
pass-1, pass-2 and pass-3 chains, the rest of each stage's files and the self-test's simulated outputs are not included. The raw inputs of the upstream stages are not all
included (the extracted DNA methylation matrices alone take 227 GB), so those stages run only where the GEO-derived
inputs are available; the historical GEO scripts in `src/historical_geo_pipeline/` are provenance copies whose
limits its README states, not a complete way to regenerate them. GEO series matrices are included up to
`!series_matrix_table_begin`, the sample annotations the verifier reads.

The per-comparison verification records (one JSON file per study, with a hash manifest and a schema description) are
the `public-export/` folder of the reference final run.

## Verify the atlas

`REFERENCE` is the reference final run shipped in the bundle (`final-2026-11-13`, recorded as `reference_final`
in `BUNDLE.json` and listed under `clockbase_data/verifier/`). The other `final-*` directories there are partial
copies and cannot serve as `REFERENCE`. `NAME` is any new directory name under that root, for example `final-rebuild`.
`SANDBOX` is a new writable directory for self-test copies and simulated outputs. With the saved stages in the
bundle, the final chain takes about 14 minutes and the planted-defect self-test, at the setting behind the reported
simulation (`--per-modality 160`, seed 20261002, 480 rows), 3 to 4 minutes.

```bash
python scripts/clockbase_verifier.py plan     --final NAME            # list the stages and whether each is reused
python scripts/clockbase_verifier.py run      --final NAME --models cache-only
python scripts/clockbase_verifier.py verify   --run clockbase_data/verifier/NAME --reference clockbase_data/verifier/REFERENCE
python scripts/clockbase_verifier.py selftest --root SANDBOX --reference clockbase_data/verifier/REFERENCE --per-modality 160
```

`run` executes the stages in dependency order and reuses every stage whose completion file exists, so with the bundle it
reruns the deterministic final chain (merge, verdicts, repairs, published-output fates, error census, atlas, ranking,
export) from the saved stage outputs. Model stages call two models independently and accept only agreement; each stage
records its model pair. `--models cache-only` never sends a request; `--models live` reruns model stages through the
APIs for stages selected for execution. Set `TOKENSWITCH_BASE_URL` (the NeoLink API base, `https://neolink.com/api`; the modules append
`/v1/chat/completions`) and `TOKENSWITCH_API_KEY` (the NeoLink key) for design reconstruction, decisive adjudication,
pending/unconfirmed adjudication and classification (the NeoLink `gp-6-as,cl-op-5` pair), and `OPENAI_API_KEY` for
supported-sample adjudication, the design audits and the individuals audit (the OpenAI `gpt-6-astra,gpt-6-sol` pair). Alternatively,
`CLOCKBASE_MODEL_LAUNCHER` supplies a command prefix with optional `{provider}` substitution. Completed stages
are reused even in live mode. Add `--force STAGE` (one or more stage names from `plan`, or `all`) to rerun them;
this neither deletes saved outputs nor bypasses a module's resume checks. Fresh upstream runs need their source
inputs, which are not all shipped in the bundle. See [the stage guide](audit_system/README.md#verifier-v3).
`verify` compares two runs row by row and reports provenance-only differences (timestamps, run paths,
model-call receipts) separately; its report is `NAME/verify.json`. `selftest` plants
known defects into copies of verified rows in a sandbox (a simulation) and reports how many each check catches.
Re-derive saved model decisions from stored raw responses without new calls:

```bash
python -m audit_system.model_decision_audit_2026_10_04 --out DIR
```

`DIR` is a new output directory; the verifier root defaults to `clockbase_data/verifier`. Records whose call
failed carry no response to re-derive.

## Release check

`scripts/release_check_2026_10_04.sh` installs a fresh environment from `requirements-verifier.txt`, then rebuilds the
final run from the code release and the bundle alone inside a bubblewrap sandbox with no network and no other
filesystem, compares the result with the reference run, runs the planted-defect self-test (a simulation), and
runs the unit tests. Run it with the unpacked code release, the bundle's parent directory, a new work directory
and the reference run name:

```bash
bash scripts/release_check_2026_10_04.sh RELEASE BUNDLE WORK REFERENCE
```

Its report is `WORK/release_check.json`; numerical results live in the run outputs. On macOS,
`scripts/release_check_macos_2026_11_12.sh` takes the same four arguments, uses the built-in sandbox in place of
bubblewrap and needs `uv`; the whole check takes about 15 minutes on a laptop.

## Companion repositories

- https://github.com/Ard-Skelling/cheap-RAG (tag 0.1.0, commit 5765880): the retrieval framework used by the scoring layer.
- https://github.com/Gladyshev-Lab/clockbase: the ClockBase web application.

## Citation and license

Citation: Ying, K. et al. Autonomous AI Agents Discover Aging Interventions from Millions of Molecular Profiles. Preprint at
bioRxiv https://doi.org/10.1101/2023.02.28.530532 (2025).

No license has been assigned to this repository yet. `minipaper_fork/` derives from
https://github.com/Ard-Skelling/minipaper and keeps the terms of that repository.
