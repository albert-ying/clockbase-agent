# Released ClockBase prompts

The `.txt` files mirror the source templates listed below. Each has one provenance line starting with `# Source:`; every byte after the first newline is the source body, including its leading and trailing whitespace. Edit the source template, then run `python -B prompts/build_prompts.py` (standard library only). Run the same command with `--check` to detect differences without writing files. Regeneration uses Python syntax parsing and does not initialize agents, contact model services, or require the saved traces.

| File | Authoritative source |
|---|---|
| `coding_agent_system.txt` | `src/coding_agent.py::DEFAULT_SYSTEM_PROMPT` |
| `task_intention_human_rna.txt` | `src/parallel_hs_rna.py::TASK_TEMPLATE` |
| `task_intention_mouse_rna.txt` | `src/parallel_mm_rna.py::TASK_TEMPLATE` |
| `task_intention_mouse_rna_all_samples_local.txt` | `src/parallel_mm_rna.py::TASK_TEMPLATE_ALL_SAMPLES_LOCAL` |
| `task_intention_mouse_rna_all_samples.txt` | `src/parallel_mm_rna.py::TASK_TEMPLATE_ALL_SAMPLES` |
| `task_intention_human_dnam.txt` | `src/parallel_hs_dnam.py::TASK_TEMPLATE` |
| `rag_query_generator.txt` | `minipaper_fork/src/nlp_process/prompt_template.py::GENERATE_QUERY` |
| `rag_relevance_classifier.txt` | `minipaper_fork/src/nlp_process/prompt_template.py::HAS_RELEVANCE` |
| `biological_interpreter_mouse_rna.txt` | `minipaper_fork/src/nlp_process/prompt_template.py::IDENTIFY_RELEVANCE` |
| `biological_interpreter_human_rna.txt` | `human_prompts_for_minipaper.py::IDENTIFY_RELEVANCE_HUMAN_RNA` |
| `biological_interpreter_human_dnam.txt` | `human_prompts_for_minipaper.py::IDENTIFY_RELEVANCE_HUMAN_DNAM` |
| `generate_paper_mouse_rna.txt` | `minipaper_fork/src/nlp_process/prompt_template.py::GENERATE_PAPER_2` |
| `generate_paper_human_rna.txt` | `human_prompts_for_minipaper.py::GENERATE_PAPER_2_HUMAN_RNA` |
| `generate_paper_human_dnam.txt` | `human_prompts_for_minipaper.py::GENERATE_PAPER_2_HUMAN_DNAM` |
| `post_audit_opus47.txt` | `prompts/post_audit_templates.py::POST_AUDIT_OPUS47` |
| `case_curator_opus47.txt` | `prompts/post_audit_templates.py::CASE_CURATOR_OPUS47` |
| `verifier_design_reconstruction.txt` | `audit_system/design_reconstruction_2026_09_28.py::PROMPT` |
| `verifier_label_adjudication.txt` | `audit_system/label_adjudication_2026_09_27.py::PROMPT` |
| `verifier_intervention_classification.txt` | `audit_system/intervention_classification_2026_09_29.py::PROMPT` |
| `verifier_relevance_rubric.txt` | `audit_system/intervention_classification_2026_09_29.py::RUBRIC` |
| `verifier_design_audit.txt` | `audit_system/design_audit_2026_09_30.py::PROMPT` |

The task templates are the per-study CodingAgent user messages. The human templates and mouse `all_samples_local`
template forbid external data fetching. The mouse `base` and `all_samples` templates direct the agent to supplied
files without an explicit prohibition. Their braces name the study ID, study-design block and runtime paths. The system prompt comes from `src/coding_agent.py`, a self-contained version of the BioAgentOS CodeAct agent that ran the historical analysis; the saved traces do not contain a system message, so this prompt cannot be checked against the historical runs. The retrieval, scoring and manuscript prompts come from the minipaper sources. Human templates use `human_prompts_for_minipaper.py` as the authoritative source. The retrieval framework's own prompts are not part of this release. The `verifier_*` files are the prompts of the verifier's model stages: design reconstruction of unrecovered rows, orientation adjudication, contrast-type classification with the relevance rubric, and the design audit of ranking candidates; each is sent to two models and only agreement is accepted. This index covers the primary analysis, minipaper, post-audit and verifier model prompts, not the update controller's or external frameworks' prompts.

The two post-audit prompts were authored as `.txt` inputs loaded by `src/autonomous_post_audit.py::load_prompt`, rather than as Python constants. Their unchanged runtime bodies are now preserved in `prompts/post_audit_templates.py`; the builder regenerates the text files consumed by that driver. The driver strips the provenance line before filling the templates.

## RNA trace recovery

The RNA task bodies were recovered from the saved CodingAgent traces under `audit_system.paths.TRACES`. Every first human message matches exactly one released template after substituting its study accession, section-2 study-design text and path fields. Leading and trailing newlines, blank lines and trailing spaces are retained. No traces are unmatched. Earlier repository templates remain in git history; the previously released human RNA template omitted the saved SONNET37 instructions, and the previously released mouse template differed even from the saved base template in whitespace and path fields.

| Mirror / driver selection | Saved traces | Distinguishing instructions | Exact earlier source after path substitution |
|---|---:|---|---|
| `task_intention_human_rna.txt` / human default | 9,546 | `using SONNET37.`, the prohibition on R/GEOquery/external fetching, the section-2 supplied-information note, `model_version: SONNET37`, `aging_analysis_sonnet37.csv`, and the final restrictions on external fetching/commands and provided inputs | Not found |
| `task_intention_mouse_rna.txt` / `base` (default) | 11,385 | No all-GSM instruction | `4c6aa975ed95fda1c2829537aa95d977b102cc99:src/parallel_prompt3.py` (also `a2a930c5d:src/parallel_prompt3.py` and `3faf1f266:src/parallel_prompt3_newapi.py`) |
| `task_intention_mouse_rna_all_samples_local.txt` / `all_samples_local` | 3,578 | `IMPORTANT: Analysis must include ALL GSM samples from this GSE study. Use ONLY the provided data files, do not retrieve additional information from online sources.` | `3faf1f2662ff67b7755adb1d77156d35dfdf3829:src/parallel_prompt3.py` |
| `task_intention_mouse_rna_all_samples.txt` / `all_samples` | 417 | `IMPORTANT: Analysis must include ALL GSM samples from this GSE study` plus `Ensure all GSM samples (individual samples/assays) belonging to this GSE (study) are included in the analysis` | Not found |

`successful_traces_human_rna.pkl` has SHA256 `b1680f27ebfe9d912a11b02ca65cd854a52b97229c87131db76be9cbece4da14`, source mtime_ns `1789833094355027000` (`2026-09-19T15:51:34.355027+00:00`). Its representative trace is `GSE212669`. `successful_traces_mouse.pkl` has SHA256 `5b1271b726458e1e4cbdef00af34d0abbf93c33f5932d9bf50a6456bd66ba24e`, source mtime_ns `1789620277602064600` (`2026-09-17T04:44:37.602065+00:00`). Representative mouse keys are `GSE212644`, `GSE208050` and `GSE210519`, respectively. The study accession is extracted from the prompt, because some mouse dictionary keys differ from the accession in their saved message. All human messages and 15,367 mouse messages occur at `result[0][0]`; the remaining 13 mouse traces start with an AI message and have their first human message at `result[0][1]`. The test locates that first human wrapper before removing it and its blank-line separator.

Human paths use `{data_dir}` and `{RESULTS_DIR}`. Mouse prompts require separate `{RESULTS_DIR}` for directory creation and `{OUTPUT_DIR}` for the saved CSV: their historical values differ. The study-design block is preserved as text before substitution, so braces, escaped newlines and accessions inside the design retain their original representation. No additional study-specific fields were found. The human target is always `Hazard2_All_All_All_scaleddiff_EN_tAge`; both mouse targets are always `Hazard_Rodents_All_All_scaleddiff_EN_tAge` and `Hazard_Mouse_All_All_scaleddiff_EN_tAge`.

The RNA drivers expose literal templates, trace hashes and per-template counts, and can be imported without model clients, credentials or network calls. `generate_task_intention(study_id, study_design)` renders the traced human template or the most common mouse template. Mouse callers can select either alternative with `template_version='all_samples_local'` or `template_version='all_samples'`. Runtime paths use `audit_system.paths.GEO`, with overviews at `GEO/gse_rna_hs.csv` and `GEO/gse_rna_mm.csv`, metadata under the corresponding `GEO_*_rna/processed_metadata_preharmonize/` directory, and clocks at `GEO/result/human_rna_age.csv` and `GEO/result/mouse_rna_age.csv`. Runtime mouse directory creation and CSV output share the configured `RESULTS_DIR`. Existing output roots remain `result_hs_rna/` and `result_mm_rna_new/`; human CSVs now use the traced SONNET37 filename. The templates preserve the historical path layout for exact trace reconstruction. No saved analyses were rerun or changed.

From the release root with the verifier environment active, run `python -B -m unittest audit_system.tests.test_trace_prompts_2026_10_05 audit_system.tests.test_prompts_release_2026_10_04` for corpus equality, counts, hashes, driver rendering and mirror checks. If a deposited trace pickle is absent, its corpus test skips with `saved agent traces not present: <path>`; mirror and runtime-rendering checks still run.

## Human DNAm recovery

`src/parallel_hs_dnam.py` is reconstructed from `src/parallel_hs_rna.py` and the saved DNAm traces. The original driver was not recovered. The recovery search checked the repository history and selected server code directories, excluding data trees, environments, caches, result directories and manuscript files, and skipping code files larger than 2 MB. The matching code inspected in that search was the trace replay extractor; this record does not establish absence from every server directory.

The task body was recovered from `successful_traces_human_dnam.pkl`, trace `GSE211704`, `result[0][0]`. The file SHA256 is `7207f7368782d23416f63a15c9937fdbfc1d824afd14c73b01e86097b8176548`; its source mtime is `2026-09-19T15:51:25.503107+00:00` (mtime_ns `1789833085503106697`). Extraction removes the saved Human Message wrapper, substitutes the study accession with `{study_id}`, the section-2 design block with `{study_design}`, and the historical data and results prefixes with `{data_dir}` and `{RESULTS_DIR}`. All remaining wording and whitespace are preserved. All 2,296 saved prompts match this 4,680-character template after those substitutions. The release test reconstructs each saved message from the template and checks exact equality when the trace file is available.

The reconstructed driver uses `audit_system.paths.GEO` for the overview and metadata, and `GEO/result/qced_human_dnam_age.csv` for the locally available clock table, matching the replay extractor. It retains the five-clock panel, SONNET37 wording and `aging_analysis_sonnet37.csv` output filename. Its dispatcher follows the RNA driver, including GEO study-design retrieval before invoking the agent. The prompt limits the agent to the provided CSVs and context. Launch selection and concurrency cannot be recovered from the traces: the reconstruction selects all study IDs in the configured clock table, skips existing SONNET37 outputs and uses the RNA driver's three-worker dispatcher. It writes under `results/hs_dnam_sonnet37/`. No historical runs were changed or re-executed.
