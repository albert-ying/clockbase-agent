"""Human DNAm CodingAgent driver reconstructed from parallel_hs_rna.py and saved traces.

A search on 2026-10-04 of this repository's history and of the server's data and
home directories (code files; data trees and caches excluded) did not find the
original driver. TASK_TEMPLATE preserves the prompt in
successful_traces_human_dnam.pkl, GSE211704, result[0][0], after removing the
Human Message wrapper and substituting the study ID, design and path fields.
All 2296 saved traces have this template. The original launch selection and
concurrency settings are unavailable; this driver follows the RNA dispatcher.
Runtime input paths are relocated through audit_system.paths; the five-clock
panel and aging_analysis_sonnet37.csv output name retain the saved wording.
"""
import os, sys, datetime
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_system import paths


llm_name = "claude-3-7-sonnet-20250219"
BASE_DIR = str(paths.REPO)
RESULTS_DIR = f"{BASE_DIR}/results/hs_dnam_sonnet37"
CLOCK_FILE = paths.GEO / 'result/qced_human_dnam_age.csv'
TRACE_SOURCE = 'successful_traces_human_dnam.pkl'
TRACE_KEY = 'GSE211704'
TRACE_SHA256 = '7207f7368782d23416f63a15c9937fdbfc1d824afd14c73b01e86097b8176548'

TASK_TEMPLATE = """
Analyze the aging effects of conditions from study {study_id} using SONNET37.

Workflow:

1. Data Acquisition and Preparation
    * Load dataset overview from: `{data_dir}/dataset_level/gse_me_hs.csv` (Filter by "accession" column to locate {study_id})
    * Extract metadata from: `{data_dir}/GEO_hs_dnam/processed_metadata_preharmonize/{study_id}.csv`
    * Import pre-calculated aging clock data from: `{data_dir}/clock_output/human_dnam_age_short.csv`
    * Create output directory if needed: `{RESULTS_DIR}/{study_id}/`
    * **IMPORTANT: DO NOT use R, GEOquery, or any external data fetching. Use ONLY the provided CSV files and context information.**

2. Study Design Information (Use this information instead of fetching from GEO):
   {study_design}

3. Clock Selection and Analysis:
    * **IMPORTANT: DO NOT compare young and old populations directly**
    * **Target Aging Clocks (DNA Methylation) - Focus on these 5 clocks only:**
    * `Horvathv1` - Primary Horvath epigenetic clock
    * `PhenoAge` - Phenotypic age predictor
    * `YingDamAge` - Ying damage-associated aging clock
    * `DunedinPACE` - Pace of aging measure
    * `GrimAgeV2` - Mortality risk predictor (version 2)

    * **For Each Condition:**
    * Identify appropriate control group(s) based on metadata
    * Determine relevant covariates (age, sex, tissue type, batch effects)
    * Focus on comparing treatment vs. control within the same age group
    * Before statistical testing, determine the experiment's trust score (0-10):
        Trust score evaluates the experimental validity and translational potential of an aging study based on intervention/perturbation type and biological plausibility, completely independent of statistical outcomes like p-values and effect sizes.

        * **0-3: Low value**
        * Typically immortalized/cancer cell lines 
        * Studies with limited connection to aging biology
        * Minimal literature precedent for mechanism

        * **4-6: Moderate value**
        * Well-designed in vitro studies or simple model organisms
        * Reasonable connection to known aging pathways
        * Some supporting evidence in aging research

        * **7-8: High value**
        * In vivo studies or promising in vitro work with strong mechanistic basis
        * Clear connection to established aging biology
        * Substantial literature support
        * Good translational potential

        * **9-10: Exceptional value**
        * Novel interventions with compelling aging-related mechanisms
        * Robust evidence across experimental systems
        * Exceptional translational potential for human applications
    * Perform statistical analysis:
        * Use appropriate tests based on data structure (t-test, ANOVA, or linear models)
        * Report raw p-values (no adjustment required)
        * Calculate standardized effect sizes

4. Deliverables:
* Standardized results table with these columns:
    * gse_id: The dataset identifier ({study_id})
    * comparison_name: Descriptive name of the comparison (e.g., "Drug X vs Vehicle")
    * condition_name: Concise intervention name (e.g., "Rapamycin")
    * condition_category: Categorization (one of: Genetic, Drug, Diet, Environment, Disease, Other)
    * covariate: Variables controlled for in analysis
    * clock: Specific aging clock used (e.g., "Horvathv1", "PhenoAge", "GrimAgeV2")
    * statistical_test: Method used (e.g., "ANOVA", "Linear Mixed Model")
    * p_value: Raw p-value (do NOT use adjusted p-values)
    * effect_size: Magnitude of effect (Cohen's d or similar)
    * direction: Age effect direction (1=accelerated aging, -1=decelerated aging)
    * trust_score: Evaluate experimental follow-up value (0-10)
    * model_version: SONNET37

Implementation Requirements:
    * Use Python with pandas, scipy, and statsmodels
    * Implement modular code structure with separate functions for:
    * Data loading and preprocessing
    * Statistical analysis
    * Results formatting
    * Include efficient documentation and inline comments
    * Handle missing data and outliers appropriately
    * Consider platform-specific effects in DNA methylation data

Output Requirements
    - Save to: {RESULTS_DIR}/{study_id}/aging_analysis_sonnet37.csv
    - DO NOT generate reports, visualizations, or any other output.
    - DO NOT provide explanations or summaries in the response.
    - DO NOT use R scripts, GEOquery, subprocess calls, or external data fetching.
    - DO NOT create temporary files or execute external commands.
    - Use ONLY the provided CSV files and study design information.
    - Return ONLY a short confirmation message when complete.
"""


def generate_task_intention(study_id, study_design):
    """
    Generate the task intention text for a specific study (Human DNAm)

    Args:
        study_id (str): The GSE ID to analyze
        study_design (dict): Metadata for the study

    Returns:
        str: Task intention text for the agent
    """
    data_dir = str(paths.GEO)
    task = TASK_TEMPLATE.format(study_id=study_id, study_design=study_design,
                                data_dir=data_dir, RESULTS_DIR=RESULTS_DIR)
    task = task.replace(f"{data_dir}/dataset_level/gse_me_hs.csv", str(paths.GEO / 'gse_me_hs.csv'))
    return task.replace(f"{data_dir}/clock_output/human_dnam_age_short.csv", str(CLOCK_FILE))


if __name__ == "__main__":
    import pandas as pd
    from src.coding_agent import CodingAgent
    from src.parallel_processing import (
        setup_signal_handlers, save_results, parallel_process_studies, print_summary
    )

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    all_results = {}
    try:
        if not CLOCK_FILE.is_file():
            raise FileNotFoundError(f"Cannot proceed - required data file not found: {CLOCK_FILE}")
        df = pd.read_csv(CLOCK_FILE, usecols=['gse_id'])
        unique_gse_ids = df['gse_id'].drop_duplicates().tolist()
        study_ids = sorted(unique_gse_ids, key=lambda gse_id: int(gse_id.replace('GSE', '')), reverse=True)
        study_ids = [gse_id for gse_id in study_ids
                     if not (Path(RESULTS_DIR) / gse_id / 'aging_analysis_sonnet37.csv').exists()]
        os.makedirs(RESULTS_DIR, exist_ok=True)
    except (OSError, ValueError) as e:
        raise SystemExit(f"Error preparing human DNAm studies: {e}") from e

    agent = CodingAgent(llm=llm_name)
    agent.configure()
    setup_signal_handlers(save_results, all_results, timestamp, RESULTS_DIR)
    print(f"Results will be periodically saved to {RESULTS_DIR}")
    print("Press Ctrl+C at any time to safely save partial results and exit")
    max_workers = 3
    print(f"Starting parallel processing with {max_workers} workers for {len(study_ids)} studies...")
    all_results = parallel_process_studies(
        study_ids, agent, generate_task_intention, RESULTS_DIR, timestamp, max_workers
    )
    print_summary(all_results)
    save_results(all_results, timestamp, RESULTS_DIR, final=True)
