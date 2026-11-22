"""Human RNA driver with the task template recovered from all 9546 saved traces.

Earlier repository templates are in git history; no exact SONNET37 version
was found there. Runtime paths are relocated through audit_system.paths.
"""
import os, sys, datetime
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_system import paths

llm_name = "claude-3-7-sonnet-20250219"

BASE_DIR = str(paths.REPO)
RESULTS_DIR = f"{BASE_DIR}/result_hs_rna"
CLOCK_FILE = paths.GEO / 'result/human_rna_age.csv'
TRACE_SOURCE = 'successful_traces_human_rna.pkl'
TRACE_SHA256 = 'b1680f27ebfe9d912a11b02ca65cd854a52b97229c87131db76be9cbece4da14'
TRACE_SOURCE_MTIME_NS = 1789833094355027000
TRACE_COUNTS = {'sonnet37': 9546}
TRACE_COMMITS = {'sonnet37': None}

TASK_TEMPLATE = """
Analyze the aging effects of conditions from study {study_id} using SONNET37.

Workflow:

1. Data Acquisition and Preparation
    * Load dataset overview from: `{data_dir}/dataset_level/gse_rna_hs.csv` (Filter by "accession" column to locate {study_id})
    * Extract metadata from: `{data_dir}/GEO_hs_rna/processed_metadata_preharmonize/{study_id}.csv`
    * Import pre-calculated aging clock data from: `{data_dir}/clock_output/human_rna_age_short.csv`
    * Create output directory if needed: `{RESULTS_DIR}/{study_id}/`
    * **IMPORTANT: DO NOT use R, GEOquery, or any external data fetching. Use ONLY the provided CSV files and context information.**

2. Study Design Information (Use this information instead of fetching from GEO):
   {study_design}

3. Clock Selection and Analysis:
    * **IMPORTANT: DO NOT compare young and old populations directly**
    * **Target Aging Clocks:**
    * `Hazard2_All_All_All_scaleddiff_EN_tAge`

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
    * clock: Specific aging clock used
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

Output Requirements
    - Save to: {RESULTS_DIR}/{study_id}/aging_analysis_sonnet37.csv
    - DO NOT generate reports, visualizations, or any other output.
    - DO NOT provide explanations or summaries in the response.
    - DO NOT use R scripts, GEOquery, subprocess calls, or external data fetching.
    - DO NOT create temporary files or execute external commands.
    - Use ONLY the provided CSV files and study design information.
    - Return ONLY a short confirmation message when complete.
"""

TASK_TEMPLATES = {'sonnet37': TASK_TEMPLATE}


def generate_task_intention(study_id, study_design):
    """
    Generate the task intention text for a specific study (Human RNA)

    Args:
        study_id (str): The GSE ID to analyze
        study_design (dict): Metadata for the study

    Returns:
        str: Task intention text for the agent
    """
    data_dir = str(paths.GEO)
    task = TASK_TEMPLATE.format(study_id=study_id, study_design=study_design,
                                data_dir=data_dir, RESULTS_DIR=RESULTS_DIR)
    task = task.replace(f"{data_dir}/dataset_level/gse_rna_hs.csv", str(paths.GEO / 'gse_rna_hs.csv'))
    return task.replace(f"{data_dir}/clock_output/human_rna_age_short.csv", str(CLOCK_FILE))


if __name__ == "__main__":
    import warnings, pandas as pd
    from src.coding_agent import CodingAgent
    from src.parallel_processing import (
        setup_signal_handlers, save_results, parallel_process_studies, print_summary
    )

    warnings.filterwarnings("ignore")
    agent = CodingAgent(llm=llm_name)
    agent.configure()

    # Create the timestamp for this run
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Initialize results directory
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    # Initialize results dictionary
    all_results = {}
    
    # Setup signal handlers for graceful interruption
    setup_signal_handlers(save_results, all_results, timestamp, RESULTS_DIR)
    
    # Inform about checkpointing
    print(f"Results will be periodically saved to {RESULTS_DIR}")
    print("Press Ctrl+C at any time to safely save partial results and exit")
    
    # Load the dataset and randomly select 100 unique GSE IDs
    try:
        # Read the dataset
        df = pd.read_csv(CLOCK_FILE)
        
        # Get unique GSE IDs
        unique_gse_ids = df["gse_id"].unique()
        
        # Convert GSE IDs to numeric values for sorting
        gse_numeric = pd.DataFrame({
            'gse_id': unique_gse_ids,
            'numeric_val': [int(gse_id.replace('GSE', '')) for gse_id in unique_gse_ids]
        })
        
        # Sort by numeric values in descending order and get top 100
        sorted_gse = gse_numeric.sort_values('numeric_val', ascending=False)
        # top_gse = sorted_gse['gse_id'].head(1000).tolist()
        top_gse = sorted_gse['gse_id'].iloc[:200].tolist()
        
        # Add the specific GSE IDs to ensure they're included
        study_ids = top_gse # + ["GSE122080", "GSE139389"]
        # Remove duplicates while preserving order
        study_ids = list(dict.fromkeys(study_ids))
        
        print(f"Selected top {len(study_ids)} GSE IDs from {len(unique_gse_ids)} available IDs")
    except Exception as e:
        print(f"Error selecting GSE IDs: {str(e)}")
        # Fallback to example list if there's an error
        study_ids = ["GSE122080", "GSE139389"]
    
    # Process all studies in parallel
    max_workers = 3
    print(f"Starting parallel processing with {max_workers} workers for {len(study_ids)} studies...")
    all_results = parallel_process_studies(
        study_ids, agent, generate_task_intention, RESULTS_DIR, timestamp, max_workers
    )
    
    # Print summary and save final results
    print_summary(all_results)
    save_results(all_results, timestamp, RESULTS_DIR, final=True)
