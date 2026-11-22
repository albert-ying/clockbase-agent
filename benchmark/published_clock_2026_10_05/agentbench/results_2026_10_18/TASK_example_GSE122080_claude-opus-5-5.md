You are running as a command-line coding agent in the directory /Entropy/kying/clock-agent-runs-ranked-claude-2026-10-18/claude-claude-opus-5-5/GSE122080.
Write Python scripts in this directory and run them with `/home/kying/anaconda3/bin/python` (pandas, scipy and statsmodels are installed).
Read only the files under /Entropy/kying/clock-agent-runs-ranked-claude-2026-10-18/claude-claude-opus-5-5/GSE122080/data. Do not read other locations and do not use the network.
Where the task below says not to execute external commands or create temporary files, it means: no R, no downloads
and no files outside this directory. Running your own Python scripts here is expected.

Analyze the aging effects of conditions from study GSE122080.

Workflow:

1. Data Acquisition and Preparation
    * Load dataset overview from: `/Entropy/kying/clock-agent-runs-ranked-claude-2026-10-18/claude-claude-opus-5-5/GSE122080/data/dataset_level/gse_rna_mm.csv` (Filter by "accession" column to locate GSE122080)
    * Extract metadata from: `/Entropy/kying/clock-agent-runs-ranked-claude-2026-10-18/claude-claude-opus-5-5/GSE122080/data/GEO_mm_rna/processed_metadata_preharmonize/GSE122080.csv`
    * Import pre-calculated aging clock data from: `/Entropy/kying/clock-agent-runs-ranked-claude-2026-10-18/claude-claude-opus-5-5/GSE122080/data/clock_output/mouse_rna_age_short.csv`
    * Create output directory if needed: `/Entropy/kying/clock-agent-runs-ranked-claude-2026-10-18/claude-claude-opus-5-5/GSE122080/results/`


2. Study Design Information:
   {'title': 'Transcriptomics analysis of the livers of old mice treated with Ouabain', 'summary': 'Methods: Liver mRNA profiles of 7 weeks old wild-type (WT) and 24 months old mice treated with Vehicle or Ouabain (1mg/Kg)', 'overall_design': 'A cohort of 16 female mice 21-22 months old were treated with Vehicle or Ouabain (1mg/Kg) for 12 weeks. Mice were treated 4 times a week on Weeks 0,3,6,9 and 12. Therefore they rested from tretament 2 weeks before a new round of treatment. Mice were culled 3 days after last round of treatment'}

3. Clock Selection and Analysis:
    * **IMPORTANT: DO NOT compare young and old populations directly**
    * **IMPORTANT: Analysis must include ALL GSM samples from this GSE study. Use ONLY the provided data files, do not retrieve additional information from online sources.**
    * **Target Aging Clocks:**
    * `Hazard_Rodents_All_All_scaleddiff_EN_tAge`
    * `Hazard_Mouse_All_All_scaleddiff_EN_tAge`

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
    * gse_id: The dataset identifier (GSE122080)
    * comparison_name: Descriptive name of the comparison (e.g., "Drug X vs Vehicle")
    * condition_name: Concise intervention name (e.g., "Rapamycin")
    * condition_category: Categorization (one of: Genetic, Drug, Diet, Environment, Disease, Other)
    * covariate: Variables controlled for in analysis
    * clock: Specific aging clock used
    * statistical_test: Method used (e.g., "ANOVA", "Linear Mixed Model")
    * p_value: Raw p-value (do NOT use adjusted p-values)
    * effect_size: Magnitude of effect (Cohen's d or similar)
    * direction: Age effect direction (1=accelerated aging, -1=decelerated aging)
    * trust_score: Evaluate experimental follow-up value (0-10):


Implementation Requirements:
    * Use Python with pandas, scipy, and statsmodels
    * Implement modular code structure with separate functions for:
    * Data loading and preprocessing
    * Statistical analysis
    * Results formatting
    * Include efficient documentation and inline comments
    * Handle missing data and outliers appropriately

Output Requirements
    - Save to: /Entropy/kying/clock-agent-runs-ranked-claude-2026-10-18/claude-claude-opus-5-5/GSE122080/results/aging_analysis.csv
    - DO NOT generate reports, visualizations, or any other output.
    - DO NOT provide explanations or summaries in the response.
    - Return ONLY a short confirmation message when complete.

Additional deliverable (required)
    - Save /Entropy/kying/clock-agent-runs-ranked-claude-2026-10-18/claude-claude-opus-5-5/GSE122080/results/comparison_samples.csv with the columns comparison_name, condition_name, arm, gsm_id.
    - One row per sample of every comparison in the results table; arm is `case` (the condition) or `control`.
    - comparison_name and condition_name must match the results table exactly.
    - For a paired or covariate-adjusted analysis, list the samples that entered the test.
