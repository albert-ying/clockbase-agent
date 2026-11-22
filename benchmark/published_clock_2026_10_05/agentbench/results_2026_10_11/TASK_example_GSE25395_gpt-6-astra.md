You are running as a command-line coding agent in the directory /Entropy/kying/clock-agent-runs-2026-10-11/codex-gpt-6-astra/GSE25395.
Write Python scripts in this directory and run them with `/home/kying/anaconda3/bin/python` (pandas, scipy and statsmodels are installed).
Read only the files under /Entropy/kying/clock-agent-runs-2026-10-11/codex-gpt-6-astra/GSE25395/data. Do not read other locations and do not use the network.
Where the task below says not to execute external commands or create temporary files, it means: no R, no downloads
and no files outside this directory. Running your own Python scripts here is expected.

Analyze the aging effects of conditions from study GSE25395 using CODEX_GPT_6_ASTRA.

Workflow:

1. Data Acquisition and Preparation
    * Load dataset overview from: `/Entropy/kying/clock-agent-runs-2026-10-11/codex-gpt-6-astra/GSE25395/data/dataset_level/gse_me_hs.csv` (Filter by "accession" column to locate GSE25395)
    * Extract metadata from: `/Entropy/kying/clock-agent-runs-2026-10-11/codex-gpt-6-astra/GSE25395/data/GEO_hs_dnam/processed_metadata_preharmonize/GSE25395.csv`
    * Import pre-calculated aging clock data from: `/Entropy/kying/clock-agent-runs-2026-10-11/codex-gpt-6-astra/GSE25395/data/clock_output/human_dnam_age_short.csv`
    * Create output directory if needed: `/Entropy/kying/clock-agent-runs-2026-10-11/codex-gpt-6-astra/GSE25395/results/`
    * **IMPORTANT: DO NOT use R, GEOquery, or any external data fetching. Use ONLY the provided CSV files and context information.**

2. Study Design Information (Use this information instead of fetching from GEO):
   {'title': 'Altered DNA Methylation in Leukocytes with Trisomy 21', 'summary': 'The primary abnormality in Down syndrome (DS), trisomy 21, is well known, but how this chromosomal gain produces the complex DS phenotype, including immune system defects, is not well understood. We profiled DNA methylation in total peripheral blood leukocytes (PBL) and T-lymphocytes from adults with DS and normal controls and found gene-specific abnormalities of CpG methylation in DS, with many of the differentially methylated genes having known or predicted roles in lymphocyte development and function.\n\nNOTE 1: the analysis in Kerkel et al., PLOS Genet, 2010 was carried out after removing probes for genes on the X and Y chromosomes.  NOTE 2: these GEO data have been analyzed with this latest version of BeadStudio; the data in the supplementary tables of Kerkel et al., PLOS Genetics, 2010, were from these same Bead chips, but extracted using an earlier version of this software.  This re-extraction does not affect the conclusions in Kerkel et al., PLOS Genetics, 2010.', 'overall_design': 'Profiling of CpG methylation status in gene promoter regions using Illumina Infinium BeadChips.'}

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
    * gse_id: The dataset identifier (GSE25395)
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
    * model_version: CODEX_GPT_6_ASTRA

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
    - Save to: /Entropy/kying/clock-agent-runs-2026-10-11/codex-gpt-6-astra/GSE25395/results/aging_analysis.csv
    - DO NOT generate reports, visualizations, or any other output.
    - DO NOT provide explanations or summaries in the response.
    - DO NOT use R scripts, GEOquery, subprocess calls, or external data fetching.
    - DO NOT create temporary files or execute external commands.
    - Use ONLY the provided CSV files and study design information.
    - Return ONLY a short confirmation message when complete.

Additional deliverable (required)
    - Save /Entropy/kying/clock-agent-runs-2026-10-11/codex-gpt-6-astra/GSE25395/results/comparison_samples.csv with the columns comparison_name, condition_name, arm, gsm_id.
    - One row per sample of every comparison in the results table; arm is `case` (the condition) or `control`.
    - comparison_name and condition_name must match the results table exactly.
    - For a paired or covariate-adjusted analysis, list the samples that entered the test.
