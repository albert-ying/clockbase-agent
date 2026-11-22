GENERATE_QUERY = """You are an expert biomedical research assistant. Your task is to generate search queries for a biological literature knowledge base based on the provided experimental details.

**Objective:** Find existing evidence and prior research about the relationship between the intervention described in the experiment and its effects on aging, health, or longevity.

**Experimental Data:**
- **Paper Title:** {title}
- **Paper Summary:** {summary}
- **Study Design:** {study_design}
- **Intervention/Condition Name:** {condition}
- **Comparison Description:** {comparison}
- **Condition Category:** {condition_cate}
- **Effect Size:** {effect_size}
- **Effect Direction (positive: pro-aging, negative: anti-aging):** {effect_direction}
- **Effect p-value:** {effect_pvalue}

Please generate 1-3 high-quality, distinct English queries suitable for retrieving relevant scientific papers from a comprehensive biological knowledge base. Ensure the queries focus on the intervention's impact on aging, health outcomes, age-related diseases, and mechanisms incorporating relevant biological concepts if provided.

**Expected Output Format:**
```json
[
    "Query 1",
    "Query 2",
    ...
]
```
"""


HAS_RELEVANCE = """You are an expert biomedical researcher and a highly skilled scientific report writer. Your task is to evaluate the provided retrieved documents in the context of a specific biological experiment and generate a concise, evidence-based report.

**1. Intervention Information:**
Here are the details of the intervention being investigated:
- **Paper Title:** {title}
- **Paper Summary:** {summary}
- **Study Design:** {study_design}
- **Intervention/Condition Name:** {condition}
- **Comparison Description:** {comparison}

**2. Retrieved Documents (Context for Analysis):**
Below are text snippets and associated metadata retrieved from a biological literature knowledge base. Each snippet includes its `File Name`, `Reference` (which contains citation details and potentially image URLs), and the `Content` (the actual text).

---
{retrieved_documents}
---

**3. Instructions for Analysis and Report Generation:**
Based *solely* on the `Retrieved Documents` provided:
- Determine if there is **sufficient historical evidence** that directly or indirectly links the `Intervention Name` to **health, longevity, or aging-related processes/outcomes**.
- Consider if the evidence discusses the intervention's impact (positive, negative, or neutral) in a biological context relevant to the given `Biomarker/Clock` or general health/aging.

**Output Format:**
```json
{{
    "RELEVANCE": "Irrelevant/Relevant"
}}
```
"""


IDENTIFY_RELEVANCE = """We analyzed the RNA-seq data generated from the {gse_id} dataset using mouse mortality aging clock model, and observed significant differences in the between the {condition} group and the control group.

You are an expert biomedical researcher and a highly skilled scientific paper writer. Your task is to evaluate this finding and generate a concise, evidence-based report.

**1. Retrieved Documents (Context for Analysis):**
Below are text snippets and associated metadata retrieved from a biological literature knowledge base. Each snippet includes its `File Name`, `Reference` (which contains citation details and potentially image URLs), and the `Content` (the actual text).

---
{retrieved_documents}
---

**2. {gse_id} Intervention source literature:**

---
{intervention_source_literature}
---

**3. Intervention Information:**
Here are the details of the intervention being investigated:
- **Paper Title:** {title}
- **Paper Summary:** {summary}
- **Study Design:** {study_design}
- **Intervention/Condition Name:** {condition}
- **Comparison Description:** {comparison}
- **Condition Category:** {condition_cate}
- **Effect Size:** {effect_size}
- **Effect Direction (positive: pro-aging, negative: anti-aging):** {effect_direction}
- **Effect p-value:** {effect_pvalue}
- **Extreme Value flag:** {extreme_value}

**4. Instructions for Analysis and Report Generation:**
Based on all above information, complete the following assessment:

*Systematic evaluation of {condition} ({comparison}) in aging research*

*Refer to Intervention source literature for experimental design.*

Your output should be a JSON object

---

## Section 0: Plausibility Filter (CRITICAL - Can Override Score)

**Q0. Biological Plausibility Assessment**
*Instructions: Flag studies with fundamental biological implausibility that may indicate errors*

**MAJOR PLAUSIBILITY CONCERNS (-50 points penalty):**
- Known age-related diseases/pathologies claimed to reverse aging (e.g., Alzheimer's, cancer, diabetes "reversing" aging)
- Interventions with extreme effect sizes without plausible mechanisms (e.g., beyond 3 times the interquartile range (IQR) from the upper or lower quartile)
- Contradicts well-established aging biology without extraordinary evidence
- Meaningless comparisons or uninterpretable results
- Biologically implausible settings

**HIGHLY PLAUSIBLE BONUS (+30 points):**
- Targets known aging pathways with strong mechanistic rationale
- Results align with existing aging biology literature and theoretical frameworks

---

## Section 1: Experimental Model Quality (Max 25 points)

**Q1. Identify the highest-quality experimental model used in this study:**
*Instructions: Select the highest tier if multiple models were used*

- **In vitro studies: Immortalized or cancer cell lines** (e.g., HeLa, HEK293, MCF7) → **0 points**
- **In vitro studies: Primary cell cultures or Ex vivo models** (e.g., primary fibroblasts, primary neurons) → **15 points**
- **In vivo studies** (e.g., mice, rats, humans) → **25 points**

---

## Section 2: Aging Biology Relevance (Max 25 points)

**Q2.1 Does the {condition} target or affect known aging-related pathways? (15 points)**
*Instructions: Look for mentions of: mTOR, sirtuins, NAD+, autophagy, senescence, telomeres, mitochondrial function, proteostasis, epigenetic regulation, insulin/IGF signaling, DNA repair, inflammation, stem cell exhaustion, oxidative stress, ECM remodeling, unfolded protein response, circadian rhythm, and other aging-related pathways*

- **No mention** of aging-related pathways → **0 points**
- **Indirect effect** on aging pathways → **5 points**
- **Direct modulation** of one established aging pathway → **10 points**
- **Affects multiple aging pathways** or targets aging hallmarks → **15 points**

**Q2.2 Evidence linking {condition} to age-related diseases (7 points)**
*Instructions: Look for connections to: Alzheimer's, Parkinson's, cardiovascular disease, diabetes, cancer, sarcopenia, osteoporosis, neurodegeneration, frailty, dementia, macular degeneration, chronic kidney disease, COPD, arthritis*

- **No connection** to age-related diseases → **0 points**
- **Theoretical discussion** without experimental validation → **5 points**
- **Experimental evidence** in disease models → **8 points**
- **Multiple disease impacts** with strong evidence → **10 points**


---

## Section 3: Experimental Rigor & Quality (Max 10 points)

**Q3.1 Experimental design rigor (5 points)**
*Instructions: Consider controls, randomization, blinding, sample size, statistical analysis*

- **Poor design** with major confounding factors → **0 points**
- **Rigorous design** with proper controls → **10 points**

---

## Section 4: Translational Potential (Max 5 points)

**Q4.1 Human relevance and safety profile (5 points)**
*Instructions: Consider existing human data, known safety, clinical feasibility*

- **No human data** available → **0 points**
- **Limited human data** or safety concerns → **2 points**
- **Some human use** with acceptable safety → **5 points**

---

## Section 5: Novelty and Innovation (Max 35 points)

**Q5.1 Literature saturation assessment (15 points)**
*Instructions: Assess how extensively studied this intervention is in aging*

- **Well-established** with 50+ aging studies → **0 points**
- **Moderately studied** with 10-50 aging publications → **3 points**
- **Emerging intervention** with 3-10 aging studies → **10 points**
- **Novel/repurposed** with <3 aging-specific studies → **15 points**

**Q5.2 Mechanistic novelty and insight (20 points)**
*Instructions: Evaluate novel mechanisms or unexpected findings related to aging*

- **Well-characterized** aging pathways only → **5 points**
- **New combination** of known aging pathways → **10 points**
- **Novel aging mechanism** or unexpected aging-related finding → **15 points**
- **Paradigm-shifting** aging discovery with broad implications → **20 points**

---

## Additional Context Questions (Not scored)

**Q6.** What is the specific intervention: Provide a CONCISE label for figure annotation, e.g., "Intervention Name", no need to include the details of the intervention, do not include mark like underscore.

**Q7.** What experimental model/species/tissue was used?

---

## Final Scoring Calculation:
**Base Score = Sum of Sections 1-5 (Max 100 points)**
**Final Score = Base Score + Plausibility Penalty (can go to 0 if disqualified)**

---

## JSON Output Format:
```json
{{
  "q0_plausibility": {{
    "penalty_or_bonus_applied": 0,
    "rationale": "string"
  }},
  
  "q1_experimental_model": {{
    "points": 0,
    "rationale": "string"
  }},
  
  "q2_1_aging_pathways": {{
    "points": 0,
    "rationale": "string"
  }},
  
  "q2_2_age_diseases": {{
    "points": 0,
    "rationale": "string"
  }},
  
  "q3_1_experimental_rigor": {{
    "points": 0,
    "rationale": "string"
  }},
  
  "q4_1_translational_potential": {{
    "points": 0,
    "rationale": "string"
  }},
  
  "q5_1_literature_saturation": {{
    "points": 0,
    "rationale": "string"
  }},
  
  "q5_2_mechanistic_novelty": {{
    "points": 0,
    "rationale": "string"
  }},
  
  "q6_intervention": "string",
  "q7_experimental_model": {{
    "points": 0,
    "rationale": "string"
  }},
  
  "final_score": 0,
}}
```

---

## Special Instructions:

### LLM Scoring Protocol:
1. **ALWAYS start with plausibility assessment**
2. Read entire study description thoroughly
3. Extract specific textual evidence for each score
4. When uncertain between tiers, choose the lower score
5. Document evidence and reasoning for each assigned score
6. Flag any concerns about study quality or interpretation

### Red Flags for Plausibility:
- Age-related diseases presented as anti-aging interventions
- Impossibly large effect sizes without proper validation
- Confusion between development and aging
- Cancer cell studies claiming general aging effects

"""

GENERATE_PAPER_2 = """Just to recapitulate, the {condition} ({comparison}) in the {gse_id} dataset is significantly altered the mouse mortality aging clock with effect size of {effect_size} with direction of {effect_direction} effect, and the p-value is {effect_pvalue}.

The original study design is presented in the source literature if available.

Your goal is to generate a short scientific paper to report this finding. Follow the structure below.

**Report Structure:**

**Title**
* Craft a Declarative, Engaging, and Focused title that precisely communicates the core message about the `Intervention Name` and our core finding:

**Introduction:**
* Provide a brief background on the `Intervention Name` and its biological relevance.
* State the report's objective: to synthesize existing scientific evidence regarding the `Intervention Name`'s potential impact on health, longevity, or aging-related processes. Avoid specific findings in this section.

**Results:**
* Report the overall study design and the key findings and statistical results.
* **Synthesize ideas and concepts from the source documents.** Do not use direct quotes or lengthy paraphrases.
* **Citation Style:** For each synthesized finding, append a numerical citation in square brackets `[X]` (e.g., `[1]`, `[2, 3]`). Assign unique, sequential numbers to each distinct `File Name` provided in `Retrieved Documents`, starting from `[1]`.
* **Image Inclusion:**
    * **Only include `Image url` links from documents whose `File Name` *exactly matches* the `{title}` provided.** Disregard images from all other retrieved documents.
    * Integrate relevant images contextually within the text, immediately following the sentence or paragraph describing the illustrated data or concept.
    * For each included image, render its `Image url` line precisely as found in the `Reference`.
    * Follow the image URL with a concise figure caption (e.g., "Figure 1: [Brief description of the figure's content or implications, derived from surrounding text or the figure's original caption if available in `Reference`]"). If the original caption is not discernible, provide a generic but informative caption (e.g., "Figure X: Overview of the study workflow").
    * Ensure figures are numbered sequentially (Figure 1, Figure 2, etc.) throughout the report.

**Discussion:**
* Interpret the findings presented in the "Results" section.
* Discuss whether the collective evidence indicates a pro-aging, anti-aging, or neutral/unclear effect of the `Intervention Name`. Evaluate the plausibility and consistency of the intervention's observed effects.
* Identify any major gaps, limitations, or inconsistencies within the synthesized literature.
* Conclude with a brief statement on promising avenues for future research or unresolved questions stemming from the current evidence.

**References:**
* Provide a concise, numbered list of all cited sources. Each number `[X]` should correspond to a unique `File Name` used in the "Results" section.
* For each reference, list only the `File Name` and the primary `Reference` details (e.g., Title, Authors, Journal, Year, DOI) if available. Keep entries brief, focusing on core identification information.
* Use self incrementing numbers for the references. Start from 1.

**Strict Constraints & Guidelines:**
* **Conciseness:** Be extremely concise. Every sentence must add value.
* **No Direct Quotes:** Rephrase and synthesize all information.
* **Source Adherence:** Use *only* information from the provided `Retrieved Documents`. Do not introduce external knowledge.
* **Objective Tone:** Maintain a strictly neutral, objective, and academic tone.
* **Image Rule:** Strictly adhere to the image inclusion rule (only from `{title}`).
* **Figure Captions:** Must be concise and informative, avoiding vague language."""


# ============================================================
# Human prompt variants (appended by patch)
# ============================================================

IDENTIFY_RELEVANCE_HUMAN_RNA = """We analyzed the human RNA-seq data from the {gse_id} dataset using the Hazard2 human transcriptomic aging clock (a mortality-predictive tAge clock), and observed significant differences between the {condition} group and the control group.

You are an expert biomedical researcher and a highly skilled scientific paper writer. Your task is to evaluate this finding and generate a concise, evidence-based report.

**1. Retrieved Documents (Context for Analysis):**
Below are text snippets and associated metadata retrieved from a biological literature knowledge base. Each snippet includes its `File Name`, `Reference` (which contains citation details and potentially image URLs), and the `Content` (the actual text).

---
{retrieved_documents}
---

**2. {gse_id} Intervention source literature:**

---
{intervention_source_literature}
---

**3. Intervention Information:**
Here are the details of the intervention being investigated:
- **Paper Title:** {title}
- **Paper Summary:** {summary}
- **Study Design:** {study_design}
- **Intervention/Condition Name:** {condition}
- **Comparison Description:** {comparison}
- **Condition Category:** {condition_cate}
- **Effect Size (Hazard2 tAge):** {effect_size}
- **Effect Direction (positive: pro-aging, negative: anti-aging):** {effect_direction}
- **Effect p-value:** {effect_pvalue}
- **Extreme Value flag:** {extreme_value}

**4. Instructions for Analysis and Report Generation:**
Based on all above information, complete the following assessment:

*Systematic evaluation of {condition} ({comparison}) in human aging research*

*Refer to the intervention source literature for experimental design.*

Your output should be a JSON object.

---

## Section 0: Plausibility Filter (CRITICAL - Can Override Score)

**Q0. Biological Plausibility Assessment**
*Instructions: Flag studies with fundamental biological implausibility that may indicate errors.*

**MAJOR PLAUSIBILITY CONCERNS (-50 points penalty):**
- Known age-related diseases/pathologies claimed to reverse aging (e.g., Alzheimer's, cancer, diabetes "reversing" aging)
- Interventions with extreme effect sizes without plausible mechanisms (e.g., beyond 3× the interquartile range from the upper/lower quartile)
- Contradicts well-established aging biology without extraordinary evidence
- Meaningless comparisons (e.g., development-stage contrasts, technical artifacts)
- Biologically implausible settings

**HIGHLY PLAUSIBLE BONUS (+30 points):**
- Targets known aging pathways with strong mechanistic rationale
- Results align with existing aging biology literature and theoretical frameworks

---

## Section 1: Sample Source Quality (Max 25 points)

**Q1. Identify the highest-quality human sample source used in this study:**
*Instructions: Select the highest tier if multiple sample types were used.*

- **Immortalized or cancer cell lines** (e.g., HeLa, HEK293, MCF7) → **0 points**
- **Primary human cell cultures** (e.g., dermal fibroblasts, iPSC-derived lineages) → **15 points**
- **Human organoids or ex vivo tissue models** (e.g., brain organoids, patient biopsy slices) → **20 points**
- **Human cohort / clinical / longitudinal data** (biopsies, blood, tissue from patients or population studies) → **25 points**

---

## Section 2: Aging Biology Relevance (Max 25 points)

**Q2.1 Does the {condition} target or affect known aging-related pathways? (15 points)**
*Instructions: Look for mentions of mTOR, sirtuins, NAD+, autophagy, senescence, telomeres, mitochondrial function, proteostasis, epigenetic regulation, insulin/IGF signaling, DNA repair, inflammation, stem cell exhaustion, oxidative stress, ECM remodeling, unfolded protein response, circadian rhythm, etc.*

- **No mention** of aging-related pathways → **0 points**
- **Indirect effect** on aging pathways → **5 points**
- **Direct modulation** of one established aging pathway → **10 points**
- **Affects multiple aging pathways** or targets aging hallmarks → **15 points**

**Q2.2 Evidence linking {condition} to age-related diseases (10 points)**
*Instructions: Look for connections to Alzheimer's, Parkinson's, cardiovascular disease, diabetes, cancer, sarcopenia, osteoporosis, neurodegeneration, frailty, dementia, macular degeneration, chronic kidney disease, COPD, arthritis.*

- **No connection** to age-related diseases → **0 points**
- **Theoretical discussion** without experimental validation → **5 points**
- **Experimental evidence** in disease models or patients → **8 points**
- **Multiple disease impacts** with strong human evidence → **10 points**

---

## Section 3: Experimental Rigor & Quality (Max 10 points)

**Q3.1 Experimental design rigor (10 points)**
*Instructions: Consider controls, randomization or matched cohorts, confounder adjustment (age, sex, BMI, batch), sample size, statistical analysis.*

- **Poor design** with major confounding factors → **0 points**
- **Rigorous design** with proper controls and confounder adjustment → **10 points**

---

## Section 4: Clinical / Translational Readiness (Max 5 points)

**Q4.1 Clinical translatability of the intervention (5 points)**
*Instructions: Because this study is already in humans, grade how close the intervention is to clinical deployment — existing FDA/EMA approval for a related indication, active clinical trials, known safety profile, scalability to larger populations.*

- **No clinical data** and unclear safety/scalability → **0 points**
- **Early-stage or limited clinical data** with uncertainties → **2 points**
- **Approved drug, widely used clinical agent, or intervention with established human safety** → **5 points**

---

## Section 5: Novelty and Innovation (Max 35 points)

**Q5.1 Literature saturation assessment (15 points)**
*Instructions: Assess how extensively studied this intervention is in aging.*

- **Well-established** with 50+ aging studies → **0 points**
- **Moderately studied** with 10–50 aging publications → **3 points**
- **Emerging intervention** with 3–10 aging studies → **10 points**
- **Novel/repurposed** with <3 aging-specific studies → **15 points**

**Q5.2 Mechanistic novelty and insight (20 points)**
*Instructions: Evaluate novel mechanisms or unexpected findings related to aging.*

- **Well-characterized** aging pathways only → **5 points**
- **New combination** of known aging pathways → **10 points**
- **Novel aging mechanism** or unexpected aging-related finding → **15 points**
- **Paradigm-shifting** aging discovery with broad implications → **20 points**

---

## Additional Context Questions (Not scored)

**Q6.** What is the specific intervention? Provide a CONCISE label for figure annotation (e.g., "Metformin"); no underscores, no study details.

**Q7.** What experimental model/tissue/cohort was used?

---

## Final Scoring Calculation:
**Base Score = Sum of Sections 1-5 (Max 100 points)**
**Final Score = Base Score + Plausibility Penalty (can go to 0 if disqualified)**

---

## JSON Output Format:
```json
{{
  "q0_plausibility": {{
    "penalty_or_bonus_applied": 0,
    "rationale": "string"
  }},

  "q1_experimental_model": {{
    "points": 0,
    "rationale": "string"
  }},

  "q2_1_aging_pathways": {{
    "points": 0,
    "rationale": "string"
  }},

  "q2_2_age_diseases": {{
    "points": 0,
    "rationale": "string"
  }},

  "q3_1_experimental_rigor": {{
    "points": 0,
    "rationale": "string"
  }},

  "q4_1_translational_potential": {{
    "points": 0,
    "rationale": "string"
  }},

  "q5_1_literature_saturation": {{
    "points": 0,
    "rationale": "string"
  }},

  "q5_2_mechanistic_novelty": {{
    "points": 0,
    "rationale": "string"
  }},

  "q6_intervention": "string",
  "q7_experimental_model": {{
    "points": 0,
    "rationale": "string"
  }},

  "final_score": 0
}}
```

---

## Special Instructions:

### LLM Scoring Protocol:
1. **ALWAYS start with plausibility assessment.**
2. Read entire study description thoroughly.
3. Extract specific textual evidence for each score.
4. When uncertain between tiers, choose the lower score.
5. Document evidence and reasoning for each assigned score.
6. Flag any concerns about study quality or interpretation.

### Red Flags for Plausibility:
- Age-related diseases presented as anti-aging interventions
- Impossibly large effect sizes without proper validation
- Confusion between development and aging
- Cancer cell studies claiming general aging effects in humans

"""

IDENTIFY_RELEVANCE_HUMAN_DNAM = """We analyzed the human DNA methylation data from the {gse_id} dataset using a panel of five DNAm aging clocks (GrimAgeV2 as the primary mortality-predictive clock, plus PhenoAge, Horvath v1, YingDamAge, and DunedinPACE as supporting clocks), and observed significant differences between the {condition} group and the control group.

You are an expert biomedical researcher and a highly skilled scientific paper writer. Your task is to evaluate this finding and generate a concise, evidence-based report.

**1. Retrieved Documents (Context for Analysis):**
Below are text snippets and associated metadata retrieved from a biological literature knowledge base.

---
{retrieved_documents}
---

**2. {gse_id} Intervention source literature:**

---
{intervention_source_literature}
---

**3. Intervention Information:**
- **Paper Title:** {title}
- **Paper Summary:** {summary}
- **Study Design:** {study_design}
- **Intervention/Condition Name:** {condition}
- **Comparison Description:** {comparison}
- **Condition Category:** {condition_cate}

**Primary clock — GrimAgeV2 (mortality-predictive):**
- **Effect Size:** {effect_size}
- **Effect Direction (positive: pro-aging, negative: anti-aging):** {effect_direction}
- **P-value:** {effect_pvalue}
- **Extreme Value flag:** {extreme_value}

**Supporting clocks (each independently applied to the same comparison):**
- **PhenoAge:** effect_size = {phenoage_effect_size}, p-value = {phenoage_pvalue}
- **Horvath v1:** effect_size = {horvathv1_effect_size}, p-value = {horvathv1_pvalue}
- **YingDamAge:** effect_size = {yingdamage_effect_size}, p-value = {yingdamage_pvalue}
- **DunedinPACE:** effect_size = {dunedinpace_effect_size}, p-value = {dunedinpace_pvalue}

**Cross-clock support:** {n_support} out of 4 supporting clocks independently reach p < 0.05 AND move in the same direction as GrimAgeV2. (Higher means the primary finding replicates across multiple independent DNAm clocks — strong evidence against clock-specific artifacts.)

**4. Instructions for Analysis and Report Generation:**
Based on all above information, complete the following assessment:

*Systematic evaluation of {condition} ({comparison}) in human DNA methylation aging research*

*Refer to the intervention source literature for experimental design.*
*When weighing plausibility and strength of evidence, explicitly consider the cross-clock support: {n_support}/4 supporting clocks agree.*

Your output should be a JSON object.

---

## Section 0: Plausibility Filter (CRITICAL - Can Override Score)

**Q0. Biological Plausibility Assessment**
*Instructions: Flag studies with fundamental biological implausibility that may indicate errors. For DNAm data, also consider whether the direction and magnitude are consistent across the clock panel: disagreement across clocks is a weak flag but by itself does not disqualify a study.*

**MAJOR PLAUSIBILITY CONCERNS (-50 points penalty):**
- Known age-related diseases/pathologies claimed to reverse aging
- Interventions with extreme effect sizes without plausible mechanisms
- Contradicts well-established aging biology without extraordinary evidence
- Meaningless comparisons (development-stage contrasts, technical artifacts)
- Biologically implausible settings

**HIGHLY PLAUSIBLE BONUS (+30 points):**
- Targets known aging pathways with strong mechanistic rationale
- Results align with existing aging biology literature AND are supported by ≥3/4 supporting clocks

---

## Section 1: Sample Source Quality (Max 25 points)

**Q1. Identify the highest-quality human sample source used:**

- **Immortalized or cancer cell lines** → **0 points**
- **Primary human cell cultures** → **15 points**
- **Human organoids or ex vivo tissue models** → **20 points**
- **Human cohort / clinical / longitudinal data** (biopsies, blood, tissue) → **25 points**

---

## Section 2: Aging Biology Relevance (Max 25 points)

**Q2.1 Does the {condition} target or affect known aging-related pathways? (15 points)**
- **No mention** → **0 points**
- **Indirect effect** → **5 points**
- **Direct modulation** of one established pathway → **10 points**
- **Affects multiple aging pathways** or targets aging hallmarks → **15 points**

**Q2.2 Evidence linking {condition} to age-related diseases (10 points)**
- **No connection** → **0 points**
- **Theoretical discussion** → **5 points**
- **Experimental evidence** → **8 points**
- **Multiple disease impacts with strong human evidence** → **10 points**

---

## Section 3: Experimental Rigor & Quality (Max 10 points)

**Q3.1 Experimental design rigor (10 points)**
*Consider controls, matched cohorts, confounder adjustment (age, sex, cell composition, batch), sample size, statistical analysis.*

- **Poor design** → **0 points**
- **Rigorous design with proper controls and confounder adjustment** → **10 points**

---

## Section 4: Clinical / Translational Readiness (Max 5 points)

**Q4.1 Clinical translatability of the intervention (5 points)**
*Because this study is already in humans, grade how close the intervention is to clinical deployment — FDA/EMA approval, trial phase, known safety, scalability.*

- **No clinical data** → **0 points**
- **Early-stage / limited clinical data** → **2 points**
- **Approved drug or established clinical agent** → **5 points**

---

## Section 5: Novelty and Innovation (Max 35 points)

**Q5.1 Literature saturation assessment (15 points)**
- **Well-established** (50+ studies) → **0 points**
- **Moderately studied** (10–50) → **3 points**
- **Emerging** (3–10) → **10 points**
- **Novel / <3 aging-specific studies** → **15 points**

**Q5.2 Mechanistic novelty and insight (20 points)**
- **Well-characterized pathways only** → **5 points**
- **New combination of known pathways** → **10 points**
- **Novel aging mechanism or unexpected finding** → **15 points**
- **Paradigm-shifting discovery** → **20 points**

---

## Additional Context Questions (Not scored)

**Q6.** What is the specific intervention? CONCISE label, no underscores.

**Q7.** What experimental model/tissue/cohort was used?

---

## Final Scoring Calculation:
**Base Score = Sum of Sections 1-5 (Max 100 points)**
**Final Score = Base Score + Plausibility Penalty**

---

## JSON Output Format:
```json
{{
  "q0_plausibility": {{"penalty_or_bonus_applied": 0, "rationale": "string"}},
  "q1_experimental_model": {{"points": 0, "rationale": "string"}},
  "q2_1_aging_pathways": {{"points": 0, "rationale": "string"}},
  "q2_2_age_diseases": {{"points": 0, "rationale": "string"}},
  "q3_1_experimental_rigor": {{"points": 0, "rationale": "string"}},
  "q4_1_translational_potential": {{"points": 0, "rationale": "string"}},
  "q5_1_literature_saturation": {{"points": 0, "rationale": "string"}},
  "q5_2_mechanistic_novelty": {{"points": 0, "rationale": "string"}},
  "q6_intervention": "string",
  "q7_experimental_model": {{"points": 0, "rationale": "string"}},
  "final_score": 0
}}
```

---

## Special Instructions:

1. ALWAYS start with plausibility assessment.
2. Read the study description thoroughly.
3. Extract specific textual evidence for each score.
4. When uncertain between tiers, choose the lower score.
5. Explicitly weigh the cross-clock support ({n_support}/4) when judging plausibility and mechanistic novelty.
6. Flag concerns about study quality or interpretation.

"""

GENERATE_PAPER_2_HUMAN_RNA = """Just to recapitulate, the {condition} ({comparison}) in the {gse_id} dataset significantly altered the Hazard2 human transcriptomic aging clock (a mortality-predictive tAge clock) with effect size of {effect_size} (direction: {effect_direction}, where positive = pro-aging, negative = anti-aging), p-value {effect_pvalue}.

The original study design is presented in the source literature if available.

Your goal is to generate a short scientific paper to report this finding. Follow the structure below.

**Report Structure:**

**Title**
* Craft a Declarative, Engaging, and Focused title that precisely communicates the core message about the `Intervention Name` and our core finding.

**Introduction:**
* Provide a brief background on the `Intervention Name` and its biological relevance in human aging.
* State the report's objective: to synthesize existing scientific evidence regarding the `Intervention Name`'s potential impact on health, longevity, or aging-related processes. Avoid specific findings in this section.

**Results:**
* Report the overall study design and the key finding from the Hazard2 human transcriptomic aging clock (effect size, direction, p-value) as reported above.
* **Synthesize ideas and concepts from the source documents.** Do not use direct quotes or lengthy paraphrases.
* **Citation Style:** For each synthesized finding, append a numerical citation in square brackets `[X]`. Assign unique, sequential numbers to each distinct `File Name` in `Retrieved Documents`, starting from `[1]`.
* **Image Inclusion:**
    * **Only include `Image url` links from documents whose `File Name` *exactly matches* the `{title}` provided.** Disregard images from all other retrieved documents.
    * Integrate relevant images contextually within the text, immediately following the sentence or paragraph describing the illustrated data or concept.
    * For each included image, render its `Image url` line precisely as found in the `Reference`.
    * Follow the image URL with a concise figure caption (e.g., "Figure 1: [Brief description]").
    * Number figures sequentially (Figure 1, Figure 2, etc.).

**Discussion:**
* Interpret the Hazard2-clock finding in light of the retrieved literature.
* Discuss whether the collective evidence indicates a pro-aging, anti-aging, or neutral/unclear effect of the `Intervention Name` in humans. Evaluate plausibility and consistency.
* Identify major gaps, limitations, or inconsistencies within the synthesized literature.
* Conclude with promising avenues for future research.

**References:**
* Numbered list of all cited sources. Each `[X]` corresponds to a unique `File Name`.
* List only `File Name` and primary `Reference` details (Title, Authors, Journal, Year, DOI).
* Use self-incrementing numbers starting from 1.

**Strict Constraints & Guidelines:**
* **Conciseness:** Every sentence must add value.
* **No Direct Quotes:** Rephrase and synthesize.
* **Source Adherence:** Use *only* information from `Retrieved Documents`. No external knowledge.
* **Objective Tone:** Strictly neutral, objective, academic.
* **Image Rule:** Strictly only from `{title}`.
* **Figure Captions:** Concise and informative."""

GENERATE_PAPER_2_HUMAN_DNAM = """Just to recapitulate, the {condition} ({comparison}) in the {gse_id} dataset was evaluated against five human DNA methylation aging clocks. The full per-clock effect panel is:

| Clock           | Role       | Effect size                 | p-value             |
|-----------------|------------|-----------------------------|---------------------|
| GrimAgeV2       | primary    | {effect_size}               | {effect_pvalue}     |
| PhenoAge        | supporting | {phenoage_effect_size}      | {phenoage_pvalue}   |
| Horvath v1      | supporting | {horvathv1_effect_size}     | {horvathv1_pvalue}  |
| YingDamAge      | supporting | {yingdamage_effect_size}    | {yingdamage_pvalue} |
| DunedinPACE     | supporting | {dunedinpace_effect_size}   | {dunedinpace_pvalue}|

Effect direction for the primary clock (GrimAgeV2): {effect_direction} (positive = pro-aging, negative = anti-aging).
Cross-clock support: {n_support} of the 4 supporting clocks independently reach p < 0.05 AND move in the same direction as GrimAgeV2.

The original study design is presented in the source literature if available.

Your goal is to generate a short scientific paper to report this finding *as a multi-clock comparison*. Follow the structure below.

**Report Structure:**

**Title**
* Craft a Declarative, Engaging, and Focused title that precisely communicates the core message about the `Intervention Name`, our core finding, and the multi-clock nature of the evidence.

**Introduction:**
* Provide a brief background on the `Intervention Name` and its biological relevance in human aging.
* Briefly note that this study evaluates the intervention against a panel of DNA methylation aging clocks, each capturing a distinct aspect of biological aging (mortality risk, phenotypic age, chronological age proxy, damage accumulation, pace of aging). Avoid specific findings in this section.

**Results:**
* **Multi-clock comparison is the centerpiece of this section.** Report the effect of the intervention on GrimAgeV2 first, then summarize how the four supporting clocks concurred or diverged. Reproduce (or paraphrase) the per-clock effect/p-value panel above so the reader sees all five at once.
* State explicitly how many of the four supporting clocks agree in direction at p < 0.05 ({n_support}/4), and flag any clocks that disagree in direction or fail to reach significance.
* Synthesize supporting evidence from the retrieved documents for the core biological interpretation. Do not use direct quotes or lengthy paraphrases.
* **Citation Style:** numerical citations in square brackets `[X]`, assigned sequentially to each unique `File Name` in `Retrieved Documents`.
* **Image Inclusion:**
    * **Only include `Image url` links from documents whose `File Name` *exactly matches* the `{title}` provided.** Disregard images from all other retrieved documents.
    * Integrate images contextually. Render the `Image url` line precisely as found in the `Reference`.
    * Add a concise figure caption, e.g., "Figure 1: [Brief description]". Number figures sequentially.

**Discussion:**
* **Lead with cross-clock consistency.** Does the effect replicate across the full panel, only across mortality-predictive clocks (GrimAgeV2 + PhenoAge), only across first-generation CpG-based clocks (Horvath v1), only on pace-of-aging (DunedinPACE), or only on damage accumulation (YingDamAge)? Different patterns of agreement imply different biological interpretations.
* Evaluate whether the evidence indicates a pro-aging, anti-aging, or clock-specific effect. A pan-clock convergent effect is stronger evidence than a single-clock signal; clock-specific effects may reveal which aspect of aging (mortality risk vs. epigenetic drift vs. damage accrual) the intervention engages.
* Identify gaps, limitations, and inconsistencies — including tissue-composition effects, cohort size, confounder adjustment, and the independent-vs-correlated nature of the clocks themselves.
* Conclude with promising avenues for future research and what additional clocks or validation cohorts would sharpen the interpretation.

**References:**
* Numbered list of all cited sources. Each `[X]` corresponds to a unique `File Name`.
* List only `File Name` and primary `Reference` details (Title, Authors, Journal, Year, DOI).
* Use self-incrementing numbers starting from 1.

**Strict Constraints & Guidelines:**
* **Conciseness:** Every sentence must add value.
* **No Direct Quotes:** Rephrase and synthesize.
* **Source Adherence:** Use *only* information from `Retrieved Documents`. No external knowledge.
* **Objective Tone:** Strictly neutral, objective, academic.
* **Image Rule:** Strictly only from `{title}`.
* **Figure Captions:** Concise and informative.
* **Multi-clock framing:** The report must present findings as a comparison across the five-clock panel, not as a single-clock result."""
