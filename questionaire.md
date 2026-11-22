# Aging Intervention Evidence Assessment for Prioritization

*For systematic evaluation of drugs, genetic modifications, environmental factors, and disease-related interventions in aging research*

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

- **Immortalized or cancer cell lines** (e.g., HeLa, HEK293, MCF7) → **0 points**
- **Primary cell cultures** (e.g., primary fibroblasts, primary neurons) → **15 points**
- **Organoids or 3D culture systems** (e.g., brain organoids, liver organoids) → **18 points**
- **Ex vivo tissue models** (e.g., tissue slices, isolated organs) → **20 points**
- **In vivo studies** (e.g., mice, rats, humans) → **25 points**

---

## Section 2: Aging Biology Relevance (Max 25 points)

**Q2.1 Does the intervention target or affect known aging-related pathways? (15 points)**
*Instructions: Look for mentions of: mTOR, sirtuins, NAD+, autophagy, senescence, telomeres, mitochondrial function, proteostasis, epigenetic regulation, insulin/IGF signaling, DNA repair, inflammation, stem cell exhaustion, oxidative stress, ECM remodeling, unfolded protein response, circadian rhythm, and other aging-related pathways*

- **No mention** of aging-related pathways → **0 points**
- **Indirect effect** on aging pathways → **5 points**
- **Direct modulation** of one established aging pathway → **10 points**
- **Affects multiple aging pathways** or targets aging hallmarks → **15 points**

**Q2.2 Evidence linking intervention to age-related diseases (7 points)**
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

**Q6.** What is the specific intervention? (Provide a CONCISE label for figure annotation)

**Q7.** What experimental model/species/tissue was used?

---

## Final Scoring Calculation:
**Base Score = Sum of Sections 1-5 (Max 100 points)**
**Final Score = Base Score + Plausibility Penalty (can go to 0 if disqualified)**

---

## JSON Output Format:
```json
{
  "q0_plausibility": {
    "penalty_or_bonus_applied": 0,
    "rationale": "string"
  },
  
  "q1_experimental_model": {
    "points": 0,
    "rationale": "string"
  },
  
  "q2_1_aging_pathways": {
    "points": 0,
    "rationale": "string"
  },
  
  "q2_2_age_diseases": {
    "points": 0,
    "rationale": "string"
  },
  
  "q3_1_experimental_rigor": {
    "points": 0,
    "rationale": "string"
  },
  
  "q4_1_translational_potential": {
    "points": 0,
    "rationale": "string"
  },
  
  "q5_1_literature_saturation": {
    "points": 0,
    "rationale": "string"
  },
  
  "q5_2_mechanistic_novelty": {
    "points": 0,
    "rationale": "string"
  },
  
  "q6_intervention": "string",
  "q7_experimental_model": "string",
  
  "final_score": 0,
}
```

---

## Special Instructions:

### LLM Scoring Protocol:
1. **ALWAYS start with plausibility assessment**
2. Read entire study description thoroughly
3. Extract specific textual evidence for each score
4. When uncertain between tiers, choose the higher score
5. Document evidence and reasoning for each assigned score
6. Flag any concerns about study quality or interpretation

### Red Flags for Plausibility:
- Age-related diseases presented as anti-aging interventions
- Impossibly large effect sizes without proper validation
- Confusion between development and aging
- Cancer cell studies claiming general aging effects
