
## read in the human DNAm data -----
human_dnam_age = dir("./GEO_hs_dnam/clock_result", full.name = T)
human_dnam_missing = dir("./GEO_hs_dnam/missing_percentage", full.name = T)
human_dnam_meta = dir("./hs_dnam_gpt-4", pattern = "dsv", full.name = T)

human_dnam_age_tb = read_csv(human_dnam_age, id = "file")
human_dnam_missing_tb = read_csv(human_dnam_missing, id = "file")
human_dnam_meta_tb = read_delim(human_dnam_meta, id = "file", col_names = F, delim = "|")

colnames(human_dnam_age_tb)
#  [1] "file"               "...1"               "Horvathv1"         
#  [4] "Hannum"             "Lin"                "PhenoAge"          
#  [7] "YingCausAge"        "YingDamAge"         "YingAdaptAge"      
# [10] "Horvathv2"          "PEDBE"              "Zhang_10"          
# [13] "DunedinPoAm38"      "DunedinPACE"        "DNAmADM...14"      
# [16] "DNAmB2M...15"       "DNAmCystatinC...16" "DNAmGDF15...17"    
# [19] "DNAmLeptin...18"    "DNAmPACKYRS...19"   "DNAmPAI1...20"     
# [22] "DNAmTIMP1...21"     "Age...22"           "Female...23"       
# [25] "DNAmGrimAge...24"   "AgeAccelGrim...25"  "DNAmADM...26"      
# [28] "DNAmB2M...27"       "DNAmCystatinC...28" "DNAmGDF15...29"    
# [31] "DNAmLeptin...30"    "DNAmPACKYRS...31"   "DNAmPAI1...32"     
# [34] "DNAmTIMP1...33"     "DNAmlogA1C"         "DNAmlogCRP"        
# [37] "Age...36"           "Female...37"        "DNAmGrimAge...38"  
# [40] "AgeAccelGrim...39"  "AlcoholMcCartney"   "BMI_McCartney"     
# [43] "DNAmTL"             "HRSInCHPhenoAge"    "Knight"            
# [46] "LeeControl"         "X"                  "Y"                 
# [49] "predicted_sex"      "LeeRefinedRobust"   "LeeRobust"         
# [52] "SmokingMcCartney"   "DownSyndrome"

human_dnam_age_tb = human_dnam_age_tb |>
  select(
    file, gsm_id = ...1, Horvathv1, Hannum, Lin, PhenoAge, YingCausAge, YingDamAge,
    YingAdaptAge, Horvathv2, PEDBE, Zhang_10, DunedinPoAm38, DunedinPACE, DNAmADM = DNAmADM...14,
    DNAmB2M = DNAmB2M...15, DNAmCystatinC = DNAmCystatinC...16, DNAmGDF15 = DNAmGDF15...17,
    DNAmLeptin = DNAmLeptin...18, DNAmPACKYRS = DNAmPACKYRS...19, DNAmPAI1 = DNAmPAI1...20,
    DNAmTIMP1 = DNAmTIMP1...21, GrimAgeV1 = DNAmGrimAge...24,
    DNAmlogA1C = DNAmlogA1C, DNAmlogCRP = DNAmlogCRP, 
    GrimAgeV2 = DNAmGrimAge...38, Alcohol = AlcoholMcCartney, BMI = BMI_McCartney,
    DNAmTL = DNAmTL, Knight = Knight, LeeRobust = LeeRobust, Smoking = SmokingMcCartney
  ) |>
  tidyr::extract(file, into = "gse_id", regex = "(GSE[0-9]+)")

human_dnam_missing_tb = human_dnam_missing_tb |>
  select(
    file, gsm_id = ...1, Horvathv1, Hannum, Lin, PhenoAge, YingCausAge, YingDamAge,
    YingAdaptAge, Horvathv2, PEDBE, Zhang_10, DunedinPoAm38, DunedinPACE,
    GrimAgeV1, GrimAgeV2, Alcohol = AlcoholMcCartney, BMI = BMI_McCartney,
    DNAmTL = DNAmTL, Knight = Knight, LeeRobust = LeeRobust, Smoking = SmokingMcCartney
  ) |>
  tidyr::extract(file, into = "gse_id", regex = "(GSE[0-9]+)")

# Mask clock values to NA if missing percentage is greater than 20%
qced_human_dnam_age_tb <- human_dnam_age_tb %>%
  inner_join(human_dnam_missing_tb, by = c("gse_id", "gsm_id"), suffix = c("", ".y")) %>%
  mutate(across(
    c(Horvathv1, Hannum, Lin, PhenoAge, YingCausAge, YingDamAge,
      YingAdaptAge, Horvathv2, PEDBE, Zhang_10, DunedinPoAm38, DunedinPACE,
      GrimAgeV1, GrimAgeV2, Alcohol, BMI, DNAmTL, Knight, LeeRobust, Smoking),
    ~ ifelse(get(paste0(cur_column(), ".y")) > 0.2, NA, .x)
  )) %>%
  select(-ends_with(".y")) |>
  # Mask DNAm* values to NA if either GrimAgeV1 or GrimAgeV2 is NA
  mutate(across(
    c(DNAmADM, DNAmB2M, DNAmCystatinC, DNAmGDF15, DNAmLeptin, DNAmPACKYRS,
      DNAmPAI1, DNAmTIMP1, DNAmlogA1C, DNAmlogCRP),
    ~ ifelse(is.na(GrimAgeV1) | is.na(GrimAgeV2), NA, .x)
  ))

clocks = c(
  "Horvathv1", "Hannum", "Lin", "PhenoAge", "YingCausAge", "YingDamAge",
  "YingAdaptAge", "Horvathv2", "PEDBE", "Zhang_10", "DunedinPoAm38", "DunedinPACE",
  "GrimAge", "GrimAge2", "DNAmTL", "Knight", "LeeRobust"
)

biomarkers = c(
  "DNAmADM", "DNAmB2M", "DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPACKYRS",
  "DNAmlogA1C", "DNAmlogCRP", "Alcohol", "BMI", "Smoking"
)

colnames(human_dnam_meta_tb) = c(
  "file",
  "gsm_id",
  "PATIENT_ID",
  "race",
  "sex",
  "age",
  "genetic_info",
  "disease",
  "tissue",
  "cell_line",
  "vivo_vitro",
  "case_control",
  "group_name",
  "treatment",
  "perturbation_category"
)
human_dnam_meta_tb = human_dnam_meta_tb |>
  tidyr::extract(file, into = "gse_id", regex = "(GSE[0-9]+)")

write_csv(qced_human_dnam_age_tb, "result/qced_human_dnam_age.csv")
write_csv(human_dnam_meta_tb, "result/gpt4_human_dnam_metadata.csv")

# ## Make annotation for the umap plot ------
# human_dnam_meta_all = read_csv("result/gpt4_human_dnam_metadata.csv")
# make_annotation = function(value_vec) {
#   annot_str = "</br>"
#   value_vec = value_vec[!duplicated(value_vec)]
#   for (i in 1:length(value_vec)) {
#     if (!is.na(value_vec[[i]])) {
#       annot_str = paste0(
#         annot_str,
#         names(value_vec)[i] |> str_remove("harmonized."),
#         ": ",
#         value_vec[[i]],
#         "</br>"
#       )
#     }
#   }
#   annot_str
# }

# all_annot = human_dnam_meta_all |>
#   # select(1:6, harmonized.disease, harmonized.treatment) |>
#   as.matrix() |>
#   apply(1, make_annotation)

# lite_annot = human_dnam_meta_all |>
#   select(1:2) |>
#   as.matrix() |>
#   apply(1, make_annotation)

# annot_df = human_dnam_meta_all |>
#   select(1:2) |>
#   mutate(annotation = all_annot) |>
#   mutate(annotation_lite = lite_annot)

# write_csv(annot_df, "data/human_dnam/annotation.csv")

# annot_df = read_csv("data/human_dnam/annotation.csv")
