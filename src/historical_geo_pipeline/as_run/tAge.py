## Script for aging clock application
print("Loading libraries...")
import pandas as pd
import numpy as np
import joblib
from scipy.special import logit, expit
import os
import sys
import re
import multiprocessing
np.random.seed(42)

print("Loading models...")
## get parameters from command line
data_dir = sys.argv[1]
out_dir = sys.argv[2]
species = sys.argv[3]

# data_dir = "./GEO_hs_rna/processed_data"
# species = "human"
# data = "./GEO_hs_rna/processed_data/GSE100866_scaled.csv"
# mouse aging clocks
# gse_id = extract_GSE_ID(data)
print("Loading models...")
clock1 = joblib.load(open("rna_clocks/tAge_models/Chronoage_liver_All_Liver_WT_Age_Relative.Species_scaled_loglog.pkl", 'rb'))
clock2 = joblib.load(open("rna_clocks/tAge_models/Chronoage_brain_Rodents_Brain_WT_Age_Relative.Species_scaled_loglog.pkl", 'rb'))
clock3 = joblib.load(open("rna_clocks/tAge_models/Chronoage_Rodents_Muscle_WT_Age_Relative.Species_scaled_loglog.pkl", 'rb'))
clock4 = joblib.load(open("rna_clocks/tAge_models/Chronoage_Mouse_All_WT_Age_Relative.Species_scaled_loglog.pkl", 'rb'))
clock5 = joblib.load(open("rna_clocks/tAge_models/Chronoage_Rodents_All_WT_Age_Relative.Species_scaled_loglog.pkl", 'rb'))
clock6 = joblib.load(open("rna_clocks/tAge_models/Chronoagediff_Rodents_Liver_WT_Difference_ChronoAge_scaled.pkl", 'rb'))
clock7 = joblib.load(open("rna_clocks/tAge_models/Chronoagediff_brain_Rodents_Brain_WT_Difference_ChronoAge_scaled.pkl", 'rb'))
clock8 = joblib.load(open("rna_clocks/tAge_models/Chronoagediff_Rodents_Muscle_WT_Difference_ChronoAge_scaled.pkl", 'rb'))
clock9 = joblib.load(open("rna_clocks/tAge_models/Chronoagediff_Mouse_All_WT_Difference_ChronoAge_scaled.pkl", 'rb'))
clock10 = joblib.load(open("rna_clocks/tAge_models/Chronoagediff_Rodents_All_WT_Difference_ChronoAge_scaled.pkl", 'rb'))
clock11 = joblib.load(open("rna_clocks/tAge_models/INTS_relage_best_All_Liver_All_Age_Relative.Max.Best_scaled_loglog.pkl", 'rb'))
clock12 = joblib.load(open("rna_clocks/tAge_models/INTS_relage_best_Mouse_All_All_Age_Relative.Max.Best_scaled_loglog.pkl", 'rb'))
clock13 = joblib.load(open("rna_clocks/tAge_models/INTS_relage_best_Rodents_All_All_Age_Relative.Max.Best_scaled_loglog.pkl", 'rb'))
clock14 = joblib.load(open("rna_clocks/tAge_models/INTS_relage_diff_best_All_Liver_All_Difference_RelAge.Max.Best_scaled.pkl", 'rb'))
clock15 = joblib.load(open("rna_clocks/tAge_models/INTS_relage_diff_best_Mouse_All_All_Difference_RelAge.Max.Best_scaled.pkl", 'rb'))
clock16 = joblib.load(open("rna_clocks/tAge_models/INTS_relage_diff_best_Rodents_All_All_Difference_RelAge.Max.Best_scaled.pkl", 'rb'))
clock17 = joblib.load(open("rna_clocks/tAge_models/INTS_Maxlifespan_best_All_Liver_All_Lifespan_Max.Best_scaled.pkl", 'rb'))
clock18 = joblib.load(open("rna_clocks/tAge_models/INTS_lifespan_diff_best_All_Liver_All_Difference_Lifespan.Max.Best_scaled.pkl", 'rb'))

# human aging clocks
clock19 = joblib.load(open("rna_clocks/tAge_models/Chronoage_Human_All_WT_Age_Relative.Species_scaled_loglog.pkl", 'rb'))
clock20 = joblib.load(open("rna_clocks/tAge_models/Chronoagediff_Human_All_WT_Difference_ChronoAge_scaled.pkl", 'rb'))

print("Loading genes...")
gene_req = pd.read_csv('./rna_clocks/tAge_models/Gene_list_clock.txt', header = None, dtype=str)
gene_req = gene_req[0].tolist()

def extract_GSE_ID(file):
  # extract GSE ID from file name use regex
  # format is gse888888
  gse_id = re.search('GSE[0-9]+', file)
  return gse_id.group()

def predict_tAge(data,model):
  features = model.feature_names
  features = [str(x) for x in features]
  # select the feature columns in the data, if not present, fill with NA
  data_for_calculation = data.reindex(columns=features, fill_value=np.nan)
  tAge = model.predict(data_for_calculation)
  return tAge

#==============================
def tAge_from_gse(gse_id, data_dir, out_dir, species):
  print(gse_id)
  exprs_scaled = pd.read_csv(f"{data_dir}/{gse_id}_scaled.csv",index_col=0)
  # make a dataframe with the first column be the column names of the expression matrix, the name of the column is "sample_id"
  age_df = pd.DataFrame(exprs_scaled.columns, columns=['sample_id'])
  # Transpose
  exprs_scaled = exprs_scaled.T
  # Change column names to string
  exprs_scaled.columns = [str(x) for x in exprs_scaled.columns]
  if species == 'mouse':
    aging_models = mouse_clock_models
    aging_models_files = mouse_clock_models_files
    aging_models_names = mouse_clock_models_names
  elif species == 'human':
    aging_models = human_clock_models
    aging_models_files = human_clock_models_files
    aging_models_names = human_clock_models_names

  for i in range(len(aging_models)):
    temp_clock = eval(aging_models[i])
    temp_tAge = predict_tAge(exprs_scaled,temp_clock)
    if (aging_models_files[i].split('_')[-1] == 'loglog'):
        temp_tAge = 10**(10**temp_tAge)
    age_df.loc[:,aging_models_names[i]] = temp_tAge 
  # Save data annotations
  age_df.to_csv(f'{out_dir}/{gse_id}.csv')

#==============================

mouse_clock_models = [
  "clock1", "clock2", "clock3", "clock4", "clock5", "clock6", "clock7", "clock8",
  "clock9", "clock10", "clock11", "clock12", "clock13", "clock14", "clock15",
  "clock16", "clock17", "clock18"
]
human_clock_models = ["clock19", "clock20"]
mouse_clock_models_files = [
  "Chronoage_liver_All_Liver_WT_Age_Relative.Species_scaled_loglog",
  "Chronoage_brain_Rodents_Brain_WT_Age_Relative.Species_scaled_loglog",
  "Chronoage_Rodents_Muscle_WT_Age_Relative.Species_scaled_loglog",
  "Chronoage_Mouse_All_WT_Age_Relative.Species_scaled_loglog",
  "Chronoage_Rodents_All_WT_Age_Relative.Species_scaled_loglog",
  "Chronoagediff_Rodents_Liver_WT_Difference_ChronoAge_scaled",
  "Chronoagediff_brain_Rodents_Brain_WT_Difference_ChronoAge_scaled",
  "Chronoagediff_Rodents_Muscle_WT_Difference_ChronoAge_scaled",
  "Chronoagediff_Mouse_All_WT_Difference_ChronoAge_scaled",
  "Chronoagediff_Rodents_All_WT_Difference_ChronoAge_scaled",
  "INTS_relage_best_All_Liver_All_Age_Relative.Max.Best_scaled_loglog",
  "INTS_relage_best_Mouse_All_All_Age_Relative.Max.Best_scaled_loglog",
  "INTS_relage_best_Rodents_All_All_Age_Relative.Max.Best_scaled_loglog",
  "INTS_relage_diff_best_All_Liver_All_Difference_RelAge.Max.Best_scaled",
  "INTS_relage_diff_best_Mouse_All_All_Difference_RelAge.Max.Best_scaled",
  "INTS_relage_diff_best_Rodents_All_All_Difference_RelAge.Max.Best_scaled",
  "INTS_Maxlifespan_best_All_Liver_All_Lifespan_Max.Best_scaled",
  "INTS_lifespan_diff_best_All_Liver_All_Difference_Lifespan.Max.Best_scaled"
]
human_clock_models_files = [
  "Chronoage_Human_All_WT_Age_Relative.Species_scaled_loglog",
  "Chronoagediff_Human_All_WT_Difference_ChronoAge_scaled"
]
mouse_clock_models_names = [
    "Chrono_tAge_liver",
    "Chrono_tAge_brain",
    "Chrono_tAge_muscle",
    "Chrono_tAge_mouse",
    "Chrono_tAge_rodents",
    "Chrono_tAge_diff_liver",
    "Chrono_tAge_diff_brain",
    "Chrono_tAge_diff_muscle",
    "Chrono_tAge_diff_mouse",
    "Chrono_tAge_diff_rodents",
    "Biological_tAge_liver",
    "Biological_tAge_mouse",
    "Biological_tAge_rodents",
    "Biological_tAge_diff_liver",
    "Biological_tAge_diff_mouse",
    "Biological_tAge_diff_rodents",
    "Lifespan_clock_liver",
    "Lifespan_clock_diff_liver"
]
human_clock_models_names = [
    "Chrono_tAge_human",
    "Chrono_tAge_diff human"
]

geo_files = os.listdir(data_dir)
print("Number of files: ", len(geo_files))
print("Files: ", geo_files)
# Subset geo_files to only include META
geo_files = [x for x in geo_files if 'scaled.csv' in x]
gse_ids = [extract_GSE_ID(x) for x in geo_files]

print("Starting to calculate tAge")

# for gse_id in gse_ids:
#   try:
#     print(gse_id)
#     tAge_from_gse(gse_id, data_dir, out_dir, species)
#   except:
#     pass


def process_gse_id(gse_id):
    if f"{gse_id}.csv" in os.listdir(out_dir):
        return
    else:
        try:
            print(gse_id)
            tAge_from_gse(gse_id, data_dir, out_dir, species)
        except Exception as e:
            print(f"Error processing {gse_id}: {e}")

with multiprocessing.Pool(processes=30) as pool:
    pool.map(process_gse_id, gse_ids)
