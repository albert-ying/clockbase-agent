import glob
import os
import re
import traceback
from importlib import import_module
from multiprocessing import Pool

import numpy as np
import pandas as pd
from biolearn.model_gallery import ModelGallery
from tqdm import tqdm

tmp = import_module('src.4_compute_age')

# Set the directories
data_dir = "./GEO_hs_dnam/series_matrix/"
out_folder = "./GEO_hs_dnam/missing_percentage"

# Create the output directory if it doesn't exist
if not os.path.exists(out_folder):
   os.makedirs(out_folder)

# List all files in the data directory
files = glob.glob(f"{data_dir}/*")

# Initialize the model gallery
clock_models = ModelGallery()

def calculate_missing_percentage(args):
   file, clock_names = args
   try:
       gse_id = re.search(r"GSE\d+", file).group(0)

       out_file = f"{out_folder}/{gse_id}.csv"

       if os.path.exists(out_file):
           print(f"{gse_id} already exists")
       else:
           # Read the data
           test_data = tmp.read_gse_matrix(file)

           missing_percentages = pd.DataFrame()

           for clock_name in clock_names:
               print(f"Processing {clock_name}")
               clock_model = clock_models.get(clock_name)
               clock_sites = clock_model.methylation_sites()

               missing_count = test_data.reindex(clock_sites).isnull().sum(axis=0)
               total_sites = len(clock_sites)
               missing_percentage = missing_count / total_sites

               missing_percentages[clock_name] = missing_percentage

           # Save the results
           missing_percentages.to_csv(out_file)
           print(f"Processed {gse_id}")

   except Exception as e:
       print(f"Error processing {file}: {e}")
       traceback.print_exc()

# Specify the path to a single dataset file
if __name__ == "__main__":

  model_definitions = list(clock_models.model_definitions.keys())
  files_clock_definitions = [(file, model_definitions) for file in files]
  with Pool(12) as p:
      list(
          tqdm(
              p.imap_unordered(calculate_missing_percentage, files_clock_definitions),
              total=len(files_clock_definitions)
          )
      )
