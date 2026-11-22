import glob
import gzip
import os
import re
import traceback
from io import StringIO
from multiprocessing import Pool

import biolearn
import numpy as np
import pandas as pd
import requests
from biolearn.data_library import DataLibrary, GeoData
from biolearn.imputation import impute_from_standard
from biolearn.model_gallery import ModelGallery
from tqdm import tqdm

## Configure the variables
reference_dataset = DataLibrary().get("GSE40279").load()
reference_averages = reference_dataset.dnam.mean(axis=1)
clock_models = ModelGallery()
clock_list = list(clock_models.model_definitions.keys())

# Set the directories
data_dir = "./GEO_hs_dnam/series_matrix/"
out_folder = "./GEO_hs_dnam/biolearn_output"

# Create the output directory if it doesn't exist
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

# List all files in the data directory
files = glob.glob(f"{data_dir}/*")


def read_gse_matrix(filename, data=True):
    def download_file(url, local_filename):
        with requests.get(url, stream=True) as r:
            with open(local_filename, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        return local_filename

    def read_until_marker(marker, file):
        lines = []
        for line in file:
            if line.strip() == marker:
                break
            lines.append(line.strip())
        return lines

    # Check if filename is a URL
    if filename.startswith("http://") or filename.startswith("https://"):
        filename = download_file(filename, os.path.basename(filename))
    # Check if file is gzipped
    is_gzipped = filename.endswith(".gz")
    # Open the file accordingly
    open_func = gzip.open if is_gzipped else open
    mode = "rt" if is_gzipped else "r"
    with open_func(filename, mode) as file:
        lines = read_until_marker("!series_matrix_table_begin", file)
        # matrix_start = lines number + 2
        matrix_start = len(lines) + 1
    data = pd.read_csv(filename, sep="\t", header=0, skiprows=matrix_start)
    # Process the data
    data.index = data.iloc[:, 0]
    data.index.name = "id"
    data = data.drop(data.columns[0], axis=1)
    # remove the row with !series_matrix_table_end
    data = data.iloc[:-1, :]
    return data


def run_human_clock(args):
    file, clock_list = args
    try:
        gse_id = re.search(r"GSE\d+", file).group(0)

        out_file = f"{out_folder}/{gse_id}.csv"

        if os.path.exists(out_file):
            print(f"{gse_id} already exists")
        else:
            # Read the data (adapt this part according to your data format)
            test_data = read_gse_matrix(file)
            all_id = reference_averages.index.union(test_data.index)
            test_data = test_data.reindex(all_id)

            imputed_data = pd.DataFrame()
            imputed_data.dnam = impute_from_standard(test_data, reference_averages)
            imputed_data.metadata = pd.DataFrame()

            # imputed_data.metadata = None
            horvath_res = clock_models.get("Horvathv1").predict(imputed_data)
            sex_res = clock_models.get(
                "SexEstimation", imputation_method="averaging"
            ).predict(imputed_data)

            imputed_data.metadata = pd.concat(
                [horvath_res, sex_res[["predicted_sex"]]], axis=1
            )
            imputed_data.metadata.columns = ["age", "sex"]

            all_results = pd.DataFrame()

            for clock_name in clock_list:
                print(f"Processing {clock_name}")
                if clock_name == "SexEstimation":
                    clock_model = clock_models.get(
                        clock_name, imputation_method="averaging"
                    )
                else:
                    clock_model = clock_models.get(clock_name)

                result = clock_model.predict(imputed_data)
                # if clock_name == "GrimAgeV1" or clock_name == "GrimAgeV2":
                if len(result.columns) > 1:
                    pass
                else:
                    result.columns = [clock_name]
                print(f"Processed {clock_name}")
                all_results = pd.concat([all_results, result], axis=1)

            # Save the results
            all_results.to_csv(out_file)
            print(f"Processed {gse_id}")

    except Exception as e:
        print(f"Error processing {file}: {e}")
        traceback.print_exc()


def main():
    # Using multiprocessing for parallel processing
    # Load the reference dataset and compute averages for imputation

    files_clock_list = [(file, clock_list) for file in files]
    with Pool(12) as p:
        # Use tqdm with p.imap or p.imap_unordered for progress bar functionality
        # imap_unordered offers better performance but does not guarantee the order of processing
        list(
            tqdm(
                p.imap_unordered(run_human_clock, files_clock_list),
                total=len(files_clock_list),
            )
        )


if __name__ == "__main__":
    main()
