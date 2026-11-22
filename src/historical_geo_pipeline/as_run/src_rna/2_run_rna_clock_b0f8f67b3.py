import os
from multiprocessing import Pool

import joblib
import numpy as np
import pandas as pd

# Load the clocks
clock_folders = [
    "/Lasagna/alex2308/tAge_clocks/BR differential models 4.6",
    "/Lasagna/alex2308/tAge_clocks/EN differential models 4.6",
    "/Lasagna/alex2308/tAge_clocks/EN absolute models 4.6",
]
clock_files = []
for folder in clock_folders:
    for file in os.listdir(folder):
        if file.endswith(".pkl"):
            clock_files.append(os.path.join(folder, file))


# Function to load a clock and get its feature names
def load_clock(file):
    clock = joblib.load(file)
    clock.feature_names = [str(x) for x in clock.feature_names]
    return clock


clocks = [load_clock(file) for file in clock_files]
clock_names = [file.strip("/Lasagna/alex2308/tAge_clocks/") for file in clock_files]


def generate_clock_name(file_path):
    # Split the file path by '/'
    path_parts = file_path.split("/")
    # Get the last part of the path (file name)
    file_name = path_parts[-1]
    # Get the model type ("differential" or "absolute")
    model_type = path_parts[-2].split(" ")[0].split("_")[-1]
    # Remove the '.pkl' extension
    file_name = file_name[:-4]
    # Split the file name by '_'
    name_parts = file_name.split("_")
    # Extract the relevant parts
    clock_type = name_parts[1]
    species = name_parts[2]
    tissue = name_parts[3]
    genotype = name_parts[4]
    data_type = name_parts[5]
    # Construct the clock name with model type
    clock_name = f"{clock_type}_{species}_{tissue}_{genotype}_{data_type}_{model_type}"
    return clock_name


clock_names = [generate_clock_name(file_path) for file_path in clock_files]


def predict_tAge(data, clock, clock_name):
    features = clock.feature_names
    data_for_calculation = data[features]

    model_type = clock_name.split("_")[-1]

    if model_type == "BR":
        tAge, tAge_std = clock.predict(data_for_calculation, return_std=True)
        return {"tAge": tAge, "tAge_std": tAge_std}
    else:
        tAge = clock.predict(data_for_calculation)
        if "scaledlog" in clock_name or "yugenelog" in clock_name:
            tAge = 10 ** (-(10 ** (-tAge)))
        return {"tAge": tAge}


# Function to apply clocks to a dataset
def apply_clocks(data, clock, clock_name, data_annotation, species):
    if species == "human" and ("Rodents" in clock_name or "Mouse" in clock_name):
        return data_annotation

    # Determine the appropriate data based on clock type
    # add breakpoint here
    temp_tAge = predict_tAge(data, clock, clock_name)
    if "tAge_std" in temp_tAge:
        data_annotation[f"{clock_name}_tAge_mean"] = temp_tAge["tAge"]
        data_annotation[f"{clock_name}_tAge_std"] = temp_tAge["tAge_std"]
    else:
        data_annotation[f"{clock_name}_tAge"] = temp_tAge["tAge"]

    return data_annotation


# Function to process a single file
def process_file(file, data_folder, output_folder, species):
    print("Processing file:", file)
    exprs_data = pd.read_csv(os.path.join(data_folder, file), index_col=0)
    data_annotation = pd.DataFrame({"sample_id": exprs_data.index})

    for clock, clock_name in zip(clocks, clock_names):
        if "scaleddiff" in clock_name:
            data_key = "scaled_diff"
        elif "yugenediff" in clock_name:
            data_key = "yugene_diff"
        elif "scaledlog" in clock_name:
            data_key = "scaled_data"
        elif "yugenelog" in clock_name:
            data_key = "yugene_data"
        if data_key in file:
            print(f"Applying clock: {clock_name} to {file}")
            data_annotation = apply_clocks(
                exprs_data, clock, clock_name, data_annotation, species
            )

    output_file = file.replace("data", "clock")
    print(f"Saving output to {output_folder}/{output_file}")
    data_annotation.to_csv(os.path.join(output_folder, output_file), index=False)


## Run clock
# Process mouse data
mouse_processed_output_folder = "GEO_mm_rna/processed_data/"
mouse_processed_files = np.sort(os.listdir(mouse_processed_output_folder))
mouse_clock_output_folder = "GEO_mm_rna/clock_result/"


# test in one dataset (GSE122080)
import fnmatch

specific_pattern = "GSE122080_*"
specific_files = fnmatch.filter(mouse_processed_files, specific_pattern)

for specific_file in specific_files:
    process_file(
        specific_file,
        mouse_processed_output_folder,
        mouse_clock_output_folder,
        "mouse",
    )


process_file(
    mouse_processed_files[1],
    mouse_processed_output_folder,
    mouse_clock_output_folder,
    "mouse",
)

with Pool() as pool:
    pool.starmap(
        process_file,
        [
            (file, mouse_processed_output_folder, mouse_clock_output_folder, "mouse")
            for file in mouse_processed_files
        ],
    )

# Process human data
human_processed_output_folder = "GEO_hs_rna/processed_data/"
human_processed_files = [
    file
    for file in os.listdir(human_processed_output_folder)
    if file.endswith(("_scaled_data.csv", "_yugene_data.csv"))
]
human_clock_output_folder = "GEO_hs_rna/clock_result/"

with Pool() as pool:
    pool.starmap(
        process_file,
        [
            (file, human_processed_output_folder, human_clock_output_folder, "human")
            for file in human_processed_files
        ],
    )
