import argparse
import csv
import os
import re
from datetime import datetime, timezone
from multiprocessing import Pool
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from tqdm import tqdm

##
# Load the clocks
clock_folders = [
    "/Lasagna/alex2308/tAge_clocks/BR differential models 4.6",
    "/Lasagna/alex2308/tAge_clocks/EN differential models 4.6",
    "/Lasagna/alex2308/tAge_clocks/EN absolute models 4.6",
]
clock_files = []
clocks = []
clock_names = []


# Function to load a clock and get its feature names
def load_clock(file):
    clock = joblib.load(file)
    clock.feature_names = [str(x) for x in clock.feature_names]
    return clock


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


def initialize_clocks():
    global clock_files, clocks, clock_names
    clock_files = []
    for folder in clock_folders:
        for file in os.listdir(folder):
            if file.endswith(".pkl"):
                clock_files.append(os.path.join(folder, file))

    clocks = [load_clock(file) for file in clock_files]
    clock_names = [generate_clock_name(file_path) for file_path in clock_files]


def get_gse_id(path):
    match = re.search(r"GSE\d+", str(path))
    return match.group(0) if match else None


##
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
def process_file_folder(file_folder, data_folder, output_folder, species):
    gse_id = file_folder
    print("Processing", gse_id)
    file_list = os.listdir(os.path.join(data_folder, gse_id))
    output_file = gse_id + ".csv"

    if os.path.exists(os.path.join(output_folder, output_file)):
        print(f"Output file {output_file} already exists, skipping...")
        return None
    output_dfs = []
    for file in file_list:
        # check if the output already exists
        exprs_data = pd.read_csv(os.path.join(data_folder, gse_id, file), index_col=0)
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
                try:
                    data_annotation = apply_clocks(
                        exprs_data, clock, clock_name, data_annotation, species
                    )
                except Exception as e:
                    print(f"Error applying clock: {clock_name} to {file}")
                    continue

        output_dfs.append(data_annotation)

    output_df = pd.concat(output_dfs, axis=1)
    output_df = output_df.loc[:, ~output_df.columns.duplicated()]
    output_df.to_csv(os.path.join(output_folder, output_file), index=False)
    return None


# def process_and_merge_files(folder):
#     clock_files = os.listdir(folder)
#     gse_dict = {}

#     for file in clock_files:
#         gse_id = file.split("_")[0]
#         if gse_id not in gse_dict:
#             gse_dict[gse_id] = []
#         gse_dict[gse_id].append(file)

#     for gse_id, files in gse_dict.items():
#         if len(files) > 1:
#             combined_df = pd.concat(
#                 [pd.read_csv(os.path.join(folder, f)) for f in files], axis=1
#             )
#             combined_df.to_csv(os.path.join(folder, f"{gse_id}.csv"), index=False)
#             # Delete original files
#             for f in files:
#                 os.remove(os.path.join(folder, f))


def safe_process_file_folder(file_folder, data_folder, output_folder, species):
    try:
        process_file_folder(file_folder, data_folder, output_folder, species)
    except Exception as e:
        print(f"Error processing file: {file_folder}")
        print(e)


def summarize_scaled_input_file(path, species):
    gse_id = get_gse_id(path)
    row = {
        "species": species,
        "GSE_ID": gse_id,
        "scaled_data_file": str(path),
    }
    try:
        with Path(path).open(newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            first_row = next(reader)

        if len(header) < 2 or len(first_row) < 2:
            raise ValueError("scaled data file has no gene columns")

        gene_values = first_row[1:]
        missing_values = {
            "",
            "NA",
            "NaN",
            "nan",
            "N/A",
            "NULL",
            "null",
        }
        missing_count = sum(value in missing_values for value in gene_values)
        missing_ratio = missing_count / len(gene_values)
        row.update(
            {
                "n_clock_gene_cols": int(len(gene_values)),
                "n_first_sample_missing_genes": int(missing_count),
                "scaled_data_first_sample_missing_ratio": float(missing_ratio),
                "scaled_data_first_sample_present_ratio": float(1 - missing_ratio),
                "input_summary_error": "",
            }
        )
    except Exception as e:
        row["input_summary_error"] = str(e)
    return row


def summarize_scaled_input_file_from_args(args):
    path, species = args
    return summarize_scaled_input_file(path, species)


def summarize_scaled_inputs(processed_dir, species, ncores):
    paths = sorted(Path(processed_dir).glob("GSE*/GSE*_scaled_data.csv"))
    if not paths:
        raise FileNotFoundError(
            f"Cannot proceed - required data file not found: {processed_dir}/GSE*/GSE*_scaled_data.csv"
        )
    tasks = [(path, species) for path in paths]
    if ncores > 1:
        with Pool(ncores) as pool:
            rows = list(
                tqdm(
                    pool.imap_unordered(summarize_scaled_input_file_from_args, tasks),
                    total=len(tasks),
                )
            )
    else:
        rows = [summarize_scaled_input_file(path, species) for path in tqdm(paths)]
    return pd.DataFrame(rows)


def summarize_age_table(path, species, near_constant_range, chunksize):
    path = Path(path)
    if not path.exists():
        return pd.DataFrame({"species": [], "GSE_ID": []})

    min_df = None
    max_df = None
    sample_counts = pd.Series(dtype="int64")

    for chunk in tqdm(pd.read_csv(path, chunksize=chunksize)):
        if "gse_id" not in chunk.columns:
            raise ValueError(f"{path} does not contain gse_id")

        prediction_cols = [
            col
            for col in chunk.columns
            if col not in {"gse_id", "gsm_id"} and not col.endswith("_tAge_std")
        ]
        numeric = chunk[prediction_cols].apply(pd.to_numeric, errors="coerce")
        group_key = chunk["gse_id"].astype(str)
        group_min = numeric.groupby(group_key).min()
        group_max = numeric.groupby(group_key).max()
        group_counts = group_key.value_counts()

        if min_df is None:
            min_df = group_min
            max_df = group_max
        else:
            min_df = pd.concat([min_df, group_min]).groupby(level=0).min()
            max_df = pd.concat([max_df, group_max]).groupby(level=0).max()

        sample_counts = sample_counts.add(group_counts, fill_value=0)

    if min_df is None:
        return pd.DataFrame({"species": [], "GSE_ID": []})

    ranges = max_df - min_df
    available_ranges = ranges.notna().sum(axis=1)
    near_constant = ranges <= near_constant_range
    near_constant_fraction = (
        near_constant.sum(axis=1) / available_ranges.replace(0, np.nan)
    )
    summary = pd.DataFrame(
        {
            "species": species,
            "GSE_ID": ranges.index,
            "raw_clock_n_samples": sample_counts.reindex(ranges.index).fillna(0).astype(int),
            "raw_clock_n_prediction_cols": available_ranges.astype(int),
            "raw_clock_n_near_constant_prediction_cols": near_constant.sum(axis=1).astype(int),
            "raw_clock_frac_near_constant_prediction_cols": near_constant_fraction,
            "raw_clock_median_abs_range": ranges.median(axis=1, skipna=True),
            "raw_clock_max_abs_range": ranges.max(axis=1, skipna=True),
        }
    )
    return summary


def read_geo_metadata(path, species):
    path = Path(path)
    if not path.exists():
        return pd.DataFrame({"species": [], "GSE_ID": []})

    metadata = pd.read_csv(path)
    if "accession" not in metadata.columns:
        return pd.DataFrame({"species": [], "GSE_ID": []})

    keep_cols = [
        col
        for col in ["title", "organism", "type", "platform", "accession", "samples"]
        if col in metadata.columns
    ]
    metadata = metadata[keep_cols].rename(
        columns={"accession": "GSE_ID", "samples": "geo_samples_reported"}
    )
    metadata["species"] = species
    return metadata.drop_duplicates(["species", "GSE_ID"])


def newest_mtime(paths):
    existing = [Path(path) for path in paths if Path(path).exists()]
    if not existing:
        return ""
    newest = max(path.stat().st_mtime for path in existing)
    return datetime.fromtimestamp(newest, timezone.utc).isoformat()


def write_source_log(output_path, args, n_all_rows, n_flagged_rows):
    output_path = Path(output_path)
    source_log = output_path.with_suffix(".sources.txt")
    hs_scaled = sorted(Path(args.human_processed_dir).glob("GSE*/GSE*_scaled_data.csv"))
    mm_scaled = sorted(Path(args.mouse_processed_dir).glob("GSE*/GSE*_scaled_data.csv"))
    lines = [
        f"generated_at_utc: {datetime.now(timezone.utc).isoformat()}",
        f"working_directory: {Path.cwd()}",
        f"output_file: {output_path}",
        f"all_rows: {n_all_rows}",
        f"flagged_rows: {n_flagged_rows}",
        "",
        "parameters:",
        f"  near_missing_threshold: {args.near_missing_threshold}",
        f"  near_constant_range: {args.near_constant_range}",
        f"  near_constant_fraction: {args.near_constant_fraction}",
        f"  age_chunksize: {args.age_chunksize}",
        f"  summary_ncores: {args.summary_ncores}",
        "",
        "input_sources:",
        f"  human_processed_dir: {args.human_processed_dir}",
        f"  human_scaled_data_file_count: {len(hs_scaled)}",
        f"  human_processed_dir_newest_mtime_utc: {newest_mtime(hs_scaled)}",
        f"  mouse_processed_dir: {args.mouse_processed_dir}",
        f"  mouse_scaled_data_file_count: {len(mm_scaled)}",
        f"  mouse_processed_dir_newest_mtime_utc: {newest_mtime(mm_scaled)}",
        f"  human_age_table: {args.human_age_table}",
        f"  human_age_table_mtime_utc: {newest_mtime([args.human_age_table])}",
        f"  mouse_age_table: {args.mouse_age_table}",
        f"  mouse_age_table_mtime_utc: {newest_mtime([args.mouse_age_table])}",
        f"  human_geo_metadata: {args.human_geo_metadata}",
        f"  human_geo_metadata_mtime_utc: {newest_mtime([args.human_geo_metadata])}",
        f"  mouse_geo_metadata: {args.mouse_geo_metadata}",
        f"  mouse_geo_metadata_mtime_utc: {newest_mtime([args.mouse_geo_metadata])}",
        "",
        "flag_logic:",
        "  scaled_data_first_sample_missing_ratio is computed from the first sample row of *_scaled_data.csv.",
        "  This is valid for this preprocessing path because absent clock genes are added as all-sample NA columns.",
        "  flag_near_all_input_missing = scaled_data_first_sample_missing_ratio >= near_missing_threshold",
        "  flag_constant_raw_clock_outputs = raw_clock_frac_near_constant_prediction_cols >= near_constant_fraction",
        "  artifact_flag = flag_near_all_input_missing",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    source_log.write_text("\n".join(lines) + "\n")


def summarize_imputation_artifacts(args):
    human_inputs = summarize_scaled_inputs(
        args.human_processed_dir, "human", args.summary_ncores
    )
    mouse_inputs = summarize_scaled_inputs(
        args.mouse_processed_dir, "mouse", args.summary_ncores
    )
    input_summary = pd.concat([human_inputs, mouse_inputs], ignore_index=True)

    human_age = summarize_age_table(
        args.human_age_table, "human", args.near_constant_range, args.age_chunksize
    )
    mouse_age = summarize_age_table(
        args.mouse_age_table, "mouse", args.near_constant_range, args.age_chunksize
    )
    age_summary = pd.concat([human_age, mouse_age], ignore_index=True)

    metadata = pd.concat(
        [
            read_geo_metadata(args.human_geo_metadata, "human"),
            read_geo_metadata(args.mouse_geo_metadata, "mouse"),
        ],
        ignore_index=True,
    )

    summary = input_summary.merge(age_summary, on=["species", "GSE_ID"], how="left")
    summary = summary.merge(metadata, on=["species", "GSE_ID"], how="left")

    summary["flag_near_all_input_missing"] = (
        summary["scaled_data_first_sample_missing_ratio"] >= args.near_missing_threshold
    )
    summary["flag_constant_raw_clock_outputs"] = (
        (summary["raw_clock_n_samples"].fillna(0) >= 2)
        & (
            summary["raw_clock_frac_near_constant_prediction_cols"].fillna(0)
            >= args.near_constant_fraction
        )
    )
    summary["artifact_flag"] = summary["flag_near_all_input_missing"]
    summary["recommended_action"] = np.where(
        summary["artifact_flag"],
        "exclude_or_flag_rna_clock_results",
        "keep",
    )

    reasons = []
    for _, row in summary.iterrows():
        row_reasons = []
        if row.get("flag_near_all_input_missing", False):
            row_reasons.append("scaled_input_missing_ratio_ge_threshold")
        if row.get("flag_constant_raw_clock_outputs", False):
            row_reasons.append("raw_clock_outputs_near_constant")
        reasons.append(";".join(row_reasons))
    summary["flag_reason"] = reasons

    first_cols = [
        "species",
        "GSE_ID",
        "artifact_flag",
        "recommended_action",
        "flag_reason",
        "flag_near_all_input_missing",
        "flag_constant_raw_clock_outputs",
        "scaled_data_first_sample_missing_ratio",
        "scaled_data_first_sample_present_ratio",
        "n_first_sample_missing_genes",
        "n_clock_gene_cols",
        "raw_clock_n_samples",
        "raw_clock_n_prediction_cols",
        "raw_clock_n_near_constant_prediction_cols",
        "raw_clock_frac_near_constant_prediction_cols",
        "raw_clock_median_abs_range",
        "raw_clock_max_abs_range",
        "title",
        "organism",
        "type",
        "platform",
        "geo_samples_reported",
        "scaled_data_file",
        "input_summary_error",
    ]
    ordered_cols = [col for col in first_cols if col in summary.columns]
    remaining_cols = [col for col in summary.columns if col not in ordered_cols]
    summary = summary[ordered_cols + remaining_cols]

    summary = summary.sort_values(
        [
            "artifact_flag",
            "scaled_data_first_sample_missing_ratio",
            "flag_constant_raw_clock_outputs",
            "species",
            "GSE_ID",
        ],
        ascending=[False, False, False, True, True],
    )
    flagged = summary[summary["artifact_flag"]].copy()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    flagged.to_csv(output_path, index=False)
    if args.all_output:
        Path(args.all_output).parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.all_output, index=False)
    write_source_log(output_path, args, len(summary), len(flagged))

    print(f"Wrote flagged GEO table: {output_path} ({len(flagged)} rows)")
    if args.all_output:
        print(f"Wrote all-GEO QC table: {args.all_output} ({len(summary)} rows)")
    print(f"Near-all-input-missing GEOs: {int(summary['flag_near_all_input_missing'].sum())}")
    print(
        "Constant raw clock output GEOs: "
        f"{int(summary['flag_constant_raw_clock_outputs'].sum())}"
    )
    print("Flagged by species:")
    print(flagged["species"].value_counts().to_string())


def run_clock_pipeline(args):
    initialize_clocks()

    # Process mouse data
    mouse_processed_output_folder = args.mouse_processed_dir
    mouse_processed_files = np.sort(os.listdir(mouse_processed_output_folder))

    # prior_list = os.listdir("GEO_mm_rna/clock_result/")
    # gse_list = [x.split(".")[0] for x in prior_list]
    # mouse_processed_files = [
    #     x for x in mouse_processed_files if x.split("_")[0] not in gse_list
    # ]

    mouse_clock_output_folder = args.mouse_clock_output_dir
    os.makedirs(mouse_clock_output_folder, exist_ok=True)

    with Pool(args.mouse_ncores) as pool:
        pool.starmap(
            safe_process_file_folder,
            [
                (file, mouse_processed_output_folder, mouse_clock_output_folder, "mouse")
                for file in mouse_processed_files
            ],
        )

    ##
    # Process human data
    human_processed_output_folder = args.human_processed_dir
    human_processed_files = np.sort(os.listdir(human_processed_output_folder))
    human_clock_output_folder = args.human_clock_output_dir
    os.makedirs(human_clock_output_folder, exist_ok=True)

    with Pool(args.human_ncores) as pool:
        pool.starmap(
            safe_process_file_folder,
            [
                (file, human_processed_output_folder, human_clock_output_folder, "human")
                for file in human_processed_files
            ],
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run RNA clocks or summarize RNA clock imputation artifacts."
    )
    parser.add_argument(
        "--summarize-artifacts",
        action="store_true",
        help="Summarize GEOs with near-total RNA clock input missingness.",
    )
    parser.add_argument("--human-processed-dir", default="GEO_hs_rna/processed_data/")
    parser.add_argument("--mouse-processed-dir", default="GEO_mm_rna/processed_data/")
    parser.add_argument("--human-clock-output-dir", default="GEO_hs_rna/clock_result/")
    parser.add_argument("--mouse-clock-output-dir", default="GEO_mm_rna/clock_result/")
    parser.add_argument("--human-age-table", default="result/human_rna_age.csv")
    parser.add_argument("--mouse-age-table", default="result/mouse_rna_age.csv")
    parser.add_argument("--human-geo-metadata", default="gse_rna_hs.csv")
    parser.add_argument("--mouse-geo-metadata", default="gse_rna_mm.csv")
    parser.add_argument(
        "--output",
        default="result/geo_rna_imputation_artifact_flags_2026-06-14.csv",
    )
    parser.add_argument(
        "--all-output",
        default="result/geo_rna_imputation_qc_all_2026-06-14.csv",
    )
    parser.add_argument("--near-missing-threshold", type=float, default=0.95)
    parser.add_argument("--near-constant-range", type=float, default=1e-8)
    parser.add_argument("--near-constant-fraction", type=float, default=0.95)
    parser.add_argument("--age-chunksize", type=int, default=100000)
    parser.add_argument("--summary-ncores", type=int, default=8)
    parser.add_argument("--mouse-ncores", type=int, default=20)
    parser.add_argument("--human-ncores", type=int, default=5)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.summarize_artifacts:
        summarize_imputation_artifacts(args)
    else:
        run_clock_pipeline(args)
