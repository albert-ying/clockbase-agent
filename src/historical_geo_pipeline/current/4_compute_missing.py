import glob
import os
import re
import traceback
import argparse
from datetime import datetime, timezone
from importlib import import_module
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

# Set the directories
data_dir = "./GEO_hs_dnam/series_matrix/"
out_folder = "./GEO_hs_dnam/missing_percentage"
clock_result_folder = "./GEO_hs_dnam/clock_result"

tmp = None
clock_models = None

# Create the output directory if it doesn't exist
if not os.path.exists(out_folder):
   os.makedirs(out_folder)

# List all files in the data directory
files = glob.glob(f"{data_dir}/*")

def load_clock_context():
   global tmp, clock_models
   if tmp is None:
       tmp = import_module("src.4_compute_age")
   if clock_models is None:
       from biolearn.model_gallery import ModelGallery

       clock_models = ModelGallery()


def get_gse_id(path):
   match = re.search(r"GSE\d+", str(path))
   return match.group(0) if match else None

def calculate_missing_percentage(args):
   file, clock_names = args
   try:
       load_clock_context()
       gse_id = get_gse_id(file)

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


def summarize_missing_file(path, all_missing_tolerance):
   gse_id = get_gse_id(path)
   row = {
       "GSE_ID": gse_id,
       "missing_percentage_file": str(path),
   }
   try:
       df = pd.read_csv(path)
       if df.shape[1] < 2:
           raise ValueError("missing percentage file has no clock columns")

       missing_ratios = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
       values = missing_ratios.to_numpy(dtype=float)
       finite = np.isfinite(values)
       if finite.sum() == 0:
           raise ValueError("missing percentage file has no numeric values")

       sample_mean = missing_ratios.mean(axis=1, skipna=True)
       sample_min = missing_ratios.min(axis=1, skipna=True)
       clock_min = missing_ratios.min(axis=0, skipna=True)

       row.update(
           {
               "n_samples_missing_file": int(len(df)),
               "n_missing_clock_cols": int(missing_ratios.shape[1]),
               "mean_missing_ratio": float(np.nanmean(values)),
               "median_missing_ratio": float(np.nanmedian(values)),
               "min_missing_ratio": float(np.nanmin(values)),
               "max_missing_ratio": float(np.nanmax(values)),
               "frac_pairs_eq_1": float(np.nanmean(values >= 1 - all_missing_tolerance)),
               "frac_pairs_ge_0_95": float(np.nanmean(values >= 0.95)),
               "frac_pairs_ge_0_99": float(np.nanmean(values >= 0.99)),
               "frac_samples_mean_ge_0_95": float(np.mean(sample_mean >= 0.95)),
               "frac_samples_all_eq_1": float(np.mean(sample_min >= 1 - all_missing_tolerance)),
               "n_clocks_all_samples_eq_1": int((clock_min >= 1 - all_missing_tolerance).sum()),
               "missing_summary_error": "",
           }
       )
   except Exception as e:
       row["missing_summary_error"] = str(e)
   return row


def summarize_clock_result_file(path, near_constant_range):
   gse_id = get_gse_id(path)
   row = {
       "GSE_ID": gse_id,
       "clock_result_file": str(path),
   }
   try:
       df = pd.read_csv(path)
       if df.shape[1] < 2:
           raise ValueError("clock result file has no clock columns")

       clock_values = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
       clock_values = clock_values.dropna(axis=1, how="all")
       if clock_values.shape[1] == 0:
           raise ValueError("clock result file has no numeric clock columns")

       ranges = clock_values.max(axis=0, skipna=True) - clock_values.min(axis=0, skipna=True)
       ranges = ranges.dropna()
       if len(ranges) == 0:
           raise ValueError("clock result file has no finite numeric ranges")

       row.update(
           {
               "raw_clock_n_samples": int(len(df)),
               "raw_clock_n_numeric_cols": int(clock_values.shape[1]),
               "raw_clock_n_near_constant_cols": int((ranges <= near_constant_range).sum()),
               "raw_clock_frac_near_constant_cols": float((ranges <= near_constant_range).mean()),
               "raw_clock_median_abs_range": float(ranges.median()),
               "raw_clock_max_abs_range": float(ranges.max()),
               "clock_summary_error": "",
           }
       )
   except Exception as e:
       row["clock_summary_error"] = str(e)
   return row


def read_geo_metadata(path):
   if not path or not Path(path).exists():
       return pd.DataFrame({"GSE_ID": []})
   metadata = pd.read_csv(path)
   if "accession" not in metadata.columns:
       return pd.DataFrame({"GSE_ID": []})
   keep_cols = [
       col
       for col in ["title", "organism", "type", "platform", "accession", "samples"]
       if col in metadata.columns
   ]
   metadata = metadata[keep_cols].rename(
       columns={"accession": "GSE_ID", "samples": "geo_samples_reported"}
   )
   return metadata.drop_duplicates("GSE_ID")


def summarize_intervention_results(path):
   if not path or not Path(path).exists():
       return pd.DataFrame({"GSE_ID": []})

   res = pd.read_csv(path)
   if "GSE_ID" not in res.columns:
       return pd.DataFrame({"GSE_ID": []})

   fdr = pd.to_numeric(res.get("fdr"), errors="coerce")
   p_value = pd.to_numeric(res.get("p.value"), errors="coerce")
   res = res.assign(
       significant_fdr_lt_0_05=fdr < 0.05,
       significant_p_lt_0_05=p_value < 0.05,
   )

   def unique_join(series):
       values = sorted({str(value) for value in series.dropna() if str(value) != ""})
       return ";".join(values)

   rows = []
   for gse_id, group in res.groupby("GSE_ID", dropna=False):
       significant = group[group["significant_fdr_lt_0_05"]]
       rows.append(
           {
               "GSE_ID": gse_id,
               "n_intervention_tests": int(len(group)),
               "n_significant_fdr_lt_0_05": int(group["significant_fdr_lt_0_05"].sum()),
               "n_significant_p_lt_0_05": int(group["significant_p_lt_0_05"].sum()),
               "significant_clocks_fdr_lt_0_05": unique_join(significant.get("clock", pd.Series(dtype=str))),
               "significant_perturbation_categories_fdr_lt_0_05": unique_join(
                   significant.get("perturbation_category_case", pd.Series(dtype=str))
               ),
           }
       )
   return pd.DataFrame(rows)


def newest_mtime(paths):
   existing = [Path(path) for path in paths if Path(path).exists()]
   if not existing:
       return ""
   newest = max(path.stat().st_mtime for path in existing)
   return datetime.fromtimestamp(newest, timezone.utc).isoformat()


def write_source_log(output_path, args, n_all_rows, n_flagged_rows):
   output_path = Path(output_path)
   source_log = output_path.with_suffix(".sources.txt")
   missing_files = sorted(Path(args.missing_dir).glob("GSE*.csv"))
   clock_files = sorted(Path(args.clock_result_dir).glob("GSE*.csv"))
   lines = [
       f"generated_at_utc: {datetime.now(timezone.utc).isoformat()}",
       f"working_directory: {Path.cwd()}",
       f"output_file: {output_path}",
       f"all_rows: {n_all_rows}",
       f"flagged_rows: {n_flagged_rows}",
       "",
       "parameters:",
       f"  near_missing_threshold: {args.near_missing_threshold}",
       f"  all_missing_tolerance: {args.all_missing_tolerance}",
       f"  near_constant_range: {args.near_constant_range}",
       f"  near_constant_fraction: {args.near_constant_fraction}",
       "",
       "input_sources:",
       f"  missing_dir: {args.missing_dir}",
       f"  missing_file_count: {len(missing_files)}",
       f"  missing_dir_newest_mtime_utc: {newest_mtime(missing_files)}",
       f"  clock_result_dir: {args.clock_result_dir}",
       f"  clock_result_file_count: {len(clock_files)}",
       f"  clock_result_dir_newest_mtime_utc: {newest_mtime(clock_files)}",
       f"  geo_metadata: {args.geo_metadata}",
       f"  geo_metadata_mtime_utc: {newest_mtime([args.geo_metadata])}",
       f"  intervention_results: {args.intervention_results}",
       f"  intervention_results_mtime_utc: {newest_mtime([args.intervention_results])}",
       "",
       "flag_logic:",
       "  flag_all_clock_sites_missing = min_missing_ratio >= 1 - all_missing_tolerance",
       "  flag_near_all_clock_sites_missing = mean_missing_ratio >= near_missing_threshold",
       "  flag_constant_raw_clock_outputs = raw_clock_frac_near_constant_cols >= near_constant_fraction",
       "  artifact_flag = flag_near_all_clock_sites_missing",
   ]
   output_path.parent.mkdir(parents=True, exist_ok=True)
   source_log.write_text("\n".join(lines) + "\n")


def summarize_imputation_artifacts(args):
   missing_paths = sorted(Path(args.missing_dir).glob("GSE*.csv"))
   if not missing_paths:
       raise FileNotFoundError(f"Cannot proceed - required data file not found: {args.missing_dir}/GSE*.csv")

   missing_summary = pd.DataFrame(
       summarize_missing_file(path, args.all_missing_tolerance) for path in tqdm(missing_paths)
   )

   clock_paths = sorted(Path(args.clock_result_dir).glob("GSE*.csv"))
   clock_summary = pd.DataFrame(
       summarize_clock_result_file(path, args.near_constant_range) for path in tqdm(clock_paths)
   )

   summary = missing_summary.merge(clock_summary, on="GSE_ID", how="left")
   summary = summary.merge(read_geo_metadata(args.geo_metadata), on="GSE_ID", how="left")
   summary = summary.merge(
       summarize_intervention_results(args.intervention_results), on="GSE_ID", how="left"
   )

   for col in ["n_intervention_tests", "n_significant_fdr_lt_0_05", "n_significant_p_lt_0_05"]:
       if col in summary.columns:
           summary[col] = summary[col].fillna(0).astype(int)

   summary["flag_all_clock_sites_missing"] = (
       summary["min_missing_ratio"] >= 1 - args.all_missing_tolerance
   )
   summary["flag_near_all_clock_sites_missing"] = (
       summary["mean_missing_ratio"] >= args.near_missing_threshold
   )
   summary["flag_constant_raw_clock_outputs"] = (
       (summary["raw_clock_n_samples"].fillna(0) >= 2)
       & (
           summary["raw_clock_frac_near_constant_cols"].fillna(0)
           >= args.near_constant_fraction
       )
   )
   summary["artifact_flag"] = summary["flag_near_all_clock_sites_missing"]
   summary["recommended_action"] = np.where(
       summary["artifact_flag"],
       "exclude_or_flag_intervention_hits",
       "keep",
   )
   summary["artifact_significant_fdr_lt_0_05"] = np.where(
       summary["artifact_flag"],
       summary["n_significant_fdr_lt_0_05"],
       0,
   )

   reasons = []
   for _, row in summary.iterrows():
       row_reasons = []
       if row.get("flag_all_clock_sites_missing", False):
           row_reasons.append("all_clock_sites_missing_ratio_1")
       elif row.get("flag_near_all_clock_sites_missing", False):
           row_reasons.append("mean_missing_ratio_ge_threshold")
       if row.get("flag_constant_raw_clock_outputs", False):
           row_reasons.append("raw_clock_outputs_near_constant")
       reasons.append(";".join(row_reasons))
   summary["flag_reason"] = reasons

   first_cols = [
       "GSE_ID",
       "artifact_flag",
       "recommended_action",
       "flag_reason",
       "flag_all_clock_sites_missing",
       "flag_near_all_clock_sites_missing",
       "flag_constant_raw_clock_outputs",
       "n_samples_missing_file",
       "n_missing_clock_cols",
       "mean_missing_ratio",
       "median_missing_ratio",
       "min_missing_ratio",
       "max_missing_ratio",
       "frac_pairs_eq_1",
       "frac_pairs_ge_0_95",
       "frac_pairs_ge_0_99",
       "frac_samples_mean_ge_0_95",
       "frac_samples_all_eq_1",
       "n_clocks_all_samples_eq_1",
       "raw_clock_n_samples",
       "raw_clock_n_numeric_cols",
       "raw_clock_frac_near_constant_cols",
       "raw_clock_median_abs_range",
       "raw_clock_max_abs_range",
       "n_intervention_tests",
       "n_significant_fdr_lt_0_05",
       "n_significant_p_lt_0_05",
       "artifact_significant_fdr_lt_0_05",
       "significant_clocks_fdr_lt_0_05",
       "significant_perturbation_categories_fdr_lt_0_05",
       "title",
       "organism",
       "type",
       "platform",
       "geo_samples_reported",
       "missing_percentage_file",
       "clock_result_file",
       "missing_summary_error",
       "clock_summary_error",
   ]
   ordered_cols = [col for col in first_cols if col in summary.columns]
   remaining_cols = [col for col in summary.columns if col not in ordered_cols]
   summary = summary[ordered_cols + remaining_cols]

   summary = summary.sort_values(
       [
           "artifact_flag",
           "flag_all_clock_sites_missing",
           "mean_missing_ratio",
           "artifact_significant_fdr_lt_0_05",
           "GSE_ID",
       ],
       ascending=[False, False, False, False, True],
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
   print(f"Strict all-missing GEOs: {int(summary['flag_all_clock_sites_missing'].sum())}")
   print(f"Near-all-missing GEOs: {int(summary['flag_near_all_clock_sites_missing'].sum())}")
   print(f"Constant raw clock output GEOs: {int(summary['flag_constant_raw_clock_outputs'].sum())}")
   print(
       "Artifact significant FDR<0.05 rows: "
       f"{int(flagged['artifact_significant_fdr_lt_0_05'].sum())}"
   )


def parse_args():
   parser = argparse.ArgumentParser(
       description="Compute DNAm clock missing percentages or summarize imputation artifacts."
   )
   parser.add_argument(
       "--summarize-artifacts",
       action="store_true",
       help="Summarize GEOs with near-total clock-site imputation and near-constant clock outputs.",
   )
   parser.add_argument("--missing-dir", default=out_folder)
   parser.add_argument("--clock-result-dir", default=clock_result_folder)
   parser.add_argument("--geo-metadata", default="./gse_me_hs.csv")
   parser.add_argument("--intervention-results", default="./result/hs_dnam_core_res.csv")
   parser.add_argument(
       "--output",
       default="./result/geo_dnam_imputation_artifact_flags_2026-06-14.csv",
   )
   parser.add_argument(
       "--all-output",
       default="./result/geo_dnam_imputation_qc_all_2026-06-14.csv",
   )
   parser.add_argument("--near-missing-threshold", type=float, default=0.95)
   parser.add_argument("--all-missing-tolerance", type=float, default=1e-12)
   parser.add_argument("--near-constant-range", type=float, default=1e-8)
   parser.add_argument("--near-constant-fraction", type=float, default=0.95)
   parser.add_argument("--ncores", type=int, default=12)
   return parser.parse_args()

# Specify the path to a single dataset file
if __name__ == "__main__":
  args = parse_args()
  if args.summarize_artifacts:
      summarize_imputation_artifacts(args)
      raise SystemExit(0)

  load_clock_context()
  model_definitions = list(clock_models.model_definitions.keys())
  files_clock_definitions = [(file, model_definitions) for file in files]
  with Pool(args.ncores) as p:
      list(
          tqdm(
              p.imap_unordered(calculate_missing_percentage, files_clock_definitions),
              total=len(files_clock_definitions)
          )
      )
