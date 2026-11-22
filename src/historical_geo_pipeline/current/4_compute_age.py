import csv
import glob
import gzip
import hashlib
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


def require_observed_beta_matrix(data):
    """Reject invalid study input before reindexing or population imputation.

    This is a necessary input check, not validation of per-clock CpG coverage.
    """
    if not isinstance(data, pd.DataFrame) or data.shape[0] == 0 or data.shape[1] == 0:
        raise ValueError("AUDIT_FAIL: no observed feature rows or sample columns")
    if data.index.has_duplicates or data.columns.has_duplicates:
        raise ValueError("AUDIT_FAIL: duplicate feature or sample identifiers")
    if any(pd.isna(v) or not str(v).strip() for v in [*data.index, *data.columns]):
        raise ValueError("AUDIT_FAIL: missing feature or sample identifier")
    values = data.to_numpy(dtype=float)
    finite = np.isfinite(values)
    if np.isinf(values).any() or (finite.sum(axis=0) == 0).any():
        raise ValueError("AUDIT_FAIL: infinite values or sample without observed features")
    observed = values[finite]
    if (observed < 0).any() or (observed > 1).any():
        raise ValueError("AUDIT_FAIL: input is outside the beta-value range [0,1]")
    return data

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


GEO_TABLE_BEGIN = "!series_matrix_table_begin"
GEO_TABLE_END = "!series_matrix_table_end"
GEO_SAMPLE_ACCESSION = "!Sample_geo_accession"


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _open_geo_matrix(filename):
    if str(filename).endswith(".gz"):
        return gzip.open(filename, "rt", encoding="utf-8", newline="")
    return open(filename, "r", encoding="utf-8", newline="")


def _scan_geo_matrix(file):
    """Validate one GEO matrix as logical tab-delimited records.

    The returned seek cookie points to the primary matrix header.  Reading from
    that cookie avoids converting physical preamble line counts into pandas
    ``skiprows`` counts when quoted metadata contains embedded newlines.
    """
    reader = csv.reader(iter(file.readline, ""), delimiter="\t", strict=True)
    declared_samples = None
    try:
        for record in reader:
            if record and record[0] == GEO_SAMPLE_ACCESSION:
                if declared_samples is not None:
                    raise ValueError(
                        "AUDIT_FAIL: duplicate GEO sample-accession metadata record"
                    )
                declared_samples = record[1:]
            if record == [GEO_TABLE_BEGIN]:
                table_start = file.tell()
                break
        else:
            raise ValueError("AUDIT_FAIL: GEO matrix begin marker is missing")

        header = next(reader, None)
        if header is None or len(header) < 2 or header[0] != "ID_REF":
            raise ValueError("AUDIT_FAIL: invalid GEO primary matrix header")
        sample_ids = header[1:]
        if (
            len(set(sample_ids)) != len(sample_ids)
            or any(re.fullmatch(r"GSM[0-9]+", value) is None for value in sample_ids)
        ):
            raise ValueError(
                "AUDIT_FAIL: source sample identifiers must be unique canonical GSM accessions"
            )
        if declared_samples is not None and declared_samples != sample_ids:
            raise ValueError(
                "AUDIT_FAIL: GEO sample-accession metadata differs from the primary header"
            )

        feature_ids = []
        seen_features = set()
        expected_width = len(header)
        for record_number, record in enumerate(reader, start=2):
            if record == [GEO_TABLE_END]:
                if next(reader, None) is not None:
                    raise ValueError(
                        "AUDIT_FAIL: records follow the GEO matrix end marker"
                    )
                break
            if record and record[0] == GEO_TABLE_END:
                raise ValueError(
                    "AUDIT_FAIL: GEO matrix end marker is misplaced"
                )
            if len(record) != expected_width:
                raise ValueError(
                    "AUDIT_FAIL: GEO matrix record "
                    f"{record_number} has {len(record)} fields; expected {expected_width}"
                )
            feature_id = record[0]
            if (
                not feature_id.strip()
                or feature_id in {GEO_TABLE_BEGIN, GEO_TABLE_END}
                or feature_id in seen_features
            ):
                raise ValueError(
                    "AUDIT_FAIL: missing or duplicate source feature identifier"
                )
            seen_features.add(feature_id)
            feature_ids.append(feature_id)
        else:
            raise ValueError("AUDIT_FAIL: GEO matrix end marker is missing")
    except csv.Error as exc:
        raise ValueError(f"AUDIT_FAIL: malformed quoted GEO record: {exc}") from exc

    return table_start, header, feature_ids


def read_gse_matrix(filename, data=True):
    def download_file(url, local_filename):
        with requests.get(url, stream=True) as r:
            with open(local_filename, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        return local_filename

    # Check if filename is a URL
    filename = str(filename)
    if filename.startswith("http://") or filename.startswith("https://"):
        filename = download_file(filename, os.path.basename(filename))
    source_sha256_before = _sha256_file(filename)
    parse_error = None
    parsed = None
    try:
        with _open_geo_matrix(filename) as file:
            table_start, source_header, source_features = _scan_geo_matrix(file)
            file.seek(table_start)
            parsed = pd.read_csv(
                file,
                sep="\t",
                header=0,
                nrows=len(source_features),
                converters={0: str},
                float_precision="round_trip",
            )
    except Exception as exc:
        parse_error = exc
    try:
        source_sha256_after = _sha256_file(filename)
    except OSError as exc:
        raise ValueError(
            "AUDIT_FAIL: GEO source could not be hashed after parsing"
        ) from exc
    if source_sha256_after != source_sha256_before:
        raise ValueError("AUDIT_FAIL: GEO source bytes changed during parsing")
    if parse_error is not None:
        raise parse_error

    if list(parsed.columns) != source_header or len(parsed) != len(source_features):
        raise ValueError("AUDIT_FAIL: parsed GEO table differs from source scan")
    parsed_features = parsed.iloc[:, 0].astype(str).tolist()
    if parsed_features != source_features:
        raise ValueError("AUDIT_FAIL: parsed GEO feature order differs from source scan")

    parsed.index = parsed.iloc[:, 0]
    parsed.index.name = "id"
    return parsed.drop(parsed.columns[0], axis=1)


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
            require_observed_beta_matrix(test_data)
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
