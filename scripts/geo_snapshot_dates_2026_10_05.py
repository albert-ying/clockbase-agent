#!/usr/bin/env python3
"""Summarize local GEO series dates and ARCHS4 files for the atlas studies."""
import argparse
import csv
from datetime import date, datetime, timezone
import gzip
import json
from pathlib import Path
import re
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_system import paths


ATLAS_TABLES = {
    "dnam": "2026-09-08-DNAm-atlas_WORKING_NOT_VALIDATED.csv",
    "hs_rna": "2026-09-08-human-RNA-atlas_WORKING_NOT_VALIDATED.csv",
    "mm_rna": "2026-09-08-mouse-RNA-atlas_WORKING_NOT_VALIDATED.csv",
}
DATE_FIELDS = {
    "public_date": "!Series_status",
    "submission_date": "!Series_submission_date",
    "last_update_date": "!Series_last_update_date",
}
SERIES_RE = re.compile(r"^GSE\d+$")


def atlas_studies():
    """Read the unique series identifiers in each correction table."""
    studies = {}
    for modality, filename in ATLAS_TABLES.items():
        table = paths.CORRECTIONS / filename
        with table.open(newline="", encoding="utf-8-sig") as handle:
            studies[modality] = sorted({
                row["gse_id"].strip()
                for row in csv.DictReader(handle)
                if row.get("gse_id", "").strip()
            })
    return studies


def _series_files():
    """Index local GEO series matrices by accession without reading data matrices."""
    indexed = {}
    for path in paths.GEO.rglob("*_series_matrix.txt*"):
        match = re.match(r"(GSE\d+)(?:_GPL\d+)?_series_matrix\.txt(?:\.gz)?$", path.name)
        if match:
            indexed.setdefault(match.group(1), []).append(path)
    for candidates in indexed.values():
        candidates.sort(key=lambda item: ("series_matrix" not in item.parts, str(item)))
    return indexed


def _metadata(path):
    raw_values = {name: [] for name in DATE_FIELDS}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_") or line.startswith("!series_matrix_table_begin"):
                break
            for name, tag in DATE_FIELDS.items():
                match = re.match(rf"^{re.escape(tag)}(?:\s+|\t+)(.*)$", line)
                if match:
                    raw_values[name].append(match.group(1).strip().strip('"'))
    values = {}
    date_values = {}
    for name, raw in raw_values.items():
        parsed = []
        for value in raw:
            if name == "public_date":
                status_match = re.search(r"Public on (.+)$", value)
                value = status_match.group(1) if status_match else ""
            if not value:
                continue
            for date_format in ("%b %d %Y", "%d-%b-%Y"):
                try:
                    parsed.append(datetime.strptime(value, date_format).date().isoformat())
                    break
                except ValueError:
                    pass
        date_values[name] = sorted(set(parsed))
        values[name] = date_values[name][-1] if date_values[name] else None
    conflicts = [name for name, found in date_values.items() if len(found) > 1]
    values["date_values"] = date_values
    values["date_conflicts"] = conflicts
    return values


def _percentile(values, percentile):
    """Return a linearly interpolated percentile of ISO calendar-date values."""
    ordered = sorted(date.fromisoformat(value).toordinal() for value in values)
    if not ordered:
        return None
    position = (len(ordered) - 1) * percentile / 100
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    ordinal = ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)
    return date.fromordinal(round(ordinal)).isoformat()


def _archs4_files():
    records = []
    for path in sorted(paths.GEO.rglob("*_gene_v*.h5")):
        stat = path.stat()
        records.append({
            "path": paths.portable(path),
            "version": path.name,
            "bytes": stat.st_size,
            "mtime": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        })
    return records


def build_report():
    studies = atlas_studies()
    matrix_index = _series_files()
    study_rows = []
    modalities = {}
    for modality, accessions in studies.items():
        dated_public = []
        found = 0
        for accession in accessions:
            components = re.findall(r"GSE\d+", accession)
            component_records = []
            for component in components:
                candidates = matrix_index.get(component, [])
                candidate_records = []
                for matrix in candidates:
                    stat = matrix.stat()
                    candidate_records.append({
                        "metadata_path": paths.portable(matrix),
                        "mtime": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                        **_metadata(matrix),
                    })
                fields = {}
                date_values = {}
                for field in DATE_FIELDS:
                    date_values[field] = sorted({
                        value for record in candidate_records for value in record["date_values"][field]
                    })
                    fields[field] = date_values[field][-1] if date_values[field] else None
                conflicts = sorted({
                    field for record in candidate_records for field in record["date_conflicts"]
                } | {field for field, values in date_values.items() if len(values) > 1})
                component_records.append({
                    "gse_id": component,
                    "metadata_found": bool(candidate_records),
                    "metadata_candidates": candidate_records,
                    "date_values": date_values,
                    "date_conflicts": conflicts,
                    **fields,
                })
            metadata_found = any(record["metadata_found"] for record in component_records)
            found += int(metadata_found)
            fields = {}
            date_values = {}
            for field in DATE_FIELDS:
                date_values[field] = sorted({
                    value for record in component_records for value in record["date_values"][field]
                })
                fields[field] = date_values[field][-1] if date_values[field] else None
            conflicts = sorted({field for record in component_records for field in record["date_conflicts"]} |
                               {field for field, values in date_values.items() if len(values) > 1})
            if fields["public_date"]:
                dated_public.append(fields["public_date"])
            study_rows.append({
                "modality": modality,
                "gse_id": accession,
                "metadata_found": metadata_found,
                "metadata_path": ";".join(
                    candidate["metadata_path"] for record in component_records
                    for candidate in record["metadata_candidates"]
                ),
                "series_components": json.dumps(component_records, separators=(",", ":")),
                "date_values": json.dumps(date_values, separators=(",", ":")),
                "date_conflicts": ";".join(conflicts),
                **fields,
            })
        counts = {}
        for field in DATE_FIELDS:
            counts[field] = sum(row["modality"] == modality and row[field] is not None for row in study_rows)
        modalities[modality] = {
            "correction_table": paths.portable(paths.CORRECTIONS / ATLAS_TABLES[modality]),
            "studies": len(accessions),
            "metadata_found": found,
            "studies_with_public_date": len(dated_public),
            "date_field_counts": counts,
            "date_field_missing_counts": {field: len(accessions) - counts[field] for field in DATE_FIELDS},
            "studies_with_conflicting_dates": sum(
                row["modality"] == modality and bool(row["date_conflicts"]) for row in study_rows
            ),
            "latest_public_date": max(dated_public) if dated_public else None,
            "public_date_distribution": {
                "p50": _percentile(dated_public, 50),
                "p90": _percentile(dated_public, 90),
                "p99": _percentile(dated_public, 99),
                "method": "linear interpolation of ordinal calendar dates",
            },
        }
    return {
        "schema": "clockbase-geo-snapshot-dates-v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "geo_root": paths.portable(paths.GEO),
        "limitations": [
            "Series public dates and local file modification times do not establish when ClockBase acquired its inputs.",
            "For composite GSE identifiers, the row date is the latest recorded date among its component series; component dates are retained.",
        ],
        "correction_tables": [
            {"path": paths.portable(paths.CORRECTIONS / filename),
             "mtime": datetime.fromtimestamp((paths.CORRECTIONS / filename).stat().st_mtime,
                                             timezone.utc).isoformat()}
            for filename in ATLAS_TABLES.values()
        ],
        "modalities": modalities,
        "archs4_files": _archs4_files(),
    }, study_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path, help="new directory for JSON and per-study CSV")
    args = parser.parse_args()
    try:
        args.output_dir.mkdir(parents=True, exist_ok=False)
        report, rows = build_report()
        json_path = args.output_dir / "geo_snapshot_dates.json"
        csv_path = args.output_dir / "geo_snapshot_studies.csv"
        json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=("modality", "gse_id", "metadata_found", "metadata_path",
                                                        "series_components", "date_values", "date_conflicts",
                                                        "public_date", "submission_date", "last_update_date"))
            writer.writeheader()
            writer.writerows(rows)
    except (OSError, ValueError) as exc:
        parser.error(f"could not write GEO snapshot report: {exc}")
    print(json.dumps({"json": str(json_path), "csv": str(csv_path), "modalities": report["modalities"],
                      "archs4_files": report["archs4_files"]}, indent=2))


if __name__ == "__main__":
    main()
