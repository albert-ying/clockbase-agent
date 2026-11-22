#!/usr/bin/env python3
"""Apply bounded supplementary worksheet corrections.

The default operation withdraws S8. The opt-in S1 operation reconciles exact
accession identities from pinned correction tables, retaining unresolved fields.
Only the selected worksheet XML changes; all other XLSX members are preserved.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
import posixpath
import re
import shutil
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
DOC_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"


def locate_sheet_member(xlsx: Path, sheet_name: str) -> str:
    with ZipFile(xlsx) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))

    sheet = workbook.find(
        f".//{{{MAIN_NS}}}sheet[@name='{sheet_name}']",
    )
    if sheet is None:
        raise ValueError(f"Worksheet not found: {sheet_name}")
    relationship_id = sheet.attrib[f"{{{DOC_REL_NS}}}id"]

    relationship = relationships.find(
        f".//{{{PKG_REL_NS}}}Relationship[@Id='{relationship_id}']",
    )
    if relationship is None:
        raise ValueError(f"Worksheet relationship not found: {relationship_id}")
    target = relationship.attrib["Target"].lstrip("/")
    if target.startswith("xl/"):
        return target
    return posixpath.normpath(posixpath.join("xl", target))


def replacement_xml() -> bytes:
    title = (
        "Supplementary Table S8. Candidate-level novelty classification "
        "withdrawn during revision"
    )
    note = (
        "The candidate table was removed because the 106-, 96-, and 92-record "
        "versions could not be traced to one reproducible cohort definition. "
        "The manuscript and response no longer report the classification "
        "percentages or rank statistics."
    )
    root = ET.Element(f"{{{MAIN_NS}}}worksheet")
    ET.SubElement(root, f"{{{MAIN_NS}}}dimension", {"ref": "A1:A3"})
    sheet_views = ET.SubElement(root, f"{{{MAIN_NS}}}sheetViews")
    ET.SubElement(
        sheet_views,
        f"{{{MAIN_NS}}}sheetView",
        {"showGridLines": "0", "workbookViewId": "0"},
    )
    ET.SubElement(root, f"{{{MAIN_NS}}}sheetFormatPr", {"defaultRowHeight": "15"})
    columns = ET.SubElement(root, f"{{{MAIN_NS}}}cols")
    ET.SubElement(
        columns,
        f"{{{MAIN_NS}}}col",
        {"min": "1", "max": "1", "width": "120", "customWidth": "1"},
    )
    data = ET.SubElement(root, f"{{{MAIN_NS}}}sheetData")
    for row_number, text, style, height in ((1, title, "1", "24"), (3, note, "0", "45")):
        row = ET.SubElement(
            data,
            f"{{{MAIN_NS}}}row",
            {"r": str(row_number), "ht": height, "customHeight": "1"},
        )
        cell = ET.SubElement(
            row,
            f"{{{MAIN_NS}}}c",
            {"r": f"A{row_number}", "s": style, "t": "inlineStr"},
        )
        inline = ET.SubElement(cell, f"{{{MAIN_NS}}}is")
        text_node = ET.SubElement(inline, f"{{{MAIN_NS}}}t")
        text_node.text = text
    ET.SubElement(
        root,
        f"{{{MAIN_NS}}}pageMargins",
        {"left": "0.7", "right": "0.7", "top": "0.75", "bottom": "0.75", "header": "0.3", "footer": "0.3"},
    )
    ET.register_namespace("", MAIN_NS)
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def replace_sheet(input_path: Path, output_path: Path, sheet_name: str = "S8", payload: bytes | None = None) -> None:
    member = locate_sheet_member(input_path, sheet_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=output_path.parent,
        prefix=f".{output_path.stem}.",
        suffix=".xlsx",
        delete=False,
    ) as handle:
        temporary_path = Path(handle.name)
    try:
        with ZipFile(input_path) as source, ZipFile(temporary_path, "w", ZIP_DEFLATED) as destination:
            for item in source.infolist():
                data = (replacement_xml() if payload is None else payload) if item.filename == member else source.read(item.filename)
                destination.writestr(item, data)
        os.replace(temporary_path, output_path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def pinned_file(item: dict) -> Path:
    path = Path(item["path"]).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Cannot proceed - required data file not found: {path}")
    if path.stat().st_size != item["size_bytes"] or digest(path) != item["sha256"]:
        raise ValueError(f"Pinned source differs: {path}")
    return path


def source_roster(path: Path, modality: str) -> tuple[dict, list[dict], int]:
    """Require structural CSV identity before interpreting any gse_id field."""
    groups: dict[tuple[str, str], list[int]] = {}
    unresolved = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader)
        if len(set(header)) != len(header) or header.count("gse_id") != 1:
            raise ValueError(f"Unique source headers and one gse_id required: {path}")
        column = header.index("gse_id")
        comparison = header.index("comparison_name") if "comparison_name" in header else None
        count = 0
        for record_number, row in enumerate(reader, 2):
            count += 1
            if len(row) != len(header):
                raise ValueError(f"CSV record {record_number} has {len(row)} fields; expected {len(header)}: {path}")
            accession = row[column]
            if not re.fullmatch(r"GSE[0-9]+", accession):
                unresolved.append({"modality": modality, "source_file": path.name,
                    "source_record_number": record_number, "source_identifier": accession,
                    "comparison_name": row[comparison] if comparison is not None else "",
                    "status": "CANNOT VERIFY", "reason": "Missing or non-single study identifier; no assignment inferred"})
                continue
            groups.setdefault((accession, modality), []).append(record_number)
    if count == 0 or not groups:
        raise ValueError(f"Source roster has no named study records: {path}")
    return groups, unresolved, count


def read_s1(archive: ZipFile, member: str) -> tuple[ET.Element, list[dict]]:
    ns = "{" + MAIN_NS + "}"
    root = ET.fromstring(archive.read(member))
    strings = []
    if "xl/sharedStrings.xml" in archive.namelist():
        strings = ["".join(s.itertext()) for s in ET.fromstring(archive.read("xl/sharedStrings.xml"))]
    def value(cell):
        if cell.get("t") == "inlineStr":
            return "".join(node.text or "" for node in cell.iter(ns + "t"))
        node = cell.find(ns + "v")
        if node is None:
            return ""
        return strings[int(node.text)] if cell.get("t") == "s" else (node.text or "")
    rows = []
    for row in root.find(ns + "sheetData"):
        number = int(row.attrib["r"])
        cells = {re.sub(r"[0-9]+$", "", c.attrib["r"]): value(c) for c in row}
        if number == 3 and [cells.get(k) for k in "ABC"] != ["GEO_accession", "Modality", "GEO_title"]:
            raise ValueError("Unexpected S1 header; no automatic worksheet rewrite")
        if number >= 4 and cells.get("A", ""):
            if any(cells.get(k, "") for k in "DEFGH"):
                raise ValueError("Existing S1 extra fields would be overwritten")
            rows.append({"source_workbook_row": number, "accession": cells["A"],
                         "modality": cells.get("B", ""), "title": cells.get("C", "")})
    if not rows:
        raise ValueError("S1 contains no source rows")
    pairs = [(r["accession"], r["modality"]) for r in rows]
    if len(pairs) != len(set(pairs)):
        raise ValueError("Duplicate S1 identities require explicit reconciliation")
    return root, rows


def roster_xml(root: ET.Element, roster: list[dict], unresolved: list[dict]) -> bytes:
    ns = "{" + MAIN_NS + "}"
    data = root.find(ns + "sheetData")
    data.clear()
    columns = root.find(ns + "cols")
    if columns is None:
        columns = ET.Element(ns + "cols")
        root.insert(list(root).index(data), columns)
    columns.clear()
    for number, width in enumerate((18, 28, 80, 4, 44, 28, 42, 62), 1):
        ET.SubElement(columns, ns + "col", {"min": str(number), "max": str(number),
                      "width": str(width), "customWidth": "1"})
    def row(number, values, style=None):
        node = ET.SubElement(data, ns + "row", {"r": str(number)})
        for column, text in values.items():
            attrs = {"r": f"{column}{number}", "t": "inlineStr"}
            if style is not None:
                attrs["s"] = style
            cell = ET.SubElement(node, ns + "c", attrs)
            ET.SubElement(ET.SubElement(cell, ns + "is"), ns + "t").text = str(text)
    row(1, {"A": "S1: Source-table accession roster; exact single-GSE identifiers"}, "1")
    row(2, {"A": "Reconciled correction snapshot; source lineage and extension eligibility remain unverified. Exceptions are retained in columns E–H."})
    row(3, {"A": "GEO_accession", "B": "Modality", "C": "GEO_title", "E": "Unresolved source record",
            "F": "Source identifier", "G": "Comparison name", "H": "Disposition"}, "2")
    for index in range(max(len(roster), len(unresolved))):
        fields = {}
        if index < len(roster):
            item = roster[index]
            fields.update({"A": item["accession"], "B": item["modality"], "C": item["title"]})
        if index < len(unresolved):
            item = unresolved[index]
            fields.update({"E": f'{item["modality"]}; CSV record {item["source_record_number"]}',
                "F": item["source_identifier"] or "[blank]", "G": item["comparison_name"] or "[blank]",
                "H": "CANNOT VERIFY; no study assignment inferred"})
        row(index + 4, fields)
    root.find(ns + "dimension").set("ref", f"A1:H{max(len(roster), len(unresolved)) + 3}")
    merges = root.find(ns + "mergeCells")
    if merges is not None:
        ET.SubElement(merges, ns + "mergeCell", {"ref": "A2:H2"})
        merges.set("count", str(len(merges)))
    ET.register_namespace("", MAIN_NS)
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def reconcile_s1(plan_path: Path, input_path: Path, output_path: Path, receipt_path: Path) -> dict:
    started = time.monotonic()
    plan = json.loads(plan_path.read_text())
    if plan.get("schema") != "clockbase-supplementary-roster-plan-v1":
        raise ValueError("Unknown source-roster plan")
    if output_path.exists() or receipt_path.exists() or output_path.resolve() == input_path.resolve():
        raise FileExistsError("Source-roster outputs must be new and distinct from the original workbook")
    inputs = [plan["input_workbook"], *plan["correction_sources"], *plan["identity_witnesses"]]
    for item in inputs:
        pinned_file(item)
    if input_path.resolve() != Path(plan["input_workbook"]["path"]).resolve():
        raise ValueError("Workbook differs from the frozen plan")
    expected_modalities = {"Human RNA-seq", "Mouse RNA-seq", "Human DNA methylation"}
    if {x["modality"] for x in plan["correction_sources"]} != expected_modalities or len(plan["correction_sources"]) != 3:
        raise ValueError("Exactly one source per supported modality is required")
    groups = {}; unresolved = []; source_rows = {}
    for item in plan["correction_sources"]:
        records, exceptions, count = source_roster(Path(item["path"]), item["modality"])
        groups.update(records); unresolved.extend(exceptions); source_rows[item["modality"]] = count
    member = locate_sheet_member(input_path, "S1")
    with ZipFile(input_path) as archive:
        sheet, original = read_s1(archive, member)
    old = {(r["accession"], r["modality"]): r for r in original}
    roster = [{"accession": accession, "modality": modality, "title": old.get((accession, modality), {}).get("title", ""),
               "source_record_numbers": records} for (accession, modality), records in sorted(groups.items(), key=lambda x: (x[0][1], x[0][0]))]
    removed = [r for r in original if (r["accession"], r["modality"]) not in groups]
    added = [r for r in roster if (r["accession"], r["modality"]) not in old]
    if any(re.fullmatch(r"GSE[0-9]+", r["accession"]) for r in removed):
        raise ValueError("A previously named single-GSE row would disappear; requires a separate reviewed plan")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    backup = output_path.parent / "original_workbook.xlsx"
    with backup.open("xb") as destination, input_path.open("rb") as source:
        shutil.copyfileobj(source, destination)
    replace_sheet(input_path, output_path, "S1", roster_xml(sheet, roster, unresolved))
    unchanged = []
    with ZipFile(input_path) as before, ZipFile(output_path) as after:
        if before.namelist() != after.namelist() or after.testzip() is not None:
            raise ValueError("Workbook package membership or ZIP integrity changed")
        for name in before.namelist():
            if name == member:
                continue
            if before.read(name) != after.read(name):
                raise ValueError(f"Non-S1 package member changed: {name}")
            unchanged.append(name)
    for item in inputs:
        pinned_file(item)
    result = {"schema": "clockbase-supplementary-roster-receipt-v1", "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS", "scientific_status": "CANNOT VERIFY", "plan_sha256": digest(plan_path),
        "implementation_sha256": digest(Path(__file__)), "inputs": inputs,
        "output": {"path": str(output_path), "sha256": digest(output_path)},
        "original_copy": {"path": str(backup), "sha256": digest(backup)}, "source_comparison_rows": source_rows,
        "original_roster_rows": len(original), "single_gse_roster_rows": len(roster),
        "unique_accessions_across_modalities": len({r["accession"] for r in roster}),
        "by_modality": dict(Counter(r["modality"] for r in roster)), "removed_original_entries": removed,
        "added_entries": added, "unresolved_source_records": unresolved, "roster": roster,
        "unchanged_workbook_members": unchanged, "scope": plan["rules"]["scope_note"],
        "elapsed_seconds": time.monotonic() - started, "network_requests": 0, "model_calls": 0, "clock_executions": 0}
    with receipt_path.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2); handle.write("\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reconcile-s1-plan", type=Path)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    if args.reconcile_s1_plan:
        if args.receipt is None:
            parser.error("S1 reconciliation requires --receipt")
        result = reconcile_s1(args.reconcile_s1_plan.resolve(), args.input.resolve(), args.output.resolve(), args.receipt.resolve())
        print(json.dumps({k: result[k] for k in ("status", "single_gse_roster_rows", "by_modality", "scientific_status")}))
        return
    replace_sheet(args.input.resolve(), args.output.resolve())
    print(f"Replaced S8 worksheet XML in {args.output}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, csv.Error, KeyError, ET.ParseError) as error:
        print(f"Supplementary correction stopped: {type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(1) from error
