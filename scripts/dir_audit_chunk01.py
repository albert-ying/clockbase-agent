"""
Direction-plausibility audit for chunk_01.csv (~99 mouse cases).
Reads sample.json for each case and applies biological heuristics.
Outputs: result/dir_audit_results/chunk_01.json
"""

import sys
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_system import paths

import csv
import json
import os
import re

BASE = paths.DATA
CHUNK_CSV = os.path.join(BASE, "result/dir_audit_chunks/chunk_01.csv")
OUT_DIR = os.path.join(paths.REPO, "result/dir_audit_results")
os.makedirs(OUT_DIR, exist_ok=True)

def load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return {}

# Load CSV
rows = []
with open(CHUNK_CSV) as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append(r)

# Load all sample.json files
data = []
for row in rows:
    base = paths.TRACES if row["base_dir"].startswith("result_mm_rna") else BASE
    base_dir = os.path.join(base, row["base_dir"])
    case_id = row["case_id"]
    sample = load_json(os.path.join(base_dir, case_id, "sample.json"))
    relevance = load_json(os.path.join(base_dir, case_id, "identified_relevance.json"))
    data.append({
        "row": row,
        "sample": sample,
        "relevance": relevance
    })

# Print key fields for review
for d in data:
    row = d["row"]
    s = d["sample"]
    print(f"=== {row['case_id']} ===")
    print(f"  condition_category: {row['condition_category']}")
    print(f"  condition_name: {row['condition_name']}")
    print(f"  comparison_name: {row['comparison_name']}")
    print(f"  log2fc: {row['log2fc']}")
    print(f"  title: {s.get('title','')[:120]}")
    print(f"  summary: {s.get('summary','')[:200]}")
    print()
