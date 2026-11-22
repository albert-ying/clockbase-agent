"""Add the source data of Figure 7k to the source-data workbook (sheet "Fig 7k", a README row and a Statistics row).

Figure 7k is the PLX5622 meta-analysis of the verified atlas on the scale of the mouse mortality clock: the difference in
mean clock output (PLX5622 minus control) of each comparison, averaged within each project and pooled across projects by
DerSimonian-Laird random effects, as build_figure7_2026_11_08.py draws it. This script recomputes those values from the
same inputs with the same formulas and stops unless they equal the values the figure builder recorded
(evidence/figure7-rebuild-2026-11-08/numbers.json). Nothing is typed in; the workbook is otherwise unchanged.

  python add_fig7k_source_data_2026_11_11.py
"""
import json
import math
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT.parent
REV = ROOT / "paper/revision_2026-09-07"
XLSX = REV / "deliverables/FINAL_REVISION_2026-11-06/SourceData_Fig7_EDFig25_EDFig26.xlsx"
META = LOCAL / "data/clockbase-verifier-v3-2026-09-27/intervention-meta-2026-10-17/primary"
NUMBERS = REV / "evidence/figure7-rebuild-2026-11-08/numbers.json"
for p in (XLSX, META / "meta_input_effects.csv", NUMBERS):
    if not p.exists():
        raise SystemExit(f"Cannot proceed - required data file not found: {p}")


def dersimonian_laird(y, v):
    w = 1 / v
    m = (w * y).sum() / w.sum()
    q = (w * (y - m) ** 2).sum()
    tau2 = max(0.0, (q - (len(y) - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
    ws = 1 / (v + tau2)
    est = float((ws * y).sum() / ws.sum())
    se = math.sqrt(1 / ws.sum())
    return est, est - 1.959964 * se, est + 1.959964 * se, tau2


mi = pd.read_csv(META / "meta_input_effects.csv", low_memory=False)
cmp_k = mi[(mi.canonical_target == "plx5622") & (mi.condition_category == "Drug")].copy()
assert cmp_k["clock"].nunique() == 1
cmp_k["var_md"] = (cmp_k.v3_mean_diff / cmp_k.v3_cohen_d) ** 2 * (1 / cmp_k.v3_n_case + 1 / cmp_k.v3_n_control)
units = cmp_k.groupby("project").agg(md=("v3_mean_diff", "mean"), var_md=("var_md", "mean"),
                                     n_cmp=("g", "size")).reset_index()
units["lo"] = units.md - 1.959964 * np.sqrt(units.var_md)
units["hi"] = units.md + 1.959964 * np.sqrt(units.var_md)
pooled = dersimonian_laird(units.md, units.var_md)

recorded = json.loads(NUMBERS.read_text())
assert np.allclose(pooled[:3], recorded["k_pooled_mean_difference"], rtol=1e-9), (pooled, recorded["k_pooled_mean_difference"])
labels = {}
for u in recorded["k_units"]:
    row = units[units.project == u["project"]].iloc[0]
    assert np.allclose([row.md, row.lo, row.hi], [u["md"], u["lo"], u["hi"]], rtol=1e-9) and row.n_cmp == u["n_cmp"]
    labels[u["project"]] = u["label"]

wb = openpyxl.load_workbook(XLSX)
if "Fig 7k" in wb.sheetnames:
    del wb["Fig 7k"]
ws = wb.create_sheet("Fig 7k", index=wb.sheetnames.index("Fig 7j") + 1)
clock = cmp_k["clock"].iloc[0]
ws.append(["Figure 7k. PLX5622 in the verified atlas: difference in mean mouse mortality-clock output (PLX5622 minus "
           "control) by comparison, averaged within each project and pooled across projects by DerSimonian–Laird "
           "random effects"])
ws.append([f"Clock: {clock}. Variance of a mean difference = (mean difference / Cohen's d)^2 × (1/n case + 1/n control); "
           "project values are means over the project's comparisons; intervals are 95% (normal)."])
ws.append([])
ws.append(["Comparisons"])
ws.append(["GEO series", "Comparison", "Project", "n PLX5622", "n control", "Mean difference (clock output)",
           "Cohen's d", "Variance of the mean difference"])
for r in cmp_k.sort_values(["project", "gse_id", "comparison_name"]).itertuples():
    ws.append([r.gse_id, r.comparison_name, r.project, int(r.v3_n_case), int(r.v3_n_control), float(r.v3_mean_diff),
               float(r.v3_cohen_d), float(r.var_md)])
ws.append([])
ws.append(["Projects (as plotted)"])
ws.append(["Project", "Label in the figure", "Comparisons", "Mean difference", "95% CI low", "95% CI high",
           "Variance"])
for r in units.sort_values("md").itertuples():
    ws.append([r.project, labels[r.project], int(r.n_cmp), float(r.md), float(r.lo), float(r.hi), float(r.var_md)])
ws.append([])
ws.append(["Pooled (diamond)", "DerSimonian–Laird random effects", len(units), pooled[0], pooled[1], pooled[2],
           f"tau^2 = {pooled[3]:.6g}"])

readme = wb["README"]
rows = [[c.value for c in row] for row in readme.iter_rows()]
if not any(r and r[0] == "Fig 7k" for r in rows):
    at = next(i for i, r in enumerate(rows, start=1) if r and r[0] == "Fig 7j") + 1
    readme.insert_rows(at)
    for j, v in enumerate(["Fig 7k", "Fig. 7k", "PLX5622 comparisons of the verified atlas: mortality-clock mean "
                           "differences per comparison, project means and the random-effects pool"], start=1):
        readme.cell(row=at, column=j, value=v)

stats = wb["Statistics"]
srows = [[c.value for c in row] for row in stats.iter_rows()]
head = next(i for i, r in enumerate(srows) if r and "Panel" in r)
cols = {name: j for j, name in enumerate(srows[head]) if name}
if not any(r and r[cols["Panel"]] == "Fig. 7k" for r in srows):
    at = next(i for i, r in enumerate(srows, start=1) if r and r[cols["Panel"]] == "Fig. 7n")
    stats.insert_rows(at)
    values = {"Panel": "Fig. 7k", "Measure": "Mortality-clock output, PLX5622 minus control",
              "Comparison": "pooled across projects", "Test": "DerSimonian–Laird random effects",
              "n": f"{len(units)} projects ({len(cmp_k)} comparisons)", "Estimate": pooled[0],
              "95% CI": f"{pooled[1]:.2f} to {pooled[2]:.2f}", "Note": "estimate = difference in mean clock output"}
    for name, v in values.items():
        stats.cell(row=at, column=cols[name] + 1, value=v)
wb.save(XLSX)
print(json.dumps({"pooled": pooled[:3], "projects": len(units), "comparisons": len(cmp_k)}, indent=1))
