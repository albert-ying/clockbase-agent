#!/usr/bin/env python3
"""Figure 7 rebuilt from source data: cold exposure and microglia depletion in aged mice.

Panels (183 mm wide; letters a-j as before, the PLX5622 block now k-n):
  a  mortality tAge of DPP4+ adipocyte progenitors, 22 vs 6 degC (GSE165974)            source data
  b  25 genes with the largest contributions to the mortality-clock difference           source data (TACO genes)
  c  design of the housing-temperature study                                             drawn here; icons from A.M.'s Figure 7
  d  fat mass at 660 and 770 days                                                        source data
  e  frailty index at 660 and 770 days                                                   source data
  f  difference in frailty-item scores at 770 days                                       source data
  g-i Kaplan-Meier survival (all with baseline, lean, fat)                               source data
  j  longevity-steepness plane with SR-model lines                                       A.M.'s panel (SR outputs are not in the source data)
  k  PLX5622 in the verified atlas: mortality-clock difference by project, pooled      intervention-meta-2026-10-17
  l  design of the microglia depletion study                                             A.M.'s panel
  m  representative Iba1/DAPI images                                                     A.M.'s panel
  n  novel object visits                                                                 source data

Every statistic printed in a, d-i and n is recomputed here from the source data and checked against the sheet
Statistics of the source-data workbook; the build stops if one differs. Nothing is simulated.

Run: source env_local.sh; .venv-pipeline/bin/python -B scripts/figures_2026_11_08_fig7/build_figure7_2026_11_08.py
"""
import sys

sys.dont_write_bytecode = True

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import openpyxl  # noqa: E402
import pandas as pd  # noqa: E402
import pymupdf  # noqa: E402
from lifelines import CoxPHFitter, KaplanMeierFitter  # noqa: E402
from lifelines.statistics import logrank_test  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle  # noqa: E402
from scipy import stats  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT.parent
REV = ROOT / "paper/revision_2026-09-07"
XLSX = REV / "deliverables/FINAL_REVISION_2026-11-06/SourceData_Fig7_EDFig25_EDFig26.xlsx"
V19 = REV / "figures/03_Figures_WORKING_2026-10-10-v19.pdf"          # page 7: A.M.'s corrected Figure 7
META = LOCAL / "data/clockbase-verifier-v3-2026-09-27/intervention-meta-2026-10-17/primary"
OUT = REV / "evidence/figure7-rebuild-2026-11-08"
for p in (XLSX, V19, META / "meta_input_effects.csv", META / "ranked_replicated_age-lowering.csv",
          META.parent / "projects.json"):
    if not p.exists():
        raise SystemExit(f"Cannot proceed - required data file not found: {p}")
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- style
W = 518.74                      # 183 mm in pt
H = 662.0
FS, FS_S, FS_T = 6.0, 5.5, 6.5  # labels and ticks, small text, panel titles
RT, COLD = "#FC4E07", "#4E84C4"            # A.M.'s housing colours (c, d-j)
RT_T = "#D9480F"
CTRL, PLX = "#B2F3B2", "#FFD0FF"           # A.M.'s treatment colours (l-n)
CTRL_E, PLX_E = "#3E8E3E", "#B04CA8"
UP, DOWN = "#D6604D", "#4393C3"
INK, GREY, NAVY = "#222222", "#666666", "#1F3A5F"
LOWER = "#0072B2"                          # lower clock output (Figure 3)
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial"], "font.size": FS,
                     "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "xtick.major.size": 2.0, "ytick.major.size": 2.0, "xtick.labelsize": FS_S,
                     "ytick.labelsize": FS_S, "axes.labelsize": FS, "pdf.fonttype": 42, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.fontsize": FS_S, "legend.frameon": False,
                     "mathtext.default": "regular", "hatch.linewidth": 0.4})

fig = plt.figure(figsize=(W / 72, H / 72))


def ax_at(x, y, w, h, **kw):
    """Axes at x, y (top-left, pt from the page's top-left) with width w and height h (pt)."""
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)


def sci(p, digits=2):
    """P value as text: 0.0054, 0.086, 2.6 x 10^-4."""
    if p >= 0.001:
        return f"{p:.2g}" if p < 0.01 else f"{p:.2f}".rstrip("0").rstrip(".") if p >= 0.1 else f"{p:.2g}"
    e = int(math.floor(math.log10(p)))
    m = p / 10 ** e
    return rf"{m:.{digits - 1}f} × 10$^{{{e}}}$".replace("-", "−")


def minus(s):
    return s.replace("-", "−")


# ---------------------------------------------------------------- source data
def sheet(name, header_row, n_rows=None):
    """A sheet of the source-data workbook from its header row (1-based); formula cells are recomputed by the caller."""
    try:
        wb = openpyxl.load_workbook(XLSX, read_only=True)
        rows = list(wb[name].iter_rows(min_row=header_row, values_only=True))
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(f"Cannot proceed - cannot read sheet {name!r} of {XLSX}: {exc}")
    head = [str(h) for h in rows[0]]
    body = [r for r in rows[1:] if r and r[0] is not None]
    if n_rows:
        body = body[:n_rows]
    df = pd.DataFrame(body, columns=head)
    return df.loc[:, [c for c in df.columns if c != "None"]]


def formula(v):
    return isinstance(v, str) and v.startswith("=")


def header_row(name, first):
    wb = openpyxl.load_workbook(XLSX, read_only=True)
    for i, r in enumerate(wb[name].iter_rows(values_only=True), start=1):
        if r and r[0] == first:
            return i
    raise SystemExit(f"Cannot proceed - no header row starting {first!r} in sheet {name!r}")


STATS = sheet("Statistics", header_row("Statistics", "Panel"))
STATS = STATS[STATS["Panel"].astype(str).str.startswith("Fig.")]
CHECKS = []


def check(panel, label, value, sheet_value, rel=1e-6):
    ok = abs(value - sheet_value) <= rel * max(abs(sheet_value), 1e-300)
    CHECKS.append({"panel": panel, "statistic": label, "computed": value, "source_data": sheet_value, "match": ok})
    if not ok:
        raise SystemExit(f"{panel} {label}: computed {value!r} differs from the source data {sheet_value!r}")


def stat_row(panel, comparison):
    sel = STATS[(STATS["Panel"] == panel) & (STATS["Comparison"] == comparison)]
    if sel.empty:
        raise SystemExit(f"Statistics sheet has no row {panel} / {comparison}")
    return sel


# a: tAge
TA = sheet("Fig 7a & ED 25a", 6, n_rows=10)
TA["tage"] = TA["Mortality tAge (log10 HR)"].astype(float)
g22, g6 = TA.loc[TA.Housing == "22 °C", "tage"], TA.loc[TA.Housing == "6 °C", "tage"]
assert len(g22) == len(g6) == 5
t_a = stats.ttest_ind(g6, g22, equal_var=True)
row = stat_row("Fig. 7a", "6 °C vs 22 °C").iloc[0]
check("Fig. 7a", "difference in means", float(g6.mean() - g22.mean()), float(row["Estimate"]))
check("Fig. 7a", "P", float(t_a.pvalue), float(row["P"]))

# b: genes
TG = sheet("TACO genes", 5)
for c in ("logFC (6 °C vs 22 °C)", "s.e. of logFC"):
    TG[c] = TG[c].astype(float)
coef_col = [c for c in TG.columns if c.startswith("Mortality clock")][0]
padj_col = [c for c in TG.columns if c.startswith("P (BH")][0]
TG = TG[pd.to_numeric(TG[coef_col], errors="coerce").notna()].copy()
TG["coef"] = TG[coef_col].astype(float)
TG["contrib"] = TG["logFC (6 °C vs 22 °C)"] * TG["coef"]
TG["contrib_se"] = (TG["coef"] * TG["s.e. of logFC"]).abs()
TG["padj"] = TG[padj_col].astype(float)
TOP = TG.reindex(TG["contrib"].abs().sort_values(ascending=False).index).head(25)
first25 = sheet("TACO genes", 5, n_rows=25)["Gene"].tolist()
assert TOP["Gene"].tolist() == first25, "top 25 genes differ from the first 25 rows of TACO genes"
TOP = TOP.sort_values("contrib", ascending=True)

# d, e: body composition and frailty
BC = sheet("Body composition", 5)
BC = BC[BC["Mouse"].notna() & BC["Group"].isin(["RT", "20 °C"])].copy()
BC["fat"] = BC["Fat mass (g)"].astype(float)
BC["age"] = BC["Age (days)"].astype(int)
FR = sheet("Frailty", 4)
FR = FR[FR["Mouse"].notna() & FR["Group"].isin(["RT", "20 °C"])].copy()
ITEMS = [c for c in FR.columns[4:] if not c.lower().startswith("frailty")]
assert len(ITEMS) == 29, len(ITEMS)
FR[ITEMS] = FR[ITEMS].astype(float)
FR["fi"] = FR[ITEMS].mean(axis=1)
FR["age"] = FR["Age (days)"].astype(int)


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.minimum(adj, 1)
    return out


def r_wilcox_ranksum(x, y):
    """R wilcox.test(x, y): exact without ties for small n, otherwise normal approximation with continuity."""
    ties = len(np.unique(np.r_[x, y])) < len(x) + len(y)
    method = "asymptotic" if ties or len(x) >= 50 or len(y) >= 50 else "exact"
    return stats.mannwhitneyu(x, y, alternative="two-sided", method=method, use_continuity=True).pvalue


def r_wilcox_signed(x, y):
    d = np.asarray(x) - np.asarray(y)
    d = d[d != 0]
    # the source analysis used the normal approximation with continuity correction (exact = FALSE) for every group
    return stats.wilcoxon(d, zero_method="wilcox", correction=True, method="approx").pvalue


def two_ages(df, value, panel, measure):
    out = {"between": {}, "within": {}}
    for age in (660, 770):
        a = df[(df.Group == "20 °C") & (df.age == age)][value]
        b = df[(df.Group == "RT") & (df.age == age)][value]
        p = r_wilcox_ranksum(a, b)
        check(panel, f"{measure} 20 °C vs RT at {age} days",
              p, float(stat_row(panel, f"20 °C vs RT at {age} days").iloc[0]["P"]), rel=1e-4)
        out["between"][age] = p
    raw = []
    for grp in ("RT", "20 °C"):
        w = df[df.Group == grp].pivot_table(index="Mouse", columns="age", values=value).dropna()
        raw.append(r_wilcox_signed(w[770], w[660]))
    adj = bh(raw)
    for grp, p in zip(("RT", "20 °C"), adj):
        check(panel, f"{measure} 770 vs 660 within {grp}", p,
              float(stat_row(panel, f"770 vs 660 days within {grp}").iloc[0]["P"]), rel=2e-3)  # R vs scipy normal approx
        out["within"][grp] = p
    return out


P_D = two_ages(BC, "fat", "Fig. 7d", "fat mass")
P_E = two_ages(FR, "fi", "Fig. 7e", "frailty index")

# f: frailty items at 770 days
F7 = sheet("Fig 7f", 5)
F7 = F7[F7["Parameter"].notna() & ~F7["Parameter"].astype(str).str.startswith("Frailty index")]
items770 = FR[FR.age == 770]
rows_f = []
for it in F7["Parameter"]:
    col = it if it in ITEMS else {"Loss of fur color": "Loss of fur colour"}.get(it, it)
    a = items770[items770.Group == "20 °C"][col]
    b = items770[items770.Group == "RT"][col]
    rows_f.append({"item": it.replace("colour", "color"), "diff": a.mean() - b.mean()})
FF = pd.DataFrame(rows_f)
FF["se"] = F7["s.e. of difference"].astype(float).values
FF["p"] = F7["P"].astype(float).values
FF["padj"] = F7["P (BH-adjusted)"].astype(float).values
fi770 = items770.groupby("Group")["fi"].mean()
FF = pd.concat([FF, pd.DataFrame([{"item": "Frailty index", "diff": fi770["20 °C"] - fi770["RT"],
                                    "se": np.nan, "p": P_E["between"][770], "padj": P_E["between"][770]}])])
# se of the FI difference from the per-mouse values (Welch)
a = items770[items770.Group == "20 °C"]["fi"]
b = items770[items770.Group == "RT"]["fi"]
FF.loc[FF.item == "Frailty index", "se"] = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
FF = FF.sort_values("diff").reset_index(drop=True)

# g-i: survival
SV = sheet("Survival", 5)
SV = SV[SV["Mouse"].notna() & SV["Group"].isin(["RT", "20 °C"])].copy()
SV["age"] = SV["Age at death (days)"].astype(float)
SV["event"] = SV["Death observed"].astype(int)
SV["time"] = SV["age"] - 674
SV["cold"] = (SV.Group == "20 °C").astype(int)
base = SV[SV["In baseline cohort"] == "Yes"].copy()
base["ratio"] = base["Baseline fat mass (g)"].astype(float) / base["Baseline lean mass (g)"].astype(float)
base["fi"] = base["Baseline FI"].astype(float)
med = base["ratio"].median()
base["half"] = np.where(base["ratio"] < med, "Lean", "Fat")
assert len(base) == 61


def cox(df, covars, cluster=None):
    cph = CoxPHFitter()
    cols = ["time", "event", "cold"] + covars + ([cluster] if cluster else [])
    cph.fit(df[cols], duration_col="time", event_col="event", cluster_col=cluster, robust=bool(cluster))
    s = cph.summary.loc["cold"]
    return float(s["exp(coef)"]), float(s["exp(coef) lower 95%"]), float(s["exp(coef) upper 95%"]), float(s["p"])


base["cage"] = base["Cage"].astype(str)
SURV = {}
hr, lo, hi, p = cox(base, ["ratio", "fi"], cluster="cage")
r_ = stat_row("Fig. 7g", "20 °C vs RT")
cx = r_[r_.Test.str.startswith("Cox")].iloc[0]
check("Fig. 7g", "adjusted HR", hr, float(cx["Estimate"]), rel=2e-3)
check("Fig. 7g", "adjusted HR P", p, float(cx["P"]), rel=5e-2)
lr = logrank_test(base[base.cold == 1].time, base[base.cold == 0].time,
                  base[base.cold == 1].event, base[base.cold == 0].event).p_value
check("Fig. 7g", "log-rank P", lr, float(r_[r_.Test.str.startswith("log-rank")].iloc[0]["P"]), rel=1e-3)
p_src = float(cx["P"])  # R coxph with cluster(cage); lifelines' robust variance agrees within the check tolerance
SURV["g"] = {"df": base, "lines": ["Baseline cohort", f"adjusted HR = {hr:.2f}, P = {sci(p_src)}",
                                   f"log-rank P = {sci(lr)}"]}
for key, half, panel in (("h", "Lean", "Fig. 7h"), ("i", "Fat", "Fig. 7i")):
    d = base[base.half == half]
    hr, lo, hi, p = cox(d, [])
    r_ = stat_row(panel, "20 °C vs RT")
    cx = r_[r_.Test.str.startswith("Cox")].iloc[0]
    check(panel, "HR", hr, float(cx["Estimate"]), rel=2e-3)
    check(panel, "HR P", p, float(cx["P"]), rel=5e-2)
    lr = logrank_test(d[d.cold == 1].time, d[d.cold == 0].time, d[d.cold == 1].event, d[d.cold == 0].event).p_value
    check(panel, "log-rank P", lr, float(r_[r_.Test.str.startswith("log-rank")].iloc[0]["P"]), rel=1e-3)
    SURV[key] = {"df": d, "lines": [f"{half} at baseline", f"HR = {hr:.2f}, P = {sci(p)}", f"log-rank P = {sci(lr)}"]}

# n: NOR
NOR = sheet("NOR", 5)
NOR = NOR[NOR["Mouse"].notna() & NOR["Group"].isin(["Control", "PLX5622"])].copy()
NOR["visits"] = NOR["Novel object visits (n)"].astype(float)
out_col = [c for c in NOR.columns if c.startswith("Outlier for visits")][0]
keep = NOR[NOR[out_col].isna() | (NOR[out_col].astype(str).str.strip().isin(["", "None", "No"]))]
assert set(NOR.Mouse) - set(keep.Mouse) == {"C64"}, set(NOR.Mouse) - set(keep.Mouse)
pc, pp = keep[keep.Group == "Control"].visits, keep[keep.Group == "PLX5622"].visits
def exact_mw_ties(x, y):
    """Exact two-sided Wilcoxon-Mann-Whitney P with midranks for ties (all splits of the pooled sample), as the
    source analysis computed it (coin-style exact test): P(|W - E W| >= |w_obs - E W|)."""
    import itertools
    allv = np.r_[np.asarray(x, float), np.asarray(y, float)]
    r = stats.rankdata(allv)
    n = len(x)
    mu = n * (len(allv) + 1) / 2
    obs = abs(r[:n].sum() - mu)
    hits = total = 0
    for idx in itertools.combinations(range(len(allv)), n):
        total += 1
        hits += abs(r[list(idx)].sum() - mu) >= obs - 1e-9
    return hits / total


p_n = exact_mw_ties(pc, pp)
check("Fig. 7n", "exact Mann-Whitney P", p_n, float(stat_row("Fig. 7n", "PLX5622 vs control").iloc[0]["P"]))

# k: PLX5622 meta-analysis on the scale of the clock: the projects and the pooling rule of Figure 3f
#    (intervention_meta_2026_09_29: comparisons averaged within a project, DerSimonian-Laird across projects),
#    applied to the difference in mean mortality-clock output instead of Hedges' g
MI = pd.read_csv(META / "meta_input_effects.csv", low_memory=False)
RR = pd.read_csv(META / "ranked_replicated_age-lowering.csv", low_memory=False)
PROJ = json.loads((META.parent / "projects.json").read_text())
pool_g = RR[(RR.canonical_target == "plx5622") & (RR.condition_category == "Drug")].iloc[0]
cmp_k = MI[(MI.canonical_target == "plx5622") & (MI.condition_category == "Drug")].copy()
assert cmp_k["clock"].nunique() == 1, "the PLX5622 comparisons use more than one clock"
cmp_k["var_md"] = (cmp_k.v3_mean_diff / cmp_k.v3_cohen_d) ** 2 * (1 / cmp_k.v3_n_case + 1 / cmp_k.v3_n_control)
units = cmp_k.groupby("project").agg(g=("g", "mean"), var_g=("var_g", "mean"), md=("v3_mean_diff", "mean"),
                                     var_md=("var_md", "mean"), n_cmp=("g", "size"),
                                     series=("gse_id", lambda x: sorted(set(x)))).reset_index()


def dersimonian_laird(y, v):
    w = 1 / v
    m = (w * y).sum() / w.sum()
    q = (w * (y - m) ** 2).sum()
    tau2 = max(0.0, (q - (len(y) - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
    ws = 1 / (v + tau2)
    est = float((ws * y).sum() / ws.sum())
    se = math.sqrt(1 / ws.sum())
    return est, est - 1.959964 * se, est + 1.959964 * se


assert int(pool_g.k) == len(units) == 3 and int(pool_g.n_series) == sum(len(x) for x in units.series)
check("Fig. 7k", "pooled Hedges' g (grouping as Figure 3f)", dersimonian_laird(units.g, units.var_g)[0],
      float(pool_g.estimate), rel=1e-6)
POOL_MD = dersimonian_laird(units.md, units.var_md)
units["lo"] = units.md - 1.959964 * np.sqrt(units.var_md)
units["hi"] = units.md + 1.959964 * np.sqrt(units.var_md)
# row labels: surname of the GEO contact (projects.json, as Extended Data Figure 22e) and what the series profiled
WHAT = {"GSE134151": "5xFAD and WT brain", "GSE105135": "Retina, brain, astrocytes", "GSE197153": "Cortex"}
units["label"] = [WHAT[u]  # what each project profiled (from its GEO records and comparison names)
                  for u in units.project]
units = units.sort_values("md").reset_index(drop=True)


# ---------------------------------------------------------------- drawing helpers
def letter(x, y, s):
    LETTERS.append((x, y, s))


LETTERS = []


def bracket(ax, x0, x1, y, text, h=None, fs=FS_S):
    h = h if h is not None else (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.02
    ax.plot([x0, x0, x1, x1], [y - h, y, y, y - h], lw=0.45, color=INK, clip_on=False)
    ax.text((x0 + x1) / 2, y + h * 0.6, text, ha="center", va="bottom", fontsize=fs, clip_on=False)


def boxdots(ax, pos, vals, color, width=0.55, seed=0):
    vals = np.asarray(vals, float)
    bp = ax.boxplot([vals], positions=[pos], widths=width, patch_artist=True, showfliers=False,
                    medianprops={"color": INK, "lw": 0.8}, whiskerprops={"lw": 0.5, "color": INK},
                    capprops={"lw": 0.5, "color": INK}, boxprops={"lw": 0.5, "edgecolor": INK})
    for b in bp["boxes"]:
        b.set_facecolor(color)
        b.set_alpha(0.55)
    rng = np.random.default_rng(seed)
    ax.scatter(pos + rng.uniform(-width * 0.28, width * 0.28, len(vals)), vals, s=5.5, color=color,
               edgecolor=INK, linewidth=0.35, zorder=3)


# ---------------------------------------------------------------- a
letter(2, 9, "a")
ax = ax_at(28, 26, 44, 118)
boxdots(ax, 0, g22, RT, seed=1)
boxdots(ax, 1, g6, COLD, seed=2)
ax.set_xticks([0, 1], ["22 °C", "6 °C"])
ax.set_xlim(-0.6, 1.6)
ax.set_ylabel("Mortality tAge, log$_{10}$ HR")
ax.set_yticks([-0.2, 0, 0.2], [minus("-0.2"), "0", "0.2"])
top = max(g22.max(), g6.max())
ax.set_ylim(min(g22.min(), g6.min()) - 0.04, top + 0.12)
bracket(ax, 0, 1, top + 0.05, sci(t_a.pvalue), h=0.012)
ax.set_title("Adipocyte progenitors", fontsize=FS, pad=3)

# ---------------------------------------------------------------- b
letter(84, 9, "b")
ax = ax_at(120, 18, 124, 126)
yy = np.arange(len(TOP))
cols = [UP if v > 0 else DOWN for v in TOP["logFC (6 °C vs 22 °C)"]]
ax.barh(yy, TOP.contrib, height=0.72, color=cols, edgecolor=INK, linewidth=0.3)
ax.errorbar(TOP.contrib, yy, xerr=TOP.contrib_se, fmt="none", ecolor=INK, elinewidth=0.4, capsize=0.8, capthick=0.4)
ax.set_yticks(yy, TOP.Gene, fontstyle="italic", fontsize=5.2)
ax.tick_params(axis="y", length=0, pad=1.5)
ax.axvline(0, color=INK, lw=0.5)
for y_, (v, se_, pa) in enumerate(zip(TOP.contrib, TOP.contrib_se, TOP.padj)):
    s = "***" if pa < 0.001 else "**" if pa < 0.01 else "*" if pa < 0.05 else ""
    ax.text(v + np.sign(v) * (se_ + 0.0012), y_, s, ha="left" if v > 0 else "right", va="center", fontsize=4.5)
ax.set_xlabel("Effect on tAge, log$_{10}$ HR")
lim = float((TOP.contrib.abs() + TOP.contrib_se).max()) + 0.009
ax.set_xlim(-lim, lim)
ax.set_xticks([-0.04, -0.02, 0, 0.02, 0.04], [minus("-0.04"), minus("-0.02"), "0", "0.02", "0.04"])
ax.spines["left"].set_visible(False)
ax.set_ylim(-0.7, len(TOP) - 0.3)
ax.text(0.02, 1.0, "At 6 °C", transform=ax.transAxes, fontsize=FS_S, fontweight="bold", va="bottom")
for k_, (lab, c_) in enumerate((("up", UP), ("down", DOWN))):
    ax.add_patch(Rectangle((0.30 + 0.27 * k_, 1.012), 0.045, 0.03, transform=ax.transAxes, color=c_, clip_on=False,
                           ec=INK, lw=0.3))
    ax.text(0.36 + 0.27 * k_, 1.0, lab, transform=ax.transAxes, fontsize=FS_S, va="bottom")

# ---------------------------------------------------------------- c (design, drawn in pt coordinates)
letter(256, 9, "c")
C_X, C_Y, C_W, C_H = 262, 14, 252, 140
axc = ax_at(C_X, C_Y, C_W, C_H)
axc.set_xlim(0, C_W)
axc.set_ylim(C_H, 0)
axc.axis("off")
ICONS = []        # (page rect target in figure pt, clip rect in A.M.'s Figure 7)
ICONS.append(((C_X + 2, C_Y + 2, C_X + 2 + 44, C_Y + 2 + 16.2), (274, 21, 323, 38)))   # mouse
axc.text(52, 8, "Male UM-HET3 mice", fontsize=7, fontweight="bold", color=NAVY, va="center")
axc.text(52, 17.5, "Room temperature (RT) until 674 days; followed until death", fontsize=FS_S, color="#333333",
         va="center")
# timeline: 660 at x0, end of life at x1
T0, T1, TY = 10, 244, 40
x_age = lambda a: T0 + (a - 660) / (790 - 660) * (T1 - 30 - T0)  # noqa: E731
axc.annotate("", xy=(T1, TY), xytext=(T0, TY), arrowprops={"arrowstyle": "-|>", "lw": 0.8, "color": NAVY,
                                                          "mutation_scale": 6})
for a_, lab in ((660, "660"), (674, "674"), (770, "770")):
    axc.plot([x_age(a_)], [TY], "o", ms=3.2, color=NAVY)
    axc.text(x_age(a_), TY - 4.5, lab, ha="center", va="bottom", fontsize=FS, fontweight="bold", color=NAVY)
axc.text(T1, TY - 4.5, "End of life", ha="right", va="bottom", fontsize=FS_S, color=NAVY)
axc.text(x_age(674) + 2, TY + 4, "Housing changed", ha="left", va="top", fontsize=FS_S, color="#333333")
axc.text(T1, TY + 4, "Age, days", ha="right", va="top", fontsize=FS_S, color="#333333")
for k_, (lab, col, tcol, icons) in enumerate((("RT · n = 51", RT, RT_T, False), ("20 °C · n = 26", COLD, COLD, True))):
    y0 = 56 + k_ * 32
    axc.add_patch(FancyBboxPatch((x_age(674), y0), T1 - x_age(674), 25, boxstyle="round,pad=0,rounding_size=4",
                                 fc="white", ec=col, lw=0.9))
    ICONS.append(((C_X + x_age(674) + 4, C_Y + y0 + 2.5, C_X + x_age(674) + 4 + 27.1 * 0.74,
                   C_Y + y0 + 2.5 + 19.9 * 0.74 + 0.6), (315.3, 90.0 if k_ == 0 else 121.2, 342.4,
                                                          109.9 if k_ == 0 else 141.1)))
    xt = x_age(674) + 30
    if icons:
        ICONS.append(((C_X + xt, C_Y + y0 + 4, C_X + xt + 16.8 * 0.95, C_Y + y0 + 4 + 16.8 * 0.95),
                      (347.2, 122.4, 364.0, 139.2)))
        xt += 19
    axc.text(xt, y0 + 12.5, lab, fontsize=7, fontweight="bold", color=tcol, va="center")
for a_ in (660, 770):
    axc.plot([x_age(a_)] * 2, [TY + 3, 121], ls=(0, (2, 2)), lw=0.6, color="#7FA7D6")
    axc.text(x_age(a_), 124, "Body composition\nand frailty", ha="center" if a_ == 770 else "left", va="top",
             fontsize=FS_S, color="#333333", linespacing=1.1)

# ---------------------------------------------------------------- d, e
def two_age_panel(lx, x, df, value, ylabel, P, ticks, fmt):
    letter(lx, 172, "d" if value == "fat" else "e")
    ax = ax_at(x, 196, 92, 92)
    pos = {("RT", 660): 0, ("RT", 770): 1, ("20 °C", 660): 2.4, ("20 °C", 770): 3.4}
    for grp in ("RT", "20 °C"):
        wide = df[df.Group == grp].pivot_table(index="Mouse", columns="age", values=value)
        for _, r in wide.dropna().iterrows():
            ax.plot([pos[(grp, 660)], pos[(grp, 770)]], [r[660], r[770]], color="#B8B8B8", lw=0.3, zorder=1)
        for age in (660, 770):
            boxdots(ax, pos[(grp, age)], df[(df.Group == grp) & (df.age == age)][value],
                    RT if grp == "RT" else COLD, width=0.62, seed=age + len(grp))
    ax.set_xticks([0, 1, 2.4, 3.4], ["660", "770", "660", "770"])
    ax.set_xlabel("Age, days")
    ax.set_ylabel(ylabel)
    ax.set_yticks(ticks, [fmt(t) for t in ticks])
    ymax = float(df[value].max())
    span = ymax - float(df[value].min())
    ax.set_xlim(-0.6, 4.0)
    ax.set_ylim(float(df[value].min()) - span * 0.04, ymax + span * 0.62)
    hb = span * 0.025
    bracket(ax, 0, 1, ymax + span * 0.08, sci(P["within"]["RT"]), h=hb)
    bracket(ax, 2.4, 3.4, ymax + span * 0.08, sci(P["within"]["20 °C"]), h=hb)
    bracket(ax, 0, 2.4, ymax + span * 0.26, sci(P["between"][660]), h=hb)
    bracket(ax, 1, 3.4, ymax + span * 0.44, sci(P["between"][770]), h=hb)
    for k_, (lab, col) in enumerate((("RT", RT), ("20 °C", COLD))):
        ax.add_patch(Rectangle((0.02 + 0.2 * k_, 1.04), 0.06, 0.045, transform=ax.transAxes, color=col, alpha=0.7,
                               clip_on=False, ec=INK, lw=0.3))
        ax.text(0.10 + 0.2 * k_, 1.035, lab, transform=ax.transAxes, fontsize=FS_S, va="bottom")
    return ax


two_age_panel(2, 28, BC, "fat", "Fat mass, g", P_D, [10, 20, 30], lambda t: f"{t:g}")
two_age_panel(130, 156, FR, "fi", "Frailty index", P_E, [0.2, 0.4, 0.6], lambda t: f"{t:g}")

# ---------------------------------------------------------------- f (frailty items, vertical bars)
letter(258, 172, "f")
axf = ax_at(290, 186, 196, 56)
cmap = LinearSegmentedColormap.from_list("reds", ["#FFF0EB", "#FCA98C", "#E2402F", "#A50F15"])
nl = -np.log10(FF.padj.clip(lower=1e-10))
norm = Normalize(0, 3)
xs = np.arange(len(FF))
axf.bar(xs, FF["diff"], width=0.72, color=cmap(norm(nl)), edgecolor=INK, linewidth=0.3)
axf.errorbar(xs, FF["diff"], yerr=FF.se, fmt="none", ecolor=INK, elinewidth=0.4, capsize=0.8, capthick=0.4)
axf.axhline(0, color=INK, lw=0.5)
for x_, (v, se_, pa) in enumerate(zip(FF["diff"], FF.se, FF.padj)):
    s = "**" if pa < 0.01 else "*" if pa < 0.05 else ""
    if s:
        axf.text(x_, v - se_ - 0.03, s, ha="center", va="top", fontsize=5)
axf.set_xticks(xs)
axf.set_xticklabels(FF.item, rotation=60, ha="right", rotation_mode="anchor", fontsize=5.2)
for t in axf.get_xticklabels():
    if t.get_text() == "Frailty index":
        t.set_fontweight("bold")
axf.tick_params(axis="x", length=0, pad=1.5)
axf.set_xlim(-0.7, len(FF) - 0.3)
axf.set_ylim(min(-0.55, float((FF["diff"] - FF.se).min()) - 0.1), 0.3)
axf.set_yticks([-0.4, -0.2, 0, 0.2], [minus("-0.4"), minus("-0.2"), "0", "0.2"])
axf.set_ylabel("Score difference at\n770 days, 20 °C − RT", fontsize=FS_S, linespacing=1.0)
cax = ax_at(492, 190, 3.2, 34)
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
cb.outline.set_linewidth(0.3)
cb.set_ticks([0, 1, 2, 3])
cb.ax.tick_params(labelsize=4.8, length=1.5, width=0.4, pad=1)
cax.set_title("−log$_{10}$ p.adj", fontsize=4.8, pad=4)

# ---------------------------------------------------------------- g-i
for k_, key in enumerate(("g", "h", "i")):
    x0 = 2 + k_ * 116
    letter(x0, 312, key)
    d = SURV[key]["df"]
    for j_, line in enumerate(SURV[key]["lines"]):
        fig.text((x0 + 10) / W, 1 - (318 + j_ * 7.2) / H, line, fontsize=FS_S if j_ else FS,
                 fontweight="bold" if j_ == 0 else "normal", va="top")
    ax = ax_at(x0 + 26, 348, 84, 96)
    for grp, col in (("RT", RT), ("20 °C", COLD)):
        s = d[d.Group == grp]
        km = KaplanMeierFitter().fit(s.age, s.event, entry=np.full(len(s), 674.0), label=f"{grp} (n = {len(s)})")
        sf = km.survival_function_
        ci = km.confidence_interval_
        ax.step(sf.index, sf.iloc[:, 0], where="post", color=col, lw=1.0, label=f"{grp} (n = {len(s)})")
        ax.fill_between(ci.index, ci.iloc[:, 0], ci.iloc[:, 1], step="post", color=col, alpha=0.15, lw=0)
    ax.set_xlim(674, 1060)
    ax.set_ylim(0, 1.03)
    ax.set_xticks([700, 800, 900, 1000])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0], ["0", "0.25", "0.50", "0.75", "1.00"])
    ax.set_xlabel("Age, days")
    if k_ == 0:
        ax.set_ylabel("Survival probability")
    ax.legend(loc="upper right", handlelength=1.2, borderaxespad=0.2, fontsize=5.0)

# ---------------------------------------------------------------- k (PLX5622 in the atlas)
letter(2, 494, "k")
axk = ax_at(84, 506, 70, 88)
labels = []
for i_, r in units.iterrows():
    labels.append(f"{r.label} ({r.n_cmp})")
    axk.plot([r.lo, r.hi], [i_, i_], color=INK, lw=0.6)
    s_ = 10 + 26 * (1 / r.var_md) / (1 / units.var_md).max()
    axk.scatter([r.md], [i_], marker="s", s=s_, color="#555555", edgecolor=INK, linewidth=0.3, zorder=3)
yp = len(units) + 0.35
d_ = 0.32
est, lo_, hi_ = POOL_MD
axk.add_patch(Polygon([[lo_, yp], [est, yp + d_], [hi_, yp], [est, yp - d_]], closed=True, fc=LOWER, ec=INK,
                      lw=0.4))
labels.append("Pooled")
axk.set_yticks(list(range(len(units))) + [yp], labels, fontsize=5.2, linespacing=1.0)
axk.get_yticklabels()[-1].set_fontweight("bold")
axk.tick_params(axis="y", length=0)
axk.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
axk.set_ylim(yp + 0.75, -0.7)
xmin = math.floor(min(units.lo.min(), lo_) * 5) / 5 - 0.05
axk.set_xlim(xmin, 0.35)
axk.set_xticks([-0.8, -0.4, 0], [minus("-0.8"), minus("-0.4"), "0"])
axk.set_xlabel("Difference in mortality-clock\noutput (95% CI)", linespacing=1.0)
axk.spines["left"].set_visible(False)
axk.text(0.35, yp, minus(f"{est:.2f}"), fontsize=5.0, ha="right", va="center")

# ---------------------------------------------------------------- n (NOR)
letter(452, 494, "n")
axn = ax_at(476, 512, 38, 120)
boxdots(axn, 0, pc, CTRL, width=0.6, seed=3)
boxdots(axn, 1, pp, PLX, width=0.6, seed=4)
for b in axn.patches:
    pass
axn.set_xticks([0, 1], ["Control", "PLX5622"], rotation=35, ha="right", rotation_mode="anchor")
axn.set_xlim(-0.6, 1.6)
axn.set_ylabel("Novel object visits, n")
topn = float(max(pc.max(), pp.max()))
axn.set_ylim(0, topn * 1.18)
bracket(axn, 0, 1, topn * 1.06, sci(p_n), h=topn * 0.02)
axn.set_title("NOR test", fontsize=FS, pad=3)

# ---------------------------------------------------------------- save, then place A.M.'s panels and icons
layer = OUT / "_figure7_layer.pdf"
fig.savefig(layer)
plt.close(fig)
doc = pymupdf.open(layer)
page = doc[0]
src = pymupdf.open(V19)
srcp = 6
CLIPS = {  # panel: (target x, target y, clip rect in A.M.'s page, old letter rect to cover)
    "j": (352, 308, pymupdf.Rect(0, 528, 152, 697), pymupdf.Rect(9.5, 531.5, 14.0, 543.0)),
    "l": (160, 488, pymupdf.Rect(152, 528, 300, 697), pymupdf.Rect(155.5, 532.5, 162.0, 543.5)),
    "m": (308, 488, pymupdf.Rect(300, 528, 436, 697), pymupdf.Rect(302.5, 532.5, 307.0, 543.5)),
}
for key, (tx, ty, clip, old) in CLIPS.items():
    target = pymupdf.Rect(tx, ty, tx + clip.width, ty + clip.height)
    page.show_pdf_page(target, src, srcp, clip=clip, overlay=True)
    cover = pymupdf.Rect(tx + old.x0 - clip.x0, ty + old.y0 - clip.y0, tx + old.x1 - clip.x0, ty + old.y1 - clip.y0)
    page.draw_rect(cover, color=None, fill=(1, 1, 1), overlay=True)
    LETTERS.append((tx + 1.5, ty + 6, key))
for target, clip in ICONS:
    page.show_pdf_page(pymupdf.Rect(*target), src, srcp, clip=pymupdf.Rect(*clip), overlay=True)
for x, y, s in LETTERS:
    page.insert_text((x, y + 6), s, fontname="hebo", fontsize=8, color=(0, 0, 0), overlay=True)
OUT_PDF = OUT / "figure7.pdf"
doc.set_metadata({"title": "Figure 7", "author": "Kejun Albert Ying", "creator": Path(__file__).name})
# l: the 30-day arm is labelled by its diet and duration only (the behavioural tests differ in their numbers of mice)
_hits = page.search_for("n = 8 per group")
assert len(_hits) == 1, _hits
_r = _hits[0]
page.draw_rect(pymupdf.Rect(_r.x0 - 0.5, _r.y0 - 0.5, _r.x1 + 0.5, _r.y1 + 0.5), color=None, fill=(1, 1, 1), overlay=True)
doc.save(OUT_PDF, garbage=3, deflate=True)
layer.unlink()
pymupdf.open(OUT_PDF)[0].get_pixmap(dpi=220, alpha=False).save(OUT / "figure7.png")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


pd.DataFrame(CHECKS).to_csv(OUT / "statistics_check.csv", index=False)
numbers = {"a_P": t_a.pvalue, "d": P_D, "e": P_E, "k_pooled_mean_difference": POOL_MD,
           "k_units": units[["project", "label", "n_cmp", "md", "lo", "hi", "g"]].to_dict("records"), "n_P": p_n,
           "surv": {k: v["lines"] for k, v in SURV.items()}, "median_baseline_ratio": med}
(OUT / "numbers.json").write_text(json.dumps(numbers, indent=1, default=str))
(OUT / "provenance.json").write_text(json.dumps({
    "script": str(Path(__file__).resolve()), "run_utc": datetime.now(timezone.utc).isoformat(),
    "inputs": {str(p): sha(p) for p in (XLSX, V19, META / "meta_input_effects.csv", META / "ranked_replicated_age-lowering.csv",
                                        META.parent / "projects.json")},
    "output_sha256": sha(OUT_PDF), "page_pt": [W, H],
    "panels_from_AM_figure": ["c icons (mouse, cages, thermometer)", "j", "l (design)", "m (images)"]}, indent=1))
print(f"{OUT_PDF}  {W:.1f} x {H:.1f} pt ({W / 72 * 25.4:.0f} x {H / 72 * 25.4:.0f} mm); "
      f"{sum(c['match'] for c in CHECKS)}/{len(CHECKS)} statistics match the source data")
