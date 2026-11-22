#!/usr/bin/env python3
"""
Cross-species validation: compare our mouse significant Genetic interventions
against GenAge Models annotations from OTHER species (C. elegans, Drosophila,
yeast) in addition to mouse.

For each scored mouse Genetic intervention with a GenAge-matched gene symbol:
  - resolve ortholog symbols in worm, fly, yeast via Alliance Genome Resources
  - check if the ortholog symbol appears in the per-species GenAge Models table
  - if so, compare its `longevity influence` to our `inferred_influence`

Reports:
  - result/cross_species_gene_comparison.csv  (per-gene per-species)
  - result/ortholog_cache.json                (cached Alliance lookups)
  - figures/cross_species_gene_waffle.{pdf,png}  (4-species summary waffle)

Ortholog resolution:
  MyGene.info (mouse symbol -> MGI ID) -> Alliance Genome Resources orthologs.
  Alliance returns symbols that match GenAge Models species-specific symbols
  directly (sir-2.1, Sir2, etc.).
"""

import json
import re
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from matplotlib.patches import Patch
from scipy.stats import binomtest

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.facecolor"] = "white"

INPUT_AGENT = "data/agent_result_genage.csv"
INPUT_DRUG_AGENT = "data/agent_result_drugage_significant.csv"
GENAGE_MODELS = "data/refdata/genage_models.csv"
DRUGAGE_CSV = "data/refdata/drugage.csv"
ORTHOLOG_CACHE = "result/ortholog_cache.json"
CURATION_FILE = "result/curation_decisions.json"


def _load_curation_drops():
    """Load drugs and genes manually flagged for removal in curation_decisions.json."""
    p = Path(CURATION_FILE)
    if not p.exists():
        return set(), set()
    j = json.loads(p.read_text())
    drug_drops = {x["drug"].lower().strip() for x in j.get("drugs", [])
                  if x.get("decision") == "DROP"}
    gene_drops = {x["gene"].upper().strip() for x in j.get("genes", [])
                  if x.get("decision") == "DROP"}
    return drug_drops, gene_drops

SPECIES_ORDER = [
    "Mus musculus",
    "Caenorhabditis elegans",
    "Drosophila melanogaster",
    "Saccharomyces cerevisiae",
]
SPECIES_SHORT = {
    "Mus musculus":             "Mouse",
    "Caenorhabditis elegans":   "C. elegans",
    "Drosophila melanogaster":  "D. melanogaster",
    "Saccharomyces cerevisiae": "S. cerevisiae",
}

# Alliance uses strain-qualified names for yeast; normalize
ALLIANCE_TAXON_MAP = {
    "Caenorhabditis elegans":         "Caenorhabditis elegans",
    "Drosophila melanogaster":        "Drosophila melanogaster",
    "Saccharomyces cerevisiae":       "Saccharomyces cerevisiae",
    "Saccharomyces cerevisiae S288C": "Saccharomyces cerevisiae",
}

COLOR_AGREE      = "#2E86AB"
COLOR_DISAGREE   = "#E63946"
COLOR_AMBIGUOUS  = "#D9C67A"   # in GenAge but non-directional (e.g. "Necessary for fitness")
COLOR_NOT_IN_DB  = "#B5B5B5"

DIRECTIONAL = {"Pro-Longevity", "Anti-Longevity"}


# ---------- ortholog resolution ----------
def load_cache():
    p = Path(ORTHOLOG_CACHE)
    if p.exists():
        return json.loads(p.read_text())
    return {}

def save_cache(cache):
    Path(ORTHOLOG_CACHE).parent.mkdir(exist_ok=True)
    Path(ORTHOLOG_CACHE).write_text(json.dumps(cache, indent=2, sort_keys=True))

def resolve_mgi(symbol):
    """Mouse symbol -> MGI ID via MyGene pantherdb field."""
    try:
        r = requests.get(
            "https://mygene.info/v3/query",
            params={"q": f"symbol:{symbol}", "species": "10090", "fields": "pantherdb"},
            timeout=20,
        )
        j = r.json()
        for hit in j.get("hits", []):
            p = hit.get("pantherdb", {}) or {}
            mgi = p.get("MGI")
            if mgi:
                return mgi if mgi.startswith("MGI:") else f"MGI:{mgi}"
    except Exception as e:
        print(f"    ! mygene lookup failed for {symbol}: {e}")
    return None

def alliance_orthologs(mgi_id):
    """Return {taxon_name: [symbol, ...]} from Alliance."""
    out = {}
    try:
        r = requests.get(
            f"https://www.alliancegenome.org/api/gene/{mgi_id}/orthologs",
            params={"limit": 100},
            timeout=30,
        )
        j = r.json()
        for res in j.get("results", []):
            obj = res["geneToGeneOrthologyGenerated"]["objectGene"]
            taxon_raw = obj["taxon"]["name"]
            sym = obj["geneSymbol"]["displayText"]
            taxon = ALLIANCE_TAXON_MAP.get(taxon_raw, taxon_raw)
            out.setdefault(taxon, []).append(sym)
    except Exception as e:
        print(f"    ! alliance lookup failed for {mgi_id}: {e}")
    return out

def get_orthologs(symbol, cache):
    """Return {species: [ortho_symbols]} using cache when possible."""
    key = symbol.upper()
    if key in cache:
        return cache[key]
    mgi = resolve_mgi(symbol)
    time.sleep(0.1)
    if not mgi:
        cache[key] = {"_mgi": None}
        return cache[key]
    orth = alliance_orthologs(mgi)
    orth["_mgi"] = mgi
    cache[key] = orth
    time.sleep(0.15)
    return orth


# ---------- GenAge lookup ----------
def build_species_lookup(df: pd.DataFrame) -> dict:
    """{organism: {SYMBOL_UPPER: longevity_influence}}"""
    out = {}
    for org, g in df.groupby("organism"):
        m = dict(zip(g["symbol"].astype(str).str.upper(),
                     g["longevity influence"].astype(str)))
        out[org] = m
    return out


# ---------- DrugAge lookup ----------
def build_drug_species_lookup(df: pd.DataFrame) -> dict:
    """{organism: {DRUG_NORMALIZED: "Pro-Longevity" | "Anti-Longevity" | "Ambiguous"}}

    A drug's direction per species is sign(mean(avg_lifespan_change_percent)) across
    all experiments — matches DrugAge convention used by create_drugage_waffle_plot_nature.py
    (which uses the curated `drugage_annotation` column, equivalent here to mean sign).
    """
    out = {}
    for org, g in df.groupby("species"):
        lookup = {}
        for drug, gg in g.groupby(g["compound_name"].astype(str).str.lower().str.strip()):
            vals = pd.to_numeric(gg["avg_lifespan_change_percent"], errors="coerce").dropna()
            if len(vals) == 0:
                lookup[drug] = "Ambiguous"
                continue
            m = vals.mean()
            if m > 0:
                lookup[drug] = "Pro-Longevity"
            elif m < 0:
                lookup[drug] = "Anti-Longevity"
            else:
                lookup[drug] = "Ambiguous"
        out[org] = lookup
    return out


def binom_p(agree, disagree):
    """Two-sided binomial p for agree vs disagree under null 0.5."""
    n = agree + disagree
    if n == 0:
        return float("nan")
    return binomtest(agree, n, p=0.5, alternative="two-sided").pvalue


# ---------- waffle renderer (shared) ----------
def render_waffle(comp, title, out_pdf, out_png, label_noun="gene"):
    """comp must have columns: <species>__status for each SPECIES_ORDER entry."""
    fig, axes = plt.subplots(1, len(SPECIES_ORDER), figsize=(5 * len(SPECIES_ORDER), 7))
    if len(SPECIES_ORDER) == 1:
        axes = [axes]
    n_total = len(comp)
    for ax, sp in zip(axes, SPECIES_ORDER):
        counts = comp[f"{sp}__status"].value_counts()
        a = int(counts.get("agree", 0))
        d = int(counts.get("disagree", 0))
        amb = int(counts.get("ambiguous", 0))
        notin = n_total - a - d - amb
        values = [a, d, amb, notin]
        colors = [COLOR_AGREE, COLOR_DISAGREE, COLOR_AMBIGUOUS, COLOR_NOT_IN_DB]

        grid_rows, grid_cols = 20, 5
        total_cells = grid_rows * grid_cols
        scaled = [int(round(v * total_cells / max(1, n_total))) for v in values]
        while sum(scaled) < total_cells: scaled[scaled.index(max(scaled))] += 1
        while sum(scaled) > total_cells: scaled[scaled.index(max(scaled))] -= 1

        cell_colors = []
        for v, c in zip(scaled, colors):
            cell_colors.extend([c] * v)
        for i in range(total_cells):
            r = i // grid_cols
            col = i % grid_cols
            y = grid_rows - 1 - r
            ax.add_patch(plt.Rectangle((col, y), 0.92, 0.92, color=cell_colors[i]))
        ax.set_xlim(-0.1, grid_cols + 0.1)
        ax.set_ylim(-0.1, grid_rows + 0.1)
        ax.set_aspect("equal")
        ax.axis("off")
        rate = 100 * a / max(1, a + d)
        p = binom_p(a, d)
        p_txt = ("n/a" if np.isnan(p) else
                 f"p = {p:.2g}" if p >= 0.001 else f"p = {p:.1e}")
        ax.set_title(
            f"{SPECIES_SHORT[sp]}\n"
            f"{a+d}/{n_total} directional  |  agreement {rate:.0f}%\n"
            f"{p_txt} (two-sided binomial)",
            fontsize=11, fontweight="bold", pad=10,
        )

    db_name = "GenAge" if label_noun == "gene" else "DrugAge"
    handles = [
        Patch(facecolor=COLOR_AGREE,     label=f"Agree with {db_name}"),
        Patch(facecolor=COLOR_DISAGREE,  label=f"Disagree with {db_name}"),
        Patch(facecolor=COLOR_AMBIGUOUS, label=f"Non-directional in {db_name}"),
        Patch(facecolor=COLOR_NOT_IN_DB, label=f"No ortholog in {db_name}" if label_noun == "gene"
                                               else f"Not in {db_name}"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4,
               bbox_to_anchor=(0.5, 0.02), frameon=False, fontsize=10)
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0.04, 1, 0.94])
    Path("figures").mkdir(exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ {out_pdf}")
    print(f"✓ {out_png}")


# ---------- gene analysis ----------
def run_gene_analysis():
    print(f"\n=== GENE (GenAge Models) cross-species ===")
    print(f"loading mouse interventions: {INPUT_AGENT}")
    agent = pd.read_csv(INPUT_AGENT)
    agent = agent.drop_duplicates()
    agent = agent[agent["trust_score"] > 80]
    # Strict filter: only interventions whose gene IS in mouse GenAge — matches
    # create_geneage_waffle_plot_nature.py
    agent = agent[agent["genage_annotated"] == True]
    agent = agent.dropna(subset=["genage_matched_genes", "inferred_influence"])
    agent = agent[agent["genage_matched_genes"].astype(str).str.strip() != ""]
    agent["gene_u"] = agent["genage_matched_genes"].astype(str).str.upper()

    single_rows = []
    for _, r in agent.iterrows():
        for g in str(r["gene_u"]).split(";"):
            g = g.strip()
            if g:
                single_rows.append({"gene": g, "inferred_influence": r["inferred_influence"]})
    single = pd.DataFrame(single_rows)
    print(f"  {len(agent)} interventions -> {len(single)} single-gene rows "
          f"({single['gene'].nunique()} unique genes)")

    gene_rows = []
    for gene, g in single.groupby("gene"):
        votes = g["inferred_influence"].value_counts()
        gene_rows.append({
            "gene": gene,
            "inferred_influence": votes.idxmax(),
            "n_interventions": len(g),
        })
    genes_df = pd.DataFrame(gene_rows)

    drug_drops, gene_drops = _load_curation_drops()
    if gene_drops:
        before = len(genes_df)
        genes_df = genes_df[~genes_df["gene"].isin(gene_drops)].reset_index(drop=True)
        print(f"  curation: dropped {before-len(genes_df)} genes ({sorted(gene_drops)})")

    print(f"\nloading GenAge Models: {GENAGE_MODELS}")
    models = pd.read_csv(GENAGE_MODELS)
    species_lookup = build_species_lookup(models)
    for sp in SPECIES_ORDER:
        n = len(species_lookup.get(sp, {}))
        print(f"  {sp:<28} {n:>5} GenAge entries")

    print("\nresolving orthologs via Alliance (cached)...")
    cache = load_cache()
    rows = []
    for _, g in genes_df.iterrows():
        gene = g["gene"]
        inferred = g["inferred_influence"]
        orth = get_orthologs(gene, cache)
        row = {"gene": gene, "inferred_influence": inferred,
               "n_interventions": g["n_interventions"], "mgi_id": orth.get("_mgi")}
        for sp in SPECIES_ORDER:
            lookup = species_lookup.get(sp, {})
            candidates = [gene] if sp == "Mus musculus" else [s.upper() for s in orth.get(sp, [])]
            ann, match_sym = None, None
            for c in candidates:
                if c in lookup:
                    ann, match_sym = lookup[c], c
                    break
            row[f"{sp}__ortholog"]   = ";".join(orth.get(sp, [])) if sp != "Mus musculus" else gene
            row[f"{sp}__match_sym"]  = match_sym or ""
            row[f"{sp}__annotation"] = ann or ""
            if not ann:
                row[f"{sp}__status"] = "not in GenAge"
            elif ann not in DIRECTIONAL or inferred not in DIRECTIONAL:
                row[f"{sp}__status"] = "ambiguous"
            elif ann == inferred:
                row[f"{sp}__status"] = "agree"
            else:
                row[f"{sp}__status"] = "disagree"
        rows.append(row)
    save_cache(cache)
    comp = pd.DataFrame(rows)

    Path("result").mkdir(exist_ok=True)
    comp.to_csv("result/cross_species_gene_comparison.csv", index=False)

    print("\nsummary (GENE):")
    for sp in SPECIES_ORDER:
        counts = comp[f"{sp}__status"].value_counts()
        a = counts.get("agree", 0); d = counts.get("disagree", 0); amb = counts.get("ambiguous", 0)
        p = binom_p(a, d)
        p_txt = "n/a" if np.isnan(p) else f"{p:.3g}"
        print(f"  {SPECIES_SHORT[sp]:<18}  directional: {a+d:>3}/{len(comp)}  "
              f"agree: {a:>3}  disagree: {d:>3}  ambiguous: {amb:>3}  "
              f"agreement: {100*a/max(1,a+d):5.1f}%  p={p_txt}")

    render_waffle(
        comp,
        f"Cross-species validation of mouse genetic interventions "
        f"({len(comp)} unique genes; orthologs via Alliance Genome Resources)",
        "figures/cross_species_gene_waffle.pdf",
        "figures/cross_species_gene_waffle.png",
        label_noun="gene",
    )


# ---------- drug analysis ----------
def _norm(s):
    s = str(s).lower().strip()
    s = re.sub(r"[^a-z0-9 ]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def match_intervention_to_drugage(iv_norm, da_compounds):
    """Return list of DrugAge compound names matched in this intervention text.

    Exact normalized match preferred; otherwise word-boundary substring (compound
    appears as a whole token within the intervention). Length >= 4 to avoid junk.
    """
    if iv_norm in da_compounds:
        return [iv_norm]
    iv_padded = " " + iv_norm + " "
    hits = []
    for c in da_compounds:
        if len(c) < 4:
            continue
        if (" " + c + " ") in iv_padded:
            hits.append(c)
    return hits


def run_drug_analysis():
    print(f"\n\n=== DRUG (DrugAge) cross-species ===")
    print(f"loading mouse drug interventions: {INPUT_DRUG_AGENT}")
    agent = pd.read_csv(INPUT_DRUG_AGENT)
    agent = agent.drop_duplicates()
    agent = agent[agent["trust_score"] > 80]
    agent = agent.dropna(subset=["intervention", "inferred_influence"])
    agent["intervention_norm"] = agent["intervention"].apply(_norm)

    # Build the DrugAge 4-species compound set (used to extract drug names from
    # intervention text — broader than mouse-only DrugAge match)
    da = pd.read_csv(DRUGAGE_CSV)
    da["species"] = da["species"].astype(str).str.strip().str.title().replace({
        "Mus Musculus": "Mus musculus",
        "Caenorhabditis Elegans": "Caenorhabditis elegans",
        "Drosophila Melanogaster": "Drosophila melanogaster",
        "Saccharomyces Cerevisiae": "Saccharomyces cerevisiae",
    })
    da["compound_norm"] = da["compound_name"].apply(_norm)
    da4 = da[da["species"].isin(SPECIES_ORDER)]
    da_compounds = set(da4["compound_norm"].dropna().unique())
    print(f"  DrugAge 4-species compound set: {len(da_compounds)}")

    # Two-tier drug discovery:
    #   curator: rows tagged by the agent annotator in `drugage_matched_drugs` —
    #            handles combos like "Angiotensin II + Quercetin" correctly.
    #   name-match: rows where the drug ≈ entire intervention text (extends to
    #               drugs the curator missed because they aren't in mouse DrugAge
    #               but are in worm/fly/yeast).
    #
    # Per-drug inferred direction: if any curator-tagged rows exist, vote among
    # those only (matches create_drugage_waffle_plot_nature.py for the original
    # 5 drugs); else fall back to name-match rows.
    curator_rows = []
    namematch_rows = []
    for _, r in agent.iterrows():
        # Curator-tagged: drugage_annotated=True (drug present in mouse DrugAge,
        # vetted by annotator). This matches the strict filter used by
        # create_drugage_waffle_plot_nature.py.
        if bool(r.get("drugage_annotated", False)):
            ct = r.get("drugage_matched_drugs", "")
            if isinstance(ct, str) and ct.strip():
                c = _norm(ct)
                if c in da_compounds:
                    curator_rows.append({"drug": c, "inferred_influence": r["inferred_influence"]})
                    continue
        iv = r["intervention_norm"]
        for c in match_intervention_to_drugage(iv, da_compounds):
            if len(c) / max(1, len(iv)) >= 0.8:
                namematch_rows.append({"drug": c, "inferred_influence": r["inferred_influence"]})
    cur = pd.DataFrame(curator_rows)
    nm  = pd.DataFrame(namematch_rows)

    drug_rows = []
    cur_drugs = set(cur["drug"].unique()) if not cur.empty else set()
    if not cur.empty:
        for drug, g in cur.groupby("drug"):
            votes = g["inferred_influence"].value_counts()
            drug_rows.append({
                "drug": drug,
                "inferred_influence": votes.idxmax(),
                "n_interventions": len(g),
                "source": "curator",
            })
    if not nm.empty:
        for drug, g in nm.groupby("drug"):
            if drug in cur_drugs:
                continue
            votes = g["inferred_influence"].value_counts()
            drug_rows.append({
                "drug": drug,
                "inferred_influence": votes.idxmax(),
                "n_interventions": len(g),
                "source": "name-match",
            })
    drugs_df = pd.DataFrame(drug_rows).sort_values("drug").reset_index(drop=True)
    print(f"  {len(cur_drugs)} curator-tagged + "
          f"{len(drugs_df)-len(cur_drugs)} name-match-only = "
          f"{len(drugs_df)} unique drugs (pre-curation)")

    drug_drops, _ = _load_curation_drops()
    if drug_drops:
        before = len(drugs_df)
        drugs_df = drugs_df[~drugs_df["drug"].isin(drug_drops)].reset_index(drop=True)
        print(f"  curation: dropped {before-len(drugs_df)} drugs ({sorted(drug_drops)})")

    drug_lookup = build_drug_species_lookup(da4)
    for sp in SPECIES_ORDER:
        n = len(drug_lookup.get(sp, {}))
        print(f"  {sp:<28} {n:>5} unique DrugAge compounds")

    rows = []
    for _, dr in drugs_df.iterrows():
        drug = dr["drug"]
        inferred = dr["inferred_influence"]
        row = {"drug": drug, "inferred_influence": inferred,
               "n_interventions": dr["n_interventions"]}
        for sp in SPECIES_ORDER:
            lookup = drug_lookup.get(sp, {})
            ann = lookup.get(drug, None)
            row[f"{sp}__annotation"] = ann or ""
            if ann is None:
                row[f"{sp}__status"] = "not in GenAge"   # reuse bucket name (= not in DrugAge)
            elif ann == "Ambiguous" or inferred not in DIRECTIONAL:
                row[f"{sp}__status"] = "ambiguous"
            elif ann == inferred:
                row[f"{sp}__status"] = "agree"
            else:
                row[f"{sp}__status"] = "disagree"
        rows.append(row)
    comp = pd.DataFrame(rows)
    comp.to_csv("result/cross_species_drug_comparison.csv", index=False)

    print("\nsummary (DRUG):")
    for sp in SPECIES_ORDER:
        counts = comp[f"{sp}__status"].value_counts()
        a = counts.get("agree", 0); d = counts.get("disagree", 0); amb = counts.get("ambiguous", 0)
        p = binom_p(a, d)
        p_txt = "n/a" if np.isnan(p) else f"{p:.3g}"
        print(f"  {SPECIES_SHORT[sp]:<18}  directional: {a+d:>3}/{len(comp)}  "
              f"agree: {a:>3}  disagree: {d:>3}  ambiguous: {amb:>3}  "
              f"agreement: {100*a/max(1,a+d):5.1f}%  p={p_txt}")

    render_waffle(
        comp,
        f"Cross-species validation of mouse drug interventions "
        f"({len(comp)} unique compounds; DrugAge direct name match)",
        "figures/cross_species_drug_waffle.pdf",
        "figures/cross_species_drug_waffle.png",
        label_noun="drug",
    )


def main():
    run_gene_analysis()
    run_drug_analysis()


if __name__ == "__main__":
    main()
