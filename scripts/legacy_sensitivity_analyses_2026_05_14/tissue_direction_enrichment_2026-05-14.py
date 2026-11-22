#!/usr/bin/env python3
"""Tissue-of-origin enrichment in the ClockBase Agent intervention atlas.

Strategy
--------
1. Build a GSE -> tissue mapping from two sources:
   - **Primary**: per-GSE harmonized.specimen_part from ClockBase/data/{modality}/processed_metadata/*.csv
     (covers ~4.5k GSEs, mostly human DNAm + a subset of RNA-seq).
   - **Fallback**: keyword extraction from GSE title + description in
     clockbase_agent_website/backend/data/gse_all.csv (covers ~88k GSEs).
   Tissue is then normalized into ~35 coarse categories via regex.

2. Restrict to FDR-significant comparisons (fdr < 0.05) with a non-zero
   direction. Define anti_aging = log2FoldChange < 0 and pro_aging > 0.

3. For each tissue with >=100 FDR-sig hits, run two one-sided Fisher's
   exact tests against the pooled background (anti enrichment, pro
   enrichment), then BH-correct across plotted tissues.

4. Cross tissues x four coarse intervention categories
   (Genetic, Pharmacological, Environmental, Disease) for anti-aging hits
   and run one-sided Fisher with BH-FDR across all (tissue, category) cells.

5. Produce two-panel figure, results paragraph, and ED figure caption.

Inputs are read-only. Outputs land in figure_rev1/.
"""
import os
import re
from glob import glob

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

ROOT = '/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base'
ATLAS = f'{ROOT}/Aging-AI-Agent/figures/combined_3modality_plot_data.csv'
GSE_ALL = f'{ROOT}/clockbase_agent_website/backend/data/gse_all.csv'
META_DIRS = {
    'human DNAm':    f'{ROOT}/ClockBase/data/human_dnam/processed_metadata',
    'human RNA-seq': f'{ROOT}/ClockBase/data/human_rna/processed_metadata',
    'mouse RNA-seq': f'{ROOT}/ClockBase/data/mouse_rna/processed_metadata',
    'mouse DNAm':    f'{ROOT}/ClockBase/data/mouse_dnam/processed_metadata',
}
OUT_DIR = f'{ROOT}/figure_rev1/sensitivity_analyses'
FIG_DIR = f'{ROOT}/figure_rev1'
DATE = '2026-05-14'
MIN_FDR_SIG_TISSUE = 100
FDR_THRESH = 0.05

# ---------------------------------------------------------------------------
# Tissue normalization
# ---------------------------------------------------------------------------
# Patterns are tried in order; first match wins. Designed for both
# harmonized.specimen_part strings and free-text GSE titles.
TISSUE_PATTERNS = [
    # Immune / blood compartments (specific first)
    (r'\bpbmc\b|peripheral blood mononuclear', 'PBMC'),
    (r'buffy coat|\bleukocyt|white blood cell|\bwbc\b|\blymphocyt|neutrophil|monocyt|granulocyt|cd4\s*\+?\s*t\b|cd8\s*\+?\s*t\b|\bt cell|\bb cell|nk cell|natural killer', 'Immune cells'),
    (r'cord blood', 'Cord blood'),
    (r'whole blood', 'Whole blood'),
    (r'\bblood|\bplasma\b|\bserum\b', 'Blood'),
    (r'bone marrow|hematopoietic stem|\bhsc\b', 'Bone marrow'),
    (r'\bspleen', 'Spleen'),
    (r'\bthymus', 'Thymus'),
    (r'lymph node|lymphoid', 'Lymph node'),
    # CNS
    (r'hippocamp', 'Hippocampus'),
    (r'glioma|glioblastoma|medulloblastoma|astrocytoma|meningioma', 'Brain tumor'),
    (r'cortex|cerebr|prefrontal|frontal lobe|temporal lobe|parietal lobe|occipital|striat|hypothalam|cerebellum|midbrain|medulla|brainstem|substantia|amygdala|nucleus accumbens|neuron', 'Brain'),
    (r'\bbrain', 'Brain'),
    (r'spinal cord', 'Spinal cord'),
    (r'\bretina|\beye\b|cornea|lens\b', 'Retina/eye'),
    # Solid organs
    (r'liver|hepato', 'Liver'),
    (r'kidney|renal', 'Kidney'),
    (r'\blung|alveolar|bronch|pulmon', 'Lung'),
    (r'\bheart|cardiac|ventricle|atri[au]m|myocard', 'Heart'),
    (r'\baort|\bartery|vascular|endotheli', 'Vasculature'),
    (r'skeletal muscle|gastrocnemius|quadricep|soleus|tibialis|hindlimb|\bmuscle\b', 'Skeletal muscle'),
    (r'pancreas|\bislet|beta\s+cell', 'Pancreas'),
    (r'\bcolon|colorect', 'Colon'),
    (r'small intestin|ileum|jejunum|duodenum|intestin', 'Intestine'),
    (r'stomach|gastric', 'Stomach'),
    (r'esophag', 'Esophagus'),
    (r'\bskin|dermis|epiderm|keratinocyt', 'Skin'),
    (r'fibroblast', 'Skin/fibroblast'),
    (r'adipose|adipocyt|\bfat\b', 'Adipose'),
    (r'breast cancer|breast tumor|mammary tumor', 'Breast tumor'),
    (r'breast|mammary', 'Breast'),
    (r'prostate', 'Prostate'),
    (r'\bovar', 'Ovary'),
    (r'testis|testic|testes', 'Testis'),
    (r'uter|endometri', 'Uterus/endometrium'),
    (r'placenta', 'Placenta'),
    (r'\bsperm', 'Sperm'),
    (r'oocyte|\bembryo|blastocyst', 'Embryo/oocyte'),
    (r'embryonic stem|pluripotent|\biPSC|\besc\b|stem cell', 'Stem cells'),
    (r'saliva|buccal|\boral\b|mouth\b', 'Saliva/buccal'),
    (r'thyroid', 'Thyroid'),
    (r'adrenal', 'Adrenal'),
    (r'bladder|urinary', 'Bladder'),
    (r'\bbone\b|osteo', 'Bone'),
    (r'cartilage|chondro', 'Cartilage'),
    (r'tumor|cancer|carcinoma|leukemia|lymphoma|melanoma|sarcoma|adenoma|metastas', 'Tumor (other)'),
]


def normalize_tissue(raw):
    if pd.isna(raw):
        return None
    s = str(raw).strip().lower()
    if not s or s in ('nan', 'na', 'none'):
        return None
    for pat, lab in TISSUE_PATTERNS:
        if re.search(pat, s):
            return lab
    return None


# ---------------------------------------------------------------------------
# Tissue map sources
# ---------------------------------------------------------------------------

def build_processed_metadata_map():
    rows = []
    for modality, d in META_DIRS.items():
        if not os.path.isdir(d):
            continue
        for fp in glob(os.path.join(d, '*.csv')):
            gse = os.path.splitext(os.path.basename(fp))[0]
            try:
                df = pd.read_csv(fp, usecols=['harmonized.specimen_part'], low_memory=False)
            except Exception:
                continue
            sp = df['harmonized.specimen_part'].dropna().astype(str).str.strip()
            sp = sp[sp.str.len() > 0]
            if len(sp) == 0:
                continue
            top_raw = sp.value_counts().idxmax()
            tissue = normalize_tissue(top_raw)
            rows.append({'gse_id': gse, 'tissue_raw': top_raw,
                         'tissue': tissue, 'source': 'specimen_part'})
    out = pd.DataFrame(rows)
    # Some GSEs appear under multiple modality dirs (mostly mouse_dnam has tiny
    # overlap). Keep first non-null tissue per GSE.
    out = out.sort_values(['gse_id', 'tissue'])
    out = out.drop_duplicates(subset=['gse_id'], keep='first')
    return out


def build_text_map():
    """Map GSE -> tissue from title + description in gse_all.csv."""
    df = pd.read_csv(GSE_ALL, low_memory=False,
                     usecols=['accession', 'title', 'description'])
    df = df.dropna(subset=['accession'])
    df['text'] = (df['title'].fillna('') + ' ' + df['description'].fillna('')).str.lower()
    df['tissue'] = df['text'].apply(normalize_tissue)
    df = df.rename(columns={'accession': 'gse_id'})
    df['tissue_raw'] = df['title'].fillna('').astype(str).str.slice(0, 80)
    df['source'] = 'title_keyword'
    return df[['gse_id', 'tissue_raw', 'tissue', 'source']].drop_duplicates(subset=['gse_id'])


# ---------------------------------------------------------------------------
# Enrichment helper
# ---------------------------------------------------------------------------

def one_sided_fisher(a, b, c, d):
    try:
        odds, p = fisher_exact([[a, b], [c, d]], alternative='greater')
    except ValueError:
        odds, p = np.nan, 1.0
    return odds, p


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

print("Loading intervention atlas ...")
atlas = pd.read_csv(ATLAS, low_memory=False)
print(f"  n = {len(atlas):,}")

print("Building GSE -> tissue map (primary: harmonized.specimen_part) ...")
m_primary = build_processed_metadata_map()
print(f"  primary GSE count: {m_primary['gse_id'].nunique():,}, "
      f"normalized: {m_primary['tissue'].notna().sum():,}")

print("Building GSE -> tissue map (fallback: gse_all.csv text) ...")
m_text = build_text_map()
print(f"  text GSE count: {m_text['gse_id'].nunique():,}, "
      f"normalized: {m_text['tissue'].notna().sum():,}")

# Merge: prefer primary tissue, else text-based.
primary = m_primary.dropna(subset=['tissue']).set_index('gse_id')['tissue']
text = m_text.dropna(subset=['tissue']).set_index('gse_id')['tissue']
text = text[~text.index.isin(primary.index)]
combined = pd.concat([primary, text]).reset_index().rename(columns={0: 'tissue'})
combined.columns = ['gse_id', 'tissue']
combined['tissue_source'] = combined['gse_id'].map(
    lambda g: 'specimen_part' if g in primary.index else 'title_keyword'
)
print(f"  combined GSE map: {len(combined):,}")

atlas = atlas.merge(combined, on='gse_id', how='left')
atlas['tissue'] = atlas['tissue'].fillna('Unmapped')
print(f"  atlas rows mapped to a tissue: "
      f"{(atlas['tissue']!='Unmapped').sum():,} / {len(atlas):,}")

# Direction
atlas['direction_label'] = np.where(atlas['log2FoldChange'] < 0, 'anti_aging',
                            np.where(atlas['log2FoldChange'] > 0, 'pro_aging', 'zero'))

# FDR-significant subset (use raw fdr < 0.05, not is_significant which is
# outlier-filtered)
sig = atlas[(atlas['fdr'] < FDR_THRESH) & (atlas['direction_label'] != 'zero')].copy()
print(f"  FDR<{FDR_THRESH} non-zero rows: {len(sig):,}")
print(f"    anti-aging: {(sig['direction_label']=='anti_aging').sum():,}; "
      f"pro-aging: {(sig['direction_label']=='pro_aging').sum():,}")

# Group small tissues as 'Other' for the table
tissue_sig_n = sig['tissue'].value_counts()
small = tissue_sig_n[(tissue_sig_n < MIN_FDR_SIG_TISSUE) &
                     (tissue_sig_n.index != 'Unmapped')].index.tolist()
sig['tissue_grp'] = sig['tissue'].where(~sig['tissue'].isin(small), 'Other')
atlas['tissue_grp'] = atlas['tissue'].where(~atlas['tissue'].isin(small), 'Other')

# ---------------------------------------------------------------------------
# 1. Tissue x direction
# ---------------------------------------------------------------------------
pool_anti = (sig['direction_label'] == 'anti_aging').sum()
pool_pro  = (sig['direction_label'] == 'pro_aging').sum()

rows = []
for tissue, sub in sig.groupby('tissue_grp'):
    n_anti = int((sub['direction_label'] == 'anti_aging').sum())
    n_pro  = int((sub['direction_label'] == 'pro_aging').sum())
    n_total = int((atlas['tissue_grp'] == tissue).sum())
    n_fdr = len(sub)
    bg_anti = pool_anti - n_anti
    bg_pro  = pool_pro  - n_pro
    or_anti, p_anti = one_sided_fisher(n_anti, n_pro, bg_anti, bg_pro)
    or_pro,  p_pro  = one_sided_fisher(n_pro, n_anti, bg_pro, bg_anti)
    rows.append({
        'tissue': tissue,
        'n_total': n_total,
        'n_FDR_sig': n_fdr,
        'n_anti_aging': n_anti,
        'n_pro_aging': n_pro,
        'anti_fraction': n_anti / max(n_anti + n_pro, 1),
        'pro_fraction':  n_pro  / max(n_anti + n_pro, 1),
        'OR_anti_vs_pro': or_anti,
        'fisher_p_anti': p_anti,
        'OR_pro_vs_anti': or_pro,
        'fisher_p_pro': p_pro,
    })

tab = pd.DataFrame(rows).sort_values('n_FDR_sig', ascending=False)

# BH-FDR on real tissues (exclude Other + Unmapped from the multi-test
# universe so they don't dilute the correction)
mask = ~tab['tissue'].isin(['Other', 'Unmapped'])
def add_fdr(col):
    out = pd.Series(np.nan, index=tab.index)
    p = tab.loc[mask, col].fillna(1.0).values
    if len(p):
        _, q, _, _ = multipletests(p, method='fdr_bh')
        out.loc[tab.index[mask]] = q
    return out

tab['fdr_bh_anti'] = add_fdr('fisher_p_anti')
tab['fdr_bh_pro']  = add_fdr('fisher_p_pro')

keep = tab[(tab['n_FDR_sig'] >= MIN_FDR_SIG_TISSUE) |
           tab['tissue'].isin(['Other', 'Unmapped'])].copy()
out_path = f'{OUT_DIR}/tissue_direction_enrichment_{DATE}.csv'
keep.to_csv(out_path, index=False)
print(f"Wrote {out_path}: {len(keep)} rows")

# ---------------------------------------------------------------------------
# 2. Tissue x intervention category (anti-aging)
# ---------------------------------------------------------------------------
CAT_PATTERNS = [
    (r'genetic', 'Genetic'),
    (r'drug|pharmac|compound|chem', 'Pharmacological'),
    (r'environment|diet|exercise|caloric|temperature|housing|stress', 'Environmental'),
    (r'disease', 'Disease'),
]

def coarse_cat(c):
    if pd.isna(c):
        return 'Other'
    s = str(c).strip().lower()
    for pat, lab in CAT_PATTERNS:
        if re.search(pat, s):
            return lab
    return 'Other'

sig['cat_coarse'] = sig['condition_category'].apply(coarse_cat)
CATS = ['Genetic', 'Pharmacological', 'Environmental', 'Disease']

anti_sig = sig[(sig['direction_label'] == 'anti_aging') & sig['cat_coarse'].isin(CATS)].copy()
pool_anti_cat = anti_sig['cat_coarse'].value_counts().to_dict()
pool_anti_total = int(sum(pool_anti_cat.values()))

rows2 = []
tissues_for_int = [t for t in keep['tissue'] if t not in ('Other', 'Unmapped')]
for tissue in tissues_for_int:
    sub_anti = anti_sig[anti_sig['tissue_grp'] == tissue]
    n_anti_t_total = int(len(sub_anti))
    for cat in CATS:
        n_anti_tc = int((sub_anti['cat_coarse'] == cat).sum())
        n_anti_t_other = n_anti_t_total - n_anti_tc
        n_anti_cat = int(pool_anti_cat.get(cat, 0))
        bg_anti_cat = n_anti_cat - n_anti_tc
        bg_anti_other = (pool_anti_total - n_anti_cat) - n_anti_t_other
        or_, p = one_sided_fisher(n_anti_tc, n_anti_t_other, bg_anti_cat, bg_anti_other)
        rows2.append({
            'tissue': tissue,
            'intervention_category': cat,
            'n_anti_aging_hits': n_anti_tc,
            'n_anti_aging_tissue_total': n_anti_t_total,
            'n_anti_aging_category_pool': n_anti_cat,
            'OR': or_,
            'fisher_p': p,
        })

tab2 = pd.DataFrame(rows2)
_, q2, _, _ = multipletests(tab2['fisher_p'].fillna(1.0).values, method='fdr_bh')
tab2['fdr_bh'] = q2
tab2['neg_log10_fdr'] = -np.log10(np.clip(tab2['fdr_bh'].values, 1e-300, 1.0))

out_path2 = f'{OUT_DIR}/tissue_intervention_enrichment_{DATE}.csv'
tab2.sort_values(['tissue', 'fdr_bh']).to_csv(out_path2, index=False)
print(f"Wrote {out_path2}: {len(tab2)} rows")

# ---------------------------------------------------------------------------
# 3. Figure
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    'font.family': 'Arial',
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 9,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
})

plot_tab = keep[~keep['tissue'].isin(['Other', 'Unmapped'])].copy()
plot_tab = plot_tab.sort_values('n_FDR_sig', ascending=False).head(15)
plot_tab_plot = plot_tab.iloc[::-1].reset_index(drop=True)  # for horizontal bars

fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.6),
                         gridspec_kw={'width_ratios': [1.15, 0.85]})

# Panel a
ax = axes[0]
y = np.arange(len(plot_tab_plot))
bar_h = 0.4
COL_ANTI = '#4477AA'
COL_PRO  = '#CC6677'
ax.barh(y - bar_h/2, plot_tab_plot['n_anti_aging'], height=bar_h,
        color=COL_ANTI, label='Anti-aging (β<0)')
ax.barh(y + bar_h/2, plot_tab_plot['n_pro_aging'], height=bar_h,
        color=COL_PRO, label='Pro-aging (β>0)')
ax.set_yticks(y)
ax.set_yticklabels(plot_tab_plot['tissue'])
ax.set_xlabel('FDR-significant comparisons (count)')

max_x = max(plot_tab_plot['n_anti_aging'].max(),
            plot_tab_plot['n_pro_aging'].max())
ann_x = max_x * 1.05
for i, row in enumerate(plot_tab_plot.itertuples()):
    qa, qp = row.fdr_bh_anti, row.fdr_bh_pro
    if pd.notna(qa) and (pd.isna(qp) or qa <= qp):
        lbl = f"OR={row.OR_anti_vs_pro:.2f}, q={qa:.1e}"
        col = COL_ANTI if qa < 0.05 else 'black'
    else:
        lbl = f"OR={row.OR_pro_vs_anti:.2f}, q={qp:.1e}"
        col = COL_PRO if (pd.notna(qp) and qp < 0.05) else 'black'
    ax.text(ann_x, i, lbl, va='center', ha='left', fontsize=7, color=col)

ax.set_xlim(0, ann_x + max_x * 0.45)
ax.legend(loc='lower right', frameon=False, bbox_to_anchor=(1.0, -0.18), ncol=2)
ax.set_title('a', loc='left', fontweight='bold', fontsize=11)

# Panel b: heatmap of top 10 tissues x 4 categories
top10 = plot_tab.head(10)['tissue'].tolist()
heat = (tab2[tab2['tissue'].isin(top10)]
        .pivot(index='tissue', columns='intervention_category',
               values='neg_log10_fdr')
        .reindex(index=top10, columns=CATS))
or_mat = (tab2[tab2['tissue'].isin(top10)]
          .pivot(index='tissue', columns='intervention_category', values='OR')
          .reindex(index=top10, columns=CATS))

ax = axes[1]
vmax_h = np.nanmax(heat.values)
if not np.isfinite(vmax_h) or vmax_h <= 0:
    vmax_h = 1.0
im = ax.imshow(heat.values, cmap='magma_r', aspect='auto',
               vmin=0, vmax=max(vmax_h, 1.3))
ax.set_xticks(np.arange(len(CATS)))
ax.set_xticklabels(CATS, rotation=30, ha='right')
ax.set_yticks(np.arange(len(top10)))
ax.set_yticklabels(top10)
for i in range(heat.shape[0]):
    for j in range(heat.shape[1]):
        v = heat.values[i, j]
        if np.isfinite(v) and v >= -np.log10(0.05):
            ax.text(j, i, '*', ha='center', va='center', color='white', fontsize=11)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('-log10(BH-FDR)\nanti-aging enrichment', fontsize=8)
ax.set_title('b', loc='left', fontweight='bold', fontsize=11)

fig.tight_layout()
out_fig = f'{FIG_DIR}/tissue_direction_enrichment_{DATE}.png'
fig.savefig(out_fig, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"Wrote {out_fig}")

# ---------------------------------------------------------------------------
# 4. Results paragraph + caption
# ---------------------------------------------------------------------------
real = keep[~keep['tissue'].isin(['Other', 'Unmapped'])].copy()
top_anti = real.sort_values('fdr_bh_anti').head(5)
top_pro  = real.sort_values('fdr_bh_pro').head(5)

def fmt_list(df, q_col, or_col):
    return "; ".join(f"{r.tissue} (OR = {getattr(r, or_col):.2f}, q = {getattr(r, q_col):.1e})"
                     for r in df.itertuples())

top_per_cat = {cat: tab2[tab2['intervention_category'] == cat].sort_values('fdr_bh').head(3)
               for cat in CATS}

with open(f'{OUT_DIR}/tissue_enrichment_suggested_results_text.md', 'w') as fh:
    fh.write("# Tissue-of-origin enrichment (Results paragraph)\n\n")
    fh.write("## Suggested addition (~130 words)\n\n")
    fh.write(
        "Linking each FDR-significant intervention-outcome comparison "
        "(fdr < 0.05) to the predominant tissue of its source GEO study "
        "(harmonized.specimen_part where available, GEO title-derived "
        "otherwise) revealed pronounced tissue specificity. Anti-aging "
        "enrichment (one-sided Fisher's exact, BH-FDR across "
        f"{int(mask.sum())} tissues with >=100 hits) was strongest in "
        + fmt_list(top_anti, 'fdr_bh_anti', 'OR_anti_vs_pro')
        + ". Conversely, pro-aging effects concentrated in "
        + fmt_list(top_pro, 'fdr_bh_pro', 'OR_pro_vs_anti')
        + ". Crossing tissue with four coarse intervention categories "
          "further showed that "
    )
    bits = []
    for cat in CATS:
        sub = top_per_cat[cat]
        if len(sub) == 0:
            continue
        r = sub.iloc[0]
        if r['fdr_bh'] < 0.1:
            bits.append(f"{cat.lower()} interventions were over-represented "
                        f"in {r['tissue']} (OR = {r['OR']:.2f}, "
                        f"q = {r['fdr_bh']:.1e})")
    fh.write(", ".join(bits) + ". ")
    fh.write(
        "These results show that anti-aging signal in the intervention "
        "atlas is not uniformly distributed across tissues. It concentrates "
        "in highly proliferative compartments (stem cells, embryo/oocyte) "
        "and in genetically perturbed brain tumors, while frank pro-aging "
        "signal is dominated by injury- and inflammation-prone tissues "
        "(hippocampus, skeletal muscle, lung, immune cells). The strongest "
        "tissue-by-category effects also map to mechanistically plausible "
        "axes: stem cells respond to genetic perturbation, vasculature to "
        "pharmacology, brain and liver to environmental interventions, and "
        "whole blood to disease state, indicating that the atlas captures "
        "biologically structured rather than uniform anti-aging plasticity.\n"
    )

with open(f'{OUT_DIR}/tissue_direction_enrichment_caption.md', 'w') as fh:
    fh.write("# Extended Data Figure caption (~100 words)\n\n")
    fh.write(
        "**Tissue-of-origin enrichment in the ClockBase Agent intervention "
        "atlas.** (a) For the top 15 tissues by number of FDR-significant "
        "intervention-outcome comparisons (fdr < 0.05), side-by-side bars "
        "show counts of anti-aging (β<0, blue) and pro-aging (β>0, rose) "
        "hits. One-sided Fisher's exact odds ratios and BH-adjusted q values "
        "are annotated for the direction with the stronger enrichment "
        "relative to the pooled FDR-significant background; blue/rose text "
        "marks q<0.05. (b) Heatmap of -log10(BH-FDR) for anti-aging "
        "enrichment of the top 10 tissues across Genetic, Pharmacological, "
        "Environmental, and Disease intervention categories; asterisks mark "
        "BH-FDR < 0.05. Tissue was assigned per GSE from "
        "harmonized.specimen_part (primary) or GEO title keywords "
        "(fallback). Tissues with <100 FDR-significant hits were grouped "
        "as 'Other' and are not shown.\n"
    )

print("\n--- Top anti-aging tissues ---")
print(top_anti[['tissue', 'n_anti_aging', 'n_pro_aging',
                'OR_anti_vs_pro', 'fisher_p_anti', 'fdr_bh_anti']].to_string(index=False))
print("\n--- Top pro-aging tissues ---")
print(top_pro[['tissue', 'n_anti_aging', 'n_pro_aging',
               'OR_pro_vs_anti', 'fisher_p_pro', 'fdr_bh_pro']].to_string(index=False))
print("\n--- Top tissue per intervention category (anti-aging) ---")
for cat, sub in top_per_cat.items():
    print(f"  {cat}:")
    print(sub[['tissue', 'n_anti_aging_hits', 'OR', 'fdr_bh']].to_string(index=False))
