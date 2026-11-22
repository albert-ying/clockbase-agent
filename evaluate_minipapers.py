#!/usr/bin/env python3
"""
Quantitative evaluation of the manuscript-agent mini-paper corpus.

Computes three families of metrics over every parsed_report.md served by the
website, cross-referencing each report against its canonical sample.json input:

  1. STRUCTURE / LENGTH   word counts (total + per section), section presence
  2. STATISTICS FIDELITY  does the effect size / p-value quoted in the prose
                          match the canonical ClockBase value the agent was
                          given? (numeric tolerance check)
  3. CITATION INTEGRITY   reference count, DOI / PMC presence, malformed-entry
                          detection, in-text [N] vs reference-list coverage

No data is fabricated: every number is parsed from real report and sample files.
Outputs a per-report CSV plus a printed summary with distributions.

Usage:
    python evaluate_minipapers.py [--limit N] [--check-doi K]

  --limit N      evaluate only the first N reports (debug)
  --check-doi K  network-validate a random sample of K DOIs via doi.org HEAD
                 (default 0 = skip network checks)
"""
import argparse
import json
import re
import statistics as stats
from pathlib import Path

REPORT_DIR = Path("../clockbase_agent_website/backend/data/interventions/report")
SAMPLE_BASES = [
    Path("minipaper_fork/output"),
    Path("data/output"),
]
SECTIONS = ["Introduction", "Results", "Discussion", "References"]

# ---- number parsing helpers -------------------------------------------------

# superscript digits / signs used in "10⁻⁷" style notation
_SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3",
        "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7",
        "⁸": "8", "⁹": "9", "⁻": "-", "⁺": "+"}
_SUP_RUN = re.compile("[" + "".join(_SUP) + "]+")


def normalize_superscripts(text):
    """Convert '10⁻⁷' / '×10⁻⁷' superscript exponents into '10^-7' so the
    numeric parser can read them. Leaves all other text intact."""
    def repl(m):
        return "^" + "".join(_SUP[c] for c in m.group(0))
    return _SUP_RUN.sub(repl, text)


def _to_float(s):
    """Parse '3.20e-05', '9.2e-4', '1.2668', '1.27 x 10^-5' style numbers."""
    s = normalize_superscripts(s.strip()).replace("−", "-")  # unicode minus
    s = s.replace("×", "x").replace("X", "x")
    m = re.match(r"^(-?\d+\.?\d*)\s*x?\s*10\^?(-?\d+)$", s)
    if m:
        return float(m.group(1)) * 10 ** int(m.group(2))
    try:
        return float(s)
    except ValueError:
        return None


NUM = r"(-?\d+\.?\d*(?:\s*[xX×]\s*10\^?-?\d+|[eE][-+]?\d+)?)"

P_PATTERNS = [
    re.compile(r"p\s*=\s*" + NUM, re.IGNORECASE),
    re.compile(r"p[-\s]?value\s*(?:of|=|:)?\s*" + NUM, re.IGNORECASE),
]
EFFECT_PATTERNS = [
    re.compile(r"log2\s*fold\s*change\s*=?\s*" + NUM, re.IGNORECASE),
    re.compile(r"log2FC\s*=?\s*" + NUM, re.IGNORECASE),
    re.compile(r"\beffect(?:\s*size)?\s*(?:was|=|of|:)?\s*" + NUM, re.IGNORECASE),
]


def first_match(patterns, text):
    for pat in patterns:
        m = pat.search(text)
        if m:
            v = _to_float(m.group(1))
            if v is not None:
                return v
    return None


# ---- section split ----------------------------------------------------------

def split_sections(md):
    """Return {section_name: body_text} using '## Header' markers."""
    out = {}
    cur = "_preamble"
    buf = []
    for line in md.splitlines():
        h = re.match(r"^#{1,3}\s+(.*)", line)
        if h:
            out[cur] = "\n".join(buf).strip()
            cur = h.group(1).strip()
            buf = []
        else:
            buf.append(line)
    out[cur] = "\n".join(buf).strip()
    return out


def wc(text):
    return len(re.findall(r"\b[\w']+\b", text))


# ---- reference / citation parsing ------------------------------------------

# Reference list entries appear in two styles across corpus vintages:
#   "[1] Title. Authors..."   (bracket style, ~28% of reports)
#   "1. Title. Authors..."     (markdown ordered-list style, ~72%)
REF_LINE = re.compile(r"^(?:\[(\d+)\]|(\d+)\.)\s+(.*)", re.MULTILINE)
DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+", re.IGNORECASE)


def clean_doi(doi):
    """Strip trailing sentence punctuation that the regex may capture."""
    return doi.rstrip(".,;)")
PMC_RE = re.compile(r"PMC\d{4,9}")
URL_RE = re.compile(r"https?://[^\s)]+", re.IGNORECASE)
INTEXT_CITE = re.compile(r"\[(\d+(?:\s*,\s*\d+)*)\]")
# A reference that is actually a figure / table fragment, not a paper
MALFORMED_RE = re.compile(r"^(figure|table|fig\.|supplementary|appendix)\b", re.IGNORECASE)


def analyze_references(ref_text, body_text):
    # each match: (bracket_id_or_'', dot_id_or_'', body)
    matches = REF_LINE.findall(ref_text)
    refs = [((b or d), t) for b, d, t in matches]
    # vintage: bracket-numbered ("[1]") vs ordered-list ("1.")
    style = "bracket" if any(b for b, d, _ in matches) else (
        "dot" if matches else "none")
    n = len(refs)
    with_doi = sum(1 for _, t in refs if DOI_RE.search(t))
    with_pmc = sum(1 for _, t in refs if PMC_RE.search(t))
    with_url = sum(1 for _, t in refs if URL_RE.search(t))
    # "linked" = reachable by SOME accessible handle: DOI, PMC id, or explicit URL
    with_link = sum(1 for _, t in refs
                    if DOI_RE.search(t) or PMC_RE.search(t) or URL_RE.search(t))
    malformed = sum(1 for _, t in refs if MALFORMED_RE.match(t.strip()))
    # title presence heuristic: a real ref has > 6 words before the doi
    titled = sum(1 for _, t in refs if wc(DOI_RE.split(t)[0]) >= 6)
    ref_ids = {int(i) for i, _ in refs}

    # in-text citations
    intext_ids = set()
    for grp in INTEXT_CITE.findall(body_text):
        for x in grp.split(","):
            x = x.strip()
            if x.isdigit():
                intext_ids.add(int(x))
    orphan_intext = sorted(intext_ids - ref_ids)      # cited but not in list
    uncited_refs = sorted(ref_ids - intext_ids)        # in list but never cited
    # Build one resolvable link per reference where possible (for network check):
    # prefer DOI -> doi.org, then PMC -> ncbi, then a raw URL.
    links = []
    for _, t in refs:
        m = DOI_RE.search(t)
        if m:
            # tag DOIs so the resolver can use the Crossref API (no bot blocking)
            links.append("doi:" + clean_doi(m.group(0))); continue
        m = PMC_RE.search(t)
        if m:
            links.append(f"https://www.ncbi.nlm.nih.gov/pmc/articles/{m.group(0)}/"); continue
        m = URL_RE.search(t)
        if m:
            links.append(m.group(0))
    return {
        "ref_style": style,
        "n_refs": n,
        "refs_with_doi": with_doi,
        "refs_with_pmc": with_pmc,
        "refs_with_url": with_url,
        "refs_with_link": with_link,
        "refs_malformed": malformed,
        "refs_titled": titled,
        "n_intext_cite_ids": len(intext_ids),
        "n_orphan_intext": len(orphan_intext),
        "orphan_intext_ids": ";".join(map(str, orphan_intext)),
        "n_uncited_refs": len(uncited_refs),
        "_links": links,
    }


# ---- canonical stats lookup -------------------------------------------------

def find_sample(case_id):
    for b in SAMPLE_BASES:
        p = b / case_id / "sample.json"
        if p.exists():
            return p
    return None


def load_canonical(case_id):
    p = find_sample(case_id)
    if not p:
        return None
    try:
        d = json.load(open(p))
    except Exception:
        return None
    eff = d.get("effect_size")
    pv = d.get("pvalue")
    return {
        "canon_effect": float(eff) if eff is not None else None,
        "canon_pvalue": float(pv) if pv is not None else None,
        "canon_category": d.get("condition_category"),
    }


def rel_match(reported, canon, rel=0.02):
    """True if reported within rel (2%) of canon, on abs value."""
    if reported is None or canon is None:
        return None
    if canon == 0:
        return abs(reported) < 1e-9
    return abs(abs(reported) - abs(canon)) / abs(canon) <= rel


# ---- main -------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--check-doi", type=int, default=0,
                    help="network-validate a random sample of K DOIs via doi.org")
    ap.add_argument("--paperclip-only", action="store_true",
                    help="restrict to paperclip-generated reports "
                         "(references carry citations.gxl.ai / PMC paperclip links)")
    ap.add_argument("--out", default="result/minipaper_eval_2026-05-21.csv")
    args = ap.parse_args()

    dirs = sorted([d for d in REPORT_DIR.iterdir() if d.is_dir()])
    if args.limit:
        dirs = dirs[: args.limit]

    # paperclip signature: the virtual-FS citation server URL
    PAPERCLIP_SIG = re.compile(r"citations\.gxl\.ai|ncbi\.nlm\.nih\.gov/pmc",
                               re.IGNORECASE)

    rows = []
    for d in dirs:
        f = d / "parsed_report.md"
        if not f.exists():
            continue
        md = f.read_text(errors="replace")
        if args.paperclip_only and not PAPERCLIP_SIG.search(md):
            continue
        secs = split_sections(md)
        body = "\n".join(v for k, v in secs.items() if k != "References")
        ref_text = secs.get("References", "")

        results_norm = normalize_superscripts(secs.get("Results", "") or md)
        rep_effect = first_match(EFFECT_PATTERNS, results_norm)
        rep_p = first_match(P_PATTERNS, results_norm)

        canon = load_canonical(d.name) or {}
        eff_ok = rel_match(rep_effect, canon.get("canon_effect"))
        p_ok = rel_match(rep_p, canon.get("canon_pvalue"), rel=0.05)

        # ---- numeric agreement (signed/relative errors) ----
        ce, cp = canon.get("canon_effect"), canon.get("canon_pvalue")
        # effect: relative error on absolute value (effect sign handled separately)
        eff_relerr = (abs(abs(rep_effect) - abs(ce)) / abs(ce)
                      if (rep_effect is not None and ce not in (None, 0)) else None)
        # p-value spans orders of magnitude -> compare on log10 scale
        import math
        p_log10_absdiff = None
        if rep_p is not None and cp not in (None, 0) and rep_p > 0 and cp > 0:
            p_log10_absdiff = abs(math.log10(rep_p) - math.log10(cp))

        ref_stats = analyze_references(ref_text, body)

        row = {
            "case_id": d.name,
            "total_words": wc(md),
            "intro_words": wc(secs.get("Introduction", "")),
            "results_words": wc(secs.get("Results", "")),
            "discussion_words": wc(secs.get("Discussion", "")),
            "has_all_sections": all(s in secs for s in SECTIONS),
            "missing_sections": ";".join(s for s in SECTIONS if s not in secs),
            "reported_effect": rep_effect,
            "reported_pvalue": rep_p,
            "canon_effect": canon.get("canon_effect"),
            "canon_pvalue": canon.get("canon_pvalue"),
            "effect_match": eff_ok,
            "pvalue_match": p_ok,
            "effect_relerr": eff_relerr,
            "pvalue_log10_absdiff": p_log10_absdiff,
            "has_canon": bool(canon),
            **ref_stats,
        }
        rows.append(row)

    # collect all resolvable links before stripping the helper column
    all_links = []
    for r in rows:
        all_links.extend(r.pop("_links", []))

    # ---- write CSV ----
    import csv
    out = Path(args.out)
    out.parent.mkdir(exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- summary ----
    n = len(rows)
    def col(name):
        return [r[name] for r in rows]
    def dist(name, fmt="{:.0f}"):
        xs = [x for x in col(name) if isinstance(x, (int, float))]
        if not xs:
            return "n/a"
        xs.sort()
        g = lambda v: fmt.format(v)
        return (f"min={g(xs[0])} p25={g(xs[len(xs)//4])} "
                f"med={g(stats.median(xs))} p75={g(xs[3*len(xs)//4])} "
                f"max={g(xs[-1])} mean={g(stats.mean(xs))}")

    def dist_vals(xs, fmt="{:.3g}"):
        xs = sorted(x for x in xs if isinstance(x, (int, float)))
        if not xs:
            return "n/a"
        g = lambda v: fmt.format(v)
        return (f"min={g(xs[0])} p25={g(xs[len(xs)//4])} "
                f"med={g(stats.median(xs))} p75={g(xs[3*len(xs)//4])} "
                f"p95={g(xs[int(0.95*len(xs))-1])} max={g(xs[-1])}")

    print("=" * 64)
    print(f"MINI-PAPER CORPUS EVALUATION  (n = {n} reports)")
    print("=" * 64)

    print("\n[1] STRUCTURE / LENGTH")
    print(f"  reports with all 4 sections : {sum(col('has_all_sections'))}/{n} "
          f"({100*sum(col('has_all_sections'))/n:.1f}%)")
    print(f"  total words   : {dist('total_words')}")
    print(f"  intro words   : {dist('intro_words')}")
    print(f"  results words : {dist('results_words')}")
    print(f"  discussion    : {dist('discussion_words')}")

    print("\n[2] STATISTICS FIDELITY (reported prose vs canonical sample.json)")
    have = [r for r in rows if r["has_canon"]]
    print(f"  reports with canonical stats     : {len(have)}/{n} ({100*len(have)/n:.0f}%)")

    # --- effect size: numeric relative-error distribution ---
    eff_err = [r["effect_relerr"] for r in have if r["effect_relerr"] is not None]
    print(f"\n  EFFECT SIZE (n={len(eff_err)} quoted & checkable)")
    print(f"    relative error |reported-canon|/|canon| : {dist_vals(eff_err)}")
    for thr in (0.001, 0.01, 0.05):
        k = sum(1 for e in eff_err if e <= thr)
        print(f"    within {thr*100:>4.1f}% : {k}/{len(eff_err)} ({100*k/len(eff_err):.1f}%)")

    # --- p-value: numeric agreement on log10 scale ---
    p_err = [r["pvalue_log10_absdiff"] for r in have
             if r["pvalue_log10_absdiff"] is not None]
    print(f"\n  P-VALUE (n={len(p_err)} quoted & checkable)")
    print(f"    |log10(reported) - log10(canon)| : {dist_vals(p_err)}")
    for thr, lbl in ((0.05, "≈ same to 2 sig figs"),
                     (0.30, "within ~2x"),
                     (1.0, "within 1 order of magnitude")):
        k = sum(1 for e in p_err if e <= thr)
        print(f"    Δlog10 ≤ {thr:<4} ({lbl:<27}): {k}/{len(p_err)} ({100*k/len(p_err):.1f}%)")
    gross = [r for r in have if r["pvalue_log10_absdiff"] is not None
             and r["pvalue_log10_absdiff"] > 1.0]
    print(f"    gross mismatches (>1 order)      : {len(gross)}")
    for r in gross[:8]:
        print(f"      {r['case_id']}: reported p={r['reported_pvalue']:.2g} "
              f"vs canon p={r['canon_pvalue']:.2g}")

    no_effect = sum(1 for r in have if r["reported_effect"] is None)
    no_p = sum(1 for r in have if r["reported_pvalue"] is None)
    print(f"\n  reports quoting NO effect size   : {no_effect}/{len(have)}")
    print(f"  reports quoting NO p-value       : {no_p}/{len(have)}")

    print("\n[3] CITATION INTEGRITY")
    print(f"  references per report : {dist('n_refs')}")
    tot_refs = sum(col("n_refs"))
    tot_doi = sum(col("refs_with_doi"))
    tot_pmc = sum(col("refs_with_pmc"))
    tot_url = sum(col("refs_with_url"))
    tot_link = sum(col("refs_with_link"))
    tot_mal = sum(col("refs_malformed"))
    tot_titled = sum(col("refs_titled"))
    print(f"  total references         : {tot_refs}")
    print(f"  with accessible link     : {tot_link} ({100*tot_link/tot_refs:.1f}%)  "
          f"<- DOI/PMC/URL (the validity criterion)")
    print(f"    - via DOI              : {tot_doi} ({100*tot_doi/tot_refs:.1f}%)")
    print(f"    - via PMC id           : {tot_pmc} ({100*tot_pmc/tot_refs:.1f}%)")
    print(f"    - via explicit URL     : {tot_url} ({100*tot_url/tot_refs:.1f}%)")
    reps_full = sum(1 for r in rows if r["n_refs"] > 0 and r["refs_with_link"] == r["n_refs"])
    reps_any = sum(1 for r in rows if r["n_refs"] > 0)
    print(f"  reports w/ ALL refs linked : {reps_full}/{reps_any} "
          f"({100*reps_full/reps_any:.1f}%)")
    print(f"  with plausible title     : {tot_titled} ({100*tot_titled/tot_refs:.1f}%)")
    print(f"  malformed (figure/table) : {tot_mal} ({100*tot_mal/tot_refs:.2f}%)")
    reps_with_mal = sum(1 for r in rows if r["refs_malformed"] > 0)
    print(f"  reports w/ >=1 malformed ref : {reps_with_mal}/{n} "
          f"({100*reps_with_mal/n:.1f}%)")
    orphan = sum(col("n_orphan_intext"))
    reps_orphan = sum(1 for r in rows if r["n_orphan_intext"] > 0)
    print(f"  in-text [N] with no matching ref : {orphan} across {reps_orphan} reports")
    uncited = sum(col("n_uncited_refs"))
    print(f"  reference-list entries never cited in text : {uncited}")

    # ---- vintage breakdown ----
    print("\n[3b] CITATION LINKING BY REPORT VINTAGE (reference style)")
    for sty in ("bracket", "dot", "none"):
        sub = [r for r in rows if r["ref_style"] == sty]
        if not sub:
            continue
        nref = sum(r["n_refs"] for r in sub)
        nlink = sum(r["refs_with_link"] for r in sub)
        ndoi = sum(r["refs_with_doi"] for r in sub)
        npmc = sum(r["refs_with_pmc"] for r in sub)
        link_pct = f"{100*nlink/nref:.0f}%" if nref else "n/a"
        doi_pct = f"{100*ndoi/nref:.0f}%" if nref else "n/a"
        pmc_pct = f"{100*npmc/nref:.0f}%" if nref else "n/a"
        print(f"  {sty:8s}: {len(sub):4d} reports, {nref:5d} refs, "
              f"linked {link_pct} (DOI {doi_pct}, PMC {pmc_pct})")

    # ---- optional network link resolution ----
    if args.check_doi:
        import random
        import urllib.request
        import urllib.error
        import urllib.parse
        random.seed(0)
        uniq = sorted(set(all_links))
        sample = random.sample(uniq, min(args.check_doi, len(uniq)))
        print(f"\n[3c] LINK RESOLUTION (network HEAD/GET, "
              f"n={len(sample)} unique sampled of {len(uniq)} total accessible links)")
        ok = bad = err = 0
        bad_list = []
        for url in sample:
            try:
                if url.startswith("doi:"):
                    # validate registration via Crossref (publishers 403 bots on
                    # doi.org HEAD; Crossref returns 200 for any registered DOI)
                    doi = url[4:]
                    api = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
                    req = urllib.request.Request(
                        api, method="GET",
                        headers={"User-Agent": "minipaper-eval/1.0 (mailto:research@example.org)"})
                    with urllib.request.urlopen(req, timeout=20) as resp:
                        ok += 1 if resp.status < 400 else 0
                        bad += 0 if resp.status < 400 else 1
                        if resp.status >= 400:
                            bad_list.append(doi)
                    continue
                req = urllib.request.Request(url, method="HEAD",
                                             headers={"User-Agent": "minipaper-eval/1.0"})
                with urllib.request.urlopen(req, timeout=20) as resp:
                    if resp.status < 400:
                        ok += 1
                    else:
                        bad += 1; bad_list.append(url)
            except urllib.error.HTTPError as e:
                # 403 = resolver reached the record but the host blocks bots
                # (the link is still valid for a human reader)
                if e.code in (301, 302, 303, 307, 308, 403):
                    ok += 1
                else:
                    bad += 1; bad_list.append(f"{url} (HTTP {e.code})")
            except Exception:
                err += 1
        print(f"  resolved OK : {ok}/{len(sample)} ({100*ok/len(sample):.0f}%)")
        print(f"  unresolved  : {bad}")
        print(f"  network err : {err}")
        if bad_list:
            print("  failed links (first 15):")
            for d in bad_list[:15]:
                print(f"    {d}")

    print(f"\nPer-report detail written to: {out}")


if __name__ == "__main__":
    main()
