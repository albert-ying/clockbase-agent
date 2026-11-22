#!/usr/bin/env python3
"""Classify the GEO inventory (geo_raw_inventory.json) and peek at candidate tables.

Peek: HTTP Range request for the first 512 KiB of each candidate processed table
(gz decompressed as a stream) to record header and whether the first numeric
values fall within [0,1] (beta-like) or not (M-values/signals). One request per
candidate table, same <=3 req/s throttle. Nothing is downloaded beyond that prefix.
"""
import csv, json, os, re, time, urllib.request, zlib, sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from audit_system import paths

OUT = os.environ.get("INVENTORY_OUT", str(paths.VERIFIER_ROOT / "empty-repair-inventory"))
SRC = os.environ.get("INVENTORY_SOURCE", str(paths.EMPTY_MATRIX_LIST))
NAME = os.environ.get("INVENTORY_NAME", "empty_dnam_geo_repair")  # output file stem
RAW = os.path.join(OUT, "raw")
PEEK = os.path.join(RAW, "peek")
os.makedirs(PEEK, exist_ok=True)
EVID = paths.REPO / "paper/revision_2026-09-07/evidence"
LINKED = json.load(open(EVID / "empty-input-impact-2026-09-18/run-001/summary.json"))[
    "linked_census_studies"]["historical_dnam_atlas"]
ATLAS_ROWS = defaultdict(int)
for row in csv.DictReader(open(EVID / "empty-input-impact-independent-2026-09-18/atlas_gse_counts.csv")):
    if row["source"] == "historical_dnam_atlas":
        ATLAS_ROWS[row["gse_id"]] += int(row["row_count"])

GPL_NAMES = {"GPL13534": "450K", "GPL16304": "450K", "GPL18809": "450K", "GPL21145": "EPIC",
             "GPL23976": "EPIC", "GPL29753": "EPIC", "GPL33022": "EPICv2", "GPL34394": "EPICv2",
             "GPL8490": "27K", "GPL9183": "27K"}

MIN_INTERVAL = 0.34
RANGE = int(os.environ.get("PEEK_RANGE", 524288))
_last = [0.0]


def throttle():
    dt = time.time() - _last[0]
    if dt < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - dt)
    _last[0] = time.time()


def stub(gse):
    return gse[:-3] + "nnn" if len(gse) > 6 else "GSEnnn"


def gpl_title(gpl):
    p = os.path.join(RAW, f"{gpl}.gpl.txt")
    if not os.path.exists(p):
        for a in range(5):
            throttle()
            try:
                url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=brief"
                with urllib.request.urlopen(url, timeout=60) as r:
                    open(p, "w").write(r.read().decode("utf-8", "replace"))
                break
            except Exception:
                time.sleep(2 ** a)
    txt = open(p).read() if os.path.exists(p) else ""
    t = re.findall(r"^!Platform_title = (.+)$", txt, re.M)
    tech = re.findall(r"^!Platform_technology = (.+)$", txt, re.M)
    return (t[0] if t else ""), (tech[0] if tech else "")


def platform_label(gpl, title, tech):
    if gpl in GPL_NAMES:
        return GPL_NAMES[gpl]
    tl = title.lower()
    if "epic" in tl and ("v2" in tl or "2.0" in tl):
        return "EPICv2"
    if "epic" in tl or "850k" in tl:
        return "EPIC"
    if "450" in tl:
        return "450K"
    if "27k" in tl or "humanmethylation27" in tl:
        return "27K"
    if "sequencing" in tech.lower() or "high-throughput" in tech.lower():
        return "SEQ"
    return "other"


SIG = re.compile(r"signal|intensit|singal|unmethylated|methylated_?signal|raw_?matrix|unnormali|_raw\.|rawdata|"
                 r"_meth\.|_unmeth\.|red_green|redgreen|detection_?p", re.I)
BETA = re.compile(r"beta|bval|b_val|processed|normali[sz]ed|norm|matrix|avg|methylation_?data|"
                  r"methylationprofile|methylation_?values|mvalue|m_value|final_?report|_data\.", re.I)
SKIP = re.compile(r"filelist\.txt$|samplesheet|sample_sheet|annotation|manifest|readme|\.pdf$|"
                  r"description|_RAW\.tar$|GPL\d+_|\.bed|\.bw$|\.bigwig|\.cov|\.tdf$|\.fastq", re.I)
SEQ_FILE = re.compile(r"\.bed(\.gz)?$|\.bw$|\.bigwig$|\.cov(\.gz)?$|\.bismark|cpg_report|\.cx_report|"
                      r"methylkit|\.bedgraph", re.I)


def peek(gse, name):
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    p = os.path.join(PEEK, safe + ".head.txt")
    if os.path.exists(p):
        return open(p, encoding="utf-8", errors="replace").read()
    if not name.lower().endswith((".gz", ".txt", ".csv", ".tsv")):
        return ""
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub(gse)}/{gse}/suppl/{name}"
    data = None
    for a in range(5):
        throttle()
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=0-%d" % (RANGE - 1),
                                                       "User-Agent": "clockbase-audit/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read(RANGE)
            break
        except Exception:
            time.sleep(2 ** a)
    if data is None:
        return ""
    if name.lower().endswith(".gz"):
        try:
            d = zlib.decompressobj(16 + zlib.MAX_WBITS)
            text = d.decompress(data).decode("utf-8", "replace")
        except Exception as e:
            text = f"__decompress_error__ {e}"
    else:
        text = data.decode("utf-8", "replace")
    lines = text.splitlines()[:60]
    text = "\n".join(l[:4000] for l in lines)
    open(p, "w").write(text)
    return text


def assess(text):
    """Return (kind, ncols, note) from table head."""
    if not text or text.startswith("__decompress_error__"):
        return "unknown", 0, "no peek"
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith("#")
             and l.replace("\t", "").replace(",", "").strip()]
    if not lines:
        return "unknown", 0, "empty head"
    def split(l):
        for d in ("\t", ",", " "):
            if l.count(d) >= 2:
                return [c.strip().strip('"') for c in l.split(d)]
        return [l.strip()]
    rows = [split(l) for l in lines]
    hdr = rows[0]
    hl = " ".join(hdr).lower()
    vals = []
    for rr in rows[1:60]:
        for c in rr[1:]:
            try:
                vals.append(float(c))
            except ValueError:
                pass
    sig_hdr = bool(re.search(r"signal|unmethylated|methylated|_meth\b|_unmeth\b|intensity", hl))
    pval_hdr = "pval" in hl or "detection" in hl or "p_val" in hl
    beta_hdr = "beta" in hl
    cg_hdr = sum(bool(re.match(r"(cg|ch\.)\d+", h)) for h in hdr)
    if not vals and cg_hdr > 0.5 * len(hdr) and len(hdr) > 100:
        return "transposed", len(hdr), "CpG-by-column (transposed) layout; values beyond peek prefix"
    if not vals:
        return "unknown", len(hdr), "no numeric values in head"
    import math
    vals = [v for v in vals if not math.isnan(v)]
    if not vals:
        return "unknown", len(hdr), "only NaN in head"
    frac01 = sum(0 <= v <= 1 for v in vals) / len(vals)
    if all(0 <= v <= 1 for v in vals) and len(vals) >= 20:
        return "beta", len(hdr), "all values in [0,1]" + ("; has pval cols" if pval_hdr else "")
    big = sum(v > 100 for v in vals) / len(vals)
    neg = sum(v < 0 for v in vals) / len(vals)
    pct = sum(0 <= v <= 100 for v in vals) / len(vals)
    if sig_hdr and not beta_hdr and big > 0.3:
        return "signals", len(hdr), f"header has signal cols; {big:.0%} values>100"
    if frac01 >= 0.95:
        return "beta", len(hdr), "values in [0,1]" + ("; has pval cols" if pval_hdr else "")
    if beta_hdr and frac01 >= 0.3:
        return "beta", len(hdr), f"beta cols mixed with others ({frac01:.0%} in [0,1])"
    if neg > 0.05 and big < 0.05:
        return "mvalue", len(hdr), f"{neg:.0%} negative values (likely M-values or log ratios)"
    if big > 0.3:
        return "signals", len(hdr), f"{big:.0%} values>100"
    if pct >= 0.95:
        return "percent", len(hdr), "values in [0,100] (possibly percent methylation)"
    return "unknown", len(hdr), f"{frac01:.0%} in [0,1]"


INFINIUM = {"450K", "EPIC", "EPICv2", "27K"}
OTHER_ILLUMINA_METH = {"GPL18246", "GPL21161", "GPL32946", "GPL19718"}
DETP = re.compile(r"detection|pval|p_val|detp|residual", re.I)


def main():
    inv = json.load(open(os.path.join(OUT, "geo_raw_inventory.json")))
    gpls = sorted({g for r in inv for g in r["platforms"]})
    gpl_info = {g: gpl_title(g) for g in gpls}
    rows, summary_rows = [], []
    for r in inv:
        gse = r["gse"]
        plats = [platform_label(g, *gpl_info[g]) for g in r["platforms"]]
        stype = "; ".join(r["series_type"])
        files = r["suppl_files"] or []
        tarlist = r["raw_tar_filelist"] or []
        idat = [f for f in tarlist if f["type"].upper() == "IDAT" or f["name"].lower().endswith((".idat", ".idat.gz"))]
        idat_bytes = sum(f["size"] or 0 for f in idat)
        tar_size = sum(s for n, s in files if n.upper().endswith("_RAW.TAR"))
        tar_other = [f for f in tarlist if f not in idat and f["kind"] == "File"
                     and not re.search(r"GPL\d+_|manifest", f["name"], re.I)]
        seq_files = [n for n, s in files if SEQ_FILE.search(n)] + [f["name"] for f in tar_other if SEQ_FILE.search(f["name"])]
        beta_c, sig_c, other_c = [], [], []
        for n, s in files:
            if n.endswith("/") or SKIP.search(n):
                continue
            if SIG.search(n):
                sig_c.append((n, s))
            elif BETA.search(n):
                beta_c.append((n, s))
            else:
                other_c.append((n, s))
        # peek at candidate beta tables and unclassified tables
        peeks = {}
        for n, s in beta_c + other_c + sig_c:
            if n.lower().endswith(".xlsx") or n.lower().endswith(".zip") or n.lower().endswith(".tar"):
                peeks[n] = ("unknown", 0, "archive/xlsx not peeked")
                continue
            peeks[n] = assess(peek(gse, n))
        allc = [(n, s) for n, s in beta_c + other_c + sig_c
                if not DETP.search(n) or (re.search("beta", n, re.I) and not re.search("residual", n, re.I))]
        beta_ok = [(n, s) for n, s in allc if peeks[n][0] == "beta"]
        mval = [(n, s) for n, s in allc if peeks[n][0] == "mvalue" and (n, s) not in sig_c]
        sig_ok = [(n, s) for n, s in allc if (peeks[n][0] == "signals" or (n, s) in sig_c) and (n, s) not in beta_ok]
        beta_unverified = [(n, s) for n, s in beta_c if peeks[n][0] in ("unknown", "percent", "transposed") and not DETP.search(n)]
        tar_perSample = [f for f in tar_other if f["type"].upper() in ("TXT", "CSV", "TSV")
                         and not SEQ_FILE.search(f["name"])]
        notes = []
        is_seq = ("sequencing" in stype.lower() and not idat) or (plats and all(p == "SEQ" for p in plats))
        non_inf = bool(plats) and not any(p in INFINIUM for p in plats) and not any(
            g in OTHER_ILLUMINA_METH for g in r["platforms"]) and not idat
        if is_seq:
            notes.append("sequencing platform")
        elif non_inf:
            notes.append("non-Infinium array platform (" + "; ".join(gpl_info[g][0][:60] for g in r["platforms"]) + ")")
        if any(g in ("GPL18246", "GPL21161") for g in r["platforms"]):
            notes.append("GoldenGate platform (~1.5k CpGs; clocks not applicable)")
        if r["soft_status"] != "ok":
            cls = "NOTHING_USABLE"; notes.append("GEO SOFT 404")
        elif is_seq or non_inf:
            cls = "NON_ARRAY"
        elif idat:
            cls = "IDAT_AVAILABLE"
        elif beta_ok or mval or beta_unverified:
            cls = "BETA_TABLE_AVAILABLE"
            if not beta_ok and mval:
                notes.append("processed table appears to be M-values (convert to beta)")
            if not beta_ok and not mval:
                notes.append("processed table not verified by peek")
        elif sig_ok:
            cls = "SIGNALS_ONLY"
        elif tar_perSample:
            cls = "BETA_TABLE_AVAILABLE"
            notes.append(f"per-sample {tar_perSample[0]['type']} files inside RAW.tar (content unverified)")
        else:
            cls = "NOTHING_USABLE"
        if files is None:
            notes.append("no suppl/ directory on GEO FTP")
        if any(re.search("residual", n, re.I) for n, _ in files):
            notes.append("residual tables present (not usable as beta)")

        # best route: smallest usable download
        routes = []
        if beta_ok:
            routes.append(("beta_table", max(s for n, s in beta_ok)))
            if len(beta_ok) > 1:
                notes.append(f"{len(beta_ok)} beta tables (may be split by subset/normalization; all together {sum(s for n, s in beta_ok)} bytes); route counts largest")
        elif mval:
            routes.append(("mvalue_table", max(s for n, s in mval)))
            if len(mval) > 1:
                notes.append(f"{len(mval)} M-value tables (all together {sum(s for n, s in mval)} bytes); route counts largest")
        elif beta_unverified and not idat:
            routes.append(("processed_table_unverified", min(s for n, s in beta_unverified)))
        if idat:
            routes.append(("idat_raw_tar", tar_size))
        if sig_ok and not routes:
            routes.append(("signals_table", min(s for n, s in sig_ok)))
        if cls == "BETA_TABLE_AVAILABLE" and not routes and tar_perSample:
            routes.append(("raw_tar_per_sample_tables", tar_size))
        best = min(routes, key=lambda x: x[1]) if routes else ("none", 0)
        # manual overrides after reading filelist.txt (see report)
        if gse == "GSE40700":
            cls = "BETA_TABLE_AVAILABLE"; best = ("raw_tar_per_sample_tables", tar_size)
            notes.append("manual: RAW.tar holds 73 ENCODE HAIB per-sample 27K BED files (score field, beta*1000 by ENCODE convention; unverified)")
        if gse == "GSE75268":
            best = ("raw_tar_per_sample_tables", tar_size)
            notes.append("manual: RAW.tar holds TCGA 27K level-1/level-2 per-sample TXT (signal-level; no level-3 beta)")
        if beta_unverified and idat and not beta_ok and not mval:
            notes.append("processed table unverified by peek; IDAT route chosen")
        if cls in ("NON_ARRAY", "NOTHING_USABLE"):
            best = ("none", 0)
        if len(set(plats)) > 1 or len(plats) > 1:
            notes.append("multi-platform series: " + ",".join(r["platforms"]))
        if not r["platforms"] and r["soft_status"] == "ok":
            notes.append("SuperSeries or no platform in header")
        if any("SuperSeries of" in x for x in r["relations"]):
            notes.append("SuperSeries")
        n_geo = r["n_samples_geo"]
        local = r["local_sample_columns"]
        if n_geo != local:
            notes.append(f"local sample columns {local} != GEO samples {n_geo}")
        if idat:
            notes.append(f"{len(idat)} IDAT files in RAW.tar ({idat_bytes} bytes)")
        if mval:
            notes.append("M-value table: " + ";".join(n for n, _ in mval))
        if best[0] != "none" and best[0] != "idat_raw_tar":
            bn = [n for n, sz in files if sz == best[1] and n in peeks]
            if bn:
                notes.append(f"best file {bn[0]} peek: {peeks[bn[0]][2]}")
        file_desc = []
        for n, s in files:
            if n == "filelist.txt":
                continue
            k = peeks.get(n, ("", 0, ""))[0]
            file_desc.append(f"{n}[{s}{',' + k if k else ''}]")
        total = sum(s for n, s in files)
        rows.append({
            "gse": gse, "platform": ";".join(f"{g}({p})" for g, p in zip(r["platforms"], plats)),
            "platform_class": "+".join(sorted(set(plats))) or "none",
            "n_samples": n_geo, "local_sample_columns": local, "class": cls,
            "best_route": best[0], "best_route_bytes": best[1],
            "idat_raw_tar_bytes": tar_size if idat else 0,
            "table_linked_32": gse in LINKED, "historical_table_rows": ATLAS_ROWS.get(gse, 0),
            "series_type": stype, "files": " | ".join(file_desc), "total_bytes": total,
            "notes": "; ".join(notes)})
    rows.sort(key=lambda x: x["gse"])
    fields = ["gse", "platform", "platform_class", "n_samples", "local_sample_columns", "class", "best_route",
              "best_route_bytes", "idat_raw_tar_bytes", "table_linked_32", "historical_table_rows",
              "series_type", "files", "total_bytes", "notes"]
    with open(os.path.join(OUT, f"{NAME}_inventory_2026-09-25.csv" if NAME == "empty_dnam_geo_repair" else f"{NAME}_inventory.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    def summarize(sub):
        return {
            "n": len(sub),
            "by_class": dict(Counter(r["class"] for r in sub)),
            "by_platform": dict(Counter(r["platform_class"] for r in sub)),
            "by_class_platform": {f"{k[0]}|{k[1]}": v for k, v in
                                  Counter((r["class"], r["platform_class"]) for r in sub).items()},
            "by_best_route": dict(Counter(r["best_route"] for r in sub)),
            "best_route_total_bytes": sum(r["best_route_bytes"] for r in sub),
            "idat_route_total_bytes": sum(r["idat_raw_tar_bytes"] for r in sub),
            "all_suppl_total_bytes": sum(r["total_bytes"] for r in sub),
            "sample_count_mismatch": [
                {"gse": r["gse"], "local": r["local_sample_columns"], "geo": r["n_samples"]}
                for r in sub if r["local_sample_columns"] != r["n_samples"]],
        }
    summ = {"created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "source_list": SRC,
            "source_list_sha256": "64f77a5f222f607c813e2e3d5554b3194f44cada815b4a70da4ee328ffe558fc",
            "linked_list_source": str(EVID / "empty-input-impact-2026-09-18/run-001/summary.json"),
            "linked_list_sha256": "7dc74f9fec68325581164d314a12b32e3b4e405589106c4d8d401e7e546618ae",
            "method": "GEO SOFT brief header, FTP LIST of suppl/, filelist.txt for RAW.tar, 512 KiB range peek of "
                      "candidate tables; best route = smallest usable download (beta table preferred, else IDAT RAW.tar, else signals)",
            "gpl_titles": {g: {"title": t, "technology": te} for g, (t, te) in gpl_info.items()},
            "all_513": summarize(rows),
            "table_linked_32": summarize([r for r in rows if r["table_linked_32"]])}
    json.dump(summ, open(os.path.join(OUT, f"{NAME}_summary_2026-09-25.json" if NAME == "empty_dnam_geo_repair" else f"{NAME}_summary.json"), "w"), indent=1)
    print(json.dumps({k: summ[k] for k in ("all_513", "table_linked_32")}, indent=1)[:6000])


if __name__ == "__main__":
    main()
