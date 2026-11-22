#!/usr/bin/env python3
"""Read-only GEO supplementary-file inventory for the 513 empty local DNAm matrices.

Source list: audit_system.paths.EMPTY_MATRIX_LIST (the empty-matrix census in the rebuild bundle).
Queries per GSE (<=3 req/s, retry with backoff):
  1. GEO SOFT brief header (acc.cgi ... form=text view=brief)
  2. FTP LIST of series suppl/ directory (exact byte sizes)
  3. suppl/filelist.txt when a RAW.tar is present (tar contents)
Raw responses cached under raw/ for provenance. No data files downloaded.
"""
import ftplib, json, os, re, sys, time, urllib.request, urllib.error
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from audit_system import paths

# INVENTORY_OUT and INVENTORY_SOURCE select another study list (e.g. studies with detached input profiles);
# the defaults reproduce the 2026-09-25 empty-matrix inventory.
OUT = os.environ.get("INVENTORY_OUT", str(paths.VERIFIER_ROOT / "empty-repair-inventory"))
RAW = os.path.join(OUT, "raw")
os.makedirs(RAW, exist_ok=True)
SRC = os.environ.get("INVENTORY_SOURCE", str(paths.EMPTY_MATRIX_LIST))
MIN_INTERVAL = 0.34
_last = [0.0]


def throttle():
    dt = time.time() - _last[0]
    if dt < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - dt)
    _last[0] = time.time()


def http_get(url, cache):
    p = os.path.join(RAW, cache)
    if os.path.exists(p):
        return open(p, encoding="utf-8", errors="replace").read()
    for attempt in range(6):
        throttle()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "clockbase-audit/1.0 (kying)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                txt = r.read().decode("utf-8", errors="replace")
            open(p, "w").write(txt)
            return txt
        except urllib.error.HTTPError as e:
            if e.code == 404:
                open(p, "w").write("__404__")
                return "__404__"
            err = e
        except Exception as e:  # network errors
            err = e
        time.sleep(2 ** attempt)
    raise RuntimeError(f"GET failed {url}: {err}")


class FTP:
    def __init__(self):
        self.f = None

    def connect(self):
        self.f = ftplib.FTP("ftp.ncbi.nlm.nih.gov", timeout=60)
        self.f.login()

    def listdir(self, path, cache):
        p = os.path.join(RAW, cache)
        if os.path.exists(p):
            return json.load(open(p))
        for attempt in range(6):
            throttle()
            try:
                if self.f is None:
                    self.connect()
                lines = []
                try:
                    self.f.retrlines("LIST " + path, lines.append)
                except (ftplib.error_perm, ftplib.error_temp) as e:
                    if "550" in str(e) or "No such file" in str(e):
                        lines = ["__missing__"]
                    else:
                        raise
                json.dump(lines, open(p, "w"))
                return lines
            except Exception as e:
                err = e
                try:
                    self.f.close()
                except Exception:
                    pass
                self.f = None
                time.sleep(2 ** attempt)
        raise RuntimeError(f"FTP LIST failed {path}: {err}")


def stub(gse):
    return gse[:-3] + "nnn" if len(gse) > 6 else "GSEnnn"


def parse_list(lines):
    files = []
    for ln in lines:
        if ln == "__missing__":
            return None
        parts = ln.split(None, 8)
        if len(parts) < 9:
            continue
        if parts[0].startswith("d"):
            files.append((parts[8] + "/", 0))
            continue
        files.append((parts[8], int(parts[4])))
    return files


def main():
    recs = [json.loads(l) for l in open(SRC)]
    gses = [r["study_id"] for r in recs]
    ftp = FTP()
    out = []
    for i, r in enumerate(recs):
        gse = r["study_id"]
        soft = http_get(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=self&form=text&view=brief",
                        f"{gse}.soft.txt")
        platforms = re.findall(r"^!Series_platform_id = (\S+)", soft, re.M)
        samples = re.findall(r"^!Series_sample_id = (\S+)", soft, re.M)
        types = re.findall(r"^!Series_type = (.+)$", soft, re.M)
        supp_urls = re.findall(r"^!Series_supplementary_file = (.+)$", soft, re.M)
        relations = re.findall(r"^!Series_relation = (.+)$", soft, re.M)
        lines = ftp.listdir(f"/geo/series/{stub(gse)}/{gse}/suppl/", f"{gse}.suppl_list.json")
        files = parse_list(lines)
        filelist = None
        if files and any(n.upper().endswith("_RAW.TAR") for n, _ in files):
            fl = http_get(f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub(gse)}/{gse}/suppl/filelist.txt",
                          f"{gse}.filelist.txt")
            if fl != "__404__":
                filelist = []
                for ln in fl.splitlines()[1:]:
                    c = ln.split("\t")
                    if len(c) >= 5:
                        filelist.append({"kind": c[0], "name": c[1], "size": int(c[3]) if c[3].isdigit() else None,
                                         "type": c[4]})
        out.append({"gse": gse, "local_sample_columns": r["sample_columns"], "local_probe_rows": r["probe_rows"],
                    "local_matrix_path": r["matrix_path"], "local_matrix_sha256": r["matrix_sha256"],
                    "platforms": platforms, "n_samples_geo": len(samples), "series_type": types,
                    "series_supplementary_file": supp_urls, "relations": relations,
                    "soft_status": "404" if soft == "__404__" else "ok",
                    "suppl_files": files, "raw_tar_filelist": filelist})
        if (i + 1) % 25 == 0:
            print(f"{i+1}/{len(recs)}", flush=True)
    json.dump(out, open(os.path.join(OUT, "geo_raw_inventory.json"), "w"), indent=1)
    print("done", len(out))


if __name__ == "__main__":
    main()
