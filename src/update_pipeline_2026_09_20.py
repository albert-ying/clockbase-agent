"""Versioned, bounded multi-study ingestion around the existing ClockBase pipeline.

Source acquisition, operational ingestion and scientific acceptance are separate.
An original-paper temporal extension remains quarantined until its corpus and
family exclusions are established elsewhere. No date is inferred from an ID,
file timestamp, download time or ARCHS4 release timestamp.
"""
import argparse
import ast
import contextlib
import csv
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import importlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "clockbase-continuous-update-plan-v1"
STAGES = ("metadata", "input", "design", "preprocess", "clocks", "agent", "audit", "report")
ROLES = {"series", "samples", "context", "matrix", "matrix_provenance", "archs4", "dnam_reference", "dnam_supplement", "dnam_idat_provenance"}
MATRIX_PROVENANCE_SCHEMA = "clockbase-archs4-derived-matrix-v1"
BATCH_HANDOFF_SCHEMA = "clockbase-archs4-update-handoff-plan-v1"
REFRESH_SCHEMA = "clockbase-future-refresh-v1"
REFRESH_VERIFIER = "src/verify_ingestion_2026_09_20.py"
# The released copy is byte-identical; retain the evidence name for saved manifests.
LEGACY_REFRESH_VERIFIER = "paper/revision_2026-09-07/evidence/continuous-audit-pipeline-2026-09-20/update-controller/verify_ingestion_2026_09_20.py"
MAX_PROVENANCE_FILE_BYTES = 16 * 1024 ** 2
MODALITIES = {"mouse_rna", "human_rna", "human_dnam"}
RNA_SPECIES = {"human_rna": "human", "mouse_rna": "mouse"}
DNAM_IDAT_ROUTE = "source_declared_idat_v2"
PIPELINE_OPTIONS = {"model", "audit_model", "max_calls", "max_cost_usd", "input_usd_per_million",
                    "output_usd_per_million", "audit_input_usd_per_million", "audit_output_usd_per_million",
                    "max_output_tokens", "max_prompt_bytes", "audit_validation_retries",
                    "coding_validation_retries", "coding_completion_policy", "rna_preprocessing", "source_contract_policy",
                    "service_tier"}
DEFAULT_OPTIONS = {"model": "gpt-5.6-luna", "audit_model": "gpt-5.6-luna", "max_calls": 8,
                   "max_cost_usd": 1.0, "input_usd_per_million": 0.5, "output_usd_per_million": 1.8,
                   "audit_input_usd_per_million": 0.5, "audit_output_usd_per_million": 1.8,
                   "max_output_tokens": 12000, "max_prompt_bytes": 393216,
                   "audit_validation_retries": 0, "coding_validation_retries": 0,
                   "rna_preprocessing": "clockbase_legacy",
                   "source_contract_policy": "strict"}
IMPLEMENTATIONS = ("src/update_pipeline_2026_09_20.py", "src/pipeline_2026_09_10.py",
                   "src/pipeline_inputs_2026_09_10.py", "src/pipeline_runtime_2026_09_10.py",
                   "src/coding_agent.py", "audit_system/source_packet.py",
                   "audit_system/report_2026_09_10.py",
                   "audit_system/agent_workflow_2026_09_10.py", "audit_system/extension_eligibility_2026_09_10.py",
                   "audit_system/statistical_rerun_2026_09_20.py", "audit_system/archs4_update_2026_09_10.py",
                   "pipeline_assets/rna_4_6_2026_09_10/manifest.json")


class UpdateIntegrityError(ValueError):
    pass


def sha(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""): result.update(block)
    return result.hexdigest()


def identity(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def implementation_hashes():
    """Bind the post-analysis auditor's declared dependencies without importing it."""
    names = set(IMPLEMENTATIONS)
    hook = ROOT / "audit_system/corpus_audit_2026_09_20.py"
    if hook.is_file():
        names.add(str(hook.relative_to(ROOT)))
        tree = ast.parse(hook.read_text())
        declarations = [node.value for node in tree.body if isinstance(node, ast.Assign)
                        and any(isinstance(target, ast.Name) and target.id == "CODE_FILES" for target in node.targets)]
        if len(declarations) != 1: raise UpdateIntegrityError("Comprehensive auditor must declare one CODE_FILES tuple")
        dependencies = ast.literal_eval(declarations[0])
        if not isinstance(dependencies, (list, tuple)) or not all(isinstance(name, str) for name in dependencies):
            raise UpdateIntegrityError("Invalid comprehensive audit dependency declaration")
        names.update(dependencies)
    for name in names:
        if Path(name).is_absolute() or not (ROOT / name).resolve().is_relative_to(ROOT.resolve()):
            raise UpdateIntegrityError("Audit dependency escapes repository: " + name)
    return {name: sha(ROOT / name) for name in sorted(names)}


def _implementation_source(name):
    return ROOT / (REFRESH_VERIFIER if name == LEGACY_REFRESH_VERIFIER else name)


def _refresh_binding(value):
    """Compare verifier locations as aliases without changing sealed historical records."""
    value = dict(value)
    declarations = dict(value["implementation"])
    if LEGACY_REFRESH_VERIFIER in declarations:
        checksum = declarations.pop(LEGACY_REFRESH_VERIFIER)
        if REFRESH_VERIFIER in declarations and declarations[REFRESH_VERIFIER] != checksum:
            raise UpdateIntegrityError("Conflicting refresh verifier aliases")
        declarations[REFRESH_VERIFIER] = checksum
    value["implementation"] = declarations
    return value


def verify_implementation(declarations):
    for name, expected in declarations.items():
        if sha(_implementation_source(name)) != expected: raise UpdateIntegrityError("Child implementation changed: " + name)


def snapshot_implementation(version_root, declarations):
    """Preserve exact declared sources even when a dependency is uncommitted."""
    directory = version_root / "implementation"
    directory.mkdir(exist_ok=False)
    for name, expected in declarations.items():
        source = _implementation_source(name)
        if sha(source) != expected: raise UpdateIntegrityError("Implementation changed before snapshot: " + name)
        target = directory / name
        target.parent.mkdir(parents=True, exist_ok=True)
        with source.open("rb") as source_file, target.open("xb") as target_file:
            shutil.copyfileobj(source_file, target_file)
        if sha(target) != expected: raise UpdateIntegrityError("Implementation changed during snapshot: " + name)


def verify_implementation_snapshot(version_root, declarations):
    for name, expected in declarations.items():
        path = version_root / "implementation" / name
        if not path.is_file() or sha(path) != expected:
            raise UpdateIntegrityError("Frozen implementation snapshot changed or is missing: " + name)


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False); handle.write("\n")


def seal(path, value):
    write_new(path, value)
    with path.with_suffix(path.suffix + ".sha256").open("x") as handle: handle.write(sha(path) + "\n")


def read_sealed(path):
    if not path.with_suffix(path.suffix + ".sha256").is_file() or sha(path) != path.with_suffix(path.suffix + ".sha256").read_text().strip():
        raise UpdateIntegrityError("Receipt changed or lacks its digest: " + str(path))
    return json.loads(path.read_text())


def number(value, name, *, integer=False, minimum=0):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < minimum:
        raise ValueError("Invalid nonnegative resource limit: " + name)
    if integer and type(value) is not int: raise ValueError("Resource limit must be an integer: " + name)
    return value


def file_source(declaration, destination, *, max_bytes, network):
    """Acquire one bounded source; local inputs require their declared digest."""
    if not isinstance(declaration, dict) or set(declaration) - {"path", "url", "sha256", "max_bytes"}:
        raise ValueError("Unsupported source declaration")
    if bool(declaration.get("path")) == bool(declaration.get("url")):
        raise ValueError("A source requires exactly one path or URL")
    cap = number(declaration.get("max_bytes", max_bytes), "source.max_bytes", integer=True, minimum=1)
    cap = min(cap, max_bytes)
    expected = declaration.get("sha256")
    if expected is not None and not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ValueError("Invalid source SHA256")
    if declaration.get("path"):
        if expected is None: raise ValueError("Local source requires its SHA256")
        path = Path(declaration["path"]).expanduser().resolve()
        if not path.is_file(): raise FileNotFoundError("Cannot proceed - required data file not found: " + str(path))
        if path.stat().st_size > cap: raise ValueError("Local input exceeds declared byte cap")
    else:
        if not network: raise ValueError("Network acquisition disabled in this update plan")
        url = declaration["url"]
        parsed = urllib.parse.urlparse(url)
        if parsed.scheme != "https" or parsed.username or parsed.password or parsed.hostname not in {
            "www.ncbi.nlm.nih.gov", "ftp.ncbi.nlm.nih.gov", "www.ebi.ac.uk", "s3.k8s.maayanlab.cloud", "archs4.org"}:
            raise ValueError("Source URL must name an allowed public-data host over HTTPS")
        suffixes = "".join(Path(parsed.path).suffixes)
        path = destination / (identity(declaration) + suffixes)
        record = path.with_name(path.name + ".source.json")
        if record.exists():
            saved = read_sealed(record)
            if saved["declaration"] != declaration or sha(path) != saved["sha256"]:
                raise UpdateIntegrityError("Cached source acquisition changed")
            return saved
        destination.mkdir(parents=True, exist_ok=True)
        if path.exists() or path.with_name(path.name + ".partial").exists():
            raise UpdateIntegrityError("Unfinished acquisition retained; no automatic overwrite")
        partial = path.with_name(path.name + ".partial")
        with urllib.request.urlopen(url, timeout=60) as response, partial.open("xb") as output:
            final_url = urllib.parse.urlparse(response.url)
            if final_url.scheme != "https" or final_url.hostname != parsed.hostname:
                raise ValueError("Unapproved source redirect")
            used = 0; started = time.monotonic()
            while True:
                block = response.read(min(1 << 20, cap - used + 1))
                if not block: break
                used += len(block)
                if used > cap: raise ValueError("Downloaded source exceeds byte cap")
                if time.monotonic() - started > 180: raise TimeoutError("Source download wall-time bound exceeded")
                output.write(block)
        partial.rename(path)
    before = path.stat()
    checksum = sha(path)
    if expected and checksum != expected: raise UpdateIntegrityError("Source SHA256 differs: " + str(path))
    if (before.st_size, before.st_mtime_ns) != (path.stat().st_size, path.stat().st_mtime_ns):
        raise UpdateIntegrityError("Source changed while hashing")
    with path.open("rb") as handle:
        if handle.read(100).startswith(b"version https://git-lfs.github.com/spec/"):
            raise ValueError("Source is an unmaterialized LFS pointer")
    result = {"path": str(path), "sha256": checksum, "bytes": before.st_size,
              "mtime_ns": before.st_mtime_ns, "declaration": declaration}
    if declaration.get("url"): seal(record, result)
    return result


def verify_sources(sources):
    for role, record in sources.items():
        if sha(record["path"]) != record["sha256"]:
            raise UpdateIntegrityError("Frozen source changed: " + role)


def read_dnam_idat_provenance(path, expected_sha256, study_id, *, sources=None, matrix_path=None):
    """Validate the DNAm producer graph and expose every file for version binding."""
    from audit_system.report_2026_09_10 import verify_dnam_idat_provenance

    path = Path(path).expanduser().resolve()
    if sha(path) != expected_sha256:
        raise UpdateIntegrityError("DNAm IDAT provenance SHA256 differs")
    if sources is not None:
        if not {"matrix", "series", "samples"} <= set(sources):
            raise ValueError("DNAm IDAT provenance requires matrix, Series and sample source bindings")
        matrix_path = Path(sources["matrix"]["path"])
    result = verify_dnam_idat_provenance(path, matrix_path=matrix_path)
    if (not isinstance(result, dict) or result.get("study_id") != study_id
            or result.get("scientific_status") != "CANNOT VERIFY"
            or not isinstance(result.get("measurement"), dict)):
        raise ValueError("DNAm IDAT provenance study or scientific scope differs")
    sample_ids = result.get("sample_ids")
    if (not isinstance(sample_ids, list) or not sample_ids
            or any(not isinstance(sample, str) or not re.fullmatch(r"GSM[0-9]+", sample) for sample in sample_ids)
            or len(sample_ids) != len(set(sample_ids))):
        raise ValueError("DNAm IDAT provenance has invalid or duplicate sample identities")
    declared = result.get("files")
    if not isinstance(declared, dict) or not {"matrix", "series", "samples"} <= set(declared):
        raise ValueError("DNAm IDAT provenance lacks its transitive source inventory")
    files = {}
    for role, record in declared.items():
        if (not isinstance(role, str) or not role or not isinstance(record, dict)
                or set(record) != {"path", "sha256"} or not isinstance(record["path"], str)
                or not Path(record["path"]).is_absolute() or not isinstance(record["sha256"], str)
                or not re.fullmatch(r"[0-9a-f]{64}", record["sha256"])):
            raise ValueError("DNAm IDAT provenance requires canonical path/SHA256 file declarations")
        source = Path(record["path"]).resolve()
        if str(source) != record["path"] or sha(source) != record["sha256"]:
            raise UpdateIntegrityError("DNAm IDAT provenance nested file differs: " + role)
        files[role] = dict(record)
    wrapper = {"path": str(path), "sha256": expected_sha256}
    if "wrapper" in files and files["wrapper"] != wrapper:
        raise ValueError("DNAm IDAT provenance wrapper declaration differs")
    files["wrapper"] = wrapper
    if sources is not None:
        for role in ("matrix", "series", "samples"):
            source = sources[role]
            if files[role] != {"path": str(Path(source["path"]).resolve()), "sha256": source["sha256"]}:
                raise ValueError("DNAm IDAT provenance source role differs: " + role)
    if sha(path) != expected_sha256:
        raise UpdateIntegrityError("DNAm IDAT provenance changed during verification")
    return {**result, "files": files}


def read_matrix_provenance(path, expected_sha256, sources, study_id, modality):
    """Bind a historical extraction/readback graph without opening its HDF5 origin.

    The returned files are small manifests and producing implementations. The
    matrix and GEO sources are verified against the supplied source bindings.
    No batch receipt or PASS label establishes clock or scientific validity.
    """
    files = {}
    path = Path(path).expanduser().resolve()

    def require(condition, message):
        if not condition:
            raise UpdateIntegrityError("Derived matrix provenance: " + message)

    def declaration(value, owner):
        require(isinstance(value, dict) and isinstance(value.get("path"), str)
                and isinstance(value.get("sha256"), str) and re.fullmatch(r"[0-9a-f]{64}", value["sha256"])
                and not value.get("url"), "requires a local file path and SHA256")
        candidate = Path(value["path"]).expanduser()
        require(not urllib.parse.urlparse(value["path"]).scheme, "URL declarations are unsupported")
        candidate = (owner / candidate).resolve() if not candidate.is_absolute() else candidate.resolve()
        return {"path": str(candidate), "sha256": value["sha256"]}

    def bound(role, value, owner):
        item = declaration(value, owner)
        source = Path(item["path"])
        require(source.is_file(), "missing nested file: " + role)
        require(source.stat().st_size <= MAX_PROVENANCE_FILE_BYTES, "nested file exceeds small-receipt cap: " + role)
        require(sha(source) == item["sha256"], "nested SHA256 differs: " + role)
        require(role not in files or files[role] == item, "conflicting nested file role: " + role)
        files[role] = item
        return source

    def document(role, value, owner):
        source = bound(role, value, owner)
        try:
            result = json.loads(source.read_text())
        except (OSError, json.JSONDecodeError) as error:
            raise UpdateIntegrityError("Derived matrix provenance: cannot read " + role + ": " + str(error)) from error
        require(isinstance(result, dict), "receipt is not an object: " + role)
        return result, source

    def selected(records, key):
        require(isinstance(records, list) and all(isinstance(r, dict) for r in records), "invalid " + key + " inventory")
        ids = [r.get("study_id") for r in records]
        require(all(isinstance(x, str) for x in ids) and len(ids) == len(set(ids)), "ambiguous " + key + " study inventory")
        matches = [r for r in records if r.get("study_id") == study_id]
        require(len(matches) == 1 and matches[0].get("modality") == modality, "foreign or absent " + key + " study")
        return matches[0]

    wrapper, _ = document("wrapper", {"path": str(path), "sha256": expected_sha256}, path.parent)
    required = {"schema", "study_id", "modality", "matrix", "series", "samples", "batch_summary", "extraction_receipt",
                "verification_receipt", "source_profile_descriptor", "source_profile_script", "source_profile_record_sha256",
                "source_hdf5", "implementation", "scientific_status", "original_paper_2025_eligibility"}
    require(set(wrapper) == required and wrapper["schema"] == MATRIX_PROVENANCE_SCHEMA, "unsupported wrapper schema or fields")
    require(all(isinstance(wrapper[role], dict) and isinstance(wrapper[role].get("path"), str)
                and Path(wrapper[role]["path"]).is_absolute() for role in (
                    "matrix", "series", "samples", "batch_summary", "extraction_receipt", "verification_receipt",
                    "source_profile_descriptor", "source_profile_script")), "wrapper file declarations must be absolute")
    require(wrapper["study_id"] == study_id and wrapper["modality"] == modality and modality in {"human_rna", "mouse_rna"},
            "study or RNA modality differs")
    require("matrix" in sources and "archs4" not in sources, "claim requires one declared matrix, without a live ARCHS4 source")
    require(wrapper["scientific_status"] == "CANNOT VERIFY" and wrapper["original_paper_2025_eligibility"] == "QUARANTINED",
            "scientific or temporal acceptance is outside this receipt scope")
    for role in ("matrix", "series", "samples"):
        require(role in sources, "missing source binding: " + role)
        item = declaration(wrapper[role], path.parent)
        expected = {"path": str(Path(sources[role]["path"]).resolve()), "sha256": sources[role]["sha256"]}
        require(item == expected and sha(item["path"]) == item["sha256"], "declared " + role + " differs from the source bytes")

    batch, batch_path = document("batch_summary", wrapper["batch_summary"], path.parent)
    extraction, extraction_path = document("extraction_receipt", wrapper["extraction_receipt"], path.parent)
    verification, verification_path = document("verification_receipt", wrapper["verification_receipt"], path.parent)
    profile, profile_path = document("source_profile_descriptor", wrapper["source_profile_descriptor"], path.parent)
    require(batch.get("schema") == "clockbase-archs4-batch-extraction-summary-v1" and batch.get("path_base") == "batch_root",
            "unsupported batch summary schema or path base")
    require(batch.get("execution_status") == "COMPLETED" and batch.get("status") in {"COMPLETED", "COMPLETED_WITH_UNRESOLVED"},
            "batch extraction is not completed")
    require(all(type(batch.get(name)) is int and batch[name] == 0 for name in ("model_api_calls", "preprocessing_runs", "clock_runs")),
            "batch scope is not extraction-only")
    require(extraction.get("schema") == "clockbase-archs4-batch-study-extraction-v1"
            and extraction.get("study_id") == study_id and extraction.get("modality") == modality
            and extraction.get("species") == {"human_rna": "human", "mouse_rna": "mouse"}[modality]
            and extraction.get("status") == "PASS" and extraction.get("input_mode") == "FRESH_ARCHS4_BATCH_EXTRACTION"
            and extraction.get("path_base") == "batch_root",
            "unsupported or foreign per-study extraction receipt")
    require(verification.get("schema") == "clockbase-archs4-independent-every-cell-v1"
            and verification.get("status") == "PASS"
            and all(type(verification.get(name)) is int and verification[name] == 0 for name in ("model_calls", "preprocessing_runs", "clock_runs")),
            "independent every-cell verification is absent, unresolved or outside numerical scope")
    require(profile.get("schema") == "clockbase-rna-source-profile-census-v1", "unsupported source profile schema")
    expected_species = RNA_SPECIES[modality]
    # Old human receipts predate the additive fields. Mouse has no such
    # historical interpretation: every linked document must declare mouse.
    require(profile.get("modality", "human_rna") == modality
            and batch.get("modality", "human_rna") == modality
            and verification.get("modality", "human_rna") == modality
            and batch.get("source_hdf5", {}).get("species", "human") == expected_species
            and verification.get("species", "human") == expected_species,
            "batch/profile/readback modality or species differs")
    require(all(doc.get("scientific_status") == "CANNOT VERIFY" for doc in (batch, extraction, verification, profile)),
            "nested scientific acceptance is outside this receipt scope")
    batch_row = selected(batch.get("results"), "batch")
    checked_row = selected(verification.get("results"), "verification")
    profile_row = selected(profile.get("studies"), "profile")
    require({r["study_id"] for r in batch["results"]} == {r["study_id"] for r in verification["results"]}
            == {r["study_id"] for r in profile["studies"]}, "batch, profile and verification study denominators differ")
    require(study_id in profile.get("supported_study_ids", []) and study_id not in profile.get("stopped_study_ids", [])
            and profile_row.get("profile_status") == batch_row.get("source_profile_status") == "SUPPORTED_SOURCE_PROFILE",
            "study is not a supported source profile")
    profile_hash = identity(profile_row)
    require(profile_hash == wrapper["source_profile_record_sha256"] == batch_row.get("source_profile_record_sha256")
            == extraction.get("source", {}).get("source_profile_record_sha256"), "source profile record differs")
    config, _ = document("batch_config", {"path": "config.json", "sha256": batch["config_sha256"]}, batch_path.parent)
    require(config.get("schema") == "clockbase-archs4-batch-extraction-config-v1"
            and config.get("plan_sha256") == batch.get("plan_sha256")
            and config.get("modality", "human_rna") == modality, "batch configuration identity differs")
    plan, plan_path = document("batch_plan", {"path": config["plan_path"], "sha256": config["plan_sha256"]}, batch_path.parent)
    require(plan.get("schema") == "clockbase-archs4-batch-extraction-plan-v1"
            and plan.get("modality", "human_rna") == modality
            and plan.get("selected_study_ids") == config.get("selected_study_ids") == profile.get("supported_study_ids"),
            "batch plan source selection differs")
    require(isinstance(plan.get("output"), str) and Path(plan["output"]).is_absolute()
            and Path(plan["output"]).expanduser().resolve() == batch_path.parent,
            "batch plan output does not identify this archived batch")
    require(declaration(batch["source_profile_descriptor"], batch_path.parent) == files["source_profile_descriptor"]
            == declaration(config["source_profile_descriptor"], batch_path.parent)
            == declaration(plan["source_profile_descriptor"], plan_path.parent)
            and verification.get("source_profile_sha256") == files["source_profile_descriptor"]["sha256"]
            == extraction["source"].get("source_profile_descriptor_sha256"), "source profile descriptor differs")
    profile_script = bound("source_profile_script", wrapper["source_profile_script"], path.parent)
    require(sha(profile_script) == profile.get("verification_script_sha256"), "source profile producing script differs")
    for role in ("series", "samples"):
        expected = declaration(wrapper[role], path.parent)
        require(declaration(batch_row[role], batch_path.parent) == expected
                and declaration(profile_row["sources"][role], profile_path.parent) == expected
                and extraction["source"].get(role + "_sha256") == expected["sha256"], "foreign " + role + " source lineage")
    require(batch_row.get("status") == "EXTRACTED" and checked_row.get("status") == "PASS", "study extraction or verification is unresolved")
    require(declaration(batch_row["matrix"], batch_path.parent) == declaration(wrapper["matrix"], path.parent)
            and declaration(extraction["matrix"], batch_path.parent) == declaration(wrapper["matrix"], path.parent)
            and checked_row.get("matrix_sha256") == sources["matrix"]["sha256"], "foreign or changed matrix lineage")
    require(declaration(batch_row["extraction_receipt"], batch_path.parent) == files["extraction_receipt"]
            and checked_row.get("extraction_receipt_sha256") == files["extraction_receipt"]["sha256"]
            and verification.get("batch_summary_sha256") == files["batch_summary"]["sha256"], "extraction/verification batch lineage differs")

    origin = wrapper["source_hdf5"]
    require(isinstance(origin, dict) and set(origin) in ({"sha256", "size_bytes", "hash_scope"},
                                                         {"sha256", "size_bytes", "hash_scope", "species"})
            and origin.get("species", "human") == expected_species
            and isinstance(origin["sha256"], str) and re.fullmatch(r"[0-9a-f]{64}", origin["sha256"])
            and type(origin["size_bytes"]) is int and origin["size_bytes"] > 0
            and origin["hash_scope"] == "historical_batch_pre_and_post", "invalid historical HDF5 origin scope")
    require(origin["sha256"] == batch["source_hdf5"].get("sha256") == extraction["source"].get("release_sha256")
            == verification.get("release_sha256") == verification.get("release_sha256_before") == verification.get("release_sha256_after")
            and origin["size_bytes"] == batch["source_hdf5"].get("size_bytes") == extraction["source"].get("release_size_bytes"),
            "historical HDF5 origin differs across receipts")
    require(isinstance(plan.get("release"), dict) and isinstance(config.get("release"), dict)
            and declaration(plan["release"], plan_path.parent) == declaration(config["release"], batch_path.parent)
            and {k: v for k, v in plan["release"].items() if k != "path"} == {k: v for k, v in config["release"].items() if k != "path"}
            and config["release"].get("sha256") == origin["sha256"] and config["release"].get("size_bytes") == origin["size_bytes"]
            and config["release"].get("species", "human") == plan["release"].get("species", "human") == expected_species
            and all(isinstance(plan.get(key), dict) and plan[key] == config.get(key) for key in ("runtime", "limits")),
            "plan/config release, runtime or limits differ from extraction")
    release_checks = []
    for role in ("release_before", "release_after"):
        released, _ = document(role, batch[role], batch_path.parent)
        require(released.get("schema") == "clockbase-archs4-batch-release-check-v1"
                and released.get("phase") == role.removeprefix("release_")
                and released.get("sha256") == released.get("expected_sha256") == origin["sha256"]
                and released.get("size_bytes") == origin["size_bytes"]
                and released.get("status") == "PASS" and released.get("hash_scope") == "FULL_HDF5",
                "historical full-hash release receipt differs: " + role)
        state = released.get("stat")
        require(isinstance(state, dict) and set(state) == {"device", "inode", "size_bytes", "mtime_ns", "ctime_ns"}
                and all(type(value) is int and value >= 0 for value in state.values())
                and state["size_bytes"] == origin["size_bytes"], "historical release stat receipt is incomplete: " + role)
        release_checks.append(released)
    require(isinstance(release_checks[0].get("path"), str) and Path(release_checks[0]["path"]).is_absolute()
            and release_checks[0]["stat"] == release_checks[1]["stat"] and release_checks[0]["path"] == release_checks[1].get("path"),
            "historical HDF5 identity/stat changed across extraction")
    execution_attempt, _ = document("execution_attempt", batch["execution_attempt"], batch_path.parent)
    require(execution_attempt.get("schema") == "clockbase-archs4-batch-execution-attempt-v1"
            and execution_attempt.get("status") == "COMPLETED"
            and execution_attempt.get("plan_sha256") == batch.get("plan_sha256")
            and execution_attempt.get("release_before") == batch["release_before"]
            and execution_attempt.get("release_after") == batch["release_after"]
            and execution_attempt.get("result_count") == len(batch["results"])
            and all(type(execution_attempt.get(name)) is int and execution_attempt[name] == 0
                    for name in ("model_api_calls", "preprocessing_runs", "clock_runs")),
            "batch execution attempt binding or extraction-only scope differs")
    stage, stage_path = document("stage_receipt", batch_row["stage_receipt"], batch_path.parent)
    output_members = stage.get("outputs", {})
    require(isinstance(output_members, dict), "extraction stage lacks its output inventory")
    for role in ("matrix", "extraction_receipt"):
        item = declaration(wrapper[role], path.parent)
        member = Path(item["path"])
        require(member.is_relative_to(stage_path.parent)
                and output_members.get(str(member.relative_to(stage_path.parent))) == item["sha256"],
                "file is absent from the completed extraction stage: " + role)
    implementations = wrapper["implementation"]
    require(isinstance(implementations, dict) and implementations and isinstance(batch.get("implementation"), dict)
            and isinstance(plan.get("implementation"), dict) and isinstance(config.get("implementation"), dict)
            and set(implementations) == set(batch["implementation"]) == set(plan["implementation"]) == set(config["implementation"]),
            "batch implementation inventory differs")
    for role, item in implementations.items():
        require(isinstance(role, str) and role and isinstance(item, dict) and isinstance(item.get("path"), str)
                and Path(item["path"]).is_absolute(), "implementation wrapper declarations must be absolute")
        bound("implementation:" + role, item, path.parent)
        require(declaration(batch["implementation"][role], batch_path.parent) == files["implementation:" + role],
                "batch implementation binding differs: " + role)
        require(declaration(plan["implementation"][role], plan_path.parent) == declaration(config["implementation"][role], batch_path.parent)
                and plan["implementation"][role]["sha256"] == files["implementation:" + role]["sha256"],
                "plan/config producing implementation differs from frozen batch: " + role)
    verifier = bound("verification_script", verification["verification_script"], verification_path.parent)
    require(sha(verifier) == verification.get("verification_script_sha256"), "independent verifier implementation differs")
    artifact_roles = {"extraction_receipt", "stage_receipt", "release_before", "release_after", "execution_attempt", "batch_config"}
    artifact_roles.update("implementation:" + role for role in implementations)
    artifacts = batch.get("artifacts")
    require(isinstance(artifacts, dict), "batch lacks its sealed artifact inventory")
    declared_matrix = Path(sources["matrix"]["path"]).resolve()
    require(declared_matrix.is_relative_to(batch_path.parent)
            and artifacts.get(str(declared_matrix.relative_to(batch_path.parent))) == sources["matrix"]["sha256"],
            "matrix is absent from batch artifact inventory")
    for role in artifact_roles:
        source = Path(files[role]["path"])
        require(source.is_relative_to(batch_path.parent)
                and artifacts.get(str(source.relative_to(batch_path.parent))) == files[role]["sha256"],
                "nested file is absent from batch artifact inventory: " + role)
    for role in ("batch_summary", "extraction_receipt", "stage_receipt", "release_before", "release_after", "execution_attempt", "batch_config"):
        # Per-study extraction.json is bound by its StageStore; it has no own
        # seal. StageStore and the batch JSON writer use different suffixes.
        if role == "extraction_receipt":
            continue
        document_path = Path(files[role]["path"])
        sidecar = document_path.with_suffix(".sha256") if role == "stage_receipt" else document_path.with_suffix(".json.sha256")
        require(sidecar.is_file() and sidecar.read_text().strip() == files[role]["sha256"], "missing or changed completion seal: " + role)
        bound(role + ":seal", {"path": str(sidecar), "sha256": sha(sidecar)}, sidecar.parent)
        if role != "batch_summary":
            require(artifacts.get(str(sidecar.relative_to(batch_path.parent))) == files[role + ":seal"]["sha256"],
                    "completion seal is absent from batch artifact inventory: " + role)

    matrix = extraction["matrix"]
    sample_ids = extraction.get("sample_ids_in_matrix_order")
    require(isinstance(sample_ids, list) and sample_ids and all(isinstance(s, str) and re.fullmatch(r"GSM[0-9]+", s) for s in sample_ids)
            and len(sample_ids) == len(set(sample_ids)) and set(sample_ids) == set(profile_row.get("sample_ids", []))
            and checked_row.get("sample_ids") == sample_ids, "matrix sample identity/order differs")
    features, samples = matrix.get("feature_count"), matrix.get("sample_count")
    require(type(features) is int and features > 0 and type(samples) is int and samples == len(sample_ids)
            and checked_row.get("feature_count") == features and checked_row.get("sample_count") == samples
            and checked_row.get("cells_compared") == features * samples
            and all(checked_row.get(name) is True for name in ("feature_ids_equal", "sample_ids_equal", "all_cells_equal")),
            "every-cell verification dimensions or scope differs")
    require(matrix.get("size_bytes") == Path(sources["matrix"]["path"]).stat().st_size, "matrix byte count differs")
    memberships = extraction.get("all_study_memberships")
    require(isinstance(memberships, dict) and set(memberships) == set(sample_ids)
            and all(isinstance(values, list) and study_id in values for values in memberships.values()),
            "ARCHS4 study membership differs")
    for role, item in files.items():
        require(sha(item["path"]) == item["sha256"], "nested source changed during readback: " + role)
    return {"files": files, "measurement": {
        "schema": MATRIX_PROVENANCE_SCHEMA, "study_id": study_id, "modality": modality,
        "matrix_sha256": sources["matrix"]["sha256"], "sample_ids": sample_ids,
        "feature_count": features, "sample_count": samples, "cells_compared": features * samples,
        "historical_cell_readback_status": "PASS", "source_hdf5_sha256": origin["sha256"],
        "source_hdf5_size_bytes": origin["size_bytes"], "hdf5_hash_scope": origin["hash_scope"],
        "hdf5_rehashed_for_this_child": False,
        "source_profile_record_sha256": profile_hash,
        "source_hashes": {role: files[role]["sha256"] for role in ("wrapper", "batch_summary", "extraction_receipt",
                        "verification_receipt", "source_profile_descriptor", "release_before", "release_after")},
        "scientific_status": "CANNOT VERIFY", "clock_method_validation": "NOT_ESTABLISHED",
        "original_paper_2025_eligibility": "QUARANTINED"}}



def build_batch_handoff(plan_path, output):
    """Make local update candidates from a completed, independently checked batch.

    Every enrolled study is retained. PASS applies only to the archived source
    cell readback; the handoff does not execute a clock, model or source-design
    assessment. The HDF5 origin is historical and is never opened here.
    """
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    inputs = {}

    def require(condition, message):
        if not condition: raise UpdateIntegrityError("Batch handoff: " + message)

    def declaration(value, owner, *, absolute=False):
        require(isinstance(value, dict) and isinstance(value.get("path"), str) and value["path"]
                and not value.get("url") and not urllib.parse.urlparse(value["path"]).scheme
                and isinstance(value.get("sha256"), str) and re.fullmatch(r"[0-9a-f]{64}", value["sha256"]),
                "requires a local path and SHA256")
        path = Path(value["path"]).expanduser()
        require(not absolute or path.is_absolute(), "plan sources must use absolute local paths")
        path = (owner / path).resolve() if not path.is_absolute() else path.resolve()
        require(path.suffix.lower() not in {".h5", ".hdf5", ".hdf"}, "HDF5 inputs are historical; handoff must not open them")
        return {"path": str(path), "sha256": value["sha256"]}

    def bind(value, owner, *, cap=MAX_PROVENANCE_FILE_BYTES, absolute=False):
        item = declaration(value, owner, absolute=absolute)
        path = Path(item["path"])
        require(not path.is_relative_to(output), "input must not be inside handoff output")
        require(path.is_file(), "missing bound input: " + str(path))
        require(path.stat().st_size <= cap, "bound input exceeds byte cap: " + str(path))
        require(sha(path) == item["sha256"], "bound input SHA256 differs: " + str(path))
        require(str(path) not in inputs or inputs[str(path)] == item["sha256"], "conflicting input identity")
        inputs[str(path)] = item["sha256"]
        return item

    def document(value, owner, *, absolute=False):
        item = bind(value, owner, absolute=absolute)
        try: record = json.loads(Path(item["path"]).read_text())
        except (OSError, json.JSONDecodeError) as error:
            raise UpdateIntegrityError("Batch handoff: cannot read bound JSON: " + str(error)) from error
        require(isinstance(record, dict), "bound JSON must be an object")
        return record, item

    def inventory(rows, name):
        require(isinstance(rows, list) and all(isinstance(row, dict) for row in rows), "invalid " + name + " inventory")
        ids = [row.get("study_id") for row in rows]
        require(all(isinstance(x, str) and re.fullmatch(r"GSE[0-9]+", x) for x in ids)
                and len(ids) == len(set(ids)), "ambiguous " + name + " study inventory")
        require(all(row.get("modality") in {"human_rna", "mouse_rna"} for row in rows), "non-RNA " + name + " inventory")
        return ids

    require(plan_path.is_file() and plan_path.stat().st_size <= MAX_PROVENANCE_FILE_BYTES, "plan is missing or too large")
    plan_hash = sha(plan_path)
    plan, _ = document({"path": str(plan_path), "sha256": plan_hash}, plan_path.parent)
    require(set(plan) in ({"schema", "sources", "limits"}, {"schema", "sources", "limits", "modality"})
            and plan["schema"] == BATCH_HANDOFF_SCHEMA, "unsupported plan schema or fields")
    require(isinstance(plan["sources"], dict) and set(plan["sources"]) == {
        "batch_summary", "verification_receipt", "source_profile_script"}, "invalid handoff source roles")
    limits = plan["limits"]
    require(isinstance(limits, dict) and set(limits) == {"max_studies", "max_total_matrix_bytes"}, "invalid handoff limits")
    max_studies = number(limits["max_studies"], "max_studies", integer=True, minimum=1)
    matrix_cap = number(limits["max_total_matrix_bytes"], "max_total_matrix_bytes", integer=True, minimum=1)
    batch, batch_source = document(plan["sources"]["batch_summary"], plan_path.parent, absolute=True)
    batch_root = Path(batch_source["path"]).parent
    require(not output.is_relative_to(batch_root), "handoff output must be outside the archived batch")
    checked, checked_source = document(plan["sources"]["verification_receipt"], plan_path.parent, absolute=True)
    profile_script = bind(plan["sources"]["source_profile_script"], plan_path.parent, absolute=True)
    require(batch.get("schema") == "clockbase-archs4-batch-extraction-summary-v1" and batch.get("path_base") == "batch_root"
            and batch.get("execution_status") == "COMPLETED" and batch.get("status") in {"COMPLETED", "COMPLETED_WITH_UNRESOLVED"},
            "extraction batch is not completed")
    require(checked.get("schema") == "clockbase-archs4-independent-every-cell-v1" and checked.get("status") == "PASS",
            "whole independent every-cell readback must PASS")
    require(checked.get("batch_summary_sha256") == batch_source["sha256"]
            and declaration(checked.get("batch_summary"), Path(checked_source["path"]).parent) == batch_source,
            "independent readback belongs to a different batch")
    profile, profile_source = document(batch["source_profile_descriptor"], batch_root)
    require(profile.get("schema") == "clockbase-rna-source-profile-census-v1" and profile.get("status") == "PASS"
            and profile.get("verification_script_sha256") == profile_script["sha256"]
            and checked.get("source_profile_sha256") == profile_source["sha256"], "source profile binding differs")
    modality = profile.get("modality")
    require(modality in RNA_SPECIES, "source profile has no supported RNA modality")
    species = RNA_SPECIES[modality]
    require(plan.get("modality", "human_rna") == batch.get("modality", "human_rna")
            == checked.get("modality", "human_rna") == modality
            and batch.get("source_hdf5", {}).get("species", "human") == checked.get("species", "human") == species,
            "handoff plan, batch, profile or readback species differs")
    for record, model_field in ((batch, "model_api_calls"), (checked, "model_calls"), (profile, "model_calls")):
        require(record.get("scientific_status") == "CANNOT VERIFY"
                and record.get("original_paper_2025_eligibility") == "QUARANTINED"
                and type(record.get(model_field)) is int and record[model_field] == 0,
                "source receipt broadens the scientific or model scope")
    require(all(type(record.get(key)) is int and record[key] == 0 for record in (batch, checked)
                for key in ("preprocessing_runs", "clock_runs")), "source receipt broadens the execution scope")
    ids = inventory(batch.get("results"), "batch")
    require(ids == inventory(profile.get("studies"), "profile") == inventory(checked.get("results"), "verification")
            and len(ids) <= max_studies, "complete study inventory differs or exceeds limit")
    require(all(row.get("modality") == modality
                for rows in (batch["results"], profile["studies"], checked["results"])
                for row in rows),
            "mixed RNA modality in the linked study inventory")
    supported, stopped = profile.get("supported_study_ids"), profile.get("stopped_study_ids")
    require(isinstance(supported, list) and isinstance(stopped, list)
            and len(supported) == len(set(supported)) and len(stopped) == len(set(stopped))
            and not set(supported) & set(stopped) and set(supported) | set(stopped) == set(ids), "source-profile partition differs")
    require(batch.get("counts", {}).get("enrolled_studies") == checked.get("enrolled_studies") == len(ids),
            "enrolled-study accounting differs")
    for source, checked_row in zip(profile["studies"], checked["results"]):
        samples = source.get("sample_ids")
        require(isinstance(samples, list) and samples and all(isinstance(sample, str) and re.fullmatch(r"GSM[0-9]+", sample)
                for sample in samples) and len(samples) == len(set(samples)), "source sample inventory is ambiguous")
        if checked_row.get("status") == "CANNOT VERIFY":
            require(checked_row.get("sample_ids") == samples, "unresolved readback sample inventory differs")
    expected_status_counts = {}
    for row in checked["results"]:
        expected_status_counts[row.get("status")] = expected_status_counts.get(row.get("status"), 0) + 1
    require(checked.get("study_status_counts") == expected_status_counts
            and set(expected_status_counts) <= {"PASS", "CANNOT VERIFY"}, "independent study-status accounting differs")
    source_hdf5 = batch.get("source_hdf5", {})
    require(isinstance(source_hdf5.get("sha256"), str) and re.fullmatch(r"[0-9a-f]{64}", source_hdf5["sha256"])
            and type(source_hdf5.get("size_bytes")) is int and source_hdf5["size_bytes"] > 0
            and all(checked.get(key) == source_hdf5["sha256"] for key in (
                "release_sha256", "release_sha256_before", "release_sha256_after")), "historical HDF5 receipt differs")
    matrix_bytes, matrices = 0, {}
    for row in batch["results"]:
        require(row.get("status") in {"EXTRACTED", "SOURCE_PROFILE_UNRESOLVED", "EXTRACTION_UNRESOLVED"}, "unknown extraction disposition")
        if row["status"] == "EXTRACTED":
            item = declaration(row.get("matrix"), batch_root)
            path = Path(item["path"])
            require(path.is_relative_to(batch_root) and path.suffix == ".csv" and str(path) not in matrices,
                    "matrix is foreign, duplicated or not a CSV")
            size = number(row["matrix"].get("size_bytes"), "matrix.size_bytes", integer=True, minimum=1)
            matrix_bytes += size
            require(path.is_file() and path.stat().st_size == size and matrix_bytes <= matrix_cap,
                    "matrix byte count differs or exceeds aggregate limit")
            matrices[str(path)] = item["sha256"]
    require(batch.get("counts", {}).get("extracted_bytes") == matrix_bytes
            and batch.get("counts", {}).get("extracted_studies") == len(matrices), "extracted matrix accounting differs")
    # Bind all batch artifacts, including stopped studies, plus every small source
    # consumed by the independent readback/profile census. No HDF5 is in this graph.
    artifacts = batch.get("artifacts")
    require(isinstance(artifacts, dict) and artifacts, "batch lacks its artifact inventory")
    for name, expected in artifacts.items():
        require(isinstance(name, str) and not Path(name).is_absolute(), "batch artifact name must be relative")
        path = (batch_root / name).resolve()
        require(path.is_relative_to(batch_root) and str(path.relative_to(batch_root)) == name, "batch artifact escapes or aliases its root")
        bind({"path": str(path), "sha256": expected}, batch_root, cap=matrix_cap if str(path) in matrices else MAX_PROVENANCE_FILE_BYTES)
    for path, expected in matrices.items():
        require(artifacts.get(str(Path(path).relative_to(batch_root))) == expected, "matrix missing from batch inventory")
    for record in (checked, profile):
        require(isinstance(record.get("source_files"), dict), "source receipt lacks its consumed-file inventory")
        for name, expected in record["source_files"].items():
            bind({"path": name, "sha256": expected}, batch_root, absolute=True,
                 cap=matrix_cap if name in matrices else MAX_PROVENANCE_FILE_BYTES)
    verifier = bind(checked.get("verification_script"), Path(checked_source["path"]).parent)
    require(verifier["sha256"] == checked.get("verification_script_sha256"), "independent verifier snapshot differs")
    batch_seal = Path(batch_source["path"]).with_suffix(".json.sha256")
    require(batch_seal.is_file() and batch_seal.read_text().strip() == batch_source["sha256"], "batch completion seal differs")
    bind({"path": str(batch_seal), "sha256": sha(batch_seal)}, batch_root)
    implementations = {role: bind(value, batch_root) for role, value in batch.get("implementation", {}).items()}
    require(bool(implementations), "producing implementation inventory is empty")
    config, _ = document({"path": str(batch_root / "config.json"), "sha256": batch.get("config_sha256")}, batch_root)
    batch_plan, _ = document({"path": config.get("plan_path"), "sha256": batch.get("plan_sha256")}, batch_root, absolute=True)
    require(config.get("schema") == "clockbase-archs4-batch-extraction-config-v1"
            and config.get("plan_sha256") == batch.get("plan_sha256")
            and config.get("modality", "human_rna") == batch_plan.get("modality", "human_rna") == modality
            and batch_plan.get("schema") == "clockbase-archs4-batch-extraction-plan-v1"
            and isinstance(batch_plan.get("output"), str) and Path(batch_plan["output"]).is_absolute()
            and Path(batch_plan["output"]).resolve() == batch_root
            and config.get("selected_study_ids") == batch_plan.get("selected_study_ids") == supported,
            "batch plan/config routing or selection differs")
    require(all(config.get(key) == batch_plan.get(key) for key in ("source_profile_descriptor", "release", "runtime", "limits", "implementation"))
            and declaration(config.get("source_profile_descriptor"), batch_root) == profile_source
            and config.get("release", {}).get("sha256") == source_hdf5["sha256"]
            and config.get("release", {}).get("size_bytes") == source_hdf5["size_bytes"]
            and config.get("release", {}).get("species", "human") == batch_plan.get("release", {}).get("species", "human") == species
            and set(config.get("implementation", {})) == set(implementations)
            and all(config["implementation"][role].get("sha256") == item["sha256"] for role, item in implementations.items()),
            "batch plan/config execution inputs differ")
    current_code = {"path": str(Path(__file__).resolve()), "sha256": sha(__file__)}
    bind(current_code, ROOT)
    run_identity = {"schema": BATCH_HANDOFF_SCHEMA, "plan_path": str(plan_path), "plan_sha256": plan_hash,
                    "output": str(output), "implementation": current_code}
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".handoff.lock").open("a") as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError("Another handoff owns this output")
        final_path = output / "handoff.json"
        if final_path.exists():
            final = read_sealed(final_path)
            require(final.get("identity") == run_identity, "completed handoff identity changed; use a new output")
            require(final.get("schema") == "clockbase-archs4-update-handoff-v1"
                    and final.get("execution_status") == "COMPLETED" and final.get("scientific_status") == "CANNOT VERIFY"
                    and final.get("original_paper_2025_eligibility") == "QUARANTINED"
                    and final.get("clock_method_validation") == "NOT_ESTABLISHED"
                    and all(type(final.get(key)) is int and final[key] == 0 for key in (
                        "model_calls", "preprocessing_runs", "clock_runs", "hdf5_files_opened")), "completed handoff scope changed")
            require(inputs == final.get("input_sha256"), "completed source inventory changed")
            members = {str(path.relative_to(output)) for path in output.rglob("*") if path.is_file()}
            require(members == set(final["output_sha256"]) | {".handoff.lock", "handoff.json", "handoff.json.sha256"},
                    "completed output inventory changed")
            for name, expected in final["output_sha256"].items():
                path = (output / name).resolve()
                require(path.is_relative_to(output) and sha(path) == expected, "completed handoff output changed: " + name)
            return final
        require(not any(path.name != ".handoff.lock" for path in output.iterdir()), "unfinished handoff retained; use a new output")
        seal(output / "config.json", {"identity": run_identity, "plan": plan})
        snapshot = output / "implementation" / Path(__file__).name
        snapshot.parent.mkdir()
        with Path(__file__).open("rb") as source, snapshot.open("xb") as target: shutil.copyfileobj(source, target)
        require(sha(snapshot) == current_code["sha256"], "implementation changed during snapshot")
        candidates, ledger, total_cells = [], [], 0
        for row, source, checked_row in zip(batch["results"], profile["studies"], checked["results"]):
            study, modality = row["study_id"], row["modality"]
            require(source["modality"] == checked_row["modality"] == modality
                    and identity(source) == row.get("source_profile_record_sha256")
                    and source.get("profile_status") == row.get("source_profile_status"), "study source-profile identity differs: " + study)
            source_items = {role: bind(row[role], batch_root) for role in ("series", "samples")}
            require(all(declaration(source["sources"][role], Path(profile_source["path"]).parent) == source_items[role]
                        for role in source_items), "GEO source lineage differs: " + study)
            entry = {"study_id": study, "modality": modality, "source_profile_status": source["profile_status"],
                "source_profile_record_sha256": identity(source), "extraction_status": row["status"],
                "verification_status": checked_row["status"], "sources": source_items,
                "sample_ids": source["sample_ids"], "sample_count": len(source["sample_ids"]),
                "limitations": source.get("limitations", []), "reason": None, "candidate": None,
                "audit_status": "NOT_YET_AUDITED", "scientific_status": "CANNOT VERIFY",
                "original_paper_2025_eligibility": "QUARANTINED"}
            if row["status"] != "EXTRACTED":
                require(checked_row["status"] == "CANNOT VERIFY" and checked_row.get("cells_compared") == 0
                        and row.get("matrix") is None and row.get("extraction_receipt") is None,
                        "unresolved study has contradictory matrix/verification state: " + study)
                require((row["status"] == "SOURCE_PROFILE_UNRESOLVED" and study in stopped)
                        or (row["status"] == "EXTRACTION_UNRESOLVED" and study in supported), "study stop differs from source partition")
                entry.update(disposition=row["status"], reason=checked_row.get("reason") or row.get("error") or source.get("limitations"))
            else:
                require(study in supported and source["profile_status"] == "SUPPORTED_SOURCE_PROFILE"
                        and checked_row["status"] == "PASS", "extracted study lacks independent PASS: " + study)
                wrapper_path = output / "wrappers" / (modality + "-" + study + ".json")
                wrapper = {"schema": MATRIX_PROVENANCE_SCHEMA, "study_id": study, "modality": modality, **source_items,
                    "matrix": declaration(row["matrix"], batch_root), "batch_summary": batch_source,
                    "extraction_receipt": declaration(row["extraction_receipt"], batch_root),
                    "verification_receipt": checked_source, "source_profile_descriptor": profile_source,
                    "source_profile_script": profile_script, "source_profile_record_sha256": identity(source),
                    "source_hdf5": {"sha256": source_hdf5["sha256"], "size_bytes": source_hdf5["size_bytes"],
                                    "species": species,
                                    "hash_scope": "historical_batch_pre_and_post"},
                    "implementation": implementations, "scientific_status": "CANNOT VERIFY",
                    "original_paper_2025_eligibility": "QUARANTINED"}
                seal(wrapper_path, wrapper)
                candidate_sources = {**source_items, "matrix": wrapper["matrix"]}
                try:
                    verified = read_matrix_provenance(wrapper_path, sha(wrapper_path), candidate_sources, study, modality)
                    for role, item in verified["files"].items():
                        if role != "wrapper": bind(item, batch_root)
                    candidate_sources["matrix_provenance"] = {"path": str(wrapper_path), "sha256": sha(wrapper_path)}
                    candidate = {"study_id": study, "modality": modality, "sources": candidate_sources}
                    candidates.append(candidate)
                    cells = verified["measurement"]["cells_compared"]
                    entry.update(disposition="READY_FOR_UPDATE", candidate=candidate, matrix=wrapper["matrix"],
                                 feature_count=verified["measurement"]["feature_count"], cells_compared=cells)
                    total_cells += cells
                except (OSError, ValueError, KeyError, TypeError) as error:
                    entry.update(disposition="PROVENANCE_UNRESOLVED", reason=type(error).__name__ + ": " + str(error),
                                 proposed_wrapper={"path": str(wrapper_path), "sha256": sha(wrapper_path)})
            ledger.append(entry)
        # The all-case ledger is complete even when source-profile stops remain.
        dispositions = {}
        for row in ledger: dispositions[row["disposition"]] = dispositions.get(row["disposition"], 0) + 1
        require(checked["cells_compared"] == sum(row.get("cells_compared", 0) for row in checked["results"]),
                "independent total cell accounting differs")
        seal(output / "candidates.json", {"schema": "clockbase-update-candidates-v1", "studies": candidates,
             "scientific_status": "CANNOT VERIFY", "original_paper_2025_eligibility": "QUARANTINED"})
        seal(output / "candidate-ledger.json", {"schema": "clockbase-archs4-update-disposition-ledger-v1", "studies": ledger,
             "batch_summary": batch_source, "verification_receipt": checked_source,
             "enrolled_studies": len(ids), "disposition_counts": dispositions})
        for name, expected in inputs.items(): require(sha(name) == expected, "input changed during handoff: " + name)
        require(sha(__file__) == current_code["sha256"], "implementation changed during handoff")
        outputs = {str(path.relative_to(output)): sha(path) for path in output.rglob("*") if path.is_file() and path.name != ".handoff.lock"}
        final = {"schema": "clockbase-archs4-update-handoff-v1", "identity": run_identity,
            "modality": modality, "species": species,
            "status": "COMPLETED_WITH_UNRESOLVED" if len(candidates) != len(ids) else "COMPLETED", "execution_status": "COMPLETED",
            "counts": {"enrolled_studies": len(ids), "candidate_studies": len(candidates), "unresolved_studies": len(ids)-len(candidates),
                       "candidate_samples": sum(len(row["sample_ids"]) for row in ledger if row["candidate"]),
                       "extracted_matrix_bytes": matrix_bytes, "candidate_cells_compared": total_cells, "dispositions": dispositions},
            "candidates": {"path": str(output / "candidates.json"), "sha256": outputs["candidates.json"]},
            "ledger": {"path": str(output / "candidate-ledger.json"), "sha256": outputs["candidate-ledger.json"]},
            "input_sha256": dict(sorted(inputs.items())), "output_sha256": dict(sorted(outputs.items())),
            "model_calls": 0, "preprocessing_runs": 0, "clock_runs": 0, "hdf5_files_opened": 0,
            "scientific_status": "CANNOT VERIFY", "clock_method_validation": "NOT_ESTABLISHED",
            "original_paper_2025_eligibility": "QUARANTINED"}
        seal(final_path, final)
        return final


def source_design_replay_binding(version_root):
    """Describe only an archived source assessment; never open producer outcomes."""
    root = Path(version_root).resolve()
    version = read_sealed(root / "version.json")
    if version["version_id"] != identity(version["identity"]):
        raise UpdateIntegrityError("Replay origin version identity changed")
    names = {"version.json", "version.json.sha256", "pipeline/config.json"}
    for stage in ("01_metadata", "03_blinded_design"):
        directory = root / "pipeline" / stage
        receipt = directory / "complete.json"
        if sha(receipt) != (directory / "complete.sha256").read_text().strip():
            raise UpdateIntegrityError("Replay stage receipt changed: " + stage)
        record = json.loads(receipt.read_text())
        if not isinstance(record.get("outputs"), dict) or not record["outputs"]:
            raise UpdateIntegrityError("Replay stage has no bound members: " + stage)
        members = set()
        names.update({str(receipt.relative_to(root)), str((directory / "complete.sha256").relative_to(root))})
        for relative, expected in record["outputs"].items():
            path = (directory / relative).resolve()
            if not path.is_relative_to(directory) or sha(path) != expected:
                raise UpdateIntegrityError("Replay stage output changed or escapes its stage")
            names.add(str(path.relative_to(root)))
            members.add(path)
        if stage == "01_metadata":
            for field in ("metadata", "study_context", "source_sample_records"):
                path = record["result"].get(field)
                if path is not None and Path(path).resolve() not in members:
                    raise UpdateIntegrityError("Replay metadata is not a completed-stage member")
    slots = sorted((root / "pipeline/api").glob("design-*"))
    if not slots or [p.name for p in slots] != [f"design-{i:03d}" for i in range(len(slots))] or len(slots) > 2:
        raise UpdateIntegrityError("Replay requires one or two complete source-model calls")
    for slot in slots:
        names.update(str((slot / name).relative_to(root)) for name in
                     ("request.json", "response.json", "reservation.json", "receipt.json"))
    files = {}
    for name in sorted(names):
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise UpdateIntegrityError("Replay source escapes its origin")
        files[name] = sha(path)
    return {"schema": "clockbase-source-design-replay-v1", "origin_version_root": str(root),
            "origin_version_id": version["version_id"], "files": files}


def read_source_design_replay(path, expected_sha256):
    path = Path(path).resolve()
    if not re.fullmatch(r"[0-9a-f]{64}", str(expected_sha256)) or sha(path) != expected_sha256:
        raise UpdateIntegrityError("Source-design replay manifest changed or lacks its SHA256")
    binding = json.loads(path.read_text())
    if (not isinstance(binding, dict) or set(binding) != {"schema", "origin_version_root", "origin_version_id", "files"}
            or binding["schema"] != "clockbase-source-design-replay-v1"):
        raise UpdateIntegrityError("Unknown source-design replay manifest")
    observed = source_design_replay_binding(binding["origin_version_root"])
    if binding != observed or sha(path) != expected_sha256:
        raise UpdateIntegrityError("Source-design replay input or inventory changed")
    return binding


def inspect_source_design(study, modality, sources, max_samples):
    from audit_system.extension_eligibility_2026_09_10 import series_record
    from src.pipeline_inputs_2026_09_10 import parse_samples
    record = series_record(study, Path(sources["series"]["path"]).read_bytes())
    samples = parse_samples(Path(sources["samples"]["path"]).read_text(), study)
    if any(not re.fullmatch(r"GSM[0-9]+", sample["GSM"]) for sample in samples):
        raise ValueError("Noncanonical GEO sample identity")
    if len(samples) > max_samples: raise ValueError("Source sample count exceeds plan limit")
    if {sample["GSM"] for sample in samples} != set(record["sample_ids"]):
        raise ValueError("Series and sample metadata disagree; no partial-study admission")
    organism = "Mus musculus" if modality == "mouse_rna" else "Homo sapiens"
    if {sample.get("organism_ch1") for sample in samples} != {organism}:
        raise ValueError("Source organism differs from declared modality")
    return record


def frozen_metadata(study, output, sources):
    """Existing GEO parsers over frozen source bytes; no extra network requests."""
    from src.pipeline_inputs_2026_09_10 import parse_samples, acquire_study_context
    from audit_system.extension_eligibility_2026_09_10 import series_record
    verify_sources(sources)
    series = output / "series.soft"; sample_source = output / "samples.soft"
    shutil.copyfile(sources["series"]["path"], series)
    shutil.copyfile(sources["samples"]["path"], sample_source)
    record = series_record(study, series.read_bytes())
    samples = parse_samples(sample_source.read_text(), study)
    if {r["GSM"] for r in samples} != set(record["sample_ids"]): raise ValueError("Frozen GEO memberships differ")
    path = output / "metadata.csv"
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted({key for row in samples for key in row}))
        writer.writeheader(); writer.writerows(samples)
    if "context" in sources:
        context = output / "study-context.json"
        shutil.copyfile(sources["context"]["path"], context)
        declaration = json.loads(context.read_text())
        if declaration.get("study_id") != study: raise ValueError("Context names a different study")
        for entry in declaration.get("sources", []):
            name = entry["path"]
            if Path(name).name != name: raise ValueError("Context source must use a local filename")
            target = output / name
            if target.exists():
                if sha(target) != entry["sha256"]: raise ValueError("Conflicting frozen context source")
            else:
                source = sources["context_source:" + name]
                shutil.copyfile(source["path"], target)
    else:
        context = acquire_study_context(study, series, output, publication_methods=False)
    record.update(metadata=str(path), source_sample_records=str(sample_source), series_source=str(series),
                  study_context=str(context), sample_metadata_view="frozen_source", acquisition_adapter="continuous_update_v1")
    write_new(output / "source-membership.json", record)
    return record


def collect_sources(entry, root, limits, input_files):
    declarations = entry.get("sources", {})
    if not isinstance(declarations, dict) or set(declarations) - ROLES or not {"series", "samples"} <= set(declarations):
        raise ValueError("Study requires frozen GEO series/sample source declarations")
    if ("matrix" in declarations) == ("archs4" in declarations):
        raise ValueError("Declare one molecular matrix or a pinned local ARCHS4 release")
    if entry["modality"] == "human_dnam" and "archs4" in declarations:
        raise ValueError("DNAm cannot use ARCHS4")
    if "dnam_idat_provenance" in declarations:
        declared_idat = declarations["dnam_idat_provenance"]
        if (entry["modality"] != "human_dnam" or "matrix" not in declarations
                or {"archs4", "matrix_provenance", "dnam_supplement"} & set(declarations)
                or not isinstance(declared_idat, dict) or set(declared_idat) != {"path", "sha256"}):
            raise ValueError("DNAm IDAT provenance requires a DNAm matrix and one local path/SHA256 declaration")
        provenance_path = Path(declared_idat["path"]).expanduser().resolve()
        if provenance_path.is_file() and provenance_path.stat().st_size > MAX_PROVENANCE_FILE_BYTES:
            raise ValueError("DNAm IDAT provenance exceeds small-receipt cap")
    if "archs4" in declarations and "path" not in declarations["archs4"]:
        raise ValueError("Acquire large ARCHS4 releases with the existing ETag-pinned updater first")
    if "matrix_provenance" in declarations:
        if ("matrix" not in declarations or "archs4" in declarations or entry["modality"] not in {"human_rna", "mouse_rna"}
                or not isinstance(declarations["matrix_provenance"], dict)
                or set(declarations["matrix_provenance"]) != {"path", "sha256"}):
            raise ValueError("Derived matrix provenance requires a declared RNA matrix and exactly a local path/SHA256")
        provenance_path = Path(declarations["matrix_provenance"]["path"]).expanduser().resolve()
        if provenance_path.is_file() and provenance_path.stat().st_size > MAX_PROVENANCE_FILE_BYTES:
            raise ValueError("Derived matrix provenance exceeds small-receipt cap")
    def acquire(value, network):
        known = input_files.get(value.get("sha256"), 0)
        remaining = limits["max_input_bytes"] - sum(input_files.values()) + known
        if remaining <= 0: raise ValueError("Aggregate source input byte cap reached before acquisition")
        source = file_source(value, root, max_bytes=remaining, network=network)
        input_files[source["sha256"]] = source["bytes"]
        return source
    sources = {role: acquire(value, limits["allow_network"]) for role, value in sorted(declarations.items())}
    if "dnam_idat_provenance" in sources:
        declared = sources["dnam_idat_provenance"]
        binding = read_dnam_idat_provenance(declared["path"], declared["sha256"], entry["study_id"], sources=sources)
        for role, item in binding["files"].items():
            if role != "wrapper":
                sources["dnam_idat_provenance:" + role] = acquire(item, False)
    if "matrix_provenance" in sources:
        declared = sources["matrix_provenance"]
        binding = read_matrix_provenance(declared["path"], declared["sha256"], sources, entry["study_id"], entry["modality"])
        for role, item in binding["files"].items():
            if role != "wrapper":
                sources["matrix_provenance:" + role] = acquire(item, False)
    if "context" in sources:
        context_path = Path(sources["context"]["path"])
        context = json.loads(context_path.read_text())
        for item in context.get("sources", []):
            name = item["path"]
            if Path(name).name != name: raise ValueError("Context path escapes its source directory")
            sources["context_source:" + name] = acquire(
                {"path": str(context_path.parent / name), "sha256": item["sha256"]}, False)
    if entry.get("source_design_replay") is not None:
        declared = entry["source_design_replay"]
        if not isinstance(declared, dict) or set(declared) != {"path", "sha256"}:
            raise ValueError("Source-design replay requires exactly a local path and SHA256")
        sources["source_design_replay"] = acquire(declared, False)
        binding = read_source_design_replay(declared["path"], declared["sha256"])
        for name, digest in binding["files"].items():
            sources["source_design_replay:" + name] = acquire(
                {"path": str(Path(binding["origin_version_root"]) / name), "sha256": digest}, False)
    return sources


def validate_replay_sources(entry, sources, options):
    if "source_design_replay" not in sources:
        return
    binding = read_source_design_replay(sources["source_design_replay"]["path"], sources["source_design_replay"]["sha256"])
    root = Path(binding["origin_version_root"])
    version = read_sealed(root / "version.json")
    current = {role: value["sha256"] for role, value in sources.items() if not role.startswith("source_design_replay")}
    original = {role: value["sha256"] for role, value in version["sources"].items()}
    config = json.loads((root / "pipeline/config.json").read_text())
    if (current != original or version["study_id"] != entry["study_id"] or version["modality"] != entry["modality"]
            or version["target"] != entry.get("target", "Source-defined intervention comparisons")
            or config.get("audit_profile") != "source_domains_v4"
            or config.get("source_binding_profile") != "literal_source_v1"
            or config.get("comparison_scope") != "source_defined_interventions"
            or config.get("source_contract_policy", "strict") != options["source_contract_policy"]
            or config.get("audit_model") != options["audit_model"]
            or ("service_tier" in config) != ("service_tier" in options)
            or config.get("service_tier") != options.get("service_tier")
            or any(("service_tier" in profile) != ("service_tier" in options)
                   or profile.get("service_tier") != options.get("service_tier")
                   for profile in (version.get("pipeline_options", {}),
                                   version.get("identity", {}).get("pipeline_options", {})))):
        raise UpdateIntegrityError("Replay source identity, target or source-review profile differs")


def reserve(root, version, phase, calls, cost, limits):
    path = root / "reservations" / (version + "-" + phase + ".json")
    expected = {"version_id": version, "phase": phase, "max_calls": calls, "max_cost_usd": cost}
    if path.exists():
        if read_sealed(path) != expected: raise UpdateIntegrityError("Budget reservation changed")
        return expected
    existing = [read_sealed(p) for p in (root / "reservations").glob("*.json")]
    if sum(r["max_calls"] for r in existing) + calls > limits["max_total_model_calls"]:
        raise RuntimeError("Aggregate model-call reservation limit reached")
    if sum(r["max_cost_usd"] for r in existing) + cost > limits["max_total_cost_usd"] + 1e-12:
        raise RuntimeError("Aggregate model-cost reservation limit reached")
    seal(path, expected)
    return expected


def validate_plan(plan):
    if plan.get("schema") != SCHEMA: raise ValueError("Unsupported continuous-update schema")
    policy=plan.get("comprehensive_audit_policy", "literal_fields_v1")
    if type(policy) is not str or policy not in {"literal_fields_v1", "source_label_pairs_v1"}:
        raise ValueError("Unknown controller-selected comprehensive source-design policy")
    limits = plan.get("limits", {})
    for name in ("max_studies", "max_input_bytes", "max_samples_per_study", "max_child_seconds"):
        number(limits.get(name), name, integer=True, minimum=1)
    number(limits.get("max_total_model_calls"), "max_total_model_calls", integer=True)
    number(limits.get("max_total_cost_usd"), "max_total_cost_usd")
    if limits.get("concurrency") != 1 or type(limits.get("allow_network")) is not bool:
        raise ValueError("Declare concurrency one and an explicit network policy")
    declared_options = plan.get("pipeline", {})
    if set(declared_options) - PIPELINE_OPTIONS: raise ValueError("Unsupported child pipeline option")
    for model_name, input_name, output_name in (
            ("model", "input_usd_per_million", "output_usd_per_million"),
            ("audit_model", "audit_input_usd_per_million", "audit_output_usd_per_million")):
        tariff_names = {input_name, output_name}
        model_changed = model_name in declared_options and declared_options[model_name] != DEFAULT_OPTIONS[model_name]
        tariffs_present = tariff_names & set(declared_options)
        if model_changed and not tariff_names <= set(declared_options):
            raise ValueError("A non-default " + model_name + " requires both explicit role-specific tariffs")
        if tariffs_present and tariffs_present != tariff_names:
            raise ValueError("Role-specific tariffs must be declared as an input/output pair for " + model_name)
    options = {**DEFAULT_OPTIONS, **declared_options}
    if "service_tier" in options and (type(options["service_tier"]) is not str or options["service_tier"] != "default"):
        raise ValueError("Explicit pipeline service_tier must be default")
    completion_policy = options.get("coding_completion_policy", "legacy")
    if type(completion_policy) is not str or completion_policy not in {"legacy", "outputs_v1"}:
        raise ValueError("Unknown coding completion policy")
    if completion_policy == "legacy":
        options.pop("coding_completion_policy", None)
    if (type(options["source_contract_policy"]) is not str or options["source_contract_policy"] not in
            {"strict", "conservative_status_v1", "literal_field_refs_v1"}):
        raise ValueError("Unknown source contract policy")
    for name in ("audit_validation_retries", "coding_validation_retries"):
        if type(options[name]) is not int or options[name] not in (0, 1):
            raise ValueError(name + " must be 0 or 1")
    for name in ("max_calls", "max_output_tokens", "max_prompt_bytes"):
        number(options[name], name, integer=True, minimum=1)
    for name in ("max_cost_usd", "input_usd_per_million", "output_usd_per_million", "audit_input_usd_per_million", "audit_output_usd_per_million"):
        number(options[name], name, minimum=1e-12)
    audit = plan.get("comprehensive_audit", {})
    required_audit = {"model", "max_calls", "max_cost_usd", "input_usd_per_million", "output_usd_per_million"}
    if not required_audit <= set(audit) or set(audit) - required_audit - {"service_tier"}:
        raise ValueError("Explicit comprehensive-audit model, tariffs and caps are required")
    if "service_tier" in audit and (type(audit["service_tier"]) is not str or audit["service_tier"] != "default"):
        raise ValueError("Explicit comprehensive-audit service_tier must be default")
    number(audit["max_calls"], "audit.max_calls", integer=True, minimum=1)
    for name in ("max_cost_usd", "input_usd_per_million", "output_usd_per_million"):
        number(audit[name], "audit." + name, minimum=1e-12)
    if not all(isinstance(v, str) and v for v in (options["model"], options["audit_model"], audit["model"])):
        raise ValueError("Explicit model names required; no fallback provider")
    if audit["model"] not in ("gpt-5.6-luna", "gpt-6-luna"):
        raise ValueError("Comprehensive audit requires explicit OpenAI gpt-5.6-luna or gpt-6-luna; unsupported model rejected before execution")
    return limits, options, audit


def run_child(request_path, log, timeout):
    with log.open("xb") as handle:
        process = subprocess.Popen([sys.executable, "-m", "src.update_pipeline_2026_09_20", "--child-request", str(request_path)],
                                   cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, start_new_session=True)
        try: return process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL); process.wait()
            raise TimeoutError("Child exceeded its declared wall-time cap; process group stopped")


def child(request):
    """Process-local acquisition adapter; restore the original pipeline function."""
    version = read_sealed(Path(request["version_path"]))
    sources = version["sources"]
    verify_sources(sources)
    verify_implementation(version["implementation"])
    if request["phase"] == "comprehensive_audit":
        module = importlib.import_module("audit_system.corpus_audit_2026_09_20")
        function = getattr(module, "audit_pipeline_result")
        result = function(request["audit_request"])
    else:
        from src import pipeline_2026_09_10 as pipeline
        argv = ["--study", version["study_id"], "--modality", version["modality"],
                "--target", version["target"], "--output", request["pipeline_root"],
                "--cache", request["cache"], "--until", request["until"],
                "--comparison-scope", "source_defined_interventions", "--source-binding-profile", "literal_source_v1",
                "--audit-profile", "source_domains_v4", "--analysis-contract", "explicit_effect_v2",
                "--no-publication-methods"]
        for name, value in version["pipeline_options"].items(): argv += ["--" + name.replace("_", "-"), str(value)]
        if "source_design_replay" in sources:
            replay = sources["source_design_replay"]
            argv += ["--source-design-replay", replay["path"], "--source-design-replay-sha256", replay["sha256"]]
        for role in ("matrix", "archs4", "dnam_reference", "dnam_supplement", "dnam_idat_provenance"):
            if role in sources: argv += ["--" + role.replace("_", "-"), sources[role]["path"]]
        original = pipeline.acquire_metadata
        pipeline.acquire_metadata = lambda study, output, **unused: frozen_metadata(study, output, sources)
        try: result = pipeline.main(argv)
        finally: pipeline.acquire_metadata = original
    verify_sources(sources)
    verify_implementation(version["implementation"])
    write_new(Path(request["result_path"]), result)
    return result


def verify_completion(path):
    record = read_sealed(path)
    for name, expected in record.get("output_sha256", {}).items():
        if sha(Path(name)) != expected: raise UpdateIntegrityError("Completed child output changed: " + name)
    return record


@contextlib.contextmanager
def source_deadline(seconds, label="Public source request"):
    """Wall-time bound for synchronous public-source reads on the Linux runner."""
    def expired(signum, frame): raise TimeoutError(label + " exceeded its wall-time bound")
    previous_handler = signal.signal(signal.SIGALRM, expired)
    previous_timer = signal.getitimer(signal.ITIMER_REAL)
    began = time.monotonic()
    signal.setitimer(signal.ITIMER_REAL, min(seconds, previous_timer[0]) if previous_timer[0] > 0 else seconds)
    try: yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)
        remaining = max(.000001, previous_timer[0] - (time.monotonic() - began)) if previous_timer[0] > 0 else 0
        signal.setitimer(signal.ITIMER_REAL, remaining, previous_timer[1])


def discover(plan_path, output, *, transport=None):
    """Bounded GEO candidate discovery; every unprocessed UID remains in the ledger.

    ESearch publication dates are a search filter. Only GEO Series SOFT Public-on
    evidence supplies a first-public date. Existing atlas IDs are an operational
    seen registry, never the submitted-paper corpus boundary.
    """
    from audit_system.extension_eligibility_2026_09_10 import series_record
    plan = json.loads(plan_path.read_text())
    if plan.get("schema") != "clockbase-update-discovery-plan-v1" or plan.get("modality") not in MODALITIES:
        raise ValueError("Invalid update discovery schema/modality")
    from datetime import date
    date_from, date_to = date.fromisoformat(plan["date_from"]), date.fromisoformat(plan["date_to"])
    if date_from > date_to: raise ValueError("Reversed source date window")
    limits = plan["limits"]
    for name in ("max_requests", "max_response_bytes", "max_total_response_bytes", "max_seconds", "max_query_ids", "max_candidates"):
        number(limits.get(name), name, integer=True, minimum=1)
    if limits["max_query_ids"] > 100000: raise ValueError("GEO discovery retmax is bounded to 100000")
    number(limits.get("minimum_interval_seconds"), "minimum_interval_seconds", minimum=0.35)
    offset = number(plan.get("candidate_offset", 0), "candidate_offset", integer=True)
    if not isinstance(plan.get("query"), str) or not plan["query"].strip(): raise ValueError("Explicit GEO query required")
    default_source = plan.get("default_molecular_source")
    if default_source is not None:
        if (plan["modality"] not in RNA_SPECIES or not isinstance(default_source, dict)
                or set(default_source) != {"archs4"} or set(default_source["archs4"]) != {"path", "sha256"}):
            raise ValueError("RNA discovery default requires a pinned local ARCHS4 source")
        source = file_source(default_source["archs4"], output / "local-release", max_bytes=plan["max_release_bytes"], network=False)
        default_source = {"archs4": {"path": source["path"], "sha256": source["sha256"]}}
    parser_path = ROOT / "audit_system/extension_eligibility_2026_09_10.py"
    config = {"plan": plan, "implementation_sha256": sha(Path(__file__)),
              "source_parser_sha256": sha(parser_path)}
    output.mkdir(parents=True, exist_ok=True)
    config_path = output / "config.json"
    if config_path.exists():
        if read_sealed(config_path) != config: raise UpdateIntegrityError("Discovery plan or implementation changed")
    else: seal(config_path, config)
    complete = output / "summary.json"
    if complete.exists():
        result = verify_completion(complete)
        for path, expected in result["input_sha256"].items():
            if sha(path) != expected: raise UpdateIntegrityError("Discovery input changed")
        return result
    started = time.monotonic(); request_count = 0; response_bytes = 0; last_request = 0.0

    def request(name, url):
        nonlocal request_count, response_bytes, last_request
        destination = output / "requests" / name
        receipt_path = destination / "receipt.json"
        if receipt_path.exists():
            receipt = read_sealed(receipt_path)
            if receipt["url"] != url: raise UpdateIntegrityError("Discovery request changed on resume")
            request_count += 1; response_bytes += receipt.get("accounted_bytes", receipt.get("bytes", 0))
            if receipt["status"] != "PASS": raise RuntimeError("Prior source request failed; retained without retry")
            if sha(destination / "response.bin") != receipt["sha256"]: raise UpdateIntegrityError("Discovery response changed")
            return (destination / "response.bin").read_bytes(), receipt
        if request_count >= limits["max_requests"] or time.monotonic() - started >= limits["max_seconds"]:
            raise RuntimeError("Discovery request/time limit reached")
        remaining = limits["max_total_response_bytes"] - response_bytes
        cap = min(limits["max_response_bytes"], remaining - 1)
        if cap <= 0: raise RuntimeError("Discovery aggregate response-byte limit reached")
        wait = limits["minimum_interval_seconds"] - (time.monotonic() - last_request)
        if wait > 0: time.sleep(wait)
        request_count += 1; last_request = time.monotonic()
        destination.mkdir(parents=True, exist_ok=False)
        receipt = {"url": url, "requested_at_utc": datetime.now(timezone.utc).isoformat(),
                   "status": "RUNNING", "byte_limit": cap}
        observed = None
        try:
            if transport is None:
                remaining_seconds = limits["max_seconds"] - (time.monotonic() - started)
                if remaining_seconds <= 0: raise TimeoutError("Discovery wall-time limit reached")
                with source_deadline(min(60, remaining_seconds)):
                    with urllib.request.urlopen(url, timeout=min(60, remaining_seconds)) as response:
                        payload = response.read(cap + 1)
            else:
                payload = transport(url, cap)
            observed = len(payload)
            if not payload or len(payload) > cap: raise ValueError("Empty or oversized source response")
            with (destination / "response.bin").open("xb") as handle: handle.write(payload)
            receipt.update(status="PASS", bytes=len(payload), sha256=sha(destination / "response.bin"),
                           response_path=str(destination / "response.bin"))
        except Exception as error:
            receipt.update(status="FAIL", error=type(error).__name__ + ": " + str(error), bytes=0)
            raise
        finally:
            receipt["accounted_bytes"] = observed if observed is not None else cap + 1
            response_bytes += receipt["accounted_bytes"]
            seal(receipt_path, receipt)
        return payload, receipt

    input_hashes = {}; seen_studies = set(); seen_samples = set(); covered_seen = set()
    if default_source is not None:
        input_hashes[default_source["archs4"]["path"]] = default_source["archs4"]["sha256"]
    for name in ("seen_worklist", "seen_sample_sources"):
        if name not in plan: continue
        source = file_source(plan[name], output / "seen", max_bytes=2000000000, network=False)
        input_hashes[source["path"]] = source["sha256"]
        with gzip.open(source["path"], "rt") if str(source["path"]).endswith(".gz") else open(source["path"]) as handle:
            if name == "seen_worklist":
                for line in handle:
                    row = json.loads(line)
                    if re.fullmatch(r"GSE[0-9]+", row.get("study_id", "")): seen_studies.add(row["study_id"])
            else:
                for row in json.load(handle):
                    if row.get("accession") in seen_studies:
                        values = {v for v in row.get("sample_ids", []) if re.fullmatch(r"GSM[0-9]+", v)}
                        if values: covered_seen.add(row["accession"]); seen_samples.update(values)
    term = f'({plan["query"]}) AND ({date_from.isoformat().replace("-", "/")}[PDAT] : {date_to.isoformat().replace("-", "/")}[PDAT]) AND gse[ETYP]'
    params = {"db": "gds", "term": term, "retmode": "json", "retmax": limits["max_query_ids"]}
    prefix = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    inherited_ledger = {}; inherited_candidates = []
    if plan.get("query_snapshot"):
        binding = file_source(plan["query_snapshot"], output / "query-source", max_bytes=20000000, network=False)
        snapshot_path = Path(binding["path"])
        previous = verify_completion(snapshot_path)
        if previous.get("schema") != "clockbase-update-discovery-result-v1" or previous.get("query") != term:
            raise ValueError("Frozen discovery snapshot belongs to a different query")
        archived_search = snapshot_path.parent / "requests/esearch/response.bin"
        if not archived_search.is_file():
            inherited = previous.get("query_snapshot_source")
            if not inherited: raise ValueError("Discovery snapshot lacks its archived ESearch source")
            archived_search = Path(inherited["response_path"])
            if sha(archived_search) != inherited["response_sha256"]:
                raise UpdateIntegrityError("Inherited query response changed")
        with (snapshot_path.parent / "candidate-ledger.jsonl").open() as handle:
            inherited_ledger = {row["geo_uid"]: row for row in map(json.loads, handle)}
        inherited_candidates = json.loads((snapshot_path.parent / "candidates.json").read_text())["studies"]
        input_hashes[binding["path"]] = binding["sha256"]
        input_hashes[str(archived_search)] = sha(archived_search)
        payload = archived_search.read_bytes()
        query_snapshot_source = {"summary_path": binding["path"], "summary_sha256": binding["sha256"],
                                 "response_path": str(archived_search), "response_sha256": sha(archived_search)}
    else:
        payload, _ = request("esearch", prefix + "esearch.fcgi?" + urllib.parse.urlencode(params))
        query_snapshot_source = None
    search = json.loads(payload)["esearchresult"]
    ids = search["idlist"]; total = int(search["count"])
    if len(ids) != len(set(ids)) or any(not re.fullmatch(r"[0-9]+", v) for v in ids): raise ValueError("Invalid query UID inventory")
    if len(ids) > total: raise ValueError("Query UID count exceeds reported count")
    selected = sorted(ids, key=int)[offset:offset + limits["max_candidates"]]
    ledger = {uid: {"geo_uid": uid, "accession": None, "status": "NOT_YET_PROCESSED", "reason": "SOURCE_VERIFICATION_CANDIDATE_CAP",
                    "first_public_date": None, "original_paper_2025_eligibility": "QUARANTINED"} for uid in ids}
    for uid in ids:
        if uid in inherited_ledger and uid not in selected: ledger[uid] = inherited_ledger[uid]
    summaries = {}
    summary_failures = {}
    summary_batches = []
    selected_batches = [selected[start:start + 100] for start in range(0, len(selected), 100)]
    for batch_number, batch_uids in enumerate(selected_batches):
        request_name = "esummary" if len(selected_batches) == 1 else f"esummary-{batch_number:04d}"
        batch_record = {"batch_number": batch_number, "request_name": request_name,
                        "requested_uids": batch_uids, "requested_uid_count": len(batch_uids),
                        "status": "RUNNING", "receipt_path": None, "receipt_sha256": None}
        try:
            summary_payload, receipt = request(request_name, prefix + "esummary.fcgi?" + urllib.parse.urlencode(
                {"db": "gds", "id": ",".join(batch_uids), "retmode": "json"}))
            result = json.loads(summary_payload).get("result")
            if not isinstance(result, dict): raise ValueError("GEO summary result is not an object")
            returned = result.get("uids")
            inventory_source = "explicit_uids"
            if returned is None:
                returned = [key for key in result if key != "uids"]
                inventory_source = "result_keys"
            if not isinstance(returned, list) or any(not isinstance(uid, str) for uid in returned):
                raise ValueError("GEO summary lacks a usable UID inventory")
            if len(returned) != len(set(returned)):
                raise ValueError("GEO summary UID inventory contains duplicates")
            expected, observed = set(batch_uids), set(returned)
            unexpected = sorted(observed - expected, key=lambda value: (not value.isascii() or not value.isdigit(), value))
            missing = sorted(expected - observed, key=int)
            keyed = {key for key in result if key != "uids"}
            unexpected_keys = sorted(keyed - expected)
            missing_keys = sorted(expected - keyed, key=int)
            if unexpected or missing or unexpected_keys or missing_keys:
                raise ValueError("GEO summary UID mismatch: " + json.dumps({
                    "missing_uids": missing, "unexpected_uids": unexpected,
                    "missing_keys": missing_keys, "unexpected_keys": unexpected_keys}, sort_keys=True))
            for uid in batch_uids:
                if not isinstance(result[uid], dict): raise ValueError("GEO summary record is not an object for UID " + uid)
                summaries[uid] = result[uid]
            batch_record.update(status="PASS", receipt_path=receipt["response_path"].replace("response.bin", "receipt.json"),
                                receipt_sha256=sha(Path(receipt["response_path"]).with_name("receipt.json")),
                                response_sha256=receipt["sha256"], returned_uids=returned,
                                uid_inventory_source=inventory_source)
        except Exception as error:
            reason = type(error).__name__ + ": " + str(error)
            receipt_path = output / "requests" / request_name / "receipt.json"
            if receipt_path.is_file():
                receipt_record = read_sealed(receipt_path)
                if receipt_record.get("status") == "FAIL" and receipt_record.get("error"):
                    reason = receipt_record["error"]
                batch_record.update(receipt_path=str(receipt_path), receipt_sha256=sha(receipt_path))
            for uid in batch_uids: summary_failures[uid] = reason
            batch_record.update(status="FAIL", error=reason)
        summary_batches.append(batch_record)
    summary_batch_record = {
        "schema": "clockbase-geo-esummary-batches-v1", "maximum_uids_per_request": 100,
        "selected_uid_count": len(selected), "batches": summary_batches}
    summary_batch_path = output / "esummary-batches.json"
    if summary_batch_path.exists():
        if read_sealed(summary_batch_path) != summary_batch_record:
            raise UpdateIntegrityError("GEO summary batch inventory changed on resume")
    else:
        seal(summary_batch_path, summary_batch_record)
    catalog_record = None
    if plan.get("inspect_archs4_catalog", False):
        try:
            payload, receipt = request("archs4-catalog", "https://archs4.org/api/versionfile")
            catalog = json.loads(payload)
            if isinstance(catalog, dict): catalog = catalog.get("versionfiles", catalog.get("data", catalog.get("files", [])))
            species = "mouse" if plan["modality"] == "mouse_rna" else "human"
            choices = [r for r in catalog if r.get("species") == species and r.get("data_level") == "gene"]
            catalog_record = {"status": "CATALOG_INSPECTED_NO_DOWNLOAD", "source_receipt": receipt,
                              "matching_release_declarations": choices, "first_public_date_evidence": False,
                              "download_entry_point": "audit_system.archs4_update_2026_09_10"}
        except Exception as error: catalog_record = {"status": "UNRESOLVED", "error": type(error).__name__ + ": " + str(error)}
    replaced_accessions = {inherited_ledger[uid]["accession"] for uid in selected if uid in inherited_ledger}
    candidates = [r for r in inherited_candidates if r["study_id"] not in replaced_accessions]
    batch_samples = {}
    for uid, record in ledger.items():
        if uid not in selected and record["status"] in {"UNRESOLVED_NO_ASSAY_INPUT", "READY_FOR_INGESTION_ATTEMPT", "SHARED_NEW_CANDIDATE_SAMPLES"}:
            for sample in record.get("sample_ids", []): batch_samples.setdefault(sample, set()).add(record["accession"])
    for uid in selected:
        record = ledger[uid]
        if uid in summary_failures:
            record.update(status="UNRESOLVED", reason="ESUMMARY_BATCH_FAILURE: " + summary_failures[uid])
            continue
        try:
            accession = summaries[uid].get("accession", "")
            if not re.fullmatch(r"GSE[0-9]+", accession): raise ValueError("GEO summary has no canonical Series accession")
            record["accession"] = accession
            if accession in seen_studies:
                record.update(status="SEEN_ACCESSION", reason="OPERATIONAL_REGISTRY_MATCH"); continue
            payload, receipt = request(accession + "-series", f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&form=text&view=full")
            source = series_record(accession, payload)
            record.update(first_public_date=source["series_first_public_date"], source_status_quote=source["source_status_quote"],
                          series_source_sha256=source["response_sha256"], sample_ids=source["sample_ids"],
                          shared_seen_samples=sorted(set(source["sample_ids"]) & seen_samples),
                          esummary_pdat=summaries[uid].get("pdat"), source_receipt=receipt)
            if not date_from.isoformat() <= source["series_first_public_date"] <= date_to.isoformat():
                record.update(status="OUTSIDE_FIRST_PUBLIC_WINDOW", reason="GEO_SOFT_FIRST_PUBLIC_DATE"); continue
            if record["shared_seen_samples"]:
                record.update(status="SHARED_SEEN_SAMPLES", reason="OPERATIONAL_SAMPLE_OVERLAP"); continue
            for sample in source["sample_ids"]: batch_samples.setdefault(sample, set()).add(accession)
            molecular = plan.get("molecular_sources", {}).get(accession, default_source)
            if not molecular:
                record.update(status="UNRESOLVED_NO_ASSAY_INPUT", reason="Candidate requires an assay-specific molecular input or pinned ARCHS4 extraction")
                continue
            if set(molecular) not in ({"matrix"}, {"matrix", "matrix_provenance"}, {"archs4"},
                                    {"matrix", "dnam_reference"}, {"matrix", "dnam_reference", "dnam_supplement"},
                                    {"matrix", "dnam_idat_provenance"},
                                    {"matrix", "dnam_reference", "dnam_idat_provenance"}):
                raise ValueError("Invalid discovered molecular-input declaration")
            _, sample_receipt = request(accession + "-samples", f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=gsm&form=text&view=brief")
            candidate = {"study_id": accession, "modality": plan["modality"], "sources": {
                "series": {"path": receipt["response_path"], "sha256": receipt["sha256"]},
                "samples": {"path": sample_receipt["response_path"], "sha256": sample_receipt["sha256"]}, **molecular}}
            candidates.append(candidate)
            record.update(status="READY_FOR_INGESTION_ATTEMPT", reason="Source date and available operational overlap checks completed; scientific design unreviewed")
        except Exception as error:
            record.update(status="UNRESOLVED", reason=type(error).__name__ + ": " + str(error))
    shared = {sample: sorted(owners) for sample, owners in batch_samples.items() if len(owners) > 1}
    related = {study for owners in shared.values() for study in owners}
    for record in ledger.values():
        if record["accession"] in related:
            record.update(status="SHARED_NEW_CANDIDATE_SAMPLES", reason="Related new candidates require a joint source-family assessment")
    candidates = [r for r in candidates if r["study_id"] not in related]
    with (output / "candidate-ledger.jsonl").open("x") as handle:
        for uid in ids: handle.write(json.dumps(ledger[uid], sort_keys=True) + "\n")
    write_new(output / "candidates.json", {"schema": "clockbase-update-candidates-v1", "studies": candidates,
                                            "scope": "Operational input candidates; all scientific judgments and original-paper eligibility remain unverified"})
    write_new(output / "shared-new-samples.json", shared)
    write_new(output / "archs4-catalog-assessment.json", catalog_record)
    for path, expected in input_hashes.items():
        if sha(path) != expected: raise UpdateIntegrityError("Operational seen source changed during discovery")
    if sha(Path(__file__)) != config["implementation_sha256"]: raise UpdateIntegrityError("Discovery code changed during execution")
    if sha(parser_path) != config["source_parser_sha256"]: raise UpdateIntegrityError("Discovery source parser changed during execution")
    counts = {}
    for record in ledger.values(): counts[record["status"]] = counts.get(record["status"], 0) + 1
    result = {"schema": "clockbase-update-discovery-result-v1", "status": "COMPLETED", "query": term,
              "query_snapshot_source": query_snapshot_source,
              "query_total_uids": total, "enumerated_uids": len(ids), "unretrieved_uid_count": total - len(ids),
              "query_fully_enumerated": len(ids) == total, "processed_uids": len(ids) - counts.get("NOT_YET_PROCESSED", 0),
              "processed_uids_in_this_execution": len(selected),
              "unprocessed_enumerated_uids": counts.get("NOT_YET_PROCESSED", 0), "candidate_status_counts": counts,
              "ready_for_ingestion_attempt": len(candidates), "requests_used": request_count,
              "response_bytes": response_bytes, "elapsed_seconds": time.monotonic() - started,
              "response_byte_accounting": "Observed bytes, with full request cap retained when failure leaves size unknown",
              "seen_registry_accessions": len(seen_studies), "seen_registry_accessions_with_sample_evidence": len(covered_seen),
              "seen_sample_ids": len(seen_samples), "input_sha256": input_hashes,
              "seen_registry_scope": "Available operational corpus, not original submitted-paper membership",
              "scientific_status": "CANNOT VERIFY", "original_paper_2025_eligibility": "QUARANTINED",
              "paid_model_calls": 0, "molecular_downloads": 0,
              "output_sha256": {str(p): sha(p) for p in output.rglob("*") if p.is_file()}}
    seal(complete, result)
    return result


def run(plan_path, output, until="report", *, executor=run_child):
    if until not in STAGES: raise ValueError("Unknown pipeline stopping stage")
    plan = json.loads(plan_path.read_text())
    limits, options, audit_options = validate_plan(plan)
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".update.lock").open("a") as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError("Another update owns this registry")
        control = {"schema": SCHEMA, "max_total_model_calls": limits["max_total_model_calls"],
                   "max_total_cost_usd": limits["max_total_cost_usd"], "concurrency": 1}
        control_path = output / "control.json"
        if control_path.exists():
            if read_sealed(control_path) != control: raise UpdateIntegrityError("Registry budget changed; use a new registry")
        else: seal(control_path, control)
        candidates = list(plan.get("studies", []))
        discovery = None
        if plan.get("discovery"):
            discovery = file_source(plan["discovery"], output / "discovery", max_bytes=limits["max_input_bytes"], network=limits["allow_network"])
            discovered = json.loads(Path(discovery["path"]).read_text())
            if discovered.get("schema") != "clockbase-update-candidates-v1":
                raise ValueError("Discovery snapshot requires source-candidate manifest schema")
            candidates += discovered["studies"]
        if not candidates or len(candidates) > limits["max_studies"]: raise ValueError("Candidate count exceeds declared bound")
        keys = [(r.get("modality"), r.get("study_id")) for r in candidates]
        if len(keys) != len(set(keys)): raise ValueError("Duplicate study/modality in candidate manifest")
        batch_key = identity({"plan": plan, "until": until, "discovery_sha256": discovery["sha256"] if discovery else None})
        batch = output / "batches" / batch_key
        batch.mkdir(parents=True, exist_ok=True)
        implementation = implementation_hashes()
        summaries, input_files = [], {}
        for entry in candidates:
            study, modality = entry.get("study_id"), entry.get("modality")
            if not isinstance(study, str) or not re.fullmatch(r"GSE[0-9]+", study) or modality not in MODALITIES:
                raise ValueError("Invalid source study/modality identity")
            record = {"study_id": study, "modality": modality, "status": "RUNNING",
                      "scientific_status": "CANNOT VERIFY", "original_paper_2025_eligibility": "QUARANTINED",
                      "eligibility_reason": "Operational ingestion does not establish original corpus boundary or family exclusions"}
            try:
                sources = collect_sources(entry, batch / "acquisition" / study, limits, input_files)
                validate_replay_sources(entry, sources, options)
                if modality == "human_dnam" and STAGES.index(until) >= STAGES.index("preprocess") and "dnam_reference" not in sources:
                    raise ValueError("Declare a pinned DNAm reference before preprocessing; implicit downloads are disabled")
                if sum(input_files.values()) > limits["max_input_bytes"]: raise ValueError("Aggregate source input byte cap exceeded")
                source_design = inspect_source_design(study, modality, sources, limits["max_samples_per_study"])
                version_identity = {"study_id": study, "modality": modality,
                    "source_sha256": {role: source["sha256"] for role, source in sorted(sources.items())},
                    "target": entry.get("target", "Source-defined intervention comparisons"),
                    "pipeline_options": options, "comprehensive_audit": audit_options,
                    "implementation": implementation, "python": sys.version,
                    "explicit_analysis_revision": entry.get("analysis_revision", "initial")}
                if plan.get("comprehensive_audit_policy", "literal_fields_v1") != "literal_fields_v1":
                    version_identity["comprehensive_audit_policy"] = plan["comprehensive_audit_policy"]
                version_id = identity(version_identity)
                study_root = output / "studies" / modality / study
                versions = sorted(p.name for p in study_root.glob("*") if (p / "version.json").is_file())
                version_root = study_root / version_id
                version_path = version_root / "version.json"
                if version_path.exists():
                    version = read_sealed(version_path)
                    if version["identity"] != version_identity: raise UpdateIntegrityError("Study version identity changed")
                    verify_implementation_snapshot(version_root, version["implementation"])
                    verify_sources(version["sources"])
                    record["version_event"] = "EXISTING_VERSION"
                else:
                    supersedes = entry.get("supersedes_version")
                    if versions and supersedes not in versions:
                        raise ValueError("Changed source/configuration requires explicit supersedes_version")
                    version = {**version_identity, "identity": version_identity, "version_id": version_id,
                               "sources": sources, "source_design": source_design, "supersedes_version": supersedes,
                               "created_at_utc": datetime.now(timezone.utc).isoformat()}
                    seal(version_path, version)
                    snapshot_implementation(version_root, implementation)
                    record["version_event"] = "NEW_VERSION" if versions else "FIRST_VERSION"
                record.update(version_id=version_id, version_root=str(version_root),
                              first_public_date=source_design["series_first_public_date"],
                              first_public_source_quote=source_design["source_status_quote"],
                              source_sample_count=len(source_design["sample_ids"]))
                pipeline_root = version_root / "pipeline"
                phase = "pipeline-" + until
                complete = version_root / (phase + ".json")
                if complete.exists():
                    outcome = verify_completion(complete)
                else:
                    for prior_stage in STAGES:
                        prior_path = version_root / ("pipeline-" + prior_stage + ".json")
                        if prior_path.exists() and verify_completion(prior_path)["status"] == "FAILED":
                            raise RuntimeError("Earlier child failed in this version; explicit new analysis revision required")
                    model_phase = STAGES.index(until) >= STAGES.index("design")
                    if model_phase:
                        reserve(output, version_id, "pipeline", options["max_calls"], options["max_cost_usd"], limits)
                    request_path = version_root / (phase + "-request.json")
                    if request_path.exists():
                        raise RuntimeError("Unfinished child retained; no automatic duplicate execution")
                    result_path = version_root / (phase + "-result.json")
                    request = {"phase": "pipeline", "version_path": str(version_path), "pipeline_root": str(pipeline_root),
                               "cache": str(output / "cache"), "until": until, "result_path": str(result_path)}
                    write_new(request_path, request)
                    try:
                        code = executor(request_path, version_root / (phase + ".log"), limits["max_child_seconds"])
                        if code or not result_path.exists(): raise RuntimeError("Existing per-study pipeline stopped; inspect preserved child log")
                        hashes = {str(p): sha(p) for p in pipeline_root.rglob("*") if p.is_file() and p.name != ".run.lock"}
                        hashes[str(result_path)] = sha(result_path)
                        outcome = {"status": "COMPLETED", "until": until, "output_sha256": hashes}
                    except Exception as error:
                        outcome = {"status": "FAILED", "until": until, "error": type(error).__name__ + ": " + str(error)}
                    seal(complete, outcome)
                record["pipeline"] = outcome
                record["status"] = "PIPELINE_COMPLETED_THROUGH_" + until.upper() if outcome["status"] == "COMPLETED" else "FAILED"
                if until == "report" and outcome["status"] == "COMPLETED":
                    audit_path = version_root / "comprehensive-audit.json"
                    if audit_path.exists():
                        audited = verify_completion(audit_path)
                    else:
                        reserve(output, version_id, "comprehensive_audit", audit_options["max_calls"], audit_options["max_cost_usd"], limits)
                        audit_request = {"study_id": study, "modality": modality, "version_id": version_id,
                            "pipeline_root": str(pipeline_root), "source_bindings": sources,
                            "output": str(version_root / "comprehensive-audit"), "budget": audit_options,
                            "scientific_status": "CANNOT VERIFY", "original_paper_2025_eligibility": "QUARANTINED"}
                        if plan.get("comprehensive_audit_policy", "literal_fields_v1") != "literal_fields_v1":
                            audit_request["comprehensive_audit_policy"] = plan["comprehensive_audit_policy"]
                        request_path = version_root / "comprehensive-audit-request.json"
                        if request_path.exists(): raise RuntimeError("Unfinished comprehensive audit retained; no duplicate call")
                        result_path = version_root / "comprehensive-audit-result.json"
                        write_new(request_path, {"phase": "comprehensive_audit", "version_path": str(version_path),
                                               "audit_request": audit_request, "result_path": str(result_path)})
                        try:
                            code = executor(request_path, version_root / "comprehensive-audit.log", limits["max_child_seconds"])
                            if code or not result_path.exists(): raise RuntimeError("Comprehensive audit unavailable or stopped; inspect retained log")
                            report = json.loads(result_path.read_text())
                            hashes = {str(p): sha(p) for p in (version_root / "comprehensive-audit").rglob("*") if p.is_file()}
                            hashes[str(result_path)] = sha(result_path)
                            audited = {"status": "COMPLETED" if report.get("execution_status") == "COMPLETED" else "UNRESOLVED",
                                       "result": report, "output_sha256": hashes}
                        except Exception as error:
                            audited = {"status": "UNRESOLVED", "error": type(error).__name__ + ": " + str(error)}
                        seal(audit_path, audited)
                    record["comprehensive_audit"] = audited
                    record["status"] = "AUDIT_EXECUTION_COMPLETED" if audited["status"] == "COMPLETED" else "PENDING_COMPREHENSIVE_AUDIT"
                verify_sources(sources)
            except Exception as error:
                record.update(status="UNRESOLVED", error=type(error).__name__ + ": " + str(error))
            summaries.append(record)
        reservation_records = [read_sealed(p) for p in (output / "reservations").glob("*.json")]
        result = {"schema": "clockbase-continuous-update-result-v1", "plan_sha256": sha(plan_path),
                  "until": until, "studies": summaries, "input_bytes": sum(input_files.values()),
                  "reserved_model_calls": sum(r["max_calls"] for r in reservation_records),
                  "reserved_cost_usd": sum(r["max_cost_usd"] for r in reservation_records),
                  "budget_accounting": "Lifetime worst-case reservations; failed or ambiguous billing is never released automatically",
                  "scientific_status": "CANNOT VERIFY", "original_paper_2025_eligibility": "QUARANTINED",
                  "status": "COMPLETED" if all(r["status"].startswith(("PIPELINE_COMPLETED", "AUDIT_EXECUTION_COMPLETED")) for r in summaries) else "COMPLETED_WITH_UNRESOLVED_CASES"}
        receipts = batch / "receipts"; receipts.mkdir(exist_ok=True)
        seal(receipts / (f"{len(list(receipts.glob('*.json'))):04d}.json"), result)
        return result


def _refresh_seen(plan, registry, directory):
    """Keep saved analyses and pipeline attempts separate from metadata catalogues."""
    members, sample_records, inputs = {}, [], {}
    for role in ("seen_worklist", "seen_sample_sources"):
        declaration = plan.get(role)
        if declaration is None: continue
        source = file_source(declaration, directory, max_bytes=2_000_000_000, network=False)
        inputs[source["path"]] = source["sha256"]
        with (gzip.open(source["path"], "rt") if source["path"].endswith(".gz") else open(source["path"])) as handle:
            if role == "seen_worklist":
                for line in handle:
                    row = json.loads(line); study = row.get("study_id", "")
                    if not re.fullmatch(r"GSE[0-9]+", study): continue
                    # A catalogue row containing only an accession is insufficient.
                    if not row.get("artifact_ids") and row.get("membership") not in {"SAVED_ANALYSIS", "PIPELINE_ATTEMPTED", "ANALYSIS_EXECUTED"}:
                        raise ValueError("Seen worklist must identify saved analyses or explicit execution membership")
                    members.setdefault(study, set()).add(row.get("membership", "SAVED_ANALYSIS"))
            else:
                sample_records = json.load(handle)
                if not isinstance(sample_records, list): raise ValueError("Seen sample source must be a record list")
    for version_path in sorted(registry.glob("studies/*/*/*/version.json")):
        version = read_sealed(version_path)
        study = version["study_id"]
        if not re.fullmatch(r"GSE[0-9]+", study): raise UpdateIntegrityError("Noncanonical study in refresh registry")
        requests = sorted(version_path.parent.glob("pipeline-*-request.json"))
        if not requests: continue  # Metadata or frozen configuration alone is not an attempt.
        inputs[str(version_path)] = sha(version_path)
        inputs[str(version_path) + ".sha256"] = sha(str(version_path) + ".sha256")
        for request in requests: inputs[str(request)] = sha(request)
        category = "PIPELINE_ATTEMPTED"
        for phase in ("agent", "audit", "report"):
            receipt = version_path.parent / ("pipeline-" + phase + ".json")
            if receipt.exists():
                record = verify_completion(receipt)
                inputs[str(receipt)] = sha(receipt)
                inputs[str(receipt) + ".sha256"] = sha(str(receipt) + ".sha256")
                inputs.update(record.get("output_sha256", {}))
                if record["status"] == "COMPLETED": category = "ANALYSIS_EXECUTED"
        members.setdefault(study, set()).add(category)
        sample_records.append({"accession": study, "sample_ids": version["source_design"]["sample_ids"],
                               "source_kind": category, "version_sha256": sha(version_path)})
    filtered = [row for row in sample_records if row.get("accession") in members]
    directory.mkdir(parents=True, exist_ok=False)
    with (directory / "studies.jsonl").open("x") as handle:
        for study, categories in sorted(members.items()):
            membership = next(label for label in ("ANALYSIS_EXECUTED", "PIPELINE_ATTEMPTED", "SAVED_ANALYSIS") if label in categories)
            handle.write(json.dumps({"study_id": study, "membership": membership,
                                     "membership_categories": sorted(categories)}, sort_keys=True) + "\n")
    write_new(directory / "sample-sources.json", filtered)
    return {"seen_worklist": {"path": str(directory / "studies.jsonl"), "sha256": sha(directory / "studies.jsonl")},
            "seen_sample_sources": {"path": str(directory / "sample-sources.json"), "sha256": sha(directory / "sample-sources.json")},
            "membership_categories": {study: sorted(values) for study, values in sorted(members.items())},
            "input_sha256": inputs, "scope": "Operational saved analyses and attempted executions; metadata-only catalogue entries do not establish membership"}


def _refresh_update_external(result):
    """Keep IDAT graph files in cached refresh checks even when run() is skipped."""
    external = {}
    for row in result["studies"]:
        if not row.get("version_root"):
            continue
        version_root = Path(row["version_root"])
        for artifact in version_root.rglob("*"):
            if artifact.is_file() and artifact.name != ".run.lock":
                external[str(artifact)] = sha(artifact)
        version = read_sealed(version_root / "version.json")
        for role, source in version["sources"].items():
            if role == "dnam_idat_provenance" or role.startswith("dnam_idat_provenance:"):
                if sha(source["path"]) != source["sha256"]:
                    raise UpdateIntegrityError("DNAm IDAT source changed before refresh completion: " + role)
                external[source["path"]] = source["sha256"]
    return external


def _refresh_phase(root, name, inputs, function):
    """Retain one terminal phase; no hidden retry after failure or interruption."""
    directory = root / "phases" / name
    receipt_path = directory / "phase.json"
    if receipt_path.exists():
        receipt = read_sealed(receipt_path)
        if receipt["input_sha256"] != inputs: raise UpdateIntegrityError("Refresh phase inputs differ: " + name)
        observed = {str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file()
                    and p.name not in {".run.lock", ".handoff.lock"} and p not in {receipt_path, Path(str(receipt_path) + ".sha256")}}
        if observed != set(receipt["output_sha256"]): raise UpdateIntegrityError("Refresh phase inventory changed: " + name)
        for relative, expected in receipt["output_sha256"].items():
            if sha(directory / relative) != expected: raise UpdateIntegrityError("Refresh phase output changed: " + name + "/" + relative)
        for path, expected in receipt.get("external_sha256", {}).items():
            if sha(path) != expected: raise UpdateIntegrityError("Refresh external output changed: " + path)
        return receipt
    directory.mkdir(parents=True, exist_ok=True)
    started = directory / "started.json"
    interrupted = started.exists()
    if not interrupted:
        seal(started, {"phase": name, "input_sha256": inputs, "started_at_utc": datetime.now(timezone.utc).isoformat()})
    elif read_sealed(started)["input_sha256"] != inputs:
        raise UpdateIntegrityError("Interrupted refresh phase belongs to another input: " + name)
    began = time.monotonic()
    try:
        if interrupted: raise RuntimeError("Interrupted phase retained; a new explicit refresh revision is required")
        value, external = function(directory)
        receipt = {"status": "COMPLETED", "result": value, "external_sha256": external}
    except Exception as error:
        receipt = {"status": "UNRESOLVED", "error": type(error).__name__ + ": " + str(error), "external_sha256": {}}
    receipt.update(schema="clockbase-refresh-phase-v1", phase=name, input_sha256=inputs,
                   elapsed_seconds=time.monotonic() - began,
                   output_sha256={str(p.relative_to(directory)): sha(p) for p in directory.rglob("*")
                                  if p.is_file() and p.name not in {".run.lock", ".handoff.lock"}},
                   scientific_status="CANNOT VERIFY")
    seal(receipt_path, receipt)
    return receipt


def _refresh_dnam_inputs(discovery_summary, reference, output, limits, *, transport=None, prior_studies=None, input_route=None):
    """Resolve source-declared human-DNAm matrices under the refresh network remainder.

    This is an input-provenance gate.  It neither infers beta semantics from a
    numeric range nor treats source samples as independent biological units.
    """
    from src.pipeline_inputs_2026_09_10 import (beta_table_header, parse_samples, propose_dnam_recovery,
                                                read_matrix, reconstruct_dnam_input)
    from audit_system.report_2026_09_10 import verify_dnam_recovery
    if input_route not in {None, DNAM_IDAT_ROUTE}:
        raise ValueError("Unknown DNAm refresh input route")

    discovery_summary, output = Path(discovery_summary).resolve(), Path(output).resolve()
    discovery = verify_completion(discovery_summary)
    source_root = discovery_summary.parent
    ledger_path = source_root / "candidate-ledger.jsonl"
    if not ledger_path.is_file() or sha(ledger_path) != discovery["output_sha256"].get(str(ledger_path)):
        raise UpdateIntegrityError("DNAm refresh requires the complete hash-bound discovery ledger")
    ledger = [json.loads(line) for line in ledger_path.read_text().splitlines() if line.strip()]
    if len(ledger) != discovery["enumerated_uids"]:
        raise UpdateIntegrityError("DNAm discovery ledger denominator differs")
    output.mkdir(parents=True, exist_ok=False)
    requests = output / "requests"; requests.mkdir()
    max_requests = number(limits["max_requests"], "DNAm remaining requests", integer=True)
    max_bytes = number(limits["max_total_response_bytes"], "DNAm remaining bytes", integer=True)
    per_response = number(limits["max_response_bytes"], "DNAm response bytes", integer=True, minimum=1)
    max_seconds = number(limits["max_seconds"], "DNAm seconds", integer=True)
    interval = number(limits["minimum_interval_seconds"], "DNAm request interval", minimum=.35)
    started, last_request = time.monotonic(), 0.0
    request_count = response_bytes = reconstructed_bytes = 0

    class OfficialSourceNoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, request_object, file_pointer, http_code, message, headers, new_url):
            raise ValueError("DNAm official source redirected; no unverified destination contacted")

    def request(name, url, *, byte_limit=None):
        nonlocal request_count, response_bytes, last_request
        if not re.fullmatch(r"GSE[0-9]+-(?:samples|matrix-listing|matrix|beta-supplement|idat-GSM[0-9]+-(?:Grn|Red))", name):
            raise ValueError("Unsafe DNAm request identity")
        parsed = urllib.parse.urlparse(url)
        if (parsed.scheme != "https" or parsed.hostname not in {"www.ncbi.nlm.nih.gov", "ftp.ncbi.nlm.nih.gov"}
                or parsed.username is not None or parsed.password is not None or parsed.fragment):
            raise ValueError("DNAm refresh accepts only bounded official HTTPS sources")
        if request_count >= max_requests or time.monotonic() - started >= max_seconds:
            raise RuntimeError("DNAm request/time remainder exhausted")
        remaining = max_bytes - response_bytes
        cap = min(per_response, remaining - 1)
        if byte_limit is not None:
            cap = min(cap, number(byte_limit, "DNAm source byte limit", integer=True, minimum=1))
        if cap <= 0: raise RuntimeError("DNAm aggregate response-byte remainder exhausted")
        destination = requests / name; destination.mkdir(exist_ok=False)
        request_path = destination / "request.json"
        receipt_path = destination / "receipt.json"
        request_record = {"url": url, "byte_limit": cap, "network_scope": "official GEO input only"}
        seal(request_path, request_record)
        wait = interval - (time.monotonic() - last_request)
        if wait > 0: time.sleep(wait)
        request_count += 1; last_request = time.monotonic()
        receipt = {**request_record, "requested_at_utc": datetime.now(timezone.utc).isoformat(), "status": "RUNNING"}
        observed = None
        try:
            if transport is None:
                remaining_seconds = max_seconds - (time.monotonic() - started)
                if remaining_seconds <= 0: raise TimeoutError("DNAm source wall-time limit reached")
                with source_deadline(min(60, remaining_seconds)):
                    with urllib.request.build_opener(OfficialSourceNoRedirect).open(
                            url, timeout=min(60, remaining_seconds)) as response:
                        payload = response.read(cap + 1)
            else:
                payload = transport(url, cap)
            observed = len(payload)
            if not payload or observed > cap: raise ValueError("Empty or oversized DNAm source response")
            response_path = destination / "response.bin"
            with response_path.open("xb") as handle: handle.write(payload)
            receipt.update(status="PASS", bytes=observed, accounted_bytes=observed,
                           sha256=sha(response_path), response_path=str(response_path))
            return response_path, receipt
        except Exception as error:
            receipt.update(status="FAIL", bytes=0, error=type(error).__name__ + ": " + str(error),
                           accounted_bytes=observed if observed is not None else cap + 1)
            raise
        finally:
            response_bytes += receipt["accounted_bytes"]
            seal(receipt_path, receipt)

    def series_fields(path):
        fields = {}
        for line in Path(path).read_text(errors="strict").splitlines():
            key, separator, value = line.partition(" = ")
            if separator and key.startswith("!Series_"):
                fields.setdefault(key[8:], []).append(value)
        return fields

    def beta_evidence(series, samples):
        evidence = []
        allowed = {"summary", "overall_design", "data_processing", "value_description"}
        for role, path in (("series", series), ("samples", samples)):
            for line in Path(path).read_text(errors="strict").splitlines():
                key, separator, value = line.partition(" = ")
                field = key.removeprefix("!Series_").removeprefix("!Sample_")
                if separator and field in allowed and re.search(r"(?i)(?<![A-Za-z])beta(?:[ -]?value)?s?(?![A-Za-z])", value):
                    evidence.append({"source": role, "field": field, "quote": value})
        return evidence

    def declaration(path):
        path = Path(path).resolve()
        return {"path": str(path), "sha256": sha(path)}

    def materialize(response_path, destination):
        destination = Path(destination)
        if destination.exists(): raise FileExistsError("Refusing to replace a DNAm source payload: " + str(destination))
        shutil.copyfile(response_path, destination)
        if sha(destination) != sha(response_path): raise UpdateIntegrityError("DNAm named source copy differs")
        return destination

    def bounded_table(path):
        """Reject compressed expansion beyond the declared source-input ceiling."""
        uncompressed = 0
        with gzip.open(path, "rb") as source:
            while block := source.read(1024 * 1024):
                uncompressed += len(block)
                if uncompressed > limits["max_reconstructed_bytes_per_study"]:
                    raise ValueError("DNAm source table exceeds decompressed input byte limit")
        return uncompressed

    def matrix_accession(path):
        identities = []
        with gzip.open(path, "rt", newline="") as source:
            for line in source:
                if line.strip() == "!series_matrix_table_begin": break
                if line.startswith("!Series_geo_accession\t"):
                    # GEO Series Matrix metadata uses a tab-separated, quoted
                    # accession. Decode that exact record with the CSV parser;
                    # stripping quote characters would accept malformed input.
                    fields = next(csv.reader([line], delimiter="\t", strict=True))
                    if len(fields) != 2:
                        raise ValueError("Official matrix study identity has an ambiguous GEO record")
                    identities.append(fields[1])
                else:
                    key, separator, value = line.partition(" = ")
                    if separator and key in {"^SERIES", "!Series_geo_accession"}:
                        identities.append(value.strip())
        if not identities or set(identities) != {study}:
            raise ValueError("Official matrix study identity differs from source Series")
        return identities

    candidates, rows, external = [], [], {}
    idat_input_files = {reference["sha256"]: reference["size_bytes"]} if input_route == DNAM_IDAT_ROUTE and reference else {}
    prior_studies = prior_studies or {}
    source_eligible = {"UNRESOLVED_NO_ASSAY_INPUT", "READY_FOR_INGESTION_ATTEMPT"}
    for source_row in ledger:
        study = source_row.get("accession")
        row = {"geo_uid": source_row["geo_uid"], "study_id": study, "discovery_status": source_row["status"],
               "status": source_row["status"], "reason": source_row.get("reason", ""),
               "sample_ids": source_row.get("sample_ids", []), "scientific_status": "CANNOT VERIFY",
               "biological_unit_status": "CANNOT VERIFY", "clock_applicability": "CANNOT VERIFY"}
        if not study or source_row["status"] not in source_eligible:
            rows.append(row); continue
        if study in prior_studies:
            row.update(status="PRIOR_DNAM_REFRESH_TERMINAL", reason="A prior sealed DNAm refresh row exists; no automatic retry or duplicate revision",
                       prior_terminal=prior_studies[study])
            rows.append(row); continue
        if reference is None:
            row.update(status="DNAM_REFERENCE_UNRESOLVED", reason="No hash-bound local DNAm reference was configured")
            rows.append(row); continue
        try:
            if request_count + (1 if input_route == DNAM_IDAT_ROUTE else 3) > max_requests:
                raise RuntimeError("Insufficient request remainder for the minimum complete DNAm source route")
            series_receipt = source_row.get("source_receipt", {})
            series_path = Path(series_receipt.get("response_path", "")).resolve()
            if (not series_path.is_file() or series_receipt.get("sha256") != source_row.get("series_source_sha256")
                    or sha(series_path) != series_receipt.get("sha256")):
                raise UpdateIntegrityError("Discovery Series source is missing or changed")
            fields = series_fields(series_path)
            if fields.get("type") != ["Methylation profiling by array"]:
                raise ValueError("DNAm refresh requires a single declared Series assay type: Methylation profiling by array; observed " + repr(fields.get("type", [])))
            sample_url = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?" +
                          urllib.parse.urlencode({"acc": study, "targ": "gsm", "form": "text", "view": "brief"}))
            samples_path, samples_receipt = request(study + "-samples", sample_url)
            samples = parse_samples(samples_path.read_text(), study)
            sample_ids = [sample["GSM"] for sample in samples]
            if (len(sample_ids) != len(set(sample_ids)) or set(sample_ids) != set(row["sample_ids"])
                    or any(not re.fullmatch(r"GSM[0-9]+", value) for value in sample_ids)):
                raise ValueError("Series/Sample GSM membership differs or is noncanonical")
            if len(sample_ids) > limits["max_samples_per_study"]:
                raise ValueError("DNAm source sample count exceeds the explicit per-study limit")
            # parse_samples retains repeated physical fields as field.1, .2,
            # etc. Check every declaration rather than selecting the first.
            def physical_values(field, *, casefold=False):
                values = set()
                for sample in samples:
                    declared = {" ".join(str(value).split()) for key, value in sample.items()
                                if re.fullmatch(re.escape(field) + r"(?:\.[0-9]+)?", key)}
                    if casefold: declared = {value.casefold() for value in declared}
                    if len(declared) != 1 or not next(iter(declared)):
                        raise ValueError("DNAm sample has absent or contradictory " + field + " declarations")
                    values.update(declared)
                return values
            if physical_values("organism_ch1", casefold=True) != {"homo sapiens"}:
                raise ValueError("DNAm refresh requires literal Homo sapiens sample declarations")
            molecules = {str(sample.get("molecule_ch1", "")).strip() for sample in samples}
            if physical_values("molecule_ch1", casefold=True) != {"genomic dna"}:
                raise ValueError("DNAm refresh requires literal genomic DNA sample declarations")
            platforms = physical_values("platform_id")
            if len(platforms) != 1 or None in platforms or not all(re.fullmatch(r"GPL[0-9]+", value) for value in platforms):
                raise ValueError("DNAm refresh requires one canonical sample platform")
            if set(fields.get("platform_id", [])) != platforms:
                raise ValueError("Series and Sample platform declarations differ")
            if input_route == DNAM_IDAT_ROUTE:
                from src.pipeline_inputs_2026_09_10 import (DNAM_IDAT_LIMITS, DNAM_IDAT_SUPPORTED_PROFILE,
                    acquire_dnam_idat, dnam_idat_source_pairs, prepare_dnam_idat_plan, process_dnam_idat)
                # Resolve the complete source-owned roster before any raw-file request.
                pairs = dnam_idat_source_pairs(study, series_path, samples_path)
                if request_count + 2 * len(pairs) > max_requests:
                    raise RuntimeError("Insufficient request remainder for every declared IDAT channel")
                output_cap = min(limits["max_reconstructed_bytes_per_study"],
                                 limits["max_total_reconstructed_bytes"] - reconstructed_bytes,
                                 limits["max_input_bytes"] - sum(idat_input_files.values())
                                 - sum(path.stat().st_size for path in (series_path, samples_path)
                                       if sha(path) not in idat_input_files))
                processing_seconds = min(DNAM_IDAT_LIMITS["processing_seconds"], limits["max_processing_seconds"])
                if output_cap <= 0: raise ValueError("DNAm IDAT reconstruction byte remainder exhausted")
                bounded = {"compressed_bytes_per_file": min(DNAM_IDAT_LIMITS["compressed_bytes_per_file"], per_response, output_cap),
                           "compressed_bytes_total": min(DNAM_IDAT_LIMITS["compressed_bytes_total"], max_bytes-response_bytes, output_cap),
                           "decompressed_bytes_per_file": min(DNAM_IDAT_LIMITS["decompressed_bytes_per_file"], output_cap),
                           "samples": min(DNAM_IDAT_LIMITS["samples"], limits["max_samples_per_study"]),
                           "processing_seconds": processing_seconds}
                study_root = output / study; study_root.mkdir(exist_ok=False)
                idat_plan = study_root / "idat-plan.json"
                prepare_dnam_idat_plan(study, series_path, samples_path, idat_plan,
                    profile_name=DNAM_IDAT_SUPPORTED_PROFILE["name"], limits=bounded)
                row.update(input_route=input_route, idat_profile=DNAM_IDAT_SUPPORTED_PROFILE["name"],
                           idat_resource_limits={"plan": bounded, "max_retained_study_bytes": output_cap,
                                                 "max_preparation_seconds": processing_seconds})
                urls = {declaration["https_url"]: study + "-idat-" + pair["sample_id"] + "-" + channel
                        for pair in pairs for channel, declaration in pair["channels"].items()}
                def retained_bytes():
                    return sum(path.stat().st_size for path in study_root.rglob("*") if path.is_file())
                def idat_fetch(url, destination, limit):
                    if url not in urls: raise ValueError("IDAT callback URL is not in the frozen source roster")
                    response_path, receipt = request(urls[url], url, byte_limit=limit)
                    materialize(response_path, destination)
                    write_new(Path(str(destination)+".source.json"), {
                        "url": url, "retrieved_unix": time.time(), "sha256": receipt["sha256"],
                        "size_bytes": receipt["bytes"]})
                    return destination
                # The enclosing preparation deadline also covers decompression and independent readback.
                # Per-request official transport has its own tighter remaining network deadline.
                with source_deadline(processing_seconds, "DNAm IDAT preparation"):
                    acquire_dnam_idat(idat_plan, study_root / "acquisition", fetch_callback=idat_fetch,
                                     max_output_bytes=output_cap-retained_bytes())
                    process_cap = output_cap-retained_bytes()
                    if process_cap <= 0: raise ValueError("DNAm IDAT processing output byte remainder exhausted")
                    processed = process_dnam_idat(idat_plan, study_root / "acquisition/acquisition.json",
                                                 study_root / "processing", max_output_bytes=process_cap)
                    candidate_sources = {"series": declaration(series_path), "samples": declaration(samples_path),
                        "matrix": declaration(processed["matrix"]["path"]),
                        "dnam_idat_provenance": dict(processed["dnam_idat_provenance"]),
                        "dnam_reference": {"path": reference["path"], "sha256": reference["sha256"]}}
                    binding = read_dnam_idat_provenance(candidate_sources["dnam_idat_provenance"]["path"],
                        candidate_sources["dnam_idat_provenance"]["sha256"], study, sources=candidate_sources)
                    if set(binding["sample_ids"]) != set(sample_ids):
                        raise ValueError("Processed IDAT/GEO sample membership differs")
                    combined_input_files = dict(idat_input_files)
                    for source in binding["files"].values():
                        combined_input_files[source["sha256"]] = Path(source["path"]).stat().st_size
                    if sum(combined_input_files.values()) > limits["max_input_bytes"]:
                        raise ValueError("DNAm IDAT source graph exceeds aggregate update input byte limit")
                    size = retained_bytes()
                    if size > output_cap: raise ValueError("DNAm IDAT retained outputs exceed reconstruction byte limits")
                    reconstructed_bytes += size
                    idat_input_files = combined_input_files
                    for source in binding["files"].values():
                        if source["path"] in external and external[source["path"]] != source["sha256"]:
                            raise UpdateIntegrityError("Conflicting IDAT source graph binding")
                        external[source["path"]] = source["sha256"]
                candidates.append({"study_id": study, "modality": "human_dnam", "sources": candidate_sources,
                                   "target": "Source-defined intervention comparisons"})
                row.update(status="READY_FOR_UPDATE", reason="Source-declared IDATs and derived matrix provenance verified",
                           sample_ids=sample_ids, platform_id=next(iter(platforms)), molecule_ch1=sorted(molecules),
                           author_method_reproduction=False, numeric_validation=binding["measurement"],
                           source_receipts={"samples": samples_receipt}, sources=candidate_sources)
                rows.append(row)
                continue
            prefix = study[:-3] + "nnn"
            listing_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{study}/matrix/"
            listing_path, listing_receipt = request(study + "-matrix-listing", listing_url)
            listing = listing_path.read_text(errors="strict")
            names = sorted(set(re.findall(r'''href=["']([A-Za-z0-9_.+-]+\.txt\.gz)["']''', listing)))
            if len(names) != 1:
                raise ValueError("Official GEO matrix listing does not contain exactly one unambiguous txt.gz matrix")
            matrix_url = listing_url + urllib.parse.quote(names[0], safe="._+-")
            matrix_response, matrix_receipt = request(study + "-matrix", matrix_url)
            study_root = output / study; study_root.mkdir(exist_ok=False)
            matrix_path = materialize(matrix_response, study_root / names[0])
            bounded_table(matrix_path)
            matrix_accession(matrix_path)
            layout = beta_table_header(matrix_path)
            original_ids = layout["header"][1:] if layout["feature_column"] else layout["header"]
            if len(original_ids) != len(set(original_ids)) or set(original_ids) != set(sample_ids):
                raise ValueError("Matrix header does not exactly cover source GSM membership")
            metadata = study_root / "metadata.csv"
            with metadata.open("x", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=sorted({key for sample in samples for key in sample}))
                writer.writeheader(); writer.writerows(samples)
            evidence = beta_evidence(series_path, samples_path)
            candidate_sources = {"series": declaration(series_path), "samples": declaration(samples_path),
                                 "matrix": declaration(matrix_path), "dnam_reference": {"path": reference["path"], "sha256": reference["sha256"]}}
            validation = None
            if layout["empty"]:
                supplement_urls = []
                expected_prefix = f"/geo/series/{prefix}/{study}/suppl/"
                for value in fields.get("supplementary_file", []):
                    parsed = urllib.parse.urlparse(value)
                    name = Path(parsed.path).name
                    if (parsed.scheme in {"ftp", "https"} and parsed.hostname == "ftp.ncbi.nlm.nih.gov"
                            and parsed.path == expected_prefix + name and not parsed.query and not parsed.fragment
                            and re.search(r"(?i)(?<![A-Za-z])betas?(?![A-Za-z])", name)
                            and re.fullmatch(r"[A-Za-z0-9_.+-]+\.(?:(?:txt|tsv|csv)(?:\.gz)?|xlsx)", name)):
                        supplement_urls.append("https://ftp.ncbi.nlm.nih.gov" + parsed.path)
                supplement_urls = sorted(set(supplement_urls))
                if len(supplement_urls) != 1:
                    raise ValueError("Empty matrix requires exactly one source-declared beta-named tabular supplement")
                if request_count >= max_requests:
                    raise RuntimeError("Insufficient request remainder for source-declared DNAm beta supplement")
                supplement_response, supplement_receipt = request(study + "-beta-supplement", supplement_urls[0])
                supplement_path = materialize(supplement_response, study_root / Path(urllib.parse.urlparse(supplement_urls[0]).path).name)
                if supplement_path.suffix == ".gz": bounded_table(supplement_path)
                elif supplement_path.suffix == ".xlsx":
                    with zipfile.ZipFile(supplement_path) as archive:
                        if sum(item.file_size for item in archive.infolist()) > limits["max_reconstructed_bytes_per_study"]:
                            raise ValueError("DNAm workbook exceeds decompressed input byte limit")
                elif supplement_path.stat().st_size > limits["max_reconstructed_bytes_per_study"]:
                    raise ValueError("DNAm source table exceeds input byte limit")
                proposal = propose_dnam_recovery(study, matrix_path, supplement_path, metadata, series_path,
                                                 sample_source=samples_path)
                validation_dir = study_root / "reconstruction-validation"; validation_dir.mkdir()
                reconstructed = reconstruct_dnam_input(proposal, validation_dir)
                reconstructed_size = Path(reconstructed["path"]).stat().st_size
                if (reconstructed_size > limits["max_reconstructed_bytes_per_study"]
                        or reconstructed_bytes + reconstructed_size > limits["max_total_reconstructed_bytes"]):
                    raise ValueError("Recovered DNAm input exceeds explicit reconstruction byte limits")
                reconstructed_bytes += reconstructed_size
                verified = verify_dnam_recovery(proposal, reconstructed)
                if verified["status"] != "PASS": raise ValueError("Source-bound DNAm reconstruction did not verify")
                candidate_sources["dnam_supplement"] = declaration(supplement_path)
                evidence = [{"source": "series", "field": "supplementary_file", "quote": supplement_urls[0],
                             "evidence_type": "source_declared_beta_named_file"}]
                validation = {"proposal_sha256": identity(proposal), "reconstruction": reconstructed,
                              "independent_verification": verified, "supplement_receipt": supplement_receipt}
            else:
                if not evidence:
                    raise ValueError("Nonempty matrix lacks a literal source beta-unit declaration")
                frame = read_matrix(matrix_path, beta=True)
                if set(frame.columns) != set(sample_ids) or len(frame.columns) != len(sample_ids):
                    raise ValueError("Parsed beta matrix membership differs from source")
                validation = {"feature_count": len(frame), "sample_count": len(frame.columns),
                              "numeric_range_consistency": "PASS", "source_beta_evidence": evidence}
            candidate = {"study_id": study, "modality": "human_dnam", "sources": candidate_sources,
                         "target": "Source-defined intervention comparisons"}
            candidates.append(candidate)
            row.update(status="READY_FOR_UPDATE", reason="Source profile, beta semantics and exact GSM input binding passed",
                       sample_ids=sample_ids, platform_id=next(iter(platforms)), molecule_ch1=sorted(molecules),
                       beta_semantics_evidence=evidence, numeric_validation=validation,
                       source_receipts={"samples": samples_receipt, "matrix_listing": listing_receipt, "matrix": matrix_receipt},
                       sources=candidate_sources)
        except Exception as error:
            row.update(status="DNAM_INPUT_UNRESOLVED", reason=type(error).__name__ + ": " + str(error))
            if input_route == DNAM_IDAT_ROUTE:
                # Failed outputs remain charged; a later study cannot reclaim their storage allowance.
                reconstructed_bytes = sum(path.stat().st_size for folder in output.iterdir()
                    if folder.is_dir() and re.fullmatch(r"GSE[0-9]+", folder.name)
                    for path in folder.rglob("*") if path.is_file())
        rows.append(row)
    candidates_path = output / "candidates.json"
    seal(candidates_path, {"schema": "clockbase-update-candidates-v1", "studies": candidates,
                           "scope": "Operational source-bound DNAm inputs; scientific acceptance and original-paper eligibility remain unresolved"})
    result = {"schema": "clockbase-dnam-refresh-input-resolution-v1", "status": "COMPLETED",
              "studies": rows, "candidates": {"path": str(candidates_path), "sha256": sha(candidates_path)},
              "ready_studies": len(candidates), "requests_used": request_count, "response_bytes": response_bytes,
              "reconstructed_bytes": reconstructed_bytes,
              "reference": reference, "source_files": {str(p): sha(p) for p in output.rglob("*") if p.is_file()},
              "scientific_status": "CANNOT VERIFY", "original_paper_2025_eligibility": "QUARANTINED",
              "limitations": ["Numeric beta-range consistency is not beta-unit evidence",
                              "GSM membership does not establish independent biological units",
                              "No preprocessing, imputation, clock applicability or scientific result is accepted"]}
    if input_route == DNAM_IDAT_ROUTE:
        result.update(input_route=input_route, external_sha256=external)
    write_new(output / "resolution.json", result)
    return result


def _refresh_source_acquisition(plan, plan_path, output, implementations, runtime, *, transport=None):
    """Resolve optional source acquisition once, with lifetime reserved resources.

    A completed immutable receipt is the sole resume path. Failed/interrupted
    attempts keep their reservations and partial objects; no transport retry.
    """
    from src.pipeline_inputs_2026_09_10 import acquire_archs4, acquire_dnam_reference, DNAM_REFERENCE_URL
    spec = plan['source_acquisition']
    if not isinstance(spec, dict) or set(spec) != {'policy', 'cache_dir', 'limits', 'budget'}:
        raise ValueError('Source acquisition requires policy, cache_dir, limits and budget')
    expected_policy = 'official_archs4_v2_latest' if plan['modality'] in RNA_SPECIES else 'geo_gse40279'
    if spec['policy'] != expected_policy or any(plan.get(k) is not None for k in ('archs4', 'dnam_reference')):
        raise ValueError('Acquisition policy must match modality and be exclusive with pinned local sources')
    cap = spec['limits']; budget = spec['budget']
    cap_keys = {'max_requests', 'max_bytes', 'max_seconds', 'max_source_bytes', 'max_metadata_bytes', 'min_free_bytes'}
    if not isinstance(cap, dict) or set(cap) != cap_keys or not isinstance(budget, dict) or set(budget) != {'max_seconds', 'max_cache_bytes'}:
        raise ValueError('Explicit source acquisition request/byte/time/storage limits required')
    for key, value in list(cap.items())+list(budget.items()): number(value, 'source acquisition '+key, integer=True, minimum=1)
    if cap['max_seconds'] > budget['max_seconds']:
        raise ValueError('Acquisition time allocation exceeds lifetime time allowance')
    if (cap['max_requests']+plan['discovery_limits']['max_requests'] > plan['network_budget']['max_requests']
            or cap['max_bytes']+plan['discovery_limits']['max_total_response_bytes'] > plan['network_budget']['max_bytes']):
        raise RuntimeError('Lifetime source acquisition network reservation cannot also cover discovery')
    if not isinstance(spec['cache_dir'], str) or not Path(spec['cache_dir']).is_absolute():
        raise ValueError('Source acquisition requires an absolute persistent cache directory')
    cache = Path(spec['cache_dir']).resolve(); output = Path(output).resolve()
    if cache == output or cache.is_relative_to(output) or output.is_relative_to(cache):
        raise ValueError('Source cache and refresh output must be separate directories')
    output.mkdir(parents=True, exist_ok=True); cache.mkdir(parents=True, exist_ok=True)
    run_id = identity(plan); directory = output/'source-acquisitions'/run_id
    terminal = directory/'result.json'
    bound = {'plan_sha256': sha(plan_path), 'implementation': implementations, 'runtime': runtime,
             'specification': spec, 'modality': plan['modality']}
    def cache_bytes():
        total = 0
        for path in cache.rglob('*'):
            if path.is_symlink(): raise UpdateIntegrityError('Symlinks are not admitted in the immutable source cache')
            if path.is_file(): total += path.stat().st_size
        return total
    with (output/'.refresh.lock').open('a') as run_lock, (cache/'.cache.lock').open('a') as cache_lock:
        for lock in (run_lock, cache_lock):
            try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError: raise RuntimeError('Another process owns source acquisition or its cache')
        controls = {
            output/'control.json': {'schema': REFRESH_SCHEMA, 'network_budget': plan['network_budget'],
                'max_total_model_calls': plan['update_limits']['max_total_model_calls'],
                'max_total_cost_usd': plan['update_limits']['max_total_cost_usd']},
            output/'acquisition-control.json': {'schema': 'clockbase-acquisition-control-v1',
                'cache_dir': str(cache), 'budget': budget},
            cache/'control.json': {'schema': 'clockbase-source-cache-control-v1', 'max_cache_bytes': budget['max_cache_bytes']}}
        for path, value in controls.items():
            if path.exists():
                if read_sealed(path) != value: raise UpdateIntegrityError('Source acquisition lifetime limits changed')
            else:
                if path.parent == cache:
                    required = len((json.dumps(value, indent=2, sort_keys=True)+'\n').encode())+65
                    if (cache_bytes()+required > budget['max_cache_bytes']
                            or shutil.disk_usage(cache).free < required+cap['min_free_bytes']):
                        raise RuntimeError('Source cache control exceeds storage/free-space allowance')
                seal(path, value)
        if cache_bytes() > budget['max_cache_bytes']: raise RuntimeError('Source cache storage cap exhausted')
        if terminal.exists():
            saved = read_sealed(terminal)
            if _refresh_binding(saved['input']) != _refresh_binding(bound): raise UpdateIntegrityError('Source acquisition input/code/runtime changed')
            observed = {str(p.relative_to(directory)) for p in directory.rglob('*') if p.is_file()
                        and p not in {terminal, Path(str(terminal)+'.sha256')}}
            if observed != set(saved['output_sha256']): raise UpdateIntegrityError('Acquisition output inventory changed')
            for name, checksum in saved['output_sha256'].items():
                if sha(directory/name) != checksum: raise UpdateIntegrityError('Acquisition output changed: '+name)
            for name, checksum in saved['external_sha256'].items():
                if sha(name) != checksum: raise UpdateIntegrityError('Cached source or acquisition binding changed: '+name)
            if saved['status'] != 'COMPLETED': raise RuntimeError('Prior source acquisition remains unresolved; no implicit retry')
            return {'path': str(terminal), 'sha256': sha(terminal), 'value': saved}
        reservations = output/'network-reservations'; reservations.mkdir(exist_ok=True)
        reservation_path = reservations/(run_id+'-source.json')
        reservation = {'refresh_id': run_id, 'scope': 'source_acquisition', 'max_requests': cap['max_requests'],
                       'max_bytes': cap['max_bytes'], 'max_seconds': cap['max_seconds']}
        if reservation_path.exists():
            if read_sealed(reservation_path) != reservation: raise UpdateIntegrityError('Acquisition reservation changed')
        else:
            prior = [read_sealed(path) for path in reservations.glob('*.json')]
            if any(sum(r[k] for r in prior)+reservation[k] > plan['network_budget'][k] for k in ('max_requests', 'max_bytes')):
                raise RuntimeError('Lifetime source acquisition network reservation limit reached')
            if sum(r.get('max_seconds', 0) for r in prior)+reservation['max_seconds'] > budget['max_seconds']:
                raise RuntimeError('Lifetime source acquisition time reservation limit reached')
            seal(reservation_path, reservation)
        directory.mkdir(parents=True, exist_ok=True)
        started_path = directory/'started.json'; interrupted = started_path.exists()
        if not interrupted: seal(started_path, {'input': bound, 'started_at_utc': datetime.now(timezone.utc).isoformat()})
        elif _refresh_binding(read_sealed(started_path)['input']) != _refresh_binding(bound): raise UpdateIntegrityError('Interrupted acquisition input differs')
        began = time.monotonic(); requests = 0; transferred = 0; last_request = 0.; external = {}

        def storage_check(extra=0):
            if cache_bytes()+extra > budget['max_cache_bytes']:
                raise RuntimeError('Source cache storage cap exhausted')
            if shutil.disk_usage(cache).free < extra+cap['min_free_bytes']:
                raise RuntimeError('Source cache free-space reserve would be violated')

        class NoRedirect(urllib.request.HTTPRedirectHandler):
            def redirect_request(self, request_object, file_pointer, http_code, message, headers, new_url):
                raise ValueError('Bounded acquisition does not follow redirects')

        def request(method, url, destination, maximum, headers):
            nonlocal requests, transferred, last_request
            allowed = {'https://archs4.org/api/versionfile', DNAM_REFERENCE_URL,
                       *{f'https://s3.k8s.maayanlab.cloud/archs4/files/{s}_gene_v2.latest.h5' for s in RNA_SPECIES.values()}}
            if method not in {'HEAD', 'GET'} or url not in allowed:
                raise ValueError('Unsupported bounded source request')
            if requests >= cap['max_requests'] or time.monotonic()-began >= cap['max_seconds']:
                raise RuntimeError('Source acquisition request/time allocation exhausted')
            remaining = cap['max_bytes']-transferred
            if remaining <= 0: raise RuntimeError('Source acquisition transfer-byte allocation exhausted')
            byte_limit = min(maximum, remaining-1) if method == 'GET' else 0
            if method == 'GET' and byte_limit < maximum:
                raise RuntimeError('Source acquisition transfer-byte allocation cannot cover request')
            if method == 'GET' and (destination is None or byte_limit <= 0): raise ValueError('Missing bounded payload destination')
            if destination is not None:
                destination = Path(destination).resolve()
                if not (destination.is_relative_to(directory) or destination.is_relative_to(cache)):
                    raise ValueError('Source acquisition destination escapes bound directories')
                destination.parent.mkdir(parents=True, exist_ok=True)
                if destination.exists(): raise FileExistsError('Source payload or partial file already exists')
                storage_check(byte_limit+2*cap['max_metadata_bytes']+1048576 if destination.is_relative_to(cache) else 0)
            requests += 1
            request_root = directory/'requests'/f'{requests:03d}'; request_root.mkdir(parents=True)
            value = {'schema': 'clockbase-bounded-source-response-v1', 'method': method, 'url': url,
                     'request_headers': headers, 'maximum_payload_bytes': byte_limit,
                     'requested_at_utc': datetime.now(timezone.utc).isoformat(), 'status': 'RUNNING',
                     'destination_at_request': str(destination) if destination else None,
                     'automatic_retries': 0, 'maximum_connections': 1}
            seal(request_root/'request.json', value)
            observed = 0; response = None; exhausted_stream = False
            try:
                wait = plan['discovery_limits']['minimum_interval_seconds']-(time.monotonic()-last_request)
                if wait > 0: time.sleep(wait)
                last_request = time.monotonic()
                timeout = min(60, max(.001, cap['max_seconds']-(time.monotonic()-began)))
                if transport is None:
                    response = urllib.request.build_opener(NoRedirect).open(
                        urllib.request.Request(url, headers=headers, method=method), timeout=timeout)
                    response_headers = dict(response.headers.items())
                else:
                    response_headers, payload = transport(method, url, headers, byte_limit)
                    response = io.BytesIO(payload)
                # Retain verification metadata only; never archive response cookies or session headers.
                response_headers = {str(k).lower(): str(v) for k, v in response_headers.items()
                                    if str(k).lower() in {'etag', 'content-length', 'content-type', 'last-modified'}}
                value['response_header_policy'] = 'source-verification-allowlist-v1'
                value['response_headers'] = response_headers
                if method == 'GET':
                    with destination.open('xb') as handle:
                        while True:
                            block = response.read(min(1048576, byte_limit-observed+1))
                            if not block: exhausted_stream = True; break
                            observed += len(block)
                            if observed > byte_limit: raise ValueError('Source payload exceeds the explicit response cap')
                            if destination.is_relative_to(cache): storage_check(len(block))
                            handle.write(block)
                    if not observed: raise ValueError('Empty bounded source payload')
                    value['payload_sha256'] = sha(destination)
                else:
                    exhausted_stream = True
                value.update(status='PASS', observed_payload_bytes=observed, accounted_bytes=observed)
            except Exception as error:
                value.update(status='FAIL', error=type(error).__name__+': '+str(error),
                             observed_payload_bytes=observed,
                             accounted_bytes=observed if exhausted_stream else max(observed, byte_limit+1))
                raise
            finally:
                if response is not None: response.close()
                transferred += value['accounted_bytes']
                value['finished_at_utc'] = datetime.now(timezone.utc).isoformat()
                seal(request_root/'receipt.json', value)
            return {'headers': response_headers, 'bytes': observed,
                    'receipt': {'path': str(request_root/'receipt.json'), 'sha256': sha(request_root/'receipt.json')}}

        saved = {'schema': 'clockbase-refresh-source-acquisition-v1', 'status': 'UNRESOLVED', 'input': bound,
                 'scientific_status': 'CANNOT VERIFY', 'original_paper_2025_eligibility': 'QUARANTINED',
                 'response_byte_accounting': 'Observed response-body bytes; the full request allowance and overflow sentinel are retained when a failed transfer leaves its size uncertain',
                 'means_computed': False, 'model_calls': 0, 'clock_runs': 0, 'preprocessing_runs': 0}
        try:
            if interrupted: raise RuntimeError('Interrupted source acquisition retained; no implicit retry')
            storage_check()
            with source_deadline(cap['max_seconds']):
                snapshot_implementation(directory, implementations)
                if plan['modality'] in RNA_SPECIES:
                    source = acquire_archs4(RNA_SPECIES[plan['modality']], cache, request=request,
                        selection_dir=directory/'selection', max_source_bytes=cap['max_source_bytes'],
                        metadata_bytes=cap['max_metadata_bytes'])
                    role = 'archs4'
                else:
                    source = acquire_dnam_reference(cache, request=request, max_source_bytes=cap['max_source_bytes'])
                    role = 'dnam_reference'
                source = source.resolve(); source_sha = sha(source)
                for path in source.parent.iterdir():
                    if path.is_file(): external[str(path)] = sha(path)
                # Bind cached objects to their original archived response receipts,
                # including cache hits from a previous explicit refresh revision.
                for path in source.parent.glob('*.json'):
                    record = json.loads(path.read_text())
                    for key in ('response_receipt', 'post_download_head_receipt'):
                        if key in record:
                            original = record[key]; receipt_path = Path(original['path'])
                            if sha(receipt_path) != original['sha256']: raise UpdateIntegrityError('Cached origin receipt changed')
                            origin = read_sealed(receipt_path)
                            if origin['status'] != 'PASS': raise UpdateIntegrityError('Cached source lacks a completed origin response')
                            expected_url = record.get('url', record.get('declaration', {}).get('url'))
                            if origin.get('url') != expected_url or origin.get('method') != ('GET' if key == 'response_receipt' else 'HEAD'):
                                raise UpdateIntegrityError('Cached acquisition origin URL or method differs')
                            if key == 'response_receipt' and (origin.get('payload_sha256') != source_sha
                                    or origin.get('observed_payload_bytes') != source.stat().st_size):
                                raise UpdateIntegrityError('Cached source differs from its archived acquisition payload')
                            request_path = receipt_path.with_name('request.json'); original_request = read_sealed(request_path)
                            if any(original_request.get(field) != origin.get(field) for field in
                                   ('method', 'url', 'request_headers', 'maximum_payload_bytes', 'requested_at_utc')):
                                raise UpdateIntegrityError('Cached origin request and response differ')
                            external[str(receipt_path)] = original['sha256']
                            external[str(receipt_path)+'.sha256'] = sha(str(receipt_path)+'.sha256')
                            external[str(request_path)] = sha(request_path)
                            external[str(request_path)+'.sha256'] = sha(str(request_path)+'.sha256')
                            origin_terminal = receipt_path.parents[2]/'result.json'
                            if origin_terminal != terminal:
                                ancestor = read_sealed(origin_terminal)
                                if (ancestor.get('schema') != 'clockbase-refresh-source-acquisition-v1'
                                        or ancestor.get('status') != 'COMPLETED'
                                        or ancestor.get('resolved_source') != {'path': str(source), 'sha256': source_sha}
                                        or ancestor.get('external_sha256', {}).get(str(path)) != sha(path)):
                                    raise UpdateIntegrityError('Cached object differs from its original completed acquisition')
                                for name, digest in ancestor['output_sha256'].items():
                                    old_path = origin_terminal.parent/name
                                    if sha(old_path) != digest: raise UpdateIntegrityError('Cached origin archive changed')
                                    external[str(old_path)] = digest
                                external[str(origin_terminal)] = sha(origin_terminal)
                                external[str(origin_terminal)+'.sha256'] = sha(str(origin_terminal)+'.sha256')
                storage_check(); verify_implementation(implementations)
                saved.update(status='COMPLETED', role=role,
                    resolved_source={'path': str(source), 'sha256': source_sha}, size_bytes=source.stat().st_size)
        except Exception as error:
            saved['error'] = type(error).__name__+': '+str(error)
        for path in [reservation_path, *controls]: external[str(path)] = sha(path)
        saved.update(requests_used=requests, response_bytes=transferred, elapsed_seconds=time.monotonic()-began,
                     finished_at_utc=datetime.now(timezone.utc).isoformat(), external_sha256=external,
                     output_sha256={str(p.relative_to(directory)): sha(p) for p in directory.rglob('*') if p.is_file()})
        seal(terminal, saved)
        if saved['status'] != 'COMPLETED': raise RuntimeError('Source acquisition unresolved: '+saved['error'])
        return {'path': str(terminal), 'sha256': sha(terminal), 'value': saved}


def refresh(plan_path, output, until="report", *, transport=None, executor=run_child, acquisition_transport=None):
    """One bounded RNA or guarded human-DNAm future refresh.

    A config revision owns immutable discovery/extraction/readback/handoff
    receipts. All revisions share network reservations and the existing update
    registry's lifetime model budget. Other modalities retain discovery evidence
    without claiming a supported assay route.
    """
    import importlib.metadata
    import platform
    import types
    from audit_system import archs4_update_2026_09_10 as intake

    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    plan = json.loads(plan_path.read_text())
    allowed = {"schema", "revision", "modality", "query", "date_from", "date_to", "candidate_offset",
               "archs4", "max_release_bytes", "dnam_reference", "seen_worklist", "seen_sample_sources", "discovery_limits",
               "network_budget", "extraction_limits", "update_limits", "pipeline", "comprehensive_audit",
               "comprehensive_audit_policy", "source_acquisition", "dnam_input_route"}
    if set(plan) - allowed or plan.get("schema") != REFRESH_SCHEMA or plan.get("modality") not in MODALITIES:
        raise ValueError("Unsupported refresh schema, fields or modality")
    if "dnam_input_route" in plan and (plan["modality"] != "human_dnam" or plan["dnam_input_route"] != DNAM_IDAT_ROUTE):
        raise ValueError("The source-declared IDAT v2 route is supported only for human_dnam")
    if until not in STAGES or not isinstance(plan.get("revision"), str) or not plan["revision"].strip():
        raise ValueError("Refresh requires an explicit revision and supported stopping stage")
    update_plan = {"schema": SCHEMA, "limits": plan["update_limits"], "pipeline": plan.get("pipeline", {}),
                   "comprehensive_audit": plan["comprehensive_audit"]}
    if "comprehensive_audit_policy" in plan:
        update_plan["comprehensive_audit_policy"] = plan["comprehensive_audit_policy"]
    limits, _, _ = validate_plan(update_plan)
    if limits["allow_network"]: raise ValueError("Refresh child inputs must be local; discovery owns the network allowance")
    discovery_limits = plan["discovery_limits"]
    expected_discovery = {"max_requests", "max_response_bytes", "max_total_response_bytes", "max_seconds",
                          "max_query_ids", "max_candidates", "minimum_interval_seconds"}
    if set(discovery_limits) != expected_discovery: raise ValueError("Explicit discovery request, byte, time and candidate limits required")
    for key in expected_discovery - {"minimum_interval_seconds"}:
        number(discovery_limits[key], "discovery." + key, integer=True, minimum=1)
    number(discovery_limits["minimum_interval_seconds"], "discovery interval", minimum=.35)
    network = plan["network_budget"]
    if set(network) != {"max_requests", "max_bytes"}: raise ValueError("Explicit lifetime network limits required")
    for key in network: number(network[key], "network." + key, integer=True, minimum=1)
    extraction = plan["extraction_limits"]
    extraction_keys = {"max_studies", "max_samples_per_study", "max_extracted_bytes_per_study", "max_total_extracted_bytes", "concurrency"}
    if set(extraction) != extraction_keys or extraction["concurrency"] != 1:
        raise ValueError("Explicit extraction limits and concurrency one required")
    for key in extraction_keys - {"concurrency"}: number(extraction[key], "extraction." + key, integer=True, minimum=1)
    if discovery_limits["max_candidates"] > min(extraction["max_studies"], limits["max_studies"]):
        raise ValueError("Discovery candidate allowance exceeds extraction/update study bounds")
    runtime = {"python": platform.python_version(), "packages": {
        name: importlib.metadata.version(name) for name in ("h5py", "numpy", "pandas")}}
    implementations = implementation_hashes()
    if plan.get("dnam_input_route") == DNAM_IDAT_ROUTE:
        worker = "src/dnam_idat_preprocess_2026_09_24.R"
        implementations[worker] = sha(ROOT / worker)
    implementations[REFRESH_VERIFIER] = sha(ROOT / REFRESH_VERIFIER)
    inputs = {str(plan_path): sha(plan_path)}
    for role in ("seen_worklist", "seen_sample_sources"):
        if plan.get(role):
            source = file_source(plan[role], output, max_bytes=2_000_000_000, network=False)
            inputs[source["path"]] = source["sha256"]
    acquisition = None
    local_archs4, local_reference = plan.get('archs4'), plan.get('dnam_reference')
    if plan['modality'] in RNA_SPECIES:
        number(plan.get('max_release_bytes'), 'max_release_bytes', integer=True, minimum=1)
    if plan.get('source_acquisition') is not None:
        acquisition_limits = plan['source_acquisition'].get('limits') if isinstance(plan['source_acquisition'], dict) else None
        source_cap = acquisition_limits.get('max_source_bytes') if isinstance(acquisition_limits, dict) else None
        number(source_cap, 'source acquisition max_source_bytes', integer=True, minimum=1)
        if source_cap > (plan['max_release_bytes'] if plan['modality'] in RNA_SPECIES else limits['max_input_bytes']):
            raise ValueError('Acquisition source-byte cap exceeds the downstream local-input allowance')
        acquisition = _refresh_source_acquisition(plan, plan_path, output, implementations, runtime, transport=acquisition_transport)
        inputs[acquisition['path']] = acquisition['sha256']
        inputs.update(acquisition['value']['external_sha256'])
        if plan['modality'] in RNA_SPECIES: local_archs4 = acquisition['value']['resolved_source']
        else: local_reference = acquisition['value']['resolved_source']
    release = dnam_reference = None
    if plan["modality"] in RNA_SPECIES:
        if not isinstance(local_archs4, dict) or set(local_archs4) != {"path", "sha256"}:
            raise ValueError("RNA refresh requires one pinned local ARCHS4 file")
        cap = number(plan.get("max_release_bytes"), "max_release_bytes", integer=True, minimum=1)
        source = file_source(local_archs4, output, max_bytes=cap, network=False)
        inputs[source["path"]] = source["sha256"]
        import h5py
        with h5py.File(source["path"], "r") as handle:
            shape = list(handle["data/expression"].shape)
            sample_count = handle["meta/samples/geo_accession"].shape[0]
            if len(shape) != 2 or not all(shape) or shape[1] != sample_count:
                raise ValueError("ARCHS4 expression and sample dimensions differ")
        release = {"path": source["path"], "sha256": source["sha256"], "size_bytes": source["bytes"],
                   "species": RNA_SPECIES[plan["modality"]], "expression_shape": shape, "sample_count": sample_count}
        if plan.get("dnam_reference") is not None:
            raise ValueError("RNA refresh cannot declare a DNAm reference")
    elif plan.get("archs4") or plan.get("max_release_bytes") is not None:
        raise ValueError("The automatic pinned-release route requires an RNA modality")
    elif plan["modality"] == "human_dnam" and local_reference is not None:
        if not isinstance(local_reference, dict) or set(local_reference) != {"path", "sha256"}:
            raise ValueError("DNAm reference requires exactly one local path and SHA256")
        source = file_source(local_reference, output, max_bytes=limits["max_input_bytes"], network=False)
        inputs[source["path"]] = source["sha256"]
        dnam_reference = {"path": source["path"], "sha256": source["sha256"], "size_bytes": source["bytes"]}
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".refresh.lock").open("a") as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError("Another refresh owns this output")
        control = {"schema": REFRESH_SCHEMA, "network_budget": network,
                   "max_total_model_calls": limits["max_total_model_calls"], "max_total_cost_usd": limits["max_total_cost_usd"]}
        control_path = output / "control.json"
        if control_path.exists():
            if read_sealed(control_path) != control: raise UpdateIntegrityError("Refresh lifetime limits changed; existing reservations must not be reset")
        else: seal(control_path, control)
        run_id = identity(plan)
        root = output / "refreshes" / run_id
        root.mkdir(parents=True, exist_ok=True)
        frozen = {"schema": REFRESH_SCHEMA, "plan": plan, "input_sha256": inputs, "implementation": implementations,
                  "runtime": runtime, "release": release, "dnam_reference": dnam_reference}
        if acquisition is not None: frozen['source_acquisition'] = {'path': acquisition['path'], 'sha256': acquisition['sha256']}
        version_path = root / "version.json"
        if version_path.exists():
            saved = read_sealed(version_path)
            if _refresh_binding(saved) != _refresh_binding(frozen): raise UpdateIntegrityError("Refresh input/code/runtime changed; use an explicit new revision")
            implementations = saved["implementation"]
            verify_implementation_snapshot(root, implementations)
        else:
            seal(version_path, frozen); snapshot_implementation(root, implementations)
        verifier_name = REFRESH_VERIFIER if REFRESH_VERIFIER in implementations else LEGACY_REFRESH_VERIFIER
        verifier_path = root / "implementation" / verifier_name
        verifier = types.ModuleType("clockbase_refresh_readback_" + run_id)
        verifier.__file__ = str(verifier_path)
        # Execute the hash-verified source bytes, without consulting a .pyc cache.
        exec(compile(verifier_path.read_bytes(), str(verifier_path), "exec"), verifier.__dict__)
        phase_inputs = {str(version_path): sha(version_path)}
        phases = {}

        def phase(name, function):
            verify_implementation(implementations)
            value = _refresh_phase(root, name, dict(phase_inputs), function)
            phases[name] = value
            path = root / "phases" / name / "phase.json"
            phase_inputs[str(path)] = sha(path)
            return value.get("result") if value["status"] == "COMPLETED" else None

        def source_binding(path):
            return {"path": str(Path(path).resolve()), "sha256": sha(path)}

        def save_result(directory, value, name="result.json"):
            path = directory / name; seal(path, value)
            return {"path": str(path), "sha256": sha(path), "value": value}, {}

        registry = output / "update-registry"
        def do_seen(directory):
            value = _refresh_seen(plan, registry, directory / "seen")
            saved, _ = save_result(directory, value)
            return saved, value["input_sha256"]
        seen = phase("seen-membership", do_seen)
        discovery = None
        if seen is not None:
            reservations = output / "network-reservations"; reservations.mkdir(exist_ok=True)
            reservation = {"refresh_id": run_id, "max_requests": discovery_limits["max_requests"],
                           "max_bytes": discovery_limits["max_total_response_bytes"]}
            reservation_path = reservations / (run_id + ".json")
            if reservation_path.exists():
                if read_sealed(reservation_path) != reservation: raise UpdateIntegrityError("Refresh network reservation changed")
            else:
                prior = [read_sealed(path) for path in reservations.glob("*.json")]
                if any(sum(r[k] for r in prior) + reservation[k] > network[k] for k in ("max_requests", "max_bytes")):
                    raise RuntimeError("Lifetime discovery reservation limit reached; prior and uncertain usage is retained")
                seal(reservation_path, reservation)
            phase_inputs[str(reservation_path)] = sha(reservation_path)

            def do_discovery(directory):
                discovery_plan = {"schema": "clockbase-update-discovery-plan-v1", "modality": plan["modality"],
                    "query": plan["query"], "date_from": plan["date_from"], "date_to": plan["date_to"],
                    "candidate_offset": plan.get("candidate_offset", 0), "limits": discovery_limits,
                    **{role: seen["value"][role] for role in ("seen_worklist", "seen_sample_sources")}}
                if release is not None:
                    discovery_plan.update(default_molecular_source={"archs4": local_archs4}, max_release_bytes=plan["max_release_bytes"])
                path = directory / "plan.json"; seal(path, discovery_plan)
                value = discover(path, directory / "source", transport=transport)
                return save_result(directory, value)

            discovery = phase("discovery", do_discovery)
        availability = profile = batch_profile = batch = readback = handoff = dnam_inputs = updated = None
        source_profile_script = verifier_path
        intake_exclusions = {}
        if discovery is not None and release is not None and discovery["value"]["ready_for_ingestion_attempt"]:
            discovery_path = root / "phases/discovery/source/summary.json"
            availability = phase("input-availability", lambda directory: save_result(directory,
                verifier.discovery_availability(discovery_path, Path(release["path"]), release["sha256"])))
            profile = phase("source-profile", lambda directory: save_result(directory, verifier.discovery_source_profiles(discovery_path), "source-profiles.json"))
        if profile is not None and availability is not None:
            availability_rows = {r["study_id"]: r for r in availability["value"]["rows"] if r["study_id"]}
            selected_rows = []
            for row in profile["value"]["studies"]:
                available = availability_rows[row["study_id"]]
                if not available["declared_as_available_input"]:
                    intake_exclusions[row["study_id"]] = {"disposition": "ASSAY_INPUT_UNRESOLVED", "reason": available["availability"]}
                elif row["sample_count"] > min(extraction["max_samples_per_study"], limits["max_samples_per_study"]):
                    intake_exclusions[row["study_id"]] = {"disposition": "SOURCE_SAMPLE_LIMIT_UNRESOLVED", "reason": "Source sample count exceeds the explicit per-study resource limit"}
                else:
                    selected_rows.append(row)
            if not intake_exclusions:
                batch_profile = profile
            elif selected_rows:
                def select_profile(directory):
                    value = {**profile["value"], "studies": selected_rows,
                        "supported_study_ids": [r["study_id"] for r in selected_rows if r["profile_status"] == "SUPPORTED_SOURCE_PROFILE"],
                        "stopped_study_ids": [r["study_id"] for r in selected_rows if r["profile_status"] != "SUPPORTED_SOURCE_PROFILE"],
                        "verification_script_sha256": implementations["src/update_pipeline_2026_09_20.py"],
                        "selection_scope": "Bounded extraction subset only; full discovery and literal source-profile receipts remain unchanged",
                        "parent_source_profile": {k: profile[k] for k in ("path", "sha256")},
                        "input_availability": {k: availability[k] for k in ("path", "sha256")}, "excluded_from_extraction": intake_exclusions}
                    molecules = {}
                    for row in selected_rows:
                        for sample in row["samples"]:
                            for molecule in sample["molecule_ch1"]: molecules[molecule] = molecules.get(molecule, 0) + 1
                    value["counts"] = {"enrolled_studies": len(selected_rows), "enrolled_samples": sum(r["sample_count"] for r in selected_rows),
                        "supported_studies": len(value["supported_study_ids"]), "stopped_studies": len(value["stopped_study_ids"]),
                        "molecule_label_occurrences": molecules}
                    value["source_files"] = {**value["source_files"], profile["path"]: profile["sha256"],
                        availability["path"]: availability["sha256"], str(verifier_path): sha(verifier_path)}
                    return save_result(directory, value, "source-profiles.json")
                batch_profile = phase("extraction-subset", select_profile)
                source_profile_script = root / "implementation/src/update_pipeline_2026_09_20.py"
        if batch_profile is not None and batch_profile["value"]["supported_study_ids"]:
            def do_batch(directory):
                batch_plan = {"schema": intake.BATCH_EXTRACTION_SCHEMA, "output": str(directory / "batch"),
                    "modality": plan["modality"],
                    "release": release, "runtime": runtime, "limits": extraction,
                    "selected_study_ids": batch_profile["value"]["supported_study_ids"],
                    "source_profile_descriptor": {**{k: batch_profile[k] for k in ("path", "sha256")},
                                                  "size_bytes": Path(batch_profile["path"]).stat().st_size},
                    "implementation": {role: {**source_binding(path), "size_bytes": Path(path).stat().st_size}
                                       for role, path in intake.BATCH_IMPLEMENTATIONS.items()}}
                path = directory / "plan.json"; seal(path, batch_plan)
                return save_result(directory, intake.run_archs4_batch_extraction(path))
            batch = phase("count-extraction", do_batch)
        if batch is not None:
            batch_path = root / "phases/count-extraction/batch/batch-summary.json"
            readback = phase("independent-source-readback", lambda directory: save_result(directory, verifier.verify_batch_extraction(batch_path)))
        if readback is not None and readback["value"]["status"] == "PASS":
            def do_handoff(directory):
                handoff_plan = {"schema": BATCH_HANDOFF_SCHEMA, "modality": plan["modality"],
                    "sources": {"batch_summary": source_binding(batch_path),
                    "verification_receipt": {k: readback[k] for k in ("path", "sha256")},
                    "source_profile_script": source_binding(source_profile_script)},
                    "limits": {"max_studies": extraction["max_studies"], "max_total_matrix_bytes": extraction["max_total_extracted_bytes"]}}
                path = directory / "plan.json"; seal(path, handoff_plan)
                return save_result(directory, build_batch_handoff(path, directory / "handoff"))
            handoff = phase("handoff", do_handoff)

        if discovery is not None and plan["modality"] == "human_dnam":
            prior_dnam, prior_hashes = {}, {}
            for prior_path in sorted((output / "refreshes").glob("*/phases/dnam-input-resolution/result.json")):
                if prior_path.is_relative_to(root): continue
                prior = read_sealed(prior_path)
                phase_path = prior_path.parent / "phase.json"
                if not phase_path.is_file(): raise UpdateIntegrityError("Prior DNAm refresh lacks its phase receipt")
                phase_receipt = read_sealed(phase_path)
                if phase_receipt.get("status") != "COMPLETED":
                    raise UpdateIntegrityError("Prior DNAm result belongs to a noncompleted phase")
                prior_wrapper = phase_receipt.get("result", {})
                if (set(prior_wrapper) != {"path", "sha256", "value"}
                        or Path(prior_wrapper["path"]).resolve() != prior_path.resolve()
                        or sha(prior_path) != prior_wrapper["sha256"] or prior_wrapper["value"] != prior):
                    raise UpdateIntegrityError("Prior DNAm refresh result binding changed")
                observed = {str(p.relative_to(prior_path.parent)) for p in prior_path.parent.rglob("*") if p.is_file()
                            and p not in {phase_path, Path(str(phase_path) + ".sha256")}}
                if observed != set(phase_receipt["output_sha256"]):
                    raise UpdateIntegrityError("Prior DNAm phase inventory changed")
                for relative, expected in phase_receipt["output_sha256"].items():
                    if Path(relative).is_absolute() or ".." in Path(relative).parts:
                        raise UpdateIntegrityError("Prior DNAm phase artifact escapes its source directory")
                    artifact = prior_path.parent / relative
                    if sha(artifact) != expected: raise UpdateIntegrityError("Prior DNAm phase output changed: " + relative)
                    prior_hashes[str(artifact)] = expected
                for path, expected in phase_receipt.get("external_sha256", {}).items():
                    if sha(path) != expected: raise UpdateIntegrityError("Prior DNAm external source changed: " + path)
                    prior_hashes[path] = expected
                prior_hashes[str(prior_path)] = sha(prior_path)
                prior_hashes[str(prior_path) + ".sha256"] = sha(str(prior_path) + ".sha256")
                prior_hashes[str(phase_path)] = sha(phase_path)
                prior_hashes[str(phase_path) + ".sha256"] = sha(str(phase_path) + ".sha256")
                for row in prior.get("studies", []):
                    study = row.get("study_id")
                    if study and study not in prior_dnam:
                        prior_dnam[study] = {"refresh_id": prior_path.parents[2].name, "status": row["status"],
                                             "result_sha256": sha(prior_path)}
            phase_inputs.update(prior_hashes)
            discovery_path = root / "phases/discovery/source/summary.json"
            remaining = {"max_requests": max(0, discovery_limits["max_requests"] - discovery["value"]["requests_used"]),
                         "max_response_bytes": discovery_limits["max_response_bytes"],
                         "max_total_response_bytes": max(0, discovery_limits["max_total_response_bytes"] - discovery["value"]["response_bytes"]),
                         "max_seconds": max(0, int(discovery_limits["max_seconds"] - discovery["value"]["elapsed_seconds"])),
                         "minimum_interval_seconds": discovery_limits["minimum_interval_seconds"],
                         "max_samples_per_study": min(extraction["max_samples_per_study"], limits["max_samples_per_study"]),
                         "max_reconstructed_bytes_per_study": extraction["max_extracted_bytes_per_study"],
                         "max_total_reconstructed_bytes": extraction["max_total_extracted_bytes"]}
            if plan.get("dnam_input_route") == DNAM_IDAT_ROUTE:
                remaining.update(max_processing_seconds=limits["max_child_seconds"], max_input_bytes=limits["max_input_bytes"])
            def do_dnam_inputs(directory):
                value = _refresh_dnam_inputs(discovery_path, dnam_reference, directory / "source", remaining,
                                             transport=transport, prior_studies=prior_dnam, input_route=plan.get("dnam_input_route"))
                saved, _ = save_result(directory, value)
                return saved, value.get("external_sha256", {})
            dnam_inputs = phase("dnam-input-resolution", do_dnam_inputs)

        candidate_manifest = None
        if handoff is not None and handoff["value"]["counts"]["candidate_studies"]:
            candidate_manifest = handoff["value"]["candidates"]
        elif dnam_inputs is not None and dnam_inputs["value"]["ready_studies"]:
            candidate_manifest = dnam_inputs["value"]["candidates"]
        if candidate_manifest is not None:
            def do_update(directory):
                path = directory / "plan.json"
                value = {**update_plan, "studies": [], "discovery": candidate_manifest}
                seal(path, value)
                result = run(path, registry, until, executor=executor)
                external = _refresh_update_external(result)
                saved, _ = save_result(directory, result)
                return saved, external
            updated = phase("update-" + until, do_update)
        ledger_path = root / "phases/discovery/source/candidate-ledger.jsonl"
        ledger = [json.loads(line) for line in ledger_path.read_text().splitlines() if line.strip()] if ledger_path.exists() else []
        if not ledger:
            search_path = root / "phases/discovery/source/requests/esearch/response.bin"
            if search_path.exists():
                try:
                    search = json.loads(search_path.read_text())["esearchresult"]
                    ledger = [{"geo_uid": uid, "accession": None, "status": "UNRESOLVED", "reason": "Discovery did not complete"} for uid in search["idlist"]]
                except (ValueError, KeyError, TypeError):
                    ledger = []  # Query denominator is unknown, not an empty discovery claim.
        profile_rows = {r["study_id"]: r for r in profile["value"]["studies"]} if profile else {}
        batch_rows = {r["study_id"]: r for r in batch["value"]["results"]} if batch else {}
        dnam_rows = {r["study_id"]: r for r in dnam_inputs["value"]["studies"] if r.get("study_id")} if dnam_inputs else {}
        update_rows = {r["study_id"]: r for r in updated["value"]["studies"]} if updated else {}
        dispositions = []
        failures = {name: record.get("error", "UNRESOLVED") for name, record in phases.items() if record["status"] != "COMPLETED"}
        for row in ledger:
            study = row.get("accession"); state, reason = row["status"], row.get("reason", "")
            if release is None and plan["modality"] != "human_dnam" and state in {"READY_FOR_INGESTION_ATTEMPT", "UNRESOLVED_NO_ASSAY_INPUT"}:
                state, reason = "UNSUPPORTED_REFRESH_ROUTE", "Automatic RNA refresh requires a pinned local ARCHS4 release; no unsupported assay route inferred"
            if study in dnam_rows:
                state, reason = dnam_rows[study]["status"], dnam_rows[study]["reason"]
            if study in profile_rows and profile_rows[study]["profile_status"] != "SUPPORTED_SOURCE_PROFILE":
                state, reason = "SOURCE_PROFILE_UNRESOLVED", "; ".join(profile_rows[study]["limitations"])
            elif study in batch_rows:
                state = batch_rows[study]["status"]
                reason = "Count extraction is not clock or scientific validation"
                if readback is None or readback["value"]["status"] != "PASS":
                    state, reason = "INDEPENDENT_SOURCE_VERIFICATION_UNRESOLVED", "Whole independent source readback did not PASS"
            if study in update_rows:
                state, reason = update_rows[study]["status"], update_rows[study].get("error", "Inspect source-bound pipeline and comprehensive audit receipts")
            elif state == "EXTRACTED":
                state, reason = "PIPELINE_UNRESOLVED", json.dumps(failures, sort_keys=True) if failures else "No completed update receipt"
            elif state in {"READY_FOR_INGESTION_ATTEMPT", "READY_FOR_UPDATE"} and failures:
                state, reason = "REFRESH_STAGE_UNRESOLVED", json.dumps(failures, sort_keys=True)
            if study in intake_exclusions:
                state, reason = intake_exclusions[study]["disposition"], intake_exclusions[study]["reason"]
            dispositions.append({"geo_uid": row["geo_uid"], "study_id": study, "discovery": row,
                "disposition": state, "reason": reason, "pipeline": update_rows.get(study),
                "scientific_status": "CANNOT VERIFY", "original_paper_2025_eligibility": "QUARANTINED"})
        verify_implementation(implementations)
        for path, expected in inputs.items():
            if sha(path) != expected: raise UpdateIntegrityError("Refresh source changed during execution: " + path)
        counts = {}
        for row in dispositions: counts[row["disposition"]] = counts.get(row["disposition"], 0) + 1
        network_reservations = [read_sealed(p) for p in (output / "network-reservations").glob("*.json")]
        model_reservations = [read_sealed(p) for p in (registry / "reservations").glob("*.json")]
        report_path = root / ("refresh-" + until + ".json")
        budget_files = {str(p): sha(p) for folder in (output / "network-reservations", registry / "reservations") for p in folder.glob("*.json")}
        if report_path.exists():
            prior_report = read_sealed(report_path)
            budget_files = prior_report["budget_snapshot_sha256"]
            for path, expected in budget_files.items():
                if sha(path) != expected: raise UpdateIntegrityError("Historical refresh budget reservation changed")
            network_reservations = [read_sealed(Path(p)) for p in budget_files if Path(p).parent == output / "network-reservations"]
            model_reservations = [read_sealed(Path(p)) for p in budget_files if Path(p).parent == registry / "reservations"]
        incomplete_query = discovery is None or not discovery["value"]["query_fully_enumerated"]
        summary = {"schema": "clockbase-future-refresh-result-v1", "refresh_id": run_id, "until": until,
            "status": "COMPLETED_WITH_UNRESOLVED_CASES" if failures or incomplete_query or any("UNRESOLVED" in r["disposition"] or r["disposition"] in {"FAILED", "NOT_YET_PROCESSED", "PENDING_COMPREHENSIVE_AUDIT", "UNSUPPORTED_REFRESH_ROUTE", "PRIOR_DNAM_REFRESH_TERMINAL"} for r in dispositions) else "COMPLETED",
            "query_total_uids": discovery["value"]["query_total_uids"] if discovery else None,
            "unretrieved_uid_count": discovery["value"]["unretrieved_uid_count"] if discovery else None,
            "enumerated_uids": len(dispositions), "dispositions": counts, "studies": dispositions,
            "phase_receipts": {name: source_binding(root / "phases" / name / "phase.json") for name in phases},
            "failures": failures, "input_sha256": inputs,
            "reserved_network_requests": sum(r["max_requests"] for r in network_reservations),
            "reserved_network_bytes": sum(r["max_bytes"] for r in network_reservations),
            "reserved_model_calls": sum(r["max_calls"] for r in model_reservations),
            "reserved_cost_usd": sum(r["max_cost_usd"] for r in model_reservations),
            "budget_snapshot_sha256": budget_files,
            "budget_accounting": "Lifetime worst-case reservations at receipt creation across refresh revisions; failed or uncertain usage is not released automatically",
            "scientific_status": "CANNOT VERIFY", "original_paper_2025_eligibility": "QUARANTINED"}
        if acquisition is not None: summary['source_acquisition'] = {'path': acquisition['path'], 'sha256': acquisition['sha256']}
        if report_path.exists():
            if read_sealed(report_path) != summary: raise UpdateIntegrityError("Completed refresh summary changed")
        else: seal(report_path, summary)
        return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--discover", type=Path, help="Bounded official GEO discovery plan; no model calls or molecular downloads")
    parser.add_argument("--batch-handoff", type=Path, help="Local sealed ARCHS4 batch/readback to update candidates; no model or HDF5 calls")
    parser.add_argument("--refresh", type=Path, help="One bounded config: GEO discovery, pinned RNA ARCHS4 or guarded DNAm input, pipeline and audit")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--until", choices=STAGES, default="report")
    parser.add_argument("--child-request", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    try:
        if sum(bool(value) for value in (args.plan, args.discover, args.batch_handoff, args.refresh, args.child_request)) != 1:
            raise ValueError("Choose exactly one of --plan, --discover, --batch-handoff, --refresh or --child-request")
        if args.refresh:
            if not args.output: raise ValueError("Refresh requires --output")
            result = refresh(args.refresh, args.output, args.until)
            print(json.dumps(result, indent=2, allow_nan=False)); return 0 if result["status"] == "COMPLETED" else 2
        if args.batch_handoff:
            if not args.output: raise ValueError("Batch handoff requires --output")
            result = build_batch_handoff(args.batch_handoff, args.output)
            print(json.dumps(result, indent=2, allow_nan=False)); return 0
        if args.child_request:
            child(json.loads(args.child_request.read_text())); return 0
        if args.discover:
            if args.plan or not args.output: raise ValueError("Discovery requires --output and excludes --plan")
            result = discover(args.discover.resolve(), args.output.resolve())
            print(json.dumps(result, indent=2)); return 0
        if not args.plan or not args.output: raise ValueError("--plan and --output are required")
        result = run(args.plan.resolve(), args.output.resolve(), args.until)
        print(json.dumps(result, indent=2, allow_nan=False))
        return 0 if result["status"] == "COMPLETED" else 2
    except Exception as error:
        print(f"Update controller stopped: {type(error).__name__}: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
