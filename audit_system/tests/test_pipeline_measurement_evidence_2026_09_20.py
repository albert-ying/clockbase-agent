from __future__ import annotations

import hashlib
import copy
import csv
import gzip
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from audit_system.pipeline_measurement_evidence_2026_09_20 import (
    SCHEMA, _recovered_dnam_source, _stage, build_measurement_evidence,
)


REPO = Path(__file__).resolve().parents[2]
UPDATE = REPO / "paper/revision_2026-09-07/evidence/continuous-audit-pipeline-2026-09-20/update-controller"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _complete(stage: Path, result: dict, files: list[Path]) -> None:
    outputs = {str(path.relative_to(stage)): _sha(path) for path in files}
    receipt = stage / "complete.json"
    _write_json(receipt, {"result": result, "outputs": outputs, "elapsed_seconds": 0.1})
    (stage / "complete.sha256").write_text(_sha(receipt) + "\n", encoding="ascii")


def _fixture(base: Path, *, archs4: bool = False) -> tuple[Path, dict, dict]:
    base.mkdir(parents=True, exist_ok=True)
    pipeline = base / "pipeline"
    samples = ["GSM1", "GSM2"]
    source = base / "source.csv"
    source.write_text(",GSM1,GSM2\nf1,1,2\nf2,3,4\n", encoding="utf-8")
    asset = base / "model.bin"
    asset.write_bytes(b"model")
    bindings = {"fixture_asset": {"path": str(asset), "sha256": _sha(asset)}}
    stage02 = pipeline / "02_molecular_input"
    raw = stage02 / "attempt-000/counts.csv"
    raw.parent.mkdir(parents=True)
    shutil.copyfile(source, raw)
    result02 = {"path": str(raw), "sha256": _sha(raw), "source": str(source)}
    stage02_files = [raw]
    if archs4:
        hdf5 = base / "archs4.h5"
        hdf5.write_bytes(b"sealed ARCHS4 fixture")
        bindings["archs4"] = {"path": str(hdf5), "sha256": _sha(hdf5)}
        result02.update({
            "feature_count": 2, "sample_count": 2, "source_hdf5": str(hdf5),
            "source_size": hdf5.stat().st_size,
            "all_study_memberships": {sample: ["GSEX"] for sample in samples},
            "source_file_sha256": None,
            "hash_scope": "Full extracted matrix hashed; HDF5 identity checked by version binding",
        })
        extraction = stage02 / "attempt-000/extraction.json"
        _write_json(extraction, result02)
        stage02_files.append(extraction)
    else:
        bindings["matrix"] = {"path": str(source), "sha256": _sha(source)}
    _complete(stage02, result02, stage02_files)

    stage04 = pipeline / "04_preprocessing"
    transformed = stage04 / "attempt-000/scaled_diff.csv"
    transformed.parent.mkdir(parents=True)
    transformed.write_text(",f1,f2\nGSM1,0.1,0.2\nGSM2,0.3,0.4\n", encoding="utf-8")
    raw_source = {"path": str(raw.resolve()), "sha256": _sha(raw)}
    asset_source = {"path": str(asset.resolve()), "sha256": _sha(asset)}
    prep = {
        "scaled": str(transformed), "scaled_diff": str(transformed), "status": "PASS",
        "sample_count": 2, "mapped_filtered_genes": 2, "clock_feature_columns": 2,
        "log_transform": "log10(count+1)", "normalization_profile": "clockbase_legacy",
        "filter_order": "mapping_then_filter", "filter_count_minimum": 10,
        "filter_sample_fraction": 0.2, "source_stratum": None,
        "measurement_sources": {"raw_counts": raw_source, "model_asset": asset_source},
        "historical_producing_version_verified": False,
        "output_sha256": {"scaled_diff.csv": _sha(transformed)},
    }
    prep_manifest = stage04 / "attempt-000/preprocessing.json"
    _write_json(prep_manifest, prep)
    _complete(stage04, prep, [transformed, prep_manifest])

    stage05 = pipeline / "05_clock_predictions"
    prediction_csv = stage05 / "attempt-000/clock_predictions.csv"
    prediction_csv.parent.mkdir(parents=True)
    prediction_csv.write_text(",Clock_tAge\nGSM1,1.1\nGSM2,1.2\n", encoding="utf-8")
    coverage_path = stage05 / "attempt-000/Clock-feature-coverage.json"
    coverage = {
        "feature_ids": ["f1", "f3", "f4"], "training_retained_ids": ["f1", "f3"],
        "training_dropped_ids": ["f4"], "selected_ids": ["f1", "f3"],
        "nonzero_coefficient_ids": ["f1"], "imputation_strategy": "median",
        "missing_ids_by_sample": {sample: ["f3", "f4"] for sample in samples},
        "imputed_ids_by_sample": {sample: ["f3"] for sample in samples},
        "observed_selected_by_sample": {sample: 1 for sample in samples},
        "observed_nonzero_by_sample": {sample: 1 for sample in samples},
    }
    _write_json(coverage_path, coverage)
    measurement_sources = {
        "raw_counts": raw_source,
        "transformed_input": {"path": str(transformed.resolve()), "sha256": _sha(transformed)},
        "model": asset_source,
        "feature_coverage": {"path": str(coverage_path.resolve()), "sha256": _sha(coverage_path)},
    }
    prediction = {
        "path": str(prediction_csv), "clocks": ["Clock_tAge"], "sample_count": 2,
        "checks": [{
            "clock": "Clock", "model_sha256": _sha(asset), "numerical_status": "PASS",
            "observed_features_min": 1, "required_features": 3,
            "max_abs_independent_difference": 0.0,
            "measurement_applicability": "CANNOT VERIFY",
            "measurement_sources": measurement_sources,
        }],
        "clock_predictions_recomputed": True, "scientific_status": "CANNOT VERIFY",
        "prediction_sha256": _sha(prediction_csv),
    }
    prediction_manifest = stage05 / "attempt-000/predictions.json"
    _write_json(prediction_manifest, prediction)
    _complete(stage05, prediction, [prediction_csv, coverage_path, prediction_manifest])
    return pipeline, {"study_id": "GSEX", "sample_ids": samples}, bindings


def _recovered_fixture(base):
    """Disclosed tiny software fixture, processed by the production recovery path."""
    from src.pipeline_inputs_2026_09_10 import propose_dnam_recovery, reconstruct_dnam_input
    from audit_system.report_2026_09_10 import verify_dnam_recovery

    study, ids = "GSE999001", ["GSM1", "GSM2"]
    root = base / "pipeline"
    metadata_dir = root / "01_metadata/attempt-000"
    metadata_dir.mkdir(parents=True)
    original = base / "empty.csv.gz"
    with gzip.open(original, "wt") as handle:
        handle.write("ID_REF,GSM1,GSM2\n")
    supplement = base / "Betas.tsv.gz"
    with gzip.open(supplement, "wt") as handle:
        handle.write("ID\tA\tB\ncg001\t0.1\t0.2\ncg002\t0.3\t0.4\n")
    metadata = metadata_dir / "metadata.csv"
    with metadata.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["GSM", "gse_id", "title"])
        writer.writeheader()
        writer.writerows({"GSM": gsm, "gse_id": study, "title": "Software fixture [" + token + "]"}
                         for gsm, token in zip(ids, ["A", "B"]))
    series = metadata_dir / "series.soft"
    series.write_text("^SERIES = " + study + "\n!Series_supplementary_file = "
                      "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE999nnn/GSE999001/suppl/Betas.tsv.gz\n"
                      "^SAMPLE = GSM1\n!Sample_title = Software fixture [A]\n"
                      "^SAMPLE = GSM2\n!Sample_title = Software fixture [B]\n")
    _complete(metadata_dir.parent, {"metadata": str(metadata), "series_source": str(series)},
              [metadata, series])
    proposal = propose_dnam_recovery(study, original, supplement, metadata, series)
    stage02 = root / "02_molecular_input/attempt-000"
    stage02.mkdir(parents=True)
    reconstruction = reconstruct_dnam_input(proposal, stage02)
    verification = verify_dnam_recovery(proposal, reconstruction)
    molecular = {**reconstruction, "repair_status": "PASS", "source": str(supplement),
                 "source_sha256": _sha(supplement), "repaired_original": str(original)}
    documents = {"recovery-proposal.json": proposal, "reconstruction.json": reconstruction,
                 "recovery-verification.json": verification, "molecular-input.json": molecular}
    for name, value in documents.items():
        _write_json(stage02 / name, value)
    _complete(stage02.parent, molecular, [Path(reconstruction["path"]),
                                         *(stage02 / name for name in documents)])
    bindings = {name: {"path": str(path), "sha256": _sha(path)}
                for name, path in (("matrix", original), ("dnam_supplement", supplement), ("series", series))}
    _, members, receipt = _stage(root, "02_molecular_input")
    return root, study, ids, molecular, bindings, members, {"02_molecular_input": receipt}


class MeasurementEvidenceSoftwareTests(unittest.TestCase):
    def test_recovered_source_chain_rechecks_numeric_values_and_primary_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            root, study, ids, molecular, bindings, members, receipts = _recovered_fixture(Path(directory))
            files = {}
            result = _recovered_dnam_source(root, study, molecular, bindings, ids, members, files, receipts)
            self.assertEqual(result["independent_input_verification"]["values_compared"], 4)
            self.assertEqual(result["independent_input_verification"]["scientific_status"], "CANNOT VERIFY")
            self.assertEqual([row["gsm_id"] for row in result["source_column_to_gsm"]], ids)
            self.assertIn("01_metadata", receipts)
            self.assertNotIn(str(root), json.dumps(result))

    def test_recovered_chain_rejects_wrong_source_roles_even_when_bytes_match(self):
        with tempfile.TemporaryDirectory() as directory:
            root, study, ids, molecular, bindings, members, receipts = _recovered_fixture(Path(directory))
            duplicate = Path(directory) / "unbound-copy.gz"
            shutil.copyfile(bindings["dnam_supplement"]["path"], duplicate)
            changed = copy.deepcopy(bindings)
            changed["dnam_supplement"]["path"] = str(duplicate)
            with self.assertRaisesRegex(ValueError, "source role differs"):
                _recovered_dnam_source(root, study, molecular, changed, ids, members, {}, receipts)
            with self.assertRaisesRegex(ValueError, "source sample membership"):
                _recovered_dnam_source(root, study, molecular, bindings, ["GSM1", "GSM3"], members, {}, receipts)

    def test_resealed_recovery_cannot_swap_columns_or_promote_scientific_acceptance(self):
        for defect in ("column_mapping", "scientific_status"):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as directory:
                root, study, ids, molecular, bindings, members, receipts = _recovered_fixture(Path(directory))
                folder = root / "02_molecular_input/attempt-000"
                if defect == "column_mapping":
                    path = folder / "recovery-proposal.json"
                    changed = json.loads(path.read_text())
                    changed["crosswalk"][0]["gsm_id"], changed["crosswalk"][1]["gsm_id"] = ids[1], ids[0]
                else:
                    path = folder / "recovery-verification.json"
                    changed = json.loads(path.read_text())
                    changed["scientific_status"] = "PASS"
                _write_json(path, changed)
                _complete(folder.parent, molecular, list(folder.iterdir()))
                _, members, receipt = _stage(root, "02_molecular_input")
                with self.assertRaisesRegex(ValueError, "no longer reproduces|does not reproduce"):
                    _recovered_dnam_source(root, study, molecular, bindings, ids, members, {},
                                           {"02_molecular_input": receipt})

    def test_compact_rna_bundle_and_mouse_route(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pipeline, index, bindings = _fixture(Path(directory))
            for modality in ("human_rna", "mouse_rna"):
                bundle = build_measurement_evidence(pipeline, "GSEX", modality, index, bindings)
                self.assertEqual(bundle["schema"], SCHEMA)
                self.assertEqual(list(bundle["evidence"]), ["M00000", "M00001", "M00002", "M00003"])
                coverage = bundle["evidence"]["M00003"]["value"]
                self.assertEqual(coverage["observed_features_by_sample"], {"GSM1": 1, "GSM2": 1})
                self.assertEqual(coverage["missing_features_by_sample"], {"GSM1": 2, "GSM2": 2})
                self.assertEqual(coverage["imputed_features_by_sample"], {"GSM1": 1, "GSM2": 1})
                self.assertEqual(coverage["missing_training_dropped_features_by_sample"], {"GSM1": 1, "GSM2": 1})
                self.assertEqual(coverage["training_dropped_feature_count"], 1)
                self.assertNotIn(str(Path(directory)), json.dumps(bundle["evidence"]))

    def test_archs4_extraction_binding(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pipeline, index, bindings = _fixture(Path(directory), archs4=True)
            bundle = build_measurement_evidence(pipeline, "GSEX", "human_rna", index, bindings)
            self.assertIn("stage02:extraction_manifest", bundle["source_files"])
            self.assertEqual(bundle["evidence"]["M00000"]["value"]["sample_count"], 2)

    def test_accepts_metadata_reorder_and_rejects_missing_identity_or_tamper(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pipeline, index, bindings = _fixture(Path(directory))
            reordered = dict(index, sample_ids=list(reversed(index["sample_ids"])))
            bundle = build_measurement_evidence(pipeline, "GSEX", "human_rna", reordered, bindings)
            self.assertFalse(bundle["evidence"]["M00000"]["value"]["raw_order_matches_source_index"])
            bad_index = dict(index, sample_ids=["GSM1", "GSM3"])
            with self.assertRaisesRegex(ValueError, "sample identity membership mismatch"):
                build_measurement_evidence(pipeline, "GSEX", "human_rna", bad_index, bindings)
            member = pipeline / "05_clock_predictions/attempt-000/clock_predictions.csv"
            member.write_text(member.read_text() + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
                build_measurement_evidence(pipeline, "GSEX", "human_rna", index, bindings)

    def test_rejects_duplicate_columns_and_unsealed_asset(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pipeline, index, bindings = _fixture(Path(directory))
            prediction_csv = pipeline / "05_clock_predictions/attempt-000/clock_predictions.csv"
            prediction_csv.write_text(",Clock_tAge,Clock_tAge\nGSM1,1,1\nGSM2,2,2\n", encoding="utf-8")
            pred_manifest = pipeline / "05_clock_predictions/attempt-000/predictions.json"
            result = json.loads(pred_manifest.read_text())
            result["prediction_sha256"] = _sha(prediction_csv)
            _write_json(pred_manifest, result)
            _complete(pipeline / "05_clock_predictions", result, [
                prediction_csv,
                pipeline / "05_clock_predictions/attempt-000/Clock-feature-coverage.json",
                pred_manifest,
            ])
            with self.assertRaisesRegex(ValueError, "duplicate or absent literal data columns"):
                build_measurement_evidence(pipeline, "GSEX", "human_rna", index, bindings)

            pipeline, index, bindings = _fixture(Path(directory) / "second")
            del bindings["fixture_asset"]
            with self.assertRaisesRegex(ValueError, "external measurement source lacks a sealed binding"):
                build_measurement_evidence(pipeline, "GSEX", "human_rna", index, bindings)


@unittest.skipUnless(os.environ.get("CLOCKBASE_EMPIRICAL_MEASUREMENT_READBACK") == "1",
                     "set CLOCKBASE_EMPIRICAL_MEASUREMENT_READBACK=1 for saved-run readbacks")
class MeasurementEvidenceSavedRunTests(unittest.TestCase):
    def test_saved_rna_and_dnam(self) -> None:
        cases = (("registry-006", "human_rna", "GSE127814", 4),
                 ("registry-007", "human_dnam", "GSE247655", 6))
        for registry, modality, study, expected in cases:
            paths = list((UPDATE / registry / "studies" / modality / study).glob("*/version.json"))
            completed = [p for p in paths if (p.parent / "pipeline/05_clock_predictions/complete.json").exists()]
            audited = [p for p in completed if (p.parent / "comprehensive-audit").exists()]
            version_path = (audited or completed)[0]
            version = json.loads(version_path.read_text())
            index_path = version_path.parent / "comprehensive-audit/study/01_source_index/attempt-000/source-index.json"
            if index_path.exists():
                index = json.loads(index_path.read_text())
            else:
                complete = json.loads((version_path.parent / "pipeline/05_clock_predictions/complete.json").read_text())
                frame = pd.read_csv(complete["result"]["path"], index_col=0)
                index = {"study_id": study, "sample_ids": frame.index.astype(str).tolist()}
            bundle = build_measurement_evidence(
                version_path.parent / "pipeline", study, modality, index, version["sources"]
            )
            self.assertEqual(len(bundle["evidence"]), expected)


class DnamIdatMeasurementInterfaceTests(unittest.TestCase):
    """Small software receipts exercise lineage; no raw assay or clocks are run."""

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.base = Path(self.directory.name)
        self.root = self.base / 'pipeline'
        self.ids = ['GSM1', 'GSM2']
        self.index = {'study_id': 'GSE1', 'sample_ids': self.ids}
        series = self.base / 'series.soft'
        series.write_text('^SERIES = GSE1\n!Series_sample_id = GSM1\n!Series_sample_id = GSM2\n')
        samples = self.base / 'samples.soft'
        samples.write_text('^SAMPLE = GSM1\n^SAMPLE = GSM2\n')
        matrix = self.base / 'beta.csv'
        matrix.write_text('CpG,GSM1,GSM2\ncg001,0.100,0.200\ncg002,,0.400\n')
        channel = self.base / 'channel.idat'
        channel.write_bytes(b'SOFTWARE FIXTURE: no decoded IDAT data')
        self.wrapper = self.base / 'idat-provenance.json'
        _write_json(self.wrapper, {'software_fixture': True})
        declaration = lambda path: {'path': str(path), 'sha256': _sha(path)}
        self.record = declaration(self.wrapper)
        files = {role: declaration(path) for role, path in
                 [('matrix', matrix), ('series', series), ('samples', samples), ('raw:GSM1:Grn', channel)]}
        self.verified = {'study_id': 'GSE1', 'sample_ids': self.ids, 'files': files,
                         'measurement': {'profile': 'SOFTWARE_FIXTURE', 'qc_policy': 'explicit_fixture_mask',
                                         'sample_count': 2, 'scientific_status': 'CANNOT VERIFY'},
                         'scientific_status': 'CANNOT VERIFY'}
        self.bindings = {role: dict(files[role]) for role in ('matrix', 'series', 'samples')}
        self.bindings['dnam_idat_provenance'] = dict(self.record)
        self.bindings.update({'dnam_idat_provenance:' + role: dict(value) for role, value in files.items()})
        reference = self.base / 'reference.csv'
        reference.write_text('CpG,reference\ncg001,0.15\ncg002,0.35\n')
        self.bindings['dnam_reference'] = declaration(reference)

        input_stage = self.root / '02_molecular_input/attempt-000'
        input_stage.mkdir(parents=True)
        self.raw = input_stage / 'beta.csv'
        self.raw.write_text('CpG,GSM1,GSM2\ncg001,0.1,0.2\ncg002,,0.4\n')
        self.molecular = {'path': str(self.raw), 'sha256': _sha(self.raw), 'source': str(matrix),
                          'source_sha256': _sha(matrix), 'feature_count': 2, 'sample_count': 2,
                          'idat_provenance': dict(self.record)}
        self.seal_input()
        prep_stage = self.root / '04_preprocessing/attempt-000'
        prep_stage.mkdir(parents=True)
        imputed = prep_stage / 'imputed_beta.csv'
        imputed.write_text('CpG,GSM1,GSM2\ncg001,0.1,0.2\ncg002,0.35,0.4\n')
        self.prep = {'path': str(imputed), 'raw': str(self.raw), 'status': 'PASS',
                     'reference': {'source': str(reference), 'source_sha256': _sha(reference),
                                   'reference_study': 'SOFTWARE_FIXTURE', 'means_sha256': '0' * 64},
                     'raw_methylation_preprocessing_recomputed': True,
                     'independent_imputation_verification': 'PASS', 'idat_processing': True,
                     'idat_provenance': dict(self.record)}
        self.seal_prep()
        clock_stage = self.root / '05_clock_predictions/attempt-000'
        clock_stage.mkdir(parents=True)
        predictions = clock_stage / 'clock_predictions.csv'
        predictions.write_text('sample_id,Hannum\nGSM1,0.8\nGSM2,1.0\n')
        coefficients = clock_stage / 'Hannum-coefficients.csv'
        coefficients.write_text('CpGmarker,CoefficientTraining\ncg001,1\ncg002,2\nintercept,0\n')
        result = {'path': str(predictions), 'clocks': ['Hannum'], 'sample_count': 2,
                  'clock_predictions_recomputed': True, 'scientific_status': 'CANNOT VERIFY',
                  'checks': [{'clock': 'Hannum', 'numerical_status': 'PASS', 'required_features': 2,
                              'observed_features_min': 1, 'max_abs_independent_difference': 0,
                              'measurement_applicability': 'CANNOT VERIFY',
                              'measurement_sources': {'raw_beta': declaration(self.raw),
                                                      'coefficients': declaration(coefficients)}}]}
        manifest = clock_stage / 'predictions.json'
        _write_json(manifest, result)
        _complete(clock_stage.parent, result, [predictions, coefficients, manifest])
        verifier = patch('audit_system.report_2026_09_10.verify_dnam_idat_provenance',
                         return_value=self.verified, create=True)
        self.verifier = verifier.start()
        self.addCleanup(verifier.stop)

    def seal_input(self):
        stage = self.raw.parent
        _write_json(stage / 'molecular-input.json', self.molecular)
        _complete(stage.parent, self.molecular, [self.raw, stage / 'molecular-input.json'])

    def seal_prep(self):
        stage = Path(self.prep['path']).parent
        _write_json(stage / 'preprocessing.json', self.prep)
        _complete(stage.parent, self.prep, [Path(self.prep['path']), stage / 'preprocessing.json'])

    def build(self):
        return build_measurement_evidence(self.root, 'GSE1', 'human_dnam', self.index, self.bindings)

    def test_idat_lineage_survives_reserialization_and_remains_separate_from_coverage(self):
        self.assertNotEqual(_sha(self.raw), self.bindings['matrix']['sha256'])
        bundle = self.build()
        records = list(bundle['evidence'].values())
        lineage = next(row for row in records if row['kind'] == 'source_bound_dnam_idat_preprocessing')
        self.assertEqual(lineage['value']['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(lineage['sample_ids'], self.ids)
        coverage = next(row for row in records if row['kind'] == 'clock_feature_coverage')['value']
        self.assertEqual(coverage['observed_features_by_sample'], {'GSM1': 1, 'GSM2': 2})
        self.assertEqual(coverage['imputed_features_by_sample'], {'GSM1': 1, 'GSM2': 0})
        self.assertNotIn(str(self.base), json.dumps(bundle['evidence']))

    def test_idat_missing_or_extra_nested_role_is_rejected(self):
        for defect in ('missing', 'extra'):
            with self.subTest(defect=defect):
                saved = copy.deepcopy(self.bindings)
                if defect == 'missing':
                    self.bindings.pop('dnam_idat_provenance:raw:GSM1:Grn')
                else:
                    self.bindings['dnam_idat_provenance:undeclared'] = dict(self.record)
                with self.assertRaisesRegex(ValueError, 'nested sources differ'):
                    self.build()
                self.bindings = saved

    def test_idat_cannot_use_same_bytes_at_an_unbound_nested_path(self):
        role = 'dnam_idat_provenance:raw:GSM1:Grn'
        duplicate = self.base / 'unbound-channel-copy.idat'
        duplicate.write_bytes(Path(self.bindings[role]['path']).read_bytes())
        self.bindings[role]['path'] = str(duplicate)
        with self.assertRaisesRegex(ValueError, 'exact sealed version binding'):
            self.build()

    def test_idat_stage02_values_cannot_be_resealed_after_change(self):
        self.raw.write_text('CpG,GSM1,GSM2\ncg001,0.9,0.2\ncg002,,0.4\n')
        self.molecular['sha256'] = _sha(self.raw)
        self.seal_input()
        with self.assertRaisesRegex(ValueError, 'values or identities differ'):
            self.build()

    def test_idat_stage04_requires_matching_lineage_and_flag(self):
        for defect in ('missing', 'changed', 'false'):
            with self.subTest(defect=defect):
                saved = copy.deepcopy(self.prep)
                if defect == 'missing':
                    self.prep.pop('idat_provenance')
                elif defect == 'changed':
                    self.prep['idat_provenance']['sha256'] = 'a' * 64
                else:
                    self.prep['idat_processing'] = False
                self.seal_prep()
                with self.assertRaisesRegex(ValueError, 'lost or changed'):
                    self.build()
                self.prep = saved

    def test_idat_claim_without_provenance_role_is_rejected(self):
        self.bindings = {role: value for role, value in self.bindings.items()
                         if not role.startswith('dnam_idat_provenance')}
        with self.assertRaisesRegex(ValueError, 'lacks a sealed provenance role'):
            self.build()

    def test_idat_compact_measurement_summary_cannot_leak_private_paths(self):
        self.verified['measurement']['unbound_private_path'] = '/Entropy/private/channel.idat'
        with self.assertRaisesRegex(ValueError, 'absolute path leaked'):
            self.build()


if __name__ == "__main__":
    unittest.main()
