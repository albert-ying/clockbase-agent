"""Software fixtures for the extraction-only ARCHS4 batch coordinator.

The tiny matrices are deterministic software fixtures, not empirical results.
"""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from audit_system import archs4_update_2026_09_10 as intake


class Archs4BatchExtractionTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        self.release = self.root / "human-release.h5"
        with h5py.File(self.release, "w") as handle:
            handle.create_dataset("data/expression", data=np.array([[1, 2, 3, 4], [5, 6, 7, 8]], dtype=np.uint32))
            handle.create_dataset("meta/genes/ensembl_gene", data=np.array([b"ENSG1", b"ENSG2"]))
            handle.create_dataset("meta/samples/geo_accession", data=np.array([b"GSM1", b"GSM2", b"GSM3", b"GSM4"]))
            handle.create_dataset("meta/samples/series_id", data=np.array([b"GSE1", b"GSE1", b"GSE2", b"GSE3"]))
            handle.create_dataset("meta/samples/organism_ch1", data=np.array([b"Homo sapiens"] * 4))
        self.sources = {}
        self.rows = []
        for study, samples, molecules, supported in (
                ("GSE1", ["GSM1", "GSM2"], ["total RNA", "total RNA"], True),
                ("GSE2", ["GSM3"], ["total RNA"], True),
                ("GSE3", ["GSM4"], ["polyA RNA"], False)):
            series = self.root / f"{study}-series.soft"
            series.write_text(
                f"^SERIES = {study}\n!Series_title = Software fixture {study}\n"
                "!Series_status = Public on Jan 01 2025\n"
                + "".join(f"!Series_sample_id = {sample}\n" for sample in samples))
            sample_source = self.root / f"{study}-samples.soft"
            sample_source.write_text("".join(
                f"^SAMPLE = {sample}\n!Sample_title = Software fixture {sample}\n"
                "!Sample_source_name_ch1 = software fixture\n"
                "!Sample_organism_ch1 = Homo sapiens\n"
                f"!Sample_molecule_ch1 = {molecule}\n!Sample_platform_id = GPL1\n"
                for sample, molecule in zip(samples, molecules)))
            declarations = {
                "series": self.declaration(series), "samples": self.declaration(sample_source)}
            self.sources[study] = declarations
            counts = {}
            for molecule in molecules:
                counts[molecule] = counts.get(molecule, 0) + 1
            self.rows.append({
                "study_id": study, "modality": "human_rna", "sample_ids": samples,
                "sample_count": len(samples), "first_public_date": "2025-01-01",
                "sources": declarations,
                "samples": [{"sample_id": sample, "organism_ch1": ["Homo sapiens"],
                             "molecule_ch1": [molecule]}
                            for sample, molecule in zip(samples, molecules)],
                "molecule_label_counts": counts,
                "profile": "existing_total_rna_intake_v21",
                "profile_status": "SUPPORTED_SOURCE_PROFILE" if supported else "UNSUPPORTED_SOURCE_PROFILE",
                "limitations": [] if supported else ["Molecule profile is outside the existing total-RNA contract"],
                "source_status_quote": "Public on Jan 01 2025",
                "scientific_status": "CANNOT VERIFY"})
        self.discovery = self.root / "discovery.json"
        self.discovery.write_text('{"software_fixture":true}\n')
        self.candidates = self.root / "candidates.json"
        self.candidates.write_text('{"software_fixture":true}\n')
        source_files = {str(path): self.sha(path) for path in
                        [self.discovery, self.candidates,
                         *[Path(value[role]["path"]) for value in self.sources.values()
                           for role in ("series", "samples")]]}
        self.descriptor = self.root / "source-profiles.json"
        self.descriptor.write_text(json.dumps({
            "schema": intake.BATCH_SOURCE_PROFILE_SCHEMA,
            "status": "PASS", "profile": "existing_total_rna_intake_v21",
            "modality": "human_rna", "discovery": self.declaration(self.discovery),
            "candidates": self.declaration(self.candidates), "studies": self.rows,
            "supported_study_ids": ["GSE1", "GSE2"], "stopped_study_ids": ["GSE3"],
            "counts": {"enrolled_studies": 3, "enrolled_samples": 4,
                       "supported_studies": 2, "stopped_studies": 1,
                       "molecule_label_occurrences": {"total RNA": 3, "polyA RNA": 1}},
            "source_files": source_files, "verification_script_sha256": "0" * 64,
            "model_calls": 0, "expression_cells_read": 0,
            "scientific_status": "CANNOT VERIFY",
            "original_paper_2025_eligibility": "QUARANTINED"}, sort_keys=True) + "\n")
        self.output = self.root / "batch"
        self.plan_path = self.root / "plan.json"
        self.write_plan(self.plan_path, self.output)

    @staticmethod
    def sha(path):
        value = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                value.update(block)
        return value.hexdigest()

    def declaration(self, path):
        path = Path(path).resolve()
        return {"path": str(path), "sha256": self.sha(path), "size_bytes": path.stat().st_size}

    def plan(self, output):
        implementation = {role: self.declaration(path) for role, path in intake.BATCH_IMPLEMENTATIONS.items()}
        return {
            "schema": intake.BATCH_EXTRACTION_SCHEMA, "output": str(Path(output).resolve()),
            "runtime": {"python": platform.python_version(), "packages": {
                name: importlib.metadata.version(name) for name in ("h5py", "numpy", "pandas")}},
            "implementation": implementation,
            "source_profile_descriptor": self.declaration(self.descriptor),
            "release": {**self.declaration(self.release), "species": "human",
                        "expression_shape": [2, 4], "sample_count": 4},
            "selected_study_ids": ["GSE1", "GSE2"],
            "limits": {"max_studies": 3, "max_samples_per_study": 2,
                       "max_extracted_bytes_per_study": 100000,
                       "max_total_extracted_bytes": 300000, "concurrency": 1}}

    def write_plan(self, path, output, mutate=None):
        plan = self.plan(output)
        if mutate:
            mutate(plan)
        path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
        return plan

    def run_basic(self, output=None, plan=None):
        output = Path(output or self.output)
        plan = Path(plan or self.plan_path)
        with patch("src.pipeline_inputs_2026_09_10.preprocess_rna",
                   side_effect=AssertionError("preprocessing must not run")), \
             patch("src.pipeline_inputs_2026_09_10.predict_rna",
                   side_effect=AssertionError("clocks must not run")):
            return intake.run_archs4_batch_extraction(plan, output=output)

    def mouse_fixture(self):
        """Change only the disclosed tiny software sources to mouse RNA."""
        mouse_release = self.root / 'mouse-release.h5'
        self.release.rename(mouse_release)
        self.release = mouse_release
        with h5py.File(self.release, 'r+') as handle:
            handle['meta/samples/organism_ch1'][:] = np.array([b'Mus musculus'] * 4)
        for row in self.rows:
            row['modality'] = 'mouse_rna'
            for sample in row['samples']:
                sample['organism_ch1'] = ['Mus musculus']
            source = Path(row['sources']['samples']['path'])
            source.write_text(source.read_text().replace('Homo sapiens', 'Mus musculus'))
            row['sources']['samples'].update(sha256=self.sha(source), size_bytes=source.stat().st_size)
        descriptor = json.loads(self.descriptor.read_text())
        descriptor['modality'] = 'mouse_rna'
        descriptor['studies'] = self.rows
        for row in self.rows:
            source = Path(row['sources']['samples']['path'])
            descriptor['source_files'][str(source)] = self.sha(source)
        self.descriptor.write_text(json.dumps(descriptor, sort_keys=True) + '\n')
        self.write_plan(self.plan_path, self.output,
                        lambda plan: plan.update(modality='mouse_rna', release={
                            **self.declaration(self.release), 'species': 'mouse',
                            'expression_shape': [2, 4], 'sample_count': 4}))

    def test_mouse_batch_binds_species_and_preserves_unsupported_profile(self):
        self.mouse_fixture()
        summary = self.run_basic()
        self.assertEqual((summary['modality'], summary['source_hdf5']['species']), ('mouse_rna', 'mouse'))
        self.assertEqual([row['modality'] for row in summary['results']], ['mouse_rna'] * 3)
        self.assertEqual(summary['results'][2]['status'], 'SOURCE_PROFILE_UNRESOLVED')
        extraction = json.loads((self.output / summary['results'][0]['extraction_receipt']['path']).read_text())
        self.assertEqual((extraction['modality'], extraction['species']), ('mouse_rna', 'mouse'))
        stop = json.loads((self.output / summary['results'][2]['stop_receipt']['path']).read_text())
        self.assertEqual((stop['modality'], stop['species']), ('mouse_rna', 'mouse'))
        self.assertEqual((summary['model_api_calls'], summary['preprocessing_runs'], summary['clock_runs']), (0, 0, 0))

    def test_mouse_profile_with_human_release_species_fails_before_extraction(self):
        self.mouse_fixture()
        plan = json.loads(self.plan_path.read_text())
        plan['release']['species'] = 'human'  # Same valid file hash and GSMs.
        self.plan_path.write_text(json.dumps(plan, sort_keys=True) + '\n')
        with self.assertRaisesRegex(ValueError, 'modality differs'):
            self.run_basic()
        self.assertFalse((self.output / 'batch-summary.json').exists())

    def test_mouse_release_human_organism_stops_only_affected_study(self):
        self.mouse_fixture()
        with h5py.File(self.release, 'r+') as handle:
            handle['meta/samples/organism_ch1'][2] = b'Homo sapiens'
        plan = json.loads(self.plan_path.read_text())
        plan['release']['sha256'] = self.sha(self.release)
        self.plan_path.write_text(json.dumps(plan, sort_keys=True) + '\n')
        summary = self.run_basic()
        by_study = {row['study_id']: row for row in summary['results']}
        self.assertEqual(by_study['GSE1']['status'], 'EXTRACTED')
        self.assertEqual(by_study['GSE2']['status'], 'EXTRACTION_UNRESOLVED')
        self.assertIsNone(by_study['GSE2']['matrix'])

    def test_extracts_supported_partition_and_retains_profile_stop(self):
        original_digest = intake.digest
        release_hash_calls = []

        def counted(path):
            if Path(path).resolve() == self.release:
                release_hash_calls.append(1)
            return original_digest(Path(path))

        original_h5 = h5py.File
        with patch.object(intake, "digest", side_effect=counted), \
             patch("h5py.File", wraps=original_h5) as opened:
            summary = self.run_basic()
        self.assertEqual(len(release_hash_calls), 2)
        self.assertEqual(opened.call_count, 1)
        self.assertEqual(summary["counts"]["enrolled_studies"], 3)
        self.assertEqual(summary["counts"]["extracted_studies"], 2)
        self.assertEqual(summary["counts"]["source_profile_stops"], 1)
        self.assertEqual(summary["counts"]["extraction_unresolved_studies"], 0)
        self.assertEqual([row["study_id"] for row in summary["results"]], ["GSE1", "GSE2", "GSE3"])
        self.assertEqual(summary["results"][2]["status"], "SOURCE_PROFILE_UNRESOLVED")
        self.assertEqual(summary["model_api_calls"], 0)
        self.assertEqual(summary["preprocessing_runs"], 0)
        self.assertEqual(summary["clock_runs"], 0)
        self.assertEqual(summary["path_base"], "batch_root")
        before = json.loads((self.output / summary["release_before"]["path"]).read_text())
        after = json.loads((self.output / summary["release_after"]["path"]).read_text())
        self.assertEqual((before["hash_scope"], after["hash_scope"]), ("FULL_HDF5", "FULL_HDF5"))
        self.assertEqual((before["size_bytes"], after["size_bytes"]),
                         (self.release.stat().st_size, self.release.stat().st_size))

        gse1 = summary["results"][0]
        observed = (self.output / gse1["matrix"]["path"]).read_text()
        self.assertEqual(observed, "ensembl_gene,GSM1,GSM2\nENSG1,1,2\nENSG2,5,6\n")
        extraction = json.loads((self.output / gse1["extraction_receipt"]["path"]).read_text())
        self.assertEqual(extraction["schema"], intake.BATCH_STUDY_EXTRACTION_SCHEMA)
        self.assertEqual(extraction["path_base"], "batch_root")
        self.assertEqual(extraction["source"]["source_profile_record_sha256"],
                         intake._intake_json_hash(self.rows[0]))

    def test_batch_matrix_bytes_match_existing_extract_rna(self):
        summary = self.run_basic()
        batch_path = self.output / summary["results"][0]["matrix"]["path"]
        direct = self.root / "direct"
        direct.mkdir()
        from src.pipeline_inputs_2026_09_10 import extract_rna
        direct_result = extract_rna(self.release, "GSE1", "human", ["GSM2", "GSM1"], direct)
        self.assertEqual(batch_path.read_bytes(), Path(direct_result["path"]).read_bytes())

    def test_completed_resume_is_read_only_and_does_not_extract(self):
        first = self.run_basic()
        before = {str(path.relative_to(self.output)): (self.sha(path), path.stat().st_mtime_ns)
                  for path in self.output.rglob("*") if path.is_file()}
        original_digest = intake.digest
        release_hash_calls = []

        def counted(path):
            if Path(path).resolve() == self.release:
                release_hash_calls.append(1)
            return original_digest(Path(path))

        with patch.object(intake, "digest", side_effect=counted), \
             patch.object(intake, "_batch_extract_from_index",
                          side_effect=AssertionError("completed resume must not extract")):
            second = self.run_basic()
        after = {str(path.relative_to(self.output)): (self.sha(path), path.stat().st_mtime_ns)
                 for path in self.output.rglob("*") if path.is_file()}
        self.assertEqual(first, second)
        self.assertEqual(before, after)
        self.assertEqual(len(release_hash_calls), 1)

    def test_completed_resume_rejects_mutated_output(self):
        summary = self.run_basic()
        matrix = self.output / summary["results"][0]["matrix"]["path"]
        matrix.write_text(matrix.read_text().replace("ENSG1,1,2", "ENSG1,9,2"))
        with self.assertRaisesRegex(ValueError, "Completed batch artifact changed"):
            self.run_basic()

    def test_plan_mutation_during_execution_is_preserved_as_failure(self):
        original = intake._batch_extract_from_index
        mutated = False

        def mutate_plan(*args, **kwargs):
            nonlocal mutated
            value = original(*args, **kwargs)
            if not mutated:
                self.plan_path.write_bytes(self.plan_path.read_bytes() + b"\n")
                mutated = True
            return value

        with patch.object(intake, "_batch_extract_from_index", side_effect=mutate_plan):
            with self.assertRaisesRegex(ValueError, "plan changed during execution"):
                self.run_basic()
        self.assertFalse((self.output / "batch-summary.json").exists())
        attempt_path = next((self.output / "execution-attempts").glob("attempt-*/attempt.json"))
        self.assertEqual(json.loads(attempt_path.read_text())["status"], "FAIL")

    def test_selected_partition_must_match_descriptor(self):
        path = self.root / "bad-plan.json"
        self.write_plan(path, self.root / "bad-output",
                        lambda plan: plan.update(selected_study_ids=["GSE1", "GSE3"]))
        with self.assertRaisesRegex(ValueError, "exactly preserve"):
            intake.run_archs4_batch_extraction(path)

    def test_membership_mismatch_is_retained_per_study(self):
        with h5py.File(self.release, "r+") as handle:
            handle["meta/samples/series_id"][2] = b"GSE9"
        path = self.root / "membership-plan.json"

        def update_release(plan):
            plan["release"]["sha256"] = self.sha(self.release)
        self.write_plan(path, self.root / "membership-output", update_release)
        summary = self.run_basic(self.root / "membership-output", path)
        by_study = {row["study_id"]: row for row in summary["results"]}
        self.assertEqual(by_study["GSE1"]["status"], "EXTRACTED")
        self.assertEqual(by_study["GSE2"]["status"], "EXTRACTION_UNRESOLVED")
        self.assertIsNone(by_study["GSE2"]["matrix"])

    def test_per_study_safe_bound_stops_without_partial_matrix(self):
        path = self.root / "cap-plan.json"
        self.write_plan(path, self.root / "cap-output",
                        lambda plan: plan["limits"].update(max_extracted_bytes_per_study=1000))
        summary = self.run_basic(self.root / "cap-output", path)
        self.assertEqual(summary["counts"]["extracted_studies"], 0)
        self.assertEqual(summary["counts"]["extraction_unresolved_studies"], 2)
        for row in summary["results"][:2]:
            self.assertIsNone(row["matrix"])

    def test_release_post_hash_mismatch_prevents_summary(self):
        output = self.root / "mutation-output"
        original_digest = intake.digest
        release_calls = 0

        def mutate_hash(path):
            nonlocal release_calls
            if Path(path).resolve() == self.release:
                release_calls += 1
                if release_calls == 2:
                    return "f" * 64
            return original_digest(Path(path))

        with patch.object(intake, "digest", side_effect=mutate_hash):
            with self.assertRaisesRegex(ValueError, "changed between full-file"):
                self.run_basic(output)
        self.assertEqual(release_calls, 2)
        self.assertFalse((output / "batch-summary.json").exists())
        attempt = json.loads(next((output / "execution-attempts").glob("attempt-*/attempt.json")).read_text())
        self.assertEqual(attempt["status"], "INTEGRITY_FAILURE")

    def test_release_pre_hash_mismatch_is_preserved_without_extraction(self):
        path = self.root / "wrong-release-plan.json"
        self.write_plan(path, self.root / "wrong-release-output",
                        lambda plan: plan["release"].update(sha256="f" * 64))
        with patch.object(intake, "_batch_extract_from_index",
                          side_effect=AssertionError("precheck failure must not extract")):
            with self.assertRaisesRegex(ValueError, "differs before"):
                intake.run_archs4_batch_extraction(path)
        attempt_path = next((self.root / "wrong-release-output" / "execution-attempts").glob("attempt-*/attempt.json"))
        attempt = json.loads(attempt_path.read_text())
        self.assertEqual(attempt["status"], "INTEGRITY_FAILURE")
        self.assertIsNone(attempt["release_after"])

    def test_fixed_identifier_bound_is_enforced_before_matrix_output(self):
        with h5py.File(self.release, "r+") as handle:
            del handle["meta/genes/ensembl_gene"]
            handle.create_dataset("meta/genes/ensembl_gene",
                                  data=np.array([b"E" * 129, b"ENSG2"], dtype="S129"))
        path = self.root / "long-id-plan.json"
        self.write_plan(path, self.root / "long-id-output",
                        lambda plan: plan["release"].update(sha256=self.sha(self.release)))
        with self.assertRaisesRegex(ValueError, "identifier exceeds 128"):
            intake.run_archs4_batch_extraction(path)
        self.assertFalse((self.root / "long-id-output" / "batch-summary.json").exists())


if __name__ == "__main__":
    unittest.main()
