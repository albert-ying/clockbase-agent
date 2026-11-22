import tempfile
import unittest
from pathlib import Path

from audit_system.dependency_graph import (
    DependencyCycleError,
    DependencyGraph,
    MissingReferenceError,
)
from audit_system.ledger import Ledger
from audit_system.models import (
    ActorRecord,
    ActorRole,
    ActorType,
    ArtifactRecord,
    ClaimRecord,
    DependencyRecord,
    DependencyRequirement,
    EvidenceClassification,
)


class DependencyGraphTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.ledger = Ledger(
            Path(self.temporary_directory.name) / "records.jsonl"
        )
        self.release_id = "release:test"
        self.created_at = "2026-09-07T00:00:00Z"
        self.ledger.append(
            ActorRecord(
                record_id="actor:custodian",
                release_id=self.release_id,
                created_at=self.created_at,
                actor_type=ActorType.PROGRAM,
                role=ActorRole.EVIDENCE_CUSTODIAN,
                organization="ClockBase audit",
                model_family=None,
                authority_scope=("artifact_registration",),
            )
        )

    def tearDown(self):
        self.temporary_directory.cleanup()

    def add_artifact(self, record_id, classification=EvidenceClassification.DERIVED):
        self.ledger.append(
            ArtifactRecord(
                record_id=record_id,
                release_id=self.release_id,
                created_at=self.created_at,
                uri=f"fixture://{record_id}",
                sha256="0" * 64,
                byte_count=1,
                classification=classification,
                created_by_actor_id="actor:custodian",
            )
        )

    def add_claim(self, record_id="claim:ouabain"):
        self.ledger.append(
            ClaimRecord(
                record_id=record_id,
                release_id=self.release_id,
                created_at=self.created_at,
                exact_text="Ouabain was selected by the reported ranking.",
                claim_type="selection",
                scope="submitted manuscript",
                required_dependency_ids=(),
                falsifier="The ranking cannot reproduce the selection.",
            )
        )

    def add_dependency(self, number, parent, child, requirement=None, release_id=None):
        self.ledger.append(
            DependencyRecord(
                record_id=f"dependency:{number}",
                release_id=release_id or self.release_id,
                created_at=self.created_at,
                parent_id=parent,
                child_id=child,
                requirement=requirement
                or DependencyRequirement.REQUIRED_DERIVED,
            )
        )

    def test_quarantine_reaches_ranking_figure_and_claim(self):
        for record_id in (
            "artifact:source",
            "artifact:effect",
            "artifact:ranking",
            "artifact:figure",
        ):
            self.add_artifact(record_id)
        self.add_claim()
        self.add_dependency(1, "artifact:source", "artifact:effect")
        self.add_dependency(2, "artifact:effect", "artifact:ranking")
        self.add_dependency(3, "artifact:ranking", "artifact:figure")
        self.add_dependency(4, "artifact:figure", "claim:ouabain")

        graph = DependencyGraph.from_snapshot(
            self.ledger.snapshot(),
            self.release_id,
        )

        self.assertEqual(
            graph.quarantine_impact({"artifact:source"}),
            (
                "artifact:effect",
                "artifact:figure",
                "artifact:ranking",
                "claim:ouabain",
            ),
        )
        self.assertEqual(
            graph.ancestors("claim:ouabain"),
            (
                "artifact:effect",
                "artifact:figure",
                "artifact:ranking",
                "artifact:source",
            ),
        )

    def test_supporting_edge_is_not_a_required_ancestor(self):
        self.add_artifact("artifact:required")
        self.add_artifact("artifact:supporting")
        self.add_claim()
        self.add_dependency(
            1,
            "artifact:required",
            "claim:ouabain",
            DependencyRequirement.REQUIRED_DERIVED,
        )
        self.add_dependency(
            2,
            "artifact:supporting",
            "claim:ouabain",
            DependencyRequirement.SUPPORTING,
        )

        graph = DependencyGraph.from_snapshot(
            self.ledger.snapshot(),
            self.release_id,
        )

        self.assertEqual(
            tuple(item.record_id for item in graph.required_ancestors("claim:ouabain")),
            ("dependency:1",),
        )

    def test_missing_endpoint_is_rejected(self):
        self.add_claim()
        self.add_dependency(1, "artifact:missing", "claim:ouabain")

        with self.assertRaisesRegex(MissingReferenceError, "artifact:missing"):
            DependencyGraph.from_snapshot(
                self.ledger.snapshot(),
                self.release_id,
            )

    def test_cycle_is_rejected(self):
        self.add_artifact("artifact:a")
        self.add_artifact("artifact:b")
        self.add_dependency(1, "artifact:a", "artifact:b")
        self.add_dependency(2, "artifact:b", "artifact:a")

        with self.assertRaises(DependencyCycleError):
            DependencyGraph.from_snapshot(
                self.ledger.snapshot(),
                self.release_id,
            )

    def test_other_release_dependencies_are_excluded(self):
        self.add_artifact("artifact:a")
        self.add_artifact("artifact:b")
        self.add_dependency(
            1,
            "artifact:a",
            "artifact:b",
            release_id="release:other",
        )

        graph = DependencyGraph.from_snapshot(
            self.ledger.snapshot(),
            self.release_id,
        )

        self.assertEqual(graph.descendants("artifact:a"), ())


if __name__ == "__main__":
    unittest.main()
