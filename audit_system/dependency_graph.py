"""Release-scoped dependency traversal and quarantine propagation."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Collection

from .ledger import LedgerSnapshot
from .models import DependencyRecord, DependencyRequirement


class DependencyGraphError(RuntimeError):
    """Base class for invalid dependency graphs."""


class DependencyCycleError(DependencyGraphError):
    """Raised when a release dependency graph contains a cycle."""


class MissingReferenceError(DependencyGraphError):
    """Raised when an edge names an unavailable release record."""


_REQUIRED = {
    DependencyRequirement.REQUIRED_EMPIRICAL,
    DependencyRequirement.REQUIRED_DERIVED,
    DependencyRequirement.REQUIRED_CHECK,
}


class DependencyGraph:
    """A validated directed acyclic graph for one release."""

    def __init__(
        self,
        release_id: str,
        records: dict[str, DependencyRecord],
    ):
        self.release_id = release_id
        self.records = records
        self._forward: dict[str, list[DependencyRecord]] = defaultdict(list)
        self._reverse: dict[str, list[DependencyRecord]] = defaultdict(list)
        for record in records.values():
            self._forward[record.parent_id].append(record)
            self._reverse[record.child_id].append(record)
        for adjacency in (self._forward, self._reverse):
            for edges in adjacency.values():
                edges.sort(key=lambda item: item.record_id)
        self._reject_cycles()

    @classmethod
    def from_snapshot(
        cls,
        snapshot: LedgerSnapshot,
        release_id: str,
    ) -> "DependencyGraph":
        release_records = {
            record_id: record
            for record_id, record in snapshot.records.items()
            if record.release_id == release_id
        }
        dependencies = {
            record_id: record
            for record_id, record in release_records.items()
            if isinstance(record, DependencyRecord)
        }
        for dependency in dependencies.values():
            for endpoint in (dependency.parent_id, dependency.child_id):
                if endpoint not in release_records:
                    raise MissingReferenceError(
                        f"Dependency {dependency.record_id} references missing "
                        f"release record {endpoint}"
                    )
        return cls(release_id, dependencies)

    def _reject_cycles(self) -> None:
        temporary: set[str] = set()
        permanent: set[str] = set()

        def visit(node: str) -> None:
            if node in permanent:
                return
            if node in temporary:
                raise DependencyCycleError(
                    f"Dependency cycle detected at {node} in {self.release_id}"
                )
            temporary.add(node)
            for edge in self._forward.get(node, ()):
                visit(edge.child_id)
            temporary.remove(node)
            permanent.add(node)

        nodes = set(self._forward) | set(self._reverse)
        for node in sorted(nodes):
            visit(node)

    def _walk(self, start: str, forward: bool) -> tuple[str, ...]:
        adjacency = self._forward if forward else self._reverse
        visited: set[str] = set()
        pending = [start]
        while pending:
            node = pending.pop()
            for edge in adjacency.get(node, ()):
                neighbor = edge.child_id if forward else edge.parent_id
                if neighbor not in visited and neighbor != start:
                    visited.add(neighbor)
                    pending.append(neighbor)
        return tuple(sorted(visited))

    def ancestors(self, subject_id: str) -> tuple[str, ...]:
        return self._walk(subject_id, forward=False)

    def descendants(self, subject_id: str) -> tuple[str, ...]:
        return self._walk(subject_id, forward=True)

    def required_ancestors(
        self,
        subject_id: str,
    ) -> tuple[DependencyRecord, ...]:
        selected: dict[str, DependencyRecord] = {}
        visited: set[str] = set()
        pending = [subject_id]
        while pending:
            node = pending.pop()
            if node in visited:
                continue
            visited.add(node)
            for edge in self._reverse.get(node, ()):
                if edge.requirement in _REQUIRED:
                    selected[edge.record_id] = edge
                    pending.append(edge.parent_id)
        return tuple(selected[key] for key in sorted(selected))

    def quarantine_impact(
        self,
        subject_ids: Collection[str],
    ) -> tuple[str, ...]:
        affected: set[str] = set()
        for subject_id in subject_ids:
            affected.update(self.descendants(subject_id))
        affected.difference_update(subject_ids)
        return tuple(sorted(affected))
