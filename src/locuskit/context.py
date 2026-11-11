"""Versioned family neighborhood records and validation."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from operator import attrgetter
from pathlib import Path

from .loci import Locus
from .tables import evalue, identifier, integer, read_table, write_table

COLUMNS = (
    "schema_version",
    "record_type",
    "focal_id",
    "species_id",
    "assembly_id",
    "seqid",
    "start",
    "end",
    "strand",
    "window",
    "observed_upstream",
    "observed_downstream",
    "status",
    "offset",
    "neighbor_id",
    "group_id",
    "group_size",
    "evalue_cutoff",
    "group_method",
)
STATUSES = {
    "observed",
    "annotation_missing",
    "locus_missing",
    "orientation_unknown",
    "legacy_unknown",
}

_ROW_FIELDS = attrgetter(*COLUMNS[1:])


@dataclass(frozen=True, slots=True)
class ContextRow:
    record_type: str
    focal_id: str
    species_id: str
    assembly_id: str = ""
    seqid: str = ""
    start: int | None = None
    end: int | None = None
    strand: str = ""
    window: int | None = None
    observed_upstream: int | None = None
    observed_downstream: int | None = None
    status: str = ""
    offset: int | None = None
    neighbor_id: str = ""
    group_id: str = ""
    group_size: int | None = None
    evalue_cutoff: float | None = None
    group_method: str = ""

    def cells(self) -> tuple[object, ...]:
        return ("1", *_ROW_FIELDS(self))


def positioned(record_type: str, focal: str, locus: Locus) -> ContextRow:
    return ContextRow(
        record_type,
        focal,
        locus.species_id,
        locus.assembly_id,
        locus.seqid,
        locus.start,
        locus.end,
        locus.strand,
    )


def has_locus(row: ContextRow) -> bool:
    bits = (
        bool(row.assembly_id),
        bool(row.seqid),
        row.start is not None,
        row.end is not None,
        bool(row.strand),
    )
    if any(bits) != all(bits):
        raise ValueError(f"{row.focal_id}: coordinates must be complete or absent")
    if all(bits):
        assert row.start is not None and row.end is not None
        if row.start < 1 or row.start > row.end or row.strand not in {"+", "-", ".", "?"}:
            raise ValueError(f"{row.focal_id}: invalid coordinates or strand")
    return all(bits)


def validate_context(rows: list[ContextRow], focals: dict[str, str] | None = None) -> None:
    by_focal: dict[str, ContextRow] = {}
    neighbors: dict[str, dict[int, ContextRow]] = defaultdict(dict)
    group_sizes: dict[str, int] = {}
    groups_by_gene: dict[tuple[str, str], str] = {}
    coordinates: dict[str, dict[str, ContextRow]] = defaultdict(dict)
    assemblies: dict[str, str] = {}
    cutoffs: set[float] = set()
    methods: set[str] = set()
    for row in rows:
        locus = has_locus(row)
        if row.record_type == "focal":
            if row.focal_id in by_focal or row.status not in STATUSES:
                raise ValueError(f"{row.focal_id}: duplicate focal row or unknown status")
            if (
                row.offset is not None
                or row.neighbor_id
                or row.group_id
                or row.group_size is not None
            ):
                raise ValueError(f"{row.focal_id}: neighbor fields on a focal row")
            counts = row.observed_upstream, row.observed_downstream
            if row.status == "observed":
                if (
                    not locus
                    or row.strand not in {"+", "-"}
                    or row.window is None
                    or not all(n is not None and 0 <= n <= row.window for n in counts)
                ):
                    raise ValueError(
                        f"{row.focal_id}: observed focal requires locus and window counts"
                    )
            elif counts != (None, None):
                raise ValueError(f"{row.focal_id}: coverage counts require an observed focal")
            if row.status == "orientation_unknown":
                if not locus or row.strand not in {".", "?"}:
                    raise ValueError(f"{row.focal_id}: unknown orientation requires a locus")
            elif row.status != "observed" and locus:
                raise ValueError(f"{row.focal_id}: coordinates conflict with focal status")
            by_focal[row.focal_id] = row
        elif row.record_type == "neighbor":
            if (
                row.status
                or row.window is not None
                or row.observed_upstream is not None
                or row.observed_downstream is not None
                or not row.neighbor_id
                or row.neighbor_id == row.focal_id
                or not row.offset
            ):
                raise ValueError(f"{row.focal_id}: invalid neighbor row")
            if row.offset in neighbors[row.focal_id]:
                raise ValueError(f"{row.focal_id}: repeated offset {row.offset}")
            neighbors[row.focal_id][row.offset] = row
            if row.group_id:
                key = row.species_id, row.neighbor_id
                if groups_by_gene.setdefault(key, row.group_id) != row.group_id:
                    raise ValueError(f"Conflicting homology groups for {key}")
            if row.group_size is not None:
                if (
                    not row.group_id
                    or group_sizes.setdefault(row.group_id, row.group_size) != row.group_size
                ):
                    raise ValueError(f"{row.focal_id}: inconsistent homology group size")
        else:
            raise ValueError(f"Unknown context record type: {row.record_type!r}")
        if locus:
            assert row.start is not None and row.end is not None
            gene_id = row.focal_id if row.record_type == "focal" else row.neighbor_id
            prior = coordinates[row.species_id].setdefault(gene_id, row)
            if prior is not row and (
                prior.assembly_id != row.assembly_id
                or prior.seqid != row.seqid
                or prior.start != row.start
                or prior.end != row.end
                or prior.strand != row.strand
            ):
                raise ValueError(f"Conflicting coordinates for {(row.species_id, gene_id)}")
            if assemblies.setdefault(row.species_id, row.assembly_id) != row.assembly_id:
                raise ValueError(f"Multiple assemblies for species {row.species_id!r}")
        if row.group_method not in {"", "provided_groups", "similarity_hits", "diamond", "legacy"}:
            raise ValueError(f"Unknown homology group method: {row.group_method!r}")
        if row.group_id and not row.group_method:
            raise ValueError(f"{row.focal_id}: group_id requires a group method")
        methods.add(row.group_method)
        if row.evalue_cutoff is not None:
            cutoffs.add(row.evalue_cutoff)
    if len(cutoffs) > 1 or len(methods) > 1:
        raise ValueError("Expected one E-value cutoff and one group method for the context file")
    if set(neighbors) - by_focal.keys():
        raise ValueError(
            f"Neighbor row without a focal row: {min(neighbors.keys() - by_focal.keys())!r}"
        )
    if focals is not None and {k: row.species_id for k, row in by_focal.items()} != focals:
        raise ValueError("Context focal IDs/species differ from the gene-to-species mapping")
    for gene, focal in by_focal.items():
        linked = neighbors.get(gene, {})
        ids = {row.neighbor_id for row in linked.values()}
        if len(ids) != len(linked):
            raise ValueError(f"{gene}: the same neighbor appears at multiple offsets")
        for row in linked.values():
            if row.species_id != focal.species_id:
                raise ValueError(f"{gene}: neighbor species differs from focal")
        if focal.status == "legacy_unknown":
            if any(has_locus(row) for row in linked.values()):
                raise ValueError(f"{gene}: imported legacy neighbors cannot have coordinates")
        elif focal.status == "observed":
            up, down = focal.observed_upstream, focal.observed_downstream
            assert up is not None and down is not None
            if len(linked) != up + down or any(not -up <= offset <= down for offset in linked):
                raise ValueError(f"{gene}: offsets do not match the stated observed coverage")
            previous: tuple[int, int, str] | None = None
            direction = 1 if focal.strand == "+" else -1
            for offset in range(-up, down + 1):
                row = linked[offset] if offset else focal
                if (
                    row.start is None
                    or row.end is None
                    or row.assembly_id != focal.assembly_id
                    or row.seqid != focal.seqid
                ):
                    raise ValueError(f"{gene}: neighbor must occupy the focal assembly and contig")
                order = (row.start, row.end, row.neighbor_id if offset else gene)
                if previous is not None and (
                    (direction == 1 and previous >= order)
                    or (direction == -1 and previous <= order)
                ):
                    raise ValueError(f"{gene}: offsets contradict coordinates or focal strand")
                previous = order
        elif linked:
            raise ValueError(f"{gene}: neighbors supplied for an unavailable focal")
    sizes: dict[str, int] = defaultdict(int)
    for group in groups_by_gene.values():
        sizes[group] += 1
    if any(group_sizes[group] < sizes[group] for group in group_sizes):
        raise ValueError("Homology group size is smaller than the observed group membership")


def parse_row(row: dict[str, str], path: Path, identifiers: dict[str, str]) -> ContextRow:
    if row["schema_version"] != "1":
        raise ValueError(f"{path}: unsupported context schema version {row['schema_version']!r}")

    def num(field: str, minimum: int | None = None) -> int | None:
        return integer(row[field], field, minimum=minimum) if row[field] else None

    def token(field: str, *, required: bool = False) -> str:
        value = row[field]
        if not value and not required:
            return ""
        return identifiers.setdefault(value, identifier(value, field))

    return ContextRow(
        identifiers.setdefault(row["record_type"], row["record_type"]),
        token("focal_id", required=True),
        token("species_id", required=True),
        token("assembly_id"),
        token("seqid"),
        num("start", 1),
        num("end", 1),
        row["strand"],
        num("window", 1),
        num("observed_upstream", 0),
        num("observed_downstream", 0),
        identifiers.setdefault(row["status"], row["status"]),
        num("offset"),
        token("neighbor_id"),
        token("group_id"),
        num("group_size", 1),
        evalue(row["evalue_cutoff"], "evalue_cutoff") if row["evalue_cutoff"] else None,
        identifiers.setdefault(row["group_method"], row["group_method"]),
    )


def read_context(path: Path, focals: dict[str, str] | None = None) -> list[ContextRow]:
    identifiers: dict[str, str] = {}
    rows = [parse_row(row, path, identifiers) for row in read_table(path, COLUMNS, exact=True)]
    identifiers.clear()
    validate_context(rows, focals)
    return rows


def write_context(path: Path, rows: list[ContextRow], inputs: list[Path]) -> None:
    validate_context(rows)
    write_table(path, COLUMNS, (row.cells() for row in rows), inputs)
