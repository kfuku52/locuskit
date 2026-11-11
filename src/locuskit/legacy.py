"""Import and export the GeneGalleon tree_plot neighborhood table."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from .context import ContextRow, read_context, write_context
from .loci import read_focals
from .tables import evalue, identifier, integer, read_table, write_table

LEGACY_COLUMNS = (
    "node_name",
    "species",
    "direction",
    "offset",
    "neighbor_gene",
    "group_id",
    "group_size",
)


def import_legacy(source: Path, target: Path, focal_path: Path | None = None) -> None:
    expected = read_focals(focal_path) if focal_path is not None else None
    focals: dict[str, str] = dict(expected or {})
    by_position: dict[tuple[str, int], ContextRow] = {}
    cutoffs: set[float] = set()
    for row in read_table(source, ("node_name", "offset", "neighbor_gene", "group_id")):
        gene = identifier(row["node_name"], "node_name")
        if expected is not None and gene not in expected:
            raise ValueError(f"{source}: unknown focal gene {gene!r}")
        species = identifier(row.get("species") or focals.get(gene, ""), "species")
        if focals.setdefault(gene, species) != species:
            raise ValueError(f"{source}: conflicting species for {gene!r}")
        offset = integer(row["offset"], "offset")
        direction = "upstream" if offset < 0 else "downstream"
        if not offset or (row.get("direction") and row["direction"] != direction):
            raise ValueError(f"{source}: invalid direction or offset for {gene!r}")
        cutoff = evalue(row["evalue_cutoff"], "evalue_cutoff") if row.get("evalue_cutoff") else None
        if cutoff is not None:
            cutoffs.add(cutoff)
        neighbor = ContextRow(
            "neighbor",
            gene,
            species,
            offset=offset,
            neighbor_id=identifier(row["neighbor_gene"], "neighbor_gene"),
            group_id=identifier(row["group_id"], "group_id"),
            group_size=integer(row["group_size"], "group_size", minimum=1)
            if row.get("group_size")
            else None,
            evalue_cutoff=cutoff,
            group_method="legacy",
        )
        if by_position.setdefault((gene, offset), neighbor) != neighbor:
            raise ValueError(f"{source}: conflicting rows for {gene!r} offset {offset}")
    if len(cutoffs) > 1:
        raise ValueError(f"{source}: mixed E-value cutoffs")
    resolved = next(iter(cutoffs), None)
    rows = [
        ContextRow(
            "focal",
            gene,
            species,
            status="legacy_unknown",
            evalue_cutoff=resolved,
            group_method="legacy",
        )
        for gene, species in sorted(focals.items())
    ]
    rows.extend(replace(row, evalue_cutoff=resolved) for _, row in sorted(by_position.items()))
    write_context(target, rows, [source, *([focal_path] if focal_path is not None else [])])


def export_legacy(source: Path, target: Path) -> None:
    rows = read_context(source)
    cutoff = next((row.evalue_cutoff for row in rows if row.evalue_cutoff is not None), None)
    data: list[tuple[object, ...]] = []
    for row in rows:
        if row.record_type != "neighbor" or not row.group_id:
            continue
        assert row.offset is not None
        direction = "upstream" if row.offset < 0 else "downstream"
        base = (
            row.focal_id,
            row.species_id,
            direction,
            row.offset,
            row.neighbor_id,
            row.group_id,
            row.group_size,
        )
        data.append(base if cutoff is None else (*base, cutoff))
    header = LEGACY_COLUMNS if cutoff is None else (*LEGACY_COLUMNS, "evalue_cutoff")
    write_table(target, header, sorted(data), [source])
