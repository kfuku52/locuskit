"""Normalize genomic locus coordinates for local neighborhood extraction."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from kffractbias.io import annotation_locus_order

from .tables import identifier, integer, read_table, write_table

LOCUS_COLUMNS = ("species_id", "assembly_id", "gene_id", "seqid", "start", "end", "strand")
FOCAL_COLUMNS = ("gene", "species")


@dataclass(frozen=True, slots=True)
class Locus:
    species_id: str
    assembly_id: str
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str


def read_focals(path: Path) -> dict[str, str]:
    focals: dict[str, str] = {}
    for row in read_table(path, FOCAL_COLUMNS):
        gene = identifier(row["gene"], "gene")
        species = identifier(row["species"], "species")
        if gene in focals:
            raise ValueError(f"{path}: duplicate focal gene {gene!r}")
        focals[gene] = species
    if not focals:
        raise ValueError(f"{path}: no focal genes")
    return focals


def read_loci(path: Path, species: set[str] | None = None) -> list[Locus]:
    loci: list[Locus] = []
    seen: set[tuple[str, str]] = set()
    assemblies: dict[str, str] = {}
    for row in read_table(path, LOCUS_COLUMNS):
        locus = Locus(
            identifier(row["species_id"], "species_id"),
            identifier(row["assembly_id"], "assembly_id"),
            identifier(row["gene_id"], "gene_id"),
            identifier(row["seqid"], "seqid"),
            integer(row["start"], "start", minimum=1),
            integer(row["end"], "end", minimum=1),
            row["strand"],
        )
        if locus.start > locus.end or locus.strand not in {"+", "-", ".", "?"}:
            raise ValueError(f"{path}: invalid coordinates or strand for {locus.gene_id!r}")
        key = locus.species_id, locus.gene_id
        if key in seen:
            raise ValueError(f"{path}: duplicate locus {key}")
        if assemblies.setdefault(locus.species_id, locus.assembly_id) != locus.assembly_id:
            raise ValueError(f"{path}: multiple assemblies for species {locus.species_id!r}")
        seen.add(key)
        if species is None or locus.species_id in species:
            loci.append(locus)
    return loci


def extract_loci(manifest: Path, output: Path, aliases: Path | None = None) -> None:
    """Export one row per genomic locus using the existing FractBias GFF reader."""
    rename: dict[tuple[str, str], str] = {}
    if aliases is not None:
        for row in read_table(aliases, ("species_id", "locus_id", "gene_id")):
            key = (
                identifier(row["species_id"], "species_id"),
                identifier(row["locus_id"], "locus_id"),
            )
            if key in rename:
                raise ValueError(f"{aliases}: duplicate locus alias {key}")
            rename[key] = identifier(row["gene_id"], "gene_id")
    inputs = [manifest, *([aliases] if aliases is not None else [])]
    genomes: list[tuple[str, str, Path, str, str]] = []
    seen: set[str] = set()
    for row in read_table(manifest, ("species_id", "assembly_id", "gff")):
        species = identifier(row["species_id"], "species_id")
        assembly = identifier(row["assembly_id"], "assembly_id")
        if species in seen:
            raise ValueError(f"{manifest}: duplicate species {species!r}")
        seen.add(species)
        if not row["gff"]:
            raise ValueError(f"{manifest}: missing GFF path for species {species!r}")
        source = manifest.parent / row["gff"]
        inputs.append(source)
        feature, attribute = row.get("feature") or "gene", row.get("attribute") or "ID"
        genomes.append((species, assembly, source, feature, attribute))
    if not genomes:
        raise ValueError(f"{manifest}: no genomes")

    def output_rows() -> Iterator[tuple[str, str, str, str, int, int, str]]:
        for species, assembly, source, feature, attribute in genomes:
            genes = annotation_locus_order(source, feature=feature, attribute=attribute)
            if not genes:
                raise ValueError(f"{source}: no loci for species {species!r}")
            identifiers: set[str] = set()
            for gene in genes:
                gene_id = identifier(gene.gene_id, "locus_id")
                gene_id = rename.pop((species, gene_id), gene_id)
                if gene_id in identifiers:
                    raise ValueError(f"{source}: locus alias collides with gene {gene_id!r}")
                identifiers.add(gene_id)
                yield (
                    species,
                    assembly,
                    gene_id,
                    identifier(gene.seqid, "seqid"),
                    gene.start + 1,
                    gene.end,
                    gene.strand,
                )
        if rename:
            raise ValueError(f"{aliases}: alias for an unknown locus {min(rename)}")

    write_table(output, LOCUS_COLUMNS, output_rows(), inputs)
