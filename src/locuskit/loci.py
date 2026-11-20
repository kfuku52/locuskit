"""Normalize genomic locus coordinates for local neighborhood extraction."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from kffractbias.io import annotation_locus_order

from .replicons import locus_span, read_replicons
from .tables import identifier, integer, read_table, write_table

LOCUS_COLUMNS = ("species_id", "assembly_id", "gene_id", "seqid", "start", "end", "strand")
REPLICON_COLUMNS = ("topology", "length", "intervals")
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
    circular_length: int | None = None


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
    replicons: dict[tuple[str, str, str], tuple[str, int]] = {}
    extended = None
    for row in read_table(path, LOCUS_COLUMNS):
        if extended is None:
            extended = any(k in row for k in REPLICON_COLUMNS)
            if extended and not all(k in row for k in REPLICON_COLUMNS):
                raise ValueError(f"{path}: incomplete replicon columns")
        topology = row.get("topology", "linear")
        length = None
        if extended:
            if topology not in {
                "linear",
                "circular",
            }:
                raise ValueError(f"{path}: invalid replicon columns")
            length = integer(row["length"], "length", minimum=1)
            pieces = [part.split("-") for part in row["intervals"].split(",")]
            if any(len(piece) != 2 for piece in pieces):
                raise ValueError(f"{path}: invalid locus intervals")
            parts = {
                (integer(a, "start", minimum=1), integer(b, "end", minimum=1)) for a, b in pieces
            }
            if any(b > length for _, b in parts):
                raise ValueError(f"{path}: interval outside replicon")
            start, end, intervals = locus_span(parts, length, topology, row["strand"])
            if (str(start), str(end), intervals) != (row["start"], row["end"], row["intervals"]):
                raise ValueError(f"{path}: inconsistent locus span or interval order")
            contig = row["species_id"], row["assembly_id"], row["seqid"]
            if replicons.setdefault(contig, (topology, length)) != (topology, length):
                raise ValueError(f"{path}: inconsistent replicon metadata for {contig}")
        locus = Locus(
            identifier(row["species_id"], "species_id"),
            identifier(row["assembly_id"], "assembly_id"),
            identifier(row["gene_id"], "gene_id"),
            identifier(row["seqid"], "seqid"),
            integer(row["start"], "start", minimum=1),
            integer(row["end"], "end", minimum=1),
            row["strand"],
            length if topology == "circular" else None,
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


def extract_loci(
    manifest: Path, output: Path, aliases: Path | None = None, replicons_tsv: Path | None = None
) -> None:
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
    replicons = read_replicons(replicons_tsv) if replicons_tsv else {}
    used_replicons = set()
    if replicons_tsv:
        inputs.append(replicons_tsv)
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

    def output_rows() -> Iterator[tuple]:
        for species, assembly, source, feature, attribute in genomes:
            parts: dict[str, set[tuple[int, int]]] | None = {} if replicons_tsv else None
            genes = annotation_locus_order(
                source, feature=feature, attribute=attribute, intervals_by_locus=parts
            )
            if not genes:
                raise ValueError(f"{source}: no loci for species {species!r}")
            identifiers: set[str] = set()
            extended_rows: list[tuple] = []
            for gene in genes:
                gene_id = identifier(gene.gene_id, "locus_id")
                gene_id = rename.pop((species, gene_id), gene_id)
                if gene_id in identifiers:
                    raise ValueError(f"{source}: locus alias collides with gene {gene_id!r}")
                identifiers.add(gene_id)
                extra: tuple = ()
                start, end = gene.start + 1, gene.end
                if parts is not None:
                    key = species, assembly, gene.seqid
                    if key not in replicons:
                        raise ValueError(f"{source}: missing replicon metadata for {key}")
                    topology, length = replicons[key]
                    used_replicons.add(key)
                    start, end, intervals = locus_span(
                        parts[gene.gene_id], length, topology, gene.strand
                    )
                    extra = (topology, length, intervals)
                row_values = (
                    species,
                    assembly,
                    gene_id,
                    identifier(gene.seqid, "seqid"),
                    start,
                    end,
                    gene.strand,
                    *extra,
                )
                if parts is None:
                    yield row_values
                else:
                    extended_rows.append(row_values)
            yield from sorted(extended_rows, key=lambda item: (item[3], item[4], item[5], item[2]))
        if replicons.keys() - used_replicons:
            raise ValueError(
                f"{replicons_tsv}: replicon without loci: {min(replicons.keys() - used_replicons)}"
            )
        if rename:
            raise ValueError(f"{aliases}: alias for an unknown locus {min(rename)}")

    write_table(
        output, LOCUS_COLUMNS + (REPLICON_COLUMNS if replicons_tsv else ()), output_rows(), inputs
    )
