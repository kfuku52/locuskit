"""Select representative proteins using explicit locus identifiers."""

from __future__ import annotations

import re
from collections.abc import Collection
from pathlib import Path
from typing import TextIO

from .tables import identifier, open_text, read_table, staged_text


def _copy_proteins(
    source: Path,
    selected: TextIO,
    wanted: Collection[str],
    seen: set[str],
    aliases: dict[str, str] | None = None,
) -> None:
    gene = ""
    has_sequence = False
    with open_text(source) as stream:
        for number, line in enumerate(stream, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if gene in wanted and not has_sequence:
                    raise ValueError(f"{source}: empty protein sequence for {gene!r}")
                header = line[1:].split()
                if not header:
                    raise ValueError(f"{source}:{number}: empty FASTA identifier")
                gene = identifier(header[0], "protein FASTA identifier")
                has_sequence = False
                if gene in wanted:
                    if gene in seen:
                        raise ValueError(f"{source}: duplicate protein sequence for {gene!r}")
                    seen.add(gene)
                    output_id = aliases[gene] if aliases is not None else gene
                    selected.write(f">{output_id}\n")
            elif not gene:
                raise ValueError(f"{source}:{number}: sequence before FASTA header")
            elif gene in wanted:
                if re.fullmatch(r"[A-Za-z*]+", line) is None:
                    raise ValueError(f"{source}:{number}: invalid protein FASTA sequence")
                selected.write(line + "\n")
                has_sequence = True
        if gene in wanted and not has_sequence:
            raise ValueError(f"{source}: empty protein sequence for {gene!r}")


def select_proteins(source: Path, output: Path, wanted: set[str]) -> set[str]:
    """Stream the protein input and retain only sequences for selected loci."""
    seen: set[str] = set()
    with output.open("w", encoding="utf-8", newline="") as selected:
        _copy_proteins(source, selected, wanted, seen)
    return seen


def map_proteins(sources: list[Path], alias_path: Path, output: Path) -> None:
    aliases: dict[str, str] = {}
    loci: set[str] = set()
    for row in read_table(alias_path, ("protein_id", "gene_id")):
        protein = identifier(row["protein_id"], "protein_id")
        gene = identifier(row["gene_id"], "gene_id")
        if protein in aliases:
            raise ValueError(f"{alias_path}: duplicate protein ID {protein!r}")
        if gene in loci:
            raise ValueError(f"{alias_path}: duplicate gene ID {gene!r}")
        aliases[protein] = gene
        loci.add(gene)
    if not aliases:
        raise ValueError(f"{alias_path}: no protein aliases")
    seen: set[str] = set()
    with staged_text(output, [alias_path, *sources]) as selected:
        for source in sources:
            _copy_proteins(source, selected, aliases, seen, aliases)
        missing = aliases.keys() - seen
        if missing:
            raise ValueError(
                f"{alias_path}: missing {len(missing)} protein IDs, first {min(missing)!r}"
            )
