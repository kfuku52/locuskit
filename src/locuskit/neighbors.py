"""Extract local neighborhoods with optional precomputed homology evidence."""

from __future__ import annotations

import re
import subprocess
import tempfile
from collections import Counter, defaultdict
from contextlib import nullcontext
from dataclasses import replace
from pathlib import Path

from .context import ContextRow, positioned, write_context
from .loci import Locus, read_focals, read_loci
from .tables import evalue, identifier, open_text, read_table

GeneKey = tuple[str, str]
GroupMap = dict[GeneKey, tuple[str, int]]


def extract(focals: dict[str, str], loci: list[Locus], window: int) -> list[ContextRow]:
    if window < 1:
        raise ValueError("Window must be a positive integer")
    contigs: dict[tuple[str, str], list[Locus]] = defaultdict(list)
    for locus in loci:
        contigs[locus.species_id, locus.seqid].append(locus)
    positions: dict[GeneKey, tuple[list[Locus], int]] = {}
    annotated_species = {locus.species_id for locus in loci}
    for contig in contigs.values():
        contig.sort(key=lambda loc: (loc.start, loc.end, loc.gene_id))
        for index, locus in enumerate(contig):
            if focals.get(locus.gene_id) == locus.species_id:
                positions[locus.species_id, locus.gene_id] = contig, index
    result: list[ContextRow] = []
    for gene, species in sorted(focals.items()):
        position = positions.get((species, gene))
        if position is None:
            status = "locus_missing" if species in annotated_species else "annotation_missing"
            result.append(ContextRow("focal", gene, species, window=window, status=status))
            continue
        contig, index = position
        focal = positioned("focal", gene, contig[index])
        if focal.strand not in {"+", "-"}:
            result.append(replace(focal, window=window, status="orientation_unknown"))
            continue
        sign = 1 if focal.strand == "+" else -1
        left = min(window, index)
        right = min(window, len(contig) - index - 1)
        up, down = (left, right) if sign == 1 else (right, left)
        result.append(
            replace(
                focal,
                window=window,
                status="observed",
                observed_upstream=up,
                observed_downstream=down,
            )
        )
        for offset in range(-up, down + 1):
            if offset:
                neighbor = contig[index + sign * offset]
                result.append(
                    replace(
                        positioned("neighbor", gene, neighbor),
                        offset=offset,
                        neighbor_id=neighbor.gene_id,
                    )
                )
    return result


def read_groups(path: Path) -> GroupMap:
    assignments: dict[GeneKey, str] = {}
    for row in read_table(path, ("species_id", "gene_id", "group_id")):
        key = identifier(row["species_id"], "species_id"), identifier(row["gene_id"], "gene_id")
        if key in assignments:
            raise ValueError(f"{path}: repeated homology assignment for {key}")
        assignments[key] = identifier(row["group_id"], "group_id")
    counts = Counter(assignments.values())
    return {key: (group, counts[group]) for key, group in assignments.items()}


def groups_from_hits(
    path: Path, neighbors: set[GeneKey], cutoff: float, *, searched: set[str] | None = None
) -> GroupMap:
    """Cluster observed neighbors using headerless DIAMOND qseqid/sseqid/evalue."""
    by_id: dict[str, GeneKey] = {}
    for key in sorted(neighbors):
        gene = key[1]
        if gene in by_id:
            raise ValueError(f"{path}: hit format needs unique sequence IDs, repeated {gene!r}")
        by_id[gene] = key
    parent = {gene: gene for gene in (searched or ())}

    def root(gene: str) -> str:
        while parent[gene] != gene:
            parent[gene] = parent[parent[gene]]
            gene = parent[gene]
        return gene

    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 3:
                raise ValueError(f"{path}:{line_number}: expected qseqid, sseqid and evalue")
            a, b, value = fields
            identifier(a, "qseqid")
            identifier(b, "sseqid")
            score = evalue(value, "hit evalue")
            if searched is not None and (a not in searched or b not in searched):
                raise ValueError(
                    f"{path}:{line_number}: DIAMOND returned an unexpected sequence ID"
                )
            for endpoint in (a, b):
                if endpoint in by_id and searched is None:
                    parent.setdefault(endpoint, endpoint)
            if score <= cutoff and a in parent and b in parent:
                aa, bb = sorted((root(a), root(b)))
                if aa != bb:
                    parent[bb] = aa
    components: dict[str, list[str]] = defaultdict(list)
    for gene in sorted(parent):
        components[root(gene)].append(gene)
    ordered = sorted(components.values(), key=lambda genes: (-len(genes), genes[0]))
    return {
        by_id[gene]: (f"SG{index:06d}", len(genes))
        for index, genes in enumerate(ordered, 1)
        for gene in genes
    }


def select_proteins(source: Path, output: Path, wanted: set[str]) -> set[str]:
    """Stream the protein input and retain only sequences for selected loci."""
    seen: set[str] = set()
    gene = ""
    has_sequence = False
    with open_text(source) as stream, output.open("w", encoding="utf-8", newline="") as selected:
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
                    selected.write(f">{gene}\n")
            elif not gene:
                raise ValueError(f"{source}:{number}: sequence before FASTA header")
            elif gene in wanted:
                if re.fullmatch(r"[A-Za-z*]+", line) is None:
                    raise ValueError(f"{source}:{number}: invalid protein FASTA sequence")
                selected.write(line + "\n")
                has_sequence = True
        if gene in wanted and not has_sequence:
            raise ValueError(f"{source}: empty protein sequence for {gene!r}")
    return seen


def groups_from_proteins(
    source: Path, neighbors: set[GeneKey], cutoff: float, threads: int
) -> GroupMap:
    with tempfile.TemporaryDirectory(prefix="locuskit-diamond-") as temp:
        directory = Path(temp)
        proteins = directory / "neighbors.faa"
        selected = select_proteins(source, proteins, {gene for _, gene in neighbors})
        hits = directory / "hits.tsv"
        if len(selected) < 2:
            hits.touch()
        else:
            database = directory / "neighbors"
            commands = [
                [
                    "diamond",
                    "makedb",
                    "--in",
                    str(proteins),
                    "--db",
                    str(database),
                    "--threads",
                    str(threads),
                ],
                [
                    "diamond",
                    "blastp",
                    "--query",
                    str(proteins),
                    "--db",
                    str(database),
                    "--out",
                    str(hits),
                    "--outfmt",
                    "6",
                    "qseqid",
                    "sseqid",
                    "evalue",
                    "--evalue",
                    str(cutoff),
                    "--max-target-seqs",
                    "2000",
                    "--threads",
                    str(threads),
                ],
            ]
            for command in commands:
                result = subprocess.run(command, capture_output=True, text=True, check=False)
                if result.returncode:
                    raise ValueError(f"DIAMOND {command[1]} failed: {result.stderr.strip()}")
        return groups_from_hits(hits, neighbors, cutoff, searched=selected)


def build_context(
    focal_path: Path,
    locus_path: Path,
    output: Path,
    window: int,
    groups: Path | None = None,
    hits: Path | None = None,
    cutoff: float | None = None,
    proteins: Path | None = None,
    threads: int = 1,
) -> None:
    sources = _homology_inputs(groups, hits, cutoff, proteins, threads)
    focals = read_focals(focal_path)
    # Raw hit endpoints omit species, so collision checking needs all loci.
    species = None if hits is not None or proteins is not None else set(focals.values())
    loci = read_loci(locus_path, species)
    rows = _family_context(focals, loci, window, groups, hits, cutoff, proteins, threads)
    write_context(output, rows, [focal_path, locus_path, *sources])


def _homology_inputs(
    groups: Path | None,
    hits: Path | None,
    cutoff: float | None,
    proteins: Path | None,
    threads: int,
) -> list[Path]:
    sources = [path for path in (groups, hits, proteins) if path is not None]
    if len(sources) > 1:
        raise ValueError("Choose one homology input")
    if cutoff is not None and not sources:
        raise ValueError("evalue_cutoff requires a homology input")
    if threads < 1:
        raise ValueError("threads must be positive")
    if proteins is not None and cutoff == 0:
        raise ValueError("DIAMOND evalue_cutoff must be positive")
    return sources


def _family_context(
    focals: dict[str, str],
    loci: list[Locus],
    window: int,
    groups: Path | None,
    hits: Path | None,
    cutoff: float | None,
    proteins: Path | None,
    threads: int,
    *,
    search_proteins: Path | None = None,
) -> list[ContextRow]:
    rows = extract(focals, loci, window)
    assignments: GroupMap = {}
    method = ""
    if groups is not None:
        assignments = read_groups(groups)
        method = "provided_groups"
    elif hits is not None or proteins is not None:
        if cutoff is None:
            cutoff = 0.01
        keys = {(row.species_id, row.neighbor_id) for row in rows if row.record_type == "neighbor"}
        selected = {gene for _, gene in keys}
        counts = Counter(loc.gene_id for loc in loci if loc.gene_id in selected)
        conflicts = sorted(gene for gene, count in counts.items() if count > 1)
        if conflicts:
            raise ValueError(
                f"{hits or proteins}: repeated sequence ID across species: {conflicts[0]!r}"
            )
        if hits is not None:
            assignments = groups_from_hits(hits, keys, cutoff)
            method = "similarity_hits"
        else:
            assert proteins is not None
            assignments = groups_from_proteins(search_proteins or proteins, keys, cutoff, threads)
            method = "diamond"
    if method:
        for i, row in enumerate(rows):
            group, size = assignments.get((row.species_id, row.neighbor_id), ("", None))
            rows[i] = replace(
                row, group_id=group, group_size=size, evalue_cutoff=cutoff, group_method=method
            )
    return rows


def _file_identity(path: Path) -> tuple[int, int] | None:
    try:
        stat = path.stat()
    except FileNotFoundError:
        return None
    return stat.st_dev, stat.st_ino


def _check_batch_outputs(outputs: list[Path], inputs: list[Path]) -> None:
    paths = {path.resolve() for path in inputs}
    identities = {identity for path in inputs if (identity := _file_identity(path)) is not None}
    for output in outputs:
        resolved, identity = output.resolve(), _file_identity(output)
        if output.is_symlink() or (identity is not None and not output.is_file()):
            raise ValueError(f"Invalid neighborhood output: {output}")
        if resolved in paths or (identity is not None and identity in identities):
            raise ValueError(f"Neighborhood output aliases another output or an input: {output}")
        paths.add(resolved)
        if identity is not None:
            identities.add(identity)
    resolved_outputs = {output.resolve() for output in outputs}
    for output in resolved_outputs:
        if any(parent in resolved_outputs for parent in output.parents):
            raise ValueError(f"Neighborhood output would also be a directory: {output}")


def _batch_proteins(
    mappings: list[dict[str, str]], loci: list[Locus], window: int, source: Path, output: Path
) -> None:
    focals: dict[str, dict[str, str]] = defaultdict(dict)
    for mapping in mappings:
        for gene, species in mapping.items():
            focals[species][gene] = species
    by_species: dict[str, list[Locus]] = defaultdict(list)
    for locus in loci:
        if locus.species_id in focals:
            by_species[locus.species_id].append(locus)
    wanted = {
        row.neighbor_id
        for species, mapping in focals.items()
        for row in extract(mapping, by_species[species], window)
        if row.record_type == "neighbor"
    }
    select_proteins(source, output, wanted)


def build_batch(
    manifest: Path,
    locus_path: Path,
    window: int,
    groups: Path | None = None,
    hits: Path | None = None,
    cutoff: float | None = None,
    proteins: Path | None = None,
    threads: int = 1,
) -> None:
    """Write independent neighborhood tables using one shared locus read."""
    sources = _homology_inputs(groups, hits, cutoff, proteins, threads)
    tasks: list[tuple[Path, Path]] = []
    for row in read_table(manifest, ("focal_tsv", "out_file")):
        if not row["focal_tsv"] or not row["out_file"]:
            raise ValueError(f"{manifest}: expected focal_tsv and out_file paths")
        tasks.append((manifest.parent / row["focal_tsv"], manifest.parent / row["out_file"]))
    if not tasks:
        raise ValueError(f"{manifest}: no neighborhood jobs")
    shared = [manifest, locus_path, *sources]
    _check_batch_outputs([output for _, output in tasks], [*shared, *(f for f, _ in tasks)])
    mappings = [read_focals(focal_path) for focal_path, _ in tasks]
    species = None
    if hits is None and proteins is None:
        species = {taxon for mapping in mappings for taxon in mapping.values()}
    loci = read_loci(locus_path, species)
    workspace = (
        tempfile.TemporaryDirectory(prefix="locuskit-batch-")
        if proteins is not None and len(tasks) > 1
        else nullcontext()
    )
    with workspace as temp:
        selected = None
        if temp is not None:
            assert proteins is not None
            shared_proteins = Path(temp) / "selected.faa"
            try:
                _batch_proteins(mappings, loci, window, proteins, shared_proteins)
            except ValueError:
                # Report malformed input at its original family, after earlier outputs.
                pass
            else:
                selected = shared_proteins
        for (focal_path, output), mapping in zip(tasks, mappings, strict=True):
            rows = _family_context(
                mapping,
                loci,
                window,
                groups,
                hits,
                cutoff,
                proteins,
                threads,
                search_proteins=selected,
            )
            write_context(output, rows, [focal_path, *shared])
