"""Cluster similarity hits for the genes observed in a family neighborhood."""

from __future__ import annotations

import sqlite3
from collections import defaultdict
from collections.abc import Iterable, Iterator
from contextlib import closing
from pathlib import Path

from .tables import evalue, identifier, open_text

GeneKey = tuple[str, str]
GroupMap = dict[GeneKey, tuple[str, int]]


def _read_hits(path: Path) -> Iterator[tuple[int, str, str, float]]:
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 3:
                raise ValueError(f"{path}:{line_number}: expected qseqid, sseqid and evalue")
            a, b, value = fields
            yield (
                line_number,
                identifier(a, "qseqid"),
                identifier(b, "sseqid"),
                evalue(value, "hit evalue"),
            )


def _genes_by_id(source: Path, neighbors: set[GeneKey]) -> dict[str, GeneKey]:
    by_id: dict[str, GeneKey] = {}
    for key in sorted(neighbors):
        gene = key[1]
        if gene in by_id:
            raise ValueError(f"{source}: hit format needs unique sequence IDs, repeated {gene!r}")
        by_id[gene] = key
    return by_id


def _components(
    by_id: dict[str, GeneKey], parent: dict[str, str], edges: Iterable[tuple[str, str]]
) -> GroupMap:
    def root(gene: str) -> str:
        while parent[gene] != gene:
            parent[gene] = parent[parent[gene]]
            gene = parent[gene]
        return gene

    for a, b in edges:
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


def groups_from_hits(
    path: Path, neighbors: set[GeneKey], cutoff: float, *, searched: set[str] | None = None
) -> GroupMap:
    """Cluster observed neighbors using headerless DIAMOND qseqid/sseqid/evalue."""
    by_id = _genes_by_id(path, neighbors)
    parent = {gene: gene for gene in (searched or ())}

    def edges() -> Iterator[tuple[str, str]]:
        for line_number, a, b, score in _read_hits(path):
            if searched is None:
                for endpoint in (a, b):
                    if endpoint in by_id:
                        parent.setdefault(endpoint, endpoint)
            elif a not in searched or b not in searched:
                raise ValueError(
                    f"{path}:{line_number}: DIAMOND returned an unexpected sequence ID"
                )
            if score <= cutoff and a in parent and b in parent:
                yield a, b

    return _components(by_id, parent, edges())


def index_hits(source: Path, index: Path, wanted: set[str], cutoff: float) -> None:
    """Validate the shared input once and store the part relevant to the batch."""
    with closing(sqlite3.connect(index)) as db, db:
        db.execute("CREATE TABLE observed (gene TEXT PRIMARY KEY) WITHOUT ROWID")
        db.execute("CREATE TABLE edge (a TEXT, b TEXT, PRIMARY KEY (a, b)) WITHOUT ROWID")
        for _, a, b, score in _read_hits(source):
            keep_a = a in wanted
            keep_b = b in wanted
            if keep_a:
                db.execute("INSERT OR IGNORE INTO observed VALUES (?)", (a,))
            if keep_b:
                db.execute("INSERT OR IGNORE INTO observed VALUES (?)", (b,))
            if keep_a and keep_b and a != b and score <= cutoff:
                db.execute("INSERT OR IGNORE INTO edge VALUES (?, ?)", tuple(sorted((a, b))))


def groups_from_index(index: Path, source: Path, neighbors: set[GeneKey]) -> GroupMap:
    by_id = _genes_by_id(source, neighbors)
    with closing(sqlite3.connect(index)) as db:
        parent = {
            gene: gene
            for gene in by_id
            if db.execute("SELECT 1 FROM observed WHERE gene = ?", (gene,)).fetchone()
        }
        edges = (
            (a, b)
            for a in by_id
            for (b,) in db.execute("SELECT b FROM edge WHERE a = ?", (a,))
            if b in by_id
        )
        return _components(by_id, parent, edges)
