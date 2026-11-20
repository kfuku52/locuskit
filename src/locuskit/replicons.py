"""Coordinates on explicitly sized linear or circular genomic replicons."""

from __future__ import annotations

from pathlib import Path

from .tables import identifier, integer, read_table


def read_replicons(path: Path) -> dict[tuple[str, str, str], tuple[str, int]]:
    result = {}
    for row in read_table(path, ("species_id", "assembly_id", "seqid", "topology", "length")):
        key = (
            identifier(row["species_id"], "species_id"),
            identifier(row["assembly_id"], "assembly_id"),
            identifier(row["seqid"], "seqid"),
        )
        topology = row["topology"]
        if key in result or topology not in {"circular", "linear"}:
            raise ValueError(f"{path}: duplicate replicon or invalid topology for {key}")
        result[key] = (topology, integer(row["length"], "length", minimum=1))
    return result


def locus_span(
    parts: set[tuple[int, int]], length: int, topology: str, strand: str
) -> tuple[int, int, str]:
    """Represent a wrapped span by its physical start and an unrolled end."""
    intervals: set[tuple[int, int]] = set()
    for start, end in parts:
        if not 1 <= start <= min(length, end) or end - start >= length:
            raise ValueError(f"Interval {start}-{end} outside replicon of length {length}")
        if end <= length:
            intervals.add((start, end))
        elif topology == "circular":
            intervals.update(((start, length), (1, end - length)))
        else:
            raise ValueError(f"Linear interval {start}-{end} crosses replicon end {length}")
    ordered = sorted(intervals)
    start, end = ordered[0][0], max(b for _, b in ordered)
    if topology == "circular" and len(ordered) > 1:
        # Put the largest uncovered gap outside the span. Ties favor the lowest start.
        gap = start + length - end - 1
        right = ordered[0][1]
        for a, b in ordered[1:]:
            if a - right - 1 > gap:
                start, end, gap = a, right + length, a - right - 1
            right = max(right, b)
        ordered.sort(key=lambda pair: (pair[0] < start, pair))
    if strand == "-":
        ordered.reverse()
    return start, end, ",".join(f"{a}-{b}" for a, b in ordered)
