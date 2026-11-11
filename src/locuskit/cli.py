"""Command line for LOCUSKIT and the inherited FractBias commands."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .context import read_context
from .legacy import export_legacy, import_legacy
from .loci import extract_loci, read_focals
from .neighbors import build_context
from .tables import evalue


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(
        prog="locuskit",
        description="LOCUSKIT: genomic neighborhoods, collinearity and fractionation analysis.",
    )
    root.add_argument("--version", action="version", version=f"LOCUSKIT {__version__}")
    commands = root.add_subparsers(dest="command", required=True)

    loci = commands.add_parser("loci", help="Export ordered loci from a genome annotation manifest")
    loci.add_argument("--manifest_tsv", required=True, type=Path)
    loci.add_argument("--aliases_tsv", type=Path)
    loci.add_argument("--out_file", required=True, type=Path)

    neighbors = commands.add_parser("neighbors", help="Build a versioned family neighborhood table")
    neighbors.add_argument("--focal_tsv", required=True, type=Path)
    neighbors.add_argument("--loci_tsv", required=True, type=Path)
    neighbors.add_argument(
        "--window", default=5, type=int, help="Maximum loci on each side of a focal gene"
    )
    homology = neighbors.add_mutually_exclusive_group()
    homology.add_argument(
        "--groups_tsv", type=Path, help="Headered species_id/gene_id/group_id table"
    )
    homology.add_argument(
        "--hits_tsv", type=Path, help="Headerless DIAMOND qseqid/sseqid/evalue output"
    )
    homology.add_argument("--protein_fasta", type=Path, help="Run DIAMOND on selected proteins")
    neighbors.add_argument("--threads", type=int, default=1, help="Threads for DIAMOND")
    neighbors.add_argument(
        "--evalue_cutoff", help="Hit cutoff (default 0.01 when hits are provided)"
    )
    neighbors.add_argument("--out_file", required=True, type=Path)

    validate = commands.add_parser("validate", help="Check a versioned family neighborhood table")
    validate.add_argument("--in_file", required=True, type=Path)
    validate.add_argument(
        "--focal_tsv", type=Path, help="Check exact focal IDs and species against this mapping"
    )

    importer = commands.add_parser(
        "import-legacy", help="Import the GeneGalleon tree_plot synteny TSV"
    )
    importer.add_argument("--in_file", required=True, type=Path)
    importer.add_argument("--focal_tsv", type=Path)
    importer.add_argument("--out_file", required=True, type=Path)

    exporter = commands.add_parser(
        "export-legacy", help="Export the GeneGalleon tree_plot synteny TSV"
    )
    exporter.add_argument("--in_file", required=True, type=Path)
    exporter.add_argument("--out_file", required=True, type=Path)

    commands.add_parser("fractbias", help="Run the kfFractBias comparison and retention commands")
    return root


def main(argv: Sequence[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == "fractbias":
        from kffractbias.cli import main as fractbias

        return fractbias(arguments[1:])
    args = parser().parse_args(arguments)
    try:
        if args.command == "loci":
            extract_loci(args.manifest_tsv, args.out_file, args.aliases_tsv)
        elif args.command == "neighbors":
            cutoff = (
                evalue(args.evalue_cutoff, "evalue_cutoff")
                if args.evalue_cutoff is not None
                else None
            )
            build_context(
                args.focal_tsv,
                args.loci_tsv,
                args.out_file,
                args.window,
                args.groups_tsv,
                args.hits_tsv,
                cutoff,
                args.protein_fasta,
                args.threads,
            )
        elif args.command == "validate":
            mapping = read_focals(args.focal_tsv) if args.focal_tsv else None
            rows = read_context(args.in_file, mapping)
            count = sum(row.record_type == "focal" for row in rows)
            print(
                f"LOCUSKIT context v1: {count} focal records, {len(rows) - count} neighbor records"
            )
        elif args.command == "import-legacy":
            import_legacy(args.in_file, args.out_file, args.focal_tsv)
        else:
            export_legacy(args.in_file, args.out_file)
    except (OSError, ValueError, EOFError) as error:
        print(f"locuskit: {error}", file=sys.stderr)
        return 2
    return 0
