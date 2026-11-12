import gzip
import subprocess
from pathlib import Path

import pytest

from locuskit.cli import main
from locuskit.context import COLUMNS, read_context
from locuskit.tables import write_table


def table(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


@pytest.mark.parametrize("second_id", ["pep-a", "pep-b", "unmapped"])
def test_map_proteins_selects_exact_ids_and_preserves_output_on_error(
    tmp_path: Path, second_id: str
) -> None:
    first = tmp_path / "first.faa.gz"
    with gzip.open(first, "wt") as stream:
        stream.write(">unselected other isoform\nMK\n>pep-a description\nMAG\nG*\n")
    second = table(tmp_path / "second.faa", f">{second_id}\nTT\n")
    aliases = table(tmp_path / "aliases.tsv", "protein_id\tgene_id\npep-b\tB\npep-a\tA\n")
    output = tmp_path / "mapped.faa.gz"
    with gzip.open(output, "wt") as stream:
        stream.write(">old\nAT\n")
    previous = output.read_bytes()

    status = main(
        [
            "map-proteins",
            "--protein_fasta",
            str(first),
            "--protein_fasta",
            str(second),
            "--aliases_tsv",
            str(aliases),
            "--out_file",
            str(output),
        ]
    )
    if second_id == "pep-b":
        assert status == 0
        assert gzip.decompress(output.read_bytes()) == b">A\nMAG\nG*\n>B\nTT\n"
    else:
        assert status == 2
        assert output.read_bytes() == previous


@pytest.mark.parametrize("group_rows", ["", "A\ta1\tgrp\nD\toutside\tgrp\n"])
def test_neighborhood_batch_matches_single_family_tables(tmp_path: Path, group_rows: str) -> None:
    loci = table(
        tmp_path / "loci.tsv",
        "species_id\tassembly_id\tgene_id\tseqid\tstart\tend\tstrand\n"
        "A\ta\ta1\tchr\t1\t3\t+\nA\ta\ta2\tchr\t4\t6\t-\nB\tb\tb1\tchr\t2\t7\t+\n",
    )
    manifest = table(
        tmp_path / "jobs.tsv",
        "focal_tsv\tout_file\nf1.tsv\tresult/one.tsv.gz\nf2.tsv\tresult/two.tsv\n",
    )
    table(tmp_path / "f1.tsv", "gene\tspecies\na2\tA\nmissing\tB\n")
    table(tmp_path / "f2.tsv", "gene\tspecies\nb1\tB\nother\tC\n")
    groups = tmp_path / "groups.tsv.gz"
    with gzip.open(groups, "wt") as stream:
        stream.write("species_id\tgene_id\tgroup_id\n" + group_rows)
    common = ["--loci_tsv", str(loci), "--groups_tsv", str(groups), "--window", "1"]
    batch = ["neighbors-batch", "--manifest_tsv", str(manifest), "--jobs", "2", *common]
    assert main(batch) == 0
    for index, name in enumerate(("one.tsv.gz", "two.tsv"), 1):
        single = tmp_path / name
        focal = tmp_path / f"f{index}.tsv"
        assert (
            main(["neighbors", "--focal_tsv", str(focal), "--out_file", str(single), *common]) == 0
        )
        assert single.read_bytes() == (tmp_path / "result" / name).read_bytes()

    first = tmp_path / "result" / "one.tsv.gz"
    (neighbor,) = (row for row in read_context(first) if row.record_type == "neighbor")
    assert neighbor.group_size == (2 if group_rows else None)
    assert neighbor.group_method == "provided_groups"

    previous = first.read_bytes()
    with gzip.open(groups, "wt") as stream:
        stream.write("species_id\tgene_id\tgroup_id\nA\ta1\tgrp\nA\ta1\tother\n")
    assert main(batch) == 2
    assert first.read_bytes() == previous


def test_batch_hits_preserves_family_groups_and_singletons(tmp_path: Path) -> None:
    loci = table(
        tmp_path / "loci.tsv",
        "species_id\tassembly_id\tgene_id\tseqid\tstart\tend\tstrand\n"
        "A\ta\ta1\tchr\t1\t3\t+\nA\ta\ta2\tchr\t4\t6\t+\n"
        "A\ta\ta3\tchr\t7\t9\t+\nA\ta\ta4\tchr\t10\t12\t+\n",
    )
    first = table(tmp_path / "first.tsv", "gene\tspecies\na2\tA\n")
    second = table(tmp_path / "second.tsv", "gene\tspecies\na3\tA\n")
    manifest = table(
        tmp_path / "jobs.tsv",
        "focal_tsv\tout_file\nfirst.tsv\tfirst-out.tsv\nsecond.tsv\tsecond-out.tsv\n",
    )
    hits = tmp_path / "hits.tsv.gz"
    with gzip.open(hits, "wt") as stream:
        stream.write("a1\ta3\t0.001\na3\ta1\t0.001\na2\toutside\t0.5\n")
    options = ["--loci_tsv", str(loci), "--hits_tsv", str(hits), "--window", "1"]
    assert main(["neighbors-batch", "--manifest_tsv", str(manifest), "--jobs", "2", *options]) == 0
    for focal in (first, second):
        single = tmp_path / f"{focal.stem}-single.tsv"
        assert (
            main(["neighbors", "--focal_tsv", str(focal), "--out_file", str(single), *options]) == 0
        )
        assert single.read_bytes() == (tmp_path / f"{focal.stem}-out.tsv").read_bytes()
    groups = {
        row.neighbor_id: row.group_size
        for row in read_context(tmp_path / "second-out.tsv")
        if row.record_type == "neighbor"
    }
    assert groups == {"a2": 1, "a4": None}


def test_require_protein_match_distinguishes_empty_and_unmatched_windows(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    loci = table(
        tmp_path / "loci.tsv",
        "species_id\tassembly_id\tgene_id\tseqid\tstart\tend\tstrand\n"
        "A\ta\ta1\tchr\t1\t3\t+\nA\ta\ta2\tchr\t4\t6\t+\nA\ta\ta3\tchr\t7\t9\t+\n"
        "B\tb\tb1\tchr\t1\t3\t+\nB\tb\tb2\tchr\t4\t6\t+\n",
    )
    proteins = table(tmp_path / "seq.faa", ">a1\nMAGG\n")
    table(tmp_path / "empty.tsv", "gene\tspecies\nunknown\tA\n")
    table(tmp_path / "matched.tsv", "gene\tspecies\na2\tA\n")
    unmatched = table(tmp_path / "unmatched.tsv", "gene\tspecies\nb2\tB\n")
    manifest = table(
        tmp_path / "jobs.tsv",
        "focal_tsv\tout_file\nempty.tsv\tempty-out.tsv\nmatched.tsv\tmatched-out.tsv\n"
        "unmatched.tsv\tunmatched-out.tsv\nmatched.tsv\tlater.tsv\n",
    )
    common = ["--loci_tsv", str(loci), "--protein_fasta", str(proteins), "--window", "1"]
    strict = [*common, "--require_protein_match"]
    assert main(["neighbors-batch", "--manifest_tsv", str(manifest), "--jobs", "2", *strict]) == 2
    assert "no FASTA IDs match 1 annotated neighbors" in capsys.readouterr().err
    assert read_context(tmp_path / "empty-out.tsv")[0].status == "locus_missing"
    observed = {
        row.neighbor_id: row.group_size
        for row in read_context(tmp_path / "matched-out.tsv")
        if row.record_type == "neighbor"
    }
    assert observed == {"a1": 1, "a3": None}
    assert not (tmp_path / "unmatched-out.tsv").exists()
    assert not (tmp_path / "later.tsv").exists()

    command = [
        "neighbors",
        "--focal_tsv",
        str(unmatched),
        "--out_file",
        str(tmp_path / "single.tsv"),
    ]
    assert main([*command, *common]) == 0
    previous = (tmp_path / "single.tsv").read_bytes()
    assert main([*command, *strict]) == 2
    assert (tmp_path / "single.tsv").read_bytes() == previous


@pytest.mark.parametrize("bad_later_sequence", [False, True])
@pytest.mark.parametrize("jobs", [1, 2])
def test_neighborhood_batch_keeps_protein_searches_local(
    tmp_path: Path, monkeypatch, bad_later_sequence: bool, jobs: int
) -> None:
    loci = table(
        tmp_path / "loci.tsv",
        "species_id\tassembly_id\tgene_id\tseqid\tstart\tend\tstrand\n"
        "A\ta\ta1\tchr\t1\t3\t+\nA\ta\ta2\tchr\t4\t6\t+\nA\ta\ta3\tchr\t7\t9\t+\n"
        "B\tb\tb1\tchr\t1\t3\t+\nB\tb\tb2\tchr\t4\t6\t-\nB\tb\tb3\tchr\t7\t9\t+\n",
    )
    table(tmp_path / "fa.tsv", "gene\tspecies\na2\tA\n")
    table(tmp_path / "fb.tsv", "gene\tspecies\nb2\tB\n")
    table(tmp_path / "fc.tsv", "gene\tspecies\nunknown\tA\n")
    manifest = table(
        tmp_path / "jobs.tsv", "focal_tsv\tout_file\nfa.tsv\ta.tsv\nfb.tsv\tb.tsv\nfc.tsv\tc.tsv\n"
    )
    proteins = tmp_path / "proteins.faa.gz"
    with gzip.open(proteins, "wt") as stream:
        second_sequence = "bad-1" if bad_later_sequence else "MAGG"
        stream.write(
            f">b3 description\nMAG\nG\n>a3\nMAAS\n>unused\nMR\n>b1\n{second_sequence}\n>a1\nMAAS\n"
        )

    makedb = []

    def diamond(command, *, capture_output, text, check):
        if command[1] == "makedb":
            makedb.append(Path(command[command.index("--in") + 1]).read_text())
        else:
            Path(command[command.index("--out") + 1]).touch()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(subprocess, "run", diamond)
    status = main(
        [
            "neighbors-batch",
            "--manifest_tsv",
            str(manifest),
            "--loci_tsv",
            str(loci),
            "--protein_fasta",
            str(proteins),
            "--window",
            "1",
            "--jobs",
            str(jobs),
        ]
    )
    assert (tmp_path / "a.tsv").is_file()
    if bad_later_sequence:
        assert status != 0
        assert not (tmp_path / "b.tsv").exists()
        assert not (tmp_path / "c.tsv").exists()
        assert makedb == [">a3\nMAAS\n>a1\nMAAS\n"]
    else:
        assert status == 0
        assert (tmp_path / "b.tsv").is_file()
        assert (tmp_path / "c.tsv").is_file()
        assert sorted(makedb) == [">a3\nMAAS\n>a1\nMAAS\n", ">b3\nMAG\nG\n>b1\nMAGG\n"]


@pytest.mark.parametrize(
    "conflict", ["focal.tsv", "linked.tsv", "symlink.tsv", "out.tsv", "out.tsv/nested.tsv"]
)
def test_neighborhood_batch_checks_all_outputs_before_writing(
    tmp_path: Path, conflict: str, capsys: pytest.CaptureFixture[str]
) -> None:
    focal = table(tmp_path / "focal.tsv", "gene\tspecies\na\tA\n")
    (tmp_path / "linked.tsv").hardlink_to(focal)
    (tmp_path / "symlink.tsv").symlink_to(tmp_path / "missing.tsv")
    loci = table(
        tmp_path / "loci.tsv",
        "species_id\tassembly_id\tgene_id\tseqid\tstart\tend\tstrand\nA\ta\ta\tchr\t1\t3\t+\n",
    )
    manifest = table(
        tmp_path / "jobs.tsv",
        f"focal_tsv\tout_file\nfocal.tsv\tout.tsv\nfocal.tsv\t{conflict}\n",
    )
    assert (
        main(
            [
                "neighbors-batch",
                "--manifest_tsv",
                str(manifest),
                "--loci_tsv",
                str(loci),
            ]
        )
        == 2
    )
    assert "neighborhood output" in capsys.readouterr().err.lower()
    assert not (tmp_path / "out.tsv").exists()
    assert focal.read_text() == "gene\tspecies\na\tA\n"


def test_annotation_to_context_with_aliases_and_missing_data(tmp_path: Path) -> None:
    table(
        tmp_path / "a.gff3",
        "##gff-version 3\n"
        "01\ta\tgene\t1\t10\t.\t+\t.\tID=left\n"
        "01\ta\tgene\t21\t30\t.\t-\t.\tID=focal\n"
        "01\ta\tgene\t41\t50\t.\t.\t.\tID=right\n"
        "02\ta\tgene\t1\t1\t.\t+\t.\tID=isolated\n",
    )
    manifest = table(tmp_path / "genomes.tsv", "species_id\tassembly_id\tgff\nspA\tA1\ta.gff3\n")
    aliases = table(tmp_path / "alias.tsv", "species_id\tlocus_id\tgene_id\nspA\tfocal\ttx1\n")
    loci = tmp_path / "loci.tsv.gz"
    assert (
        main(
            [
                "loci",
                "--manifest_tsv",
                str(manifest),
                "--aliases_tsv",
                str(aliases),
                "--out_file",
                str(loci),
            ]
        )
        == 0
    )
    mapping = table(
        tmp_path / "mapping.tsv",
        "gene\tspecies\ntx1\tspA\nright\tspA\nmissing\tspA\nnot_annotated\tspB\nisolated\tspA\n",
    )
    groups = table(tmp_path / "groups.tsv", "species_id\tgene_id\tgroup_id\nspA\tright\tfamily1\n")
    command = [
        "neighbors",
        "--focal_tsv",
        str(mapping),
        "--loci_tsv",
        str(loci),
        "--groups_tsv",
        str(groups),
        "--window",
        "2",
        "--out_file",
    ]
    output, second = tmp_path / "context.tsv.gz", tmp_path / "another.tsv.gz"
    assert main([*command, str(output)]) == 0
    assert main([*command, str(second)]) == 0
    assert output.read_bytes() == second.read_bytes()
    rows = read_context(output)
    focal = {r.focal_id: r for r in rows if r.record_type == "focal"}
    assert (focal["tx1"].seqid, focal["tx1"].start, focal["tx1"].end) == ("01", 21, 30)
    assert (focal["tx1"].observed_upstream, focal["tx1"].observed_downstream) == (1, 1)
    assert (focal["isolated"].observed_upstream, focal["isolated"].observed_downstream) == (0, 0)
    assert focal["missing"].status == "locus_missing"
    assert focal["not_annotated"].status == "annotation_missing"
    assert focal["right"].status == "orientation_unknown"
    assert [(r.offset, r.neighbor_id, r.group_id) for r in rows if r.record_type == "neighbor"] == [
        (-1, "right", "family1"),
        (1, "left", ""),
    ]
    assert main(["validate", "--in_file", str(output), "--focal_tsv", str(mapping)]) == 0


def test_precomputed_hits_record_singletons_and_unsearched_genes(tmp_path: Path) -> None:
    loci = table(
        tmp_path / "loci.tsv",
        "species_id\tassembly_id\tgene_id\tseqid\tstart\tend\tstrand\n"
        "A\ta\ta1\tchr\t1\t3\t+\nA\ta\ta2\tchr\t4\t6\t+\nA\ta\ta3\tchr\t7\t9\t+\n"
        "B\tb\tb1\tchr\t1\t3\t+\nB\tb\tb2\tchr\t4\t6\t+\nB\tb\tb3\tchr\t7\t9\t+\n"
        "B\tb\tb4\tchr\t10\t12\t+\n",
    )
    mapping = table(tmp_path / "mapping.tsv", "gene\tspecies\na2\tA\nb2\tB\n")
    hits = table(
        tmp_path / "hits.tsv",
        "a1\tb1\t1e-30\nb1\tb3\t1e-20\na3\ta3\t0\nb3\ta3\t0.02\nb1\toutside\t1e-30\n",
    )
    output = tmp_path / "ctx.tsv"
    assert (
        main(
            [
                "neighbors",
                "--focal_tsv",
                str(mapping),
                "--loci_tsv",
                str(loci),
                "--hits_tsv",
                str(hits),
                "--out_file",
                str(output),
            ]
        )
        == 0
    )
    neighbors = {r.neighbor_id: r for r in read_context(output) if r.record_type == "neighbor"}
    assert {g: (r.group_id, r.group_size) for g, r in neighbors.items()} == {
        "a1": ("SG000001", 3),
        "b1": ("SG000001", 3),
        "b3": ("SG000001", 3),
        "a3": ("SG000002", 1),
        "b4": ("", None),
    }
    assert all(
        r.group_method == "similarity_hits" and r.evalue_cutoff == 0.01 for r in neighbors.values()
    )


def test_legacy_round_trip_and_missing_focal(tmp_path: Path) -> None:
    content = (
        "node_name\tspecies\tdirection\toffset\tneighbor_gene\tgroup_id\tgroup_size\tevalue_cutoff\n"
        "leaf\tsp\tdownstream\t1\tNA\tSG000001\t2\t1e-05\n"
        "leaf\tsp\tupstream\t-2\tleft\tSG000001\t2\t1e-05\n"
    )
    original = table(tmp_path / "legacy.tsv", content)
    mapping = table(tmp_path / "mapping.tsv", "gene\tspecies\nleaf\tsp\nmissing\tother\n")
    canonical, projected = tmp_path / "context.tsv", tmp_path / "projected.tsv"
    assert (
        main(
            [
                "import-legacy",
                "--in_file",
                str(original),
                "--focal_tsv",
                str(mapping),
                "--out_file",
                str(canonical),
            ]
        )
        == 0
    )
    rows = read_context(canonical)
    assert len(rows) == 4
    assert all(r.status == "legacy_unknown" and r.observed_upstream is None for r in rows[:2])
    assert main(["export-legacy", "--in_file", str(canonical), "--out_file", str(projected)]) == 0
    assert projected.read_text() == content


def test_reject_conflicting_legacy_and_preserve_destination(tmp_path: Path) -> None:
    path = table(
        tmp_path / "bad.tsv",
        "node_name\tspecies\toffset\tneighbor_gene\tgroup_id\n"
        "a\tsp\t1\tx\tgroup1\na\tsp\t1\ty\tgroup2\n",
    )
    destination = table(tmp_path / "output.tsv", "original\n")
    assert main(["import-legacy", "--in_file", str(path), "--out_file", str(destination)]) == 2
    assert destination.read_text() == "original\n"


def test_reject_impossible_observed_coverage_and_in_place_write(tmp_path: Path) -> None:
    output = tmp_path / "context.tsv"
    record = {
        "schema_version": "1",
        "record_type": "focal",
        "focal_id": "a",
        "species_id": "sp",
        "assembly_id": "as",
        "seqid": "contig",
        "start": 1,
        "end": 1,
        "strand": "+",
        "status": "observed",
        "window": 10**12,
        "observed_upstream": 0,
        "observed_downstream": 10**12,
    }
    write_table(output, COLUMNS, [[record.get(field, "") for field in COLUMNS]], [])
    assert main(["validate", "--in_file", str(output)]) == 2
    legacy = table(
        tmp_path / "legacy.tsv",
        "node_name\tspecies\toffset\tneighbor_gene\tgroup_id\na\tsp\t1\tx\tgrp\n",
    )
    assert main(["import-legacy", "--in_file", str(legacy), "--out_file", str(legacy)]) == 2
    assert legacy.read_text().endswith("a\tsp\t1\tx\tgrp\n")


def test_legacy_dispatch_and_program_label(capsys) -> None:
    assert main(["fractbias", "version"]) == 0
    assert capsys.readouterr().out.startswith("kfFractBias ")


def test_diamond_uses_only_requested_proteins_and_keeps_missing_neighbors(
    tmp_path: Path, monkeypatch
) -> None:
    from locuskit.neighbors import groups_from_proteins

    protein = table(tmp_path / "seq.fa", ">A other words\nMAGG\n>B\nMAGG\n>unused\nMAAS\n")
    calls = []

    def diamond(command, *, capture_output, text, check):
        calls.append(command[1])
        if command[1] == "makedb":
            subset = Path(command[command.index("--in") + 1])
            assert subset.read_text() == ">A\nMAGG\n>B\nMAGG\n"
        else:
            assert command[command.index("--threads") + 1] == "2"
            table(Path(command[command.index("--out") + 1]), "A\tB\t0.001\n")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(subprocess, "run", diamond)
    assert groups_from_proteins(protein, {("x", "A"), ("y", "B"), ("x", "missing")}, 0.01, 2) == {
        ("x", "A"): ("SG000001", 2),
        ("y", "B"): ("SG000001", 2),
    }
    assert calls == ["makedb", "blastp"]


def test_context_checks_genomic_direction_and_locus_identity(tmp_path: Path, capsys) -> None:
    focal = {
        "schema_version": "1",
        "record_type": "focal",
        "focal_id": "a",
        "species_id": "sp",
        "assembly_id": "as",
        "seqid": "ctg",
        "start": 10,
        "end": 15,
        "strand": "+",
        "status": "observed",
        "window": 1,
        "observed_upstream": 0,
        "observed_downstream": 1,
    }
    neighbor = {
        "schema_version": "1",
        "record_type": "neighbor",
        "focal_id": "a",
        "species_id": "sp",
        "assembly_id": "as",
        "seqid": "ctg",
        "start": 2,
        "end": 5,
        "strand": "+",
        "offset": 1,
        "neighbor_id": "b",
    }
    output = tmp_path / "bad.tsv"

    def validate(*rows):
        write_table(
            output, COLUMNS, [[row.get(field, "") for field in COLUMNS] for row in rows], []
        )
        return main(["validate", "--in_file", str(output)])

    assert validate(focal, neighbor) == 2
    assert "offsets contradict" in capsys.readouterr().err
    other = dict(focal, focal_id="b", start=21, end=25, observed_downstream=0)
    assert validate(focal, dict(neighbor, start=20, end=25), other) == 2
    assert "Conflicting coordinates" in capsys.readouterr().err


def test_cli_reports_truncated_gzip(tmp_path: Path, capsys) -> None:
    import gzip

    output = tmp_path / "truncated.tsv.gz"
    output.write_bytes(gzip.compress(("\t".join(COLUMNS) + "\n").encode())[:-4])
    assert main(["validate", "--in_file", str(output)]) == 2
    assert capsys.readouterr().err.startswith("locuskit:")


def test_loci_rejects_multi_id_gene_and_preserves_existing_output(tmp_path: Path, capsys) -> None:
    table(tmp_path / "one.gff", "chr\t.\tgene\t1\t2\t.\t+\t.\tID=gene1\n")
    table(tmp_path / "two.gff", "chr\t.\tgene\t5\t7\t.\t+\t.\tID=gene2,gene3\n")
    manifest = table(
        tmp_path / "genomes.tsv",
        "species_id\tassembly_id\tgff\nA\tx\tone.gff\nB\ty\ttwo.gff\n",
    )
    output = table(tmp_path / "loci.tsv", "prior output\n")
    assert main(["loci", "--manifest_tsv", str(manifest), "--out_file", str(output)]) == 2
    assert "multiple gene loci" in capsys.readouterr().err
    assert output.read_text() == "prior output\n"
