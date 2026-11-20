import pytest

from locuskit.cli import main
from locuskit.context import read_context
from locuskit.replicons import locus_span


def test_physical_copies_wrap_and_validate(tmp_path):
    gff = tmp_path / "a.gff"
    gff.write_text(
        "".join(
            f"c\t.\tgene\t{a}\t{b}\t.\t{s}\t.\tID={g}\n"
            for g, a, b, s in [
                ("rpl2a", 90, 103, "+"),
                ("trnH", 4, 7, "-"),
                ("psbA", 20, 30, "+"),
                ("rpl2b", 43, 55, "-"),
            ]
        )
    )
    (tmp_path / "manifest").write_text("species_id\tassembly_id\tgff\ns\ta\ta.gff\n")
    rep = tmp_path / "rep"
    rep.write_text("species_id\tassembly_id\tseqid\tlength\ttopology\ns\ta\tc\t100\tcircular\n")
    locus = tmp_path / "locus"
    args = [
        "loci",
        "--manifest_tsv",
        str(tmp_path / "manifest"),
        "--replicons_tsv",
        str(rep),
        "--out_file",
        str(locus),
    ]
    assert main(args) == 0
    assert "90-100,1-3" in locus.read_text()
    (tmp_path / "focal").write_text("gene\tspecies\ntrnH\ts\n")
    ctx = tmp_path / "ctx"
    assert (
        main(
            [
                "neighbors",
                "--focal_tsv",
                str(tmp_path / "focal"),
                "--loci_tsv",
                str(locus),
                "--out_file",
                str(ctx),
                "--window",
                "10",
            ]
        )
        == 0
    )
    rows = read_context(ctx)
    assert [(r.offset, r.neighbor_id) for r in rows[1:]] == [
        (-1, "psbA"),
        (1, "rpl2a"),
        (2, "rpl2b"),
    ]
    assert (rows[0].observed_upstream, rows[0].observed_downstream) == (1, 2)
    assert all(r.circular_length == 100 for r in rows)
    before = locus.read_bytes()
    rep.write_text(rep.read_text().replace("circular", "linear"))
    assert main(args) == 2 and locus.read_bytes() == before


def test_span_parts_follow_strand_and_reject_extra_turn():
    assert locus_span({(90, 100), (1, 5), (9, 11)}, 100, "circular", "-") == (
        90,
        111,
        "9-11,1-5,90-100",
    )
    assert locus_span({(30, 35), (42, 45)}, 100, "linear", "+") == (30, 45, "30-35,42-45")
    with pytest.raises(ValueError):
        locus_span({(80, 180)}, 100, "circular", "+")


def test_multipart_cds_annotation_keeps_physical_segments(tmp_path):
    from kffractbias.io import annotation_locus_order

    gff = tmp_path / "split.gff"
    gff.write_text(
        "c\t.\tgene\t1\t100\t.\t+\t.\tID=A\n"
        "c\t.\tCDS\t1\t6\t.\t+\t0\tID=p;Parent=A\n"
        "c\t.\tCDS\t94\t100\t.\t+\t0\tID=p;Parent=A\n"
    )
    parts = {}
    genes = annotation_locus_order(gff, feature="CDS", attribute="ID", intervals_by_locus=parts)
    assert [g.gene_id for g in genes] == ["A"]
    assert parts == {"A": {(1, 6), (94, 100)}}
