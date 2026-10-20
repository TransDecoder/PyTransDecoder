"""FASTA reading: text before the first header is skipped (GitHub issue: FHR ';' comment lines)."""

import logging

from pytransdecoder.core.fasta import read_fasta
from pytransdecoder.longorfs import run_longorfs


RECORDS = ">tx1 first transcript\nATGAAA\nCCC TGA\n>tx2\nATGGGG\nTAA\n"


def _read(path):
    return [(r.id, r.description, str(r.seq)) for r in read_fasta(path)]


def test_read_fasta_records(tmp_path):
    fa = tmp_path / "plain.fa"
    fa.write_text(RECORDS)
    assert _read(fa) == [("tx1", "tx1 first transcript", "ATGAAACCCTGA"), ("tx2", "tx2", "ATGGGGTAA")]


def test_text_before_first_header_is_skipped(tmp_path, caplog):
    plain = tmp_path / "plain.fa"
    plain.write_text(RECORDS)
    commented = tmp_path / "commented.fa"
    commented.write_text(";FHR provenance: ATGAAACCCGGGTTTAAA\n;more ATGCCCAAAGGG\n\n" + RECORDS)
    with caplog.at_level(logging.WARNING):
        assert _read(commented) == _read(plain)
    assert "skipped text before the first '>' header" in caplog.text

    blank = tmp_path / "blank.fa"
    blank.write_text("\n\n" + RECORDS)
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        assert _read(blank) == _read(plain)
    assert "skipped" not in caplog.text


def test_longorfs_ignores_comment_block(tmp_path):
    orf = "ATG" + "GCT" * 120 + "TAA"
    body = f">tx1\nCC{orf}GG\n"
    comment = ";" + orf + "\n"          # a comment that would hold an ORF if read as sequence
    outputs = []
    for name, text in (("plain", body), ("commented", comment + body)):
        fa = tmp_path / f"{name}.fa"
        fa.write_text(text)
        run_longorfs(fa, output_dir=tmp_path / f"{name}_out", min_protein_length=50)
        workdir = tmp_path / f"{name}_out" / f"{name}.fa.transdecoder_dir"
        outputs.append({f: (workdir / f).read_text()
                        for f in ("longest_orfs.pep", "longest_orfs.gff3", "base_freqs.dat")})
    assert outputs[0]["longest_orfs.pep"].count(">") >= 1
    assert outputs[0] == outputs[1]
