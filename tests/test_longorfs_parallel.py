"""LongOrfs --cpu: output files identical to the single-process run."""

import random
from pathlib import Path

import pytest

from pytransdecoder.longorfs import run_longorfs

FILES = ["longest_orfs.gff3", "longest_orfs.cds", "longest_orfs.pep", "base_freqs.dat",
         "longorfs.params.json"]


def _random_fasta(path: Path, ids, seed=7):
    rng = random.Random(seed)
    with open(path, "w") as fh:
        for i, acc in enumerate(ids):
            n = rng.randint(100, 3000)
            # coding-like stretches (no in-frame stop) mixed with random sequence, some IUPAC / lowercase
            seq = "".join(rng.choice("ACGT") for _ in range(n))
            if i % 3 == 0:
                body = "".join(rng.choice(["GCT", "AAA", "CTG", "GAG", "TCC", "ATG"]) for _ in range(150))
                seq = seq[:50] + "ATG" + body + "TAA" + seq[50:]
            if i % 5 == 0:
                seq = seq[:30] + "N" + seq[31:]
            if i % 7 == 0:
                seq = seq.lower()
            fh.write(f">{acc} desc {i}\n")
            for k in range(0, len(seq), 70):
                fh.write(seq[k:k + 70] + "\n")


def _run(tmp_path, name, fasta, cpu, **kw):
    out = tmp_path / name
    out.mkdir()
    run_longorfs(fasta, output_dir=out, cpu=cpu, **kw)
    d = out / f"{fasta.name}.transdecoder_dir"
    return {f: (d / f).read_bytes() for f in FILES}


@pytest.mark.parametrize("strand_specific", [False, True])
@pytest.mark.parametrize("genetic_code", ["universal", "Tetrahymena"])
def test_cpu_parity(tmp_path, strand_specific, genetic_code):
    ids = [f"tx{i}" for i in range(120)] + ["A", "A.p1", "A.p1.p", "B.p2"]  # IDs that look like model IDs
    fasta = tmp_path / "t.fa"
    _random_fasta(fasta, ids)
    kw = dict(strand_specific=strand_specific, genetic_code=genetic_code, min_protein_length=30)
    one = _run(tmp_path, "c1", fasta, 1, **kw)
    assert one["longest_orfs.pep"].count(b">") > 50
    import pytransdecoder.longorfs as lo
    lo.SERIAL_BATCH, keep = 7, lo.SERIAL_BATCH
    try:
        assert _run(tmp_path, "c1b", fasta, 1, **kw) == one
    finally:
        lo.SERIAL_BATCH = keep
    for cpu in (2, 3, 8):
        assert _run(tmp_path, f"c{cpu}", fasta, cpu, **kw) == one


@pytest.mark.parametrize("batch", [2000, 3])
def test_repeated_accessions(tmp_path, monkeypatch, batch):
    """A repeated accession continues the .pN numbering of its earlier copy, at any --cpu and
    across the batches of the streamed single-process run."""
    import pytransdecoder.longorfs as lo
    monkeypatch.setattr(lo, "SERIAL_BATCH", batch)
    ids = [f"tx{i}" for i in range(40)] + ["tx3", "tx0", "tx3"]
    fasta = tmp_path / "t.fa"
    _random_fasta(fasta, ids)
    one = _run(tmp_path, "c1", fasta, 1, min_protein_length=30)
    assert _run(tmp_path, "c4", fasta, 4, min_protein_length=30) == one
    pep = one["longest_orfs.pep"].decode()
    tx3 = [l.split()[0] for l in pep.splitlines() if l.startswith(">tx3.p")]
    assert len(tx3) == len(set(tx3)) > 1


def test_gene_trans_map_and_complete_only(tmp_path):
    ids = [f"tx{i}" for i in range(60)]
    fasta = tmp_path / "t.fa"
    _random_fasta(fasta, ids)
    gtm = tmp_path / "gtm.tsv"
    gtm.write_text("".join(f"gene{i // 2}\ttx{i}\n" for i in range(0, 60, 3)))
    for kw in (dict(gene_trans_map_file=gtm), dict(complete_orfs_only=True)):
        kw["min_protein_length"] = 30
        one = _run(tmp_path, f"c1_{len(kw)}_{'g' if 'gene_trans_map_file' in kw else 'c'}", fasta, 1, **kw)
        assert _run(tmp_path, f"c3_{'g' if 'gene_trans_map_file' in kw else 'c'}", fasta, 3, **kw) == one
