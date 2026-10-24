import gzip
import random
from collections import Counter

from pytransdecoder.core import decoys, presets
from pytransdecoder.core.prng import Xoshiro128ss, fnv1a32


def test_prng_golden_values():
    """Fixed outputs, shared with the Perl implementation (t/decoys.t)."""
    assert fnv1a32("") == 0x811C9DC5
    assert fnv1a32("a") == 0xE40C292C
    g = Xoshiro128ss("42:x:0")
    assert [g.next32() for _ in range(5)] == [330423170, 325826645, 3944918289, 989894671, 3051725679]
    assert decoys.dinucleotide_shuffle("ATGCGCGTATTAGCCGATAGGCTAGCTAACG", Xoshiro128ss("42:x:0")) == \
        "ATTATAATAGCGACGCCTAGGCGCGTGCTAG"


def test_shuffle_keeps_dinucleotides_and_ends():
    r = random.Random(1)
    for i in range(300):
        s = "".join(r.choice("ACGTN") for _ in range(r.randint(3, 300)))
        t = decoys.dinucleotide_shuffle(s, Xoshiro128ss(f"7:{i}"))
        assert Counter(zip(s, s[1:])) == Counter(zip(t, t[1:])) and s[0] == t[0] and s[-1] == t[-1]


def test_strata_bins():
    assert decoys.stratum(49, 400, True) == (-1, 0, True)
    assert decoys.stratum(50, 500, False) == (0, 1, False)
    assert decoys.stratum(56, 999, True) == (1, 1, True)
    assert decoys.stratum(100, 1000, True)[1] == 2
    assert decoys.stratum(100, 4000, True)[1] == 4


def test_qvalues_running_minimum():
    targets = {f"t{i}": ("s", float(c)) for i, c in enumerate([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])}
    q = decoys.qvalues(targets, [("s", 1.5), ("s", 2.5)], 1)
    assert q["t9"] == 0.0 and q["t2"] == 0.0           # no decoys at or above 3
    assert q["t1"] == 1 / 9 and q["t0"] == 2 / 10


def test_preset_roundtrip(tmp_path):
    run = tmp_path / "run"
    (run / "__checkpoints_TDpredict").mkdir(parents=True)
    (run / "hexamer.scores").write_text("AAAAAA\t1\t1\t0.5\n")
    (run / "__checkpoints_TDpredict" / "top_training_orfs.cds").write_text(">a\nATGAAATAA\n")
    (run / "decoy_qvalues.tsv").write_text("#orf_id\tstratum\tcontrast\tq\n"
                                           "o1\t3,2,True\t10.00\t0.001\no2\t3,2,True\t4.00\t0.04\n"
                                           "o3\t3,2,True\t2.00\t0.3\n")
    out = presets.build_preset(run, tmp_path / "p", "test", 45.0, 726, "Universal", "len_tx_type")
    assert (out / "thresholds.tsv.gz").exists() and (out / "hexamer.scores.gz").exists()
    p = presets.Preset(presets.resolve(str(out)))
    assert p.threshold((3, 2, True), 0.05) == 4.0
    assert p.threshold((3, 2, True), 0.01) == 10.0
    assert p.threshold((5, 2, True), 0.5) == 2.0        # nearest ORF-length bin
    assert p.auto_accept_nt == 726
    with gzip.open(p.hexamer_scores, "rt") as fh:
        assert fh.read().startswith("AAAAAA")
