"""NumPy Markov scoring is bit-identical to the reference loop."""
import random

from pytransdecoder.core.markov_fast import compile_table
from pytransdecoder.core.parallel import score_markov
from pytransdecoder.core.sequence import reverse_complement
from pytransdecoder.core.translator import Translator


def _scores(rng, frac=0.7):
    out = {}
    for k in range(1, 7):
        for _ in range(4 ** k if k < 6 else 4096):
            kmer = "".join(rng.choice("ACGT") for _ in range(k))
            for fr in range(3):
                if rng.random() < frac:
                    out[f"{kmer}-{fr}"] = rng.choice([rng.uniform(-3, 3), -0.0, 0.0, round(rng.uniform(-3, 3), 3)])
    return out


def _seq(rng, n, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


def _same(a, b):
    return repr(a) == repr(b)


def test_bit_identical_random():
    for code in ("Standard", "Tetrahymena"):
        stops = frozenset(Translator(code).get_stop_codons())
        rng = random.Random(5)
        sc = _scores(rng)
        t = compile_table(sc, stops)
        assert t is not None
        seqs = [_seq(rng, n) for n in range(0, 41) for _ in range(5)]
        seqs += [_seq(rng, rng.randint(6, 3000)) for _ in range(300)]
        seqs += [_seq(rng, rng.randint(6, 400), "ACGTNNRY") for _ in range(200)]
        seqs += [_seq(rng, 3 * rng.randint(2, 100)) + s for s in sorted(stops) for _ in range(20)]
        seqs += [_seq(rng, 3 * rng.randint(2, 100) + k) + s for s in sorted(stops) for k in (1, 2) for _ in range(10)]
        for s in seqs:
            r = reverse_complement(s)
            ref = tuple(score_markov(x, 5, sc, stops) for x in (s, s[1:], s[2:], r, r[1:], r[2:]))
            got = t.six_frames(s, r)
            assert all(_same(a, b) for a, b in zip(ref, got)), (s, ref, got)
            assert _same(t.score(s), ref[0])


def test_uncompilable_table():
    assert compile_table({"ACGTNA-0": 1.0}, frozenset()) is None
    assert compile_table({"FRAME-0": 1.0}, frozenset()) is None
