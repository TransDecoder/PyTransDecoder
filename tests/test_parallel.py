"""Scoring with worker processes gives exactly the single-process results, in input order."""
import random

from pytransdecoder.core import parallel
from pytransdecoder.core.translator import Translator

STOPS = frozenset(Translator("Standard").get_stop_codons())


def _hex_scores(rng):
    scores = {}
    for _ in range(3000):
        k = rng.randint(1, 6)
        kmer = "".join(rng.choice("ACGT") for _ in range(k))
        scores[f"{kmer}-{rng.randint(0, 2)}"] = rng.uniform(-2, 2)
    return scores


def _seqs(rng, n):
    return ["".join(rng.choice("ACGT") for _ in range(rng.randint(3, 900))) for _ in range(n)]


def test_six_frame_scores_independent_of_cpu():
    rng = random.Random(7)
    hexs, seqs = _hex_scores(rng), _seqs(rng, 300)
    one = parallel.score_six_frames(seqs, hexs, STOPS, 1)
    assert len(one) == len(seqs)
    for cpu in (2, 3):
        assert parallel.score_six_frames(seqs, hexs, STOPS, cpu) == one


def test_score_markov_matches_single_call():
    rng = random.Random(3)
    hexs, seqs = _hex_scores(rng), _seqs(rng, 20)
    six = parallel.score_six_frames(seqs, hexs, STOPS, 2)
    for s, sc in zip(seqs, six):
        assert sc[0] == parallel.score_markov(s, 5, hexs, STOPS)
        assert sc[2] == parallel.score_markov(s[2:], 5, hexs, STOPS)


def test_decoy_scores_independent_of_cpu():
    rng = random.Random(11)
    hexs = _hex_scores(rng)
    tx = [(f"t{i}", s) for i, s in enumerate(_seqs(rng, 120))]
    finder_args = dict(min_protein_length=20, allow_5prime_partial=True, allow_3prime_partial=True,
                       forward_strand=True, reverse_strand=True, genetic_code="Standard",
                       complete_orfs_only=False)
    one = parallel.decoy_scores(tx, 2, 42, hexs, STOPS, finder_args, 1)
    assert one[0] > 0
    for cpu in (2, 4):
        assert parallel.decoy_scores(tx, 2, 42, hexs, STOPS, finder_args, cpu) == one
