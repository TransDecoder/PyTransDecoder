"""
Decoy ORFs from shuffled transcripts, for target-decoy FDR estimates of candidate ORFs (experimental).

Each transcript is shuffled preserving its dinucleotide counts (Altschul & Erickson 1985; the random
Eulerian path is drawn with Wilson's algorithm for the last-exit tree), so decoys keep the transcript's
length and base composition. Random numbers come from core/prng.py, so the Perl implementation
(PerlLib/Decoys.pm) produces the same decoys. Decoy ORFs are found with the same ORF finder, genetic code, minimum
length and strand setting as LongOrfs and scored with the run's trained Markov model.

Within strata of ORF length (log-spaced bins), transcript length (log2 bins) and ORF type (complete vs
partial), a target ORF with frame contrast c has
    FDR(c) = decoys with contrast >= c / (shuffles x targets with contrast >= c)
and its q-value is the lowest FDR of any cutoff that accepts it. Strata without decoys give q = 0.
"""
import bisect
import math
from collections import defaultdict
from typing import Dict, Hashable, Iterable, List, Tuple

from pytransdecoder.core.prng import Xoshiro128ss

LENGTH_BIN_RATIO = 1.1
LENGTH_BIN_START = 50  # codons
TX_BIN_START = 500  # nt
STRATA = ("len", "len_tx", "len_tx_type")


def dinucleotide_shuffle(seq: str, rng: Xoshiro128ss) -> str:
    """Random sequence with the same first base, last base and dinucleotide counts as seq."""
    n = len(seq)
    if n < 3:
        return seq
    edges = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        edges[a].append(b)
    end = seq[-1]
    # random last-exit tree rooted at the final base (Wilson's loop-erased random walks)
    in_tree, nxt, last = {end}, {}, {}
    for v in sorted(edges):
        u = v
        while u not in in_tree:
            nxt[u] = edges[u][rng.below(len(edges[u]))]
            u = nxt[u]
        u = v
        while u not in in_tree:
            in_tree.add(u)
            last[u] = nxt[u]
            u = nxt[u]
    order = {}
    for v in sorted(edges):
        out = edges[v][:]
        if v in last:
            out.remove(last[v])
        rng.shuffle(out)
        if v in last:
            out.append(last[v])
        order[v] = out
    pos = defaultdict(int)
    res, u = [seq[0]], seq[0]
    for _ in range(n - 1):
        w = order[u][pos[u]]
        pos[u] += 1
        res.append(w)
        u = w
    return "".join(res)


def _length_bin(codons: int) -> int:
    """Largest k with codons >= LENGTH_BIN_START * LENGTH_BIN_RATIO**k (powers by repeated
    multiplication, so the Perl port bins identically); -1 below LENGTH_BIN_START."""
    if codons < LENGTH_BIN_START:
        return -1
    k, edge = 0, LENGTH_BIN_START * LENGTH_BIN_RATIO
    while codons >= edge:
        k += 1
        edge *= LENGTH_BIN_RATIO
    return k


def stratum(codons: int, tx_len: int, complete: bool, strata: str = "len_tx_type") -> Tuple:
    """Stratum key: ORF-length bin, then (len_tx) transcript-length bin (floor(log2(tx_len / 500)) + 1,
    0 below 500 nt), then (len_tx_type) ORF type."""
    k = (_length_bin(codons),)
    if strata in ("len_tx", "len_tx_type"):
        k += ((tx_len // TX_BIN_START).bit_length() if tx_len >= TX_BIN_START else 0,)
    if strata == "len_tx_type":
        k += (bool(complete),)
    return k


def qvalues(targets: Dict[str, Tuple[Hashable, float]], decoys: Iterable[Tuple[Hashable, float]],
            n_shuffles: int) -> Dict[str, float]:
    """q-value per target ORF id from (stratum, contrast) of targets and decoys (eligible ORFs only)."""
    tby, dby = defaultdict(list), defaultdict(list)
    for orf_id, (k, c) in targets.items():
        tby[k].append(c)
    for k, c in decoys:
        dby[k].append(c)
    q_of = {}
    for k, tc in tby.items():
        tc = sorted(tc)
        dc = sorted(dby.get(k, ()))
        q, run_min = [], math.inf
        for i, c in enumerate(tc):
            fdr = (len(dc) - bisect.bisect_left(dc, c)) / n_shuffles / (len(tc) - i)
            run_min = min(run_min, fdr)
            q.append(run_min)
        q_of[k] = (tc, q)
    out = {}
    for orf_id, (k, c) in targets.items():
        tc, q = q_of[k]
        out[orf_id] = q[bisect.bisect_left(tc, c)]
    return out
