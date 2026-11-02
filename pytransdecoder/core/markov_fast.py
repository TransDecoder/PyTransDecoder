"""
Framed 5th-order Markov scoring with NumPy, bit-identical to the reference loop
(core/parallel.score_markov).

The reference adds, for i = 0 .. n-1, the score of the framed k-mer ending at i
(k = min(i, 5) + 1, frame i % 3) to a running total that starts at 0.0, and stops
before a terminal in-frame stop codon. Here the per-position values are gathered
with array lookups, and the total is formed by np.add.accumulate over
[0.0, v0, v1, ...], which performs the same IEEE additions in the same order.
Positions 0-4 (k-mers shorter than 6) use the dict, as in the reference.

Only score tables whose keys are all ^[ACGT]{1,6}-[012]$ can be compiled; for any
other table, compile() returns None and callers use the reference loop.
"""

import re
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

_KEY = re.compile(r"^[ACGT]{1,6}-[012]$")
_CODE = np.full(256, 4, dtype=np.int64)
for _i, _b in enumerate("ACGT"):
    _CODE[ord(_b)] = _i
_POW = (1024, 256, 64, 16, 4, 1)
INVALID = 4096  # index of the zero slot: 6-mers with a non-ACGT base


class MarkovTable:
    """Hexamer scores of one model, as a (3, 4097) array plus the dict for short k-mers."""

    def __init__(self, scores: Dict[str, float], stop_codons):
        self.scores = scores
        self.stop_codons = stop_codons
        tbl = np.zeros((3, INVALID + 1), dtype=np.float64)
        for key, v in scores.items():
            kmer, fr = key[:-2], int(key[-1])
            if len(kmer) == 6:
                idx = 0
                for b, p in zip(kmer, _POW):
                    idx += "ACGT".index(b) * p
                tbl[fr, idx] = v
        self.tbl = tbl

    # -- helpers -------------------------------------------------------------------
    @staticmethod
    def _hexamer_index(seq: str) -> np.ndarray:
        """idx[t] = base-4 index of seq[t:t+6], or INVALID if it contains a non-ACGT base."""
        c = _CODE[np.frombuffer(seq.encode("latin-1"), dtype=np.uint8)]
        n = len(c)
        bad = c == 4
        cc = np.where(bad, 0, c)
        m = n - 5
        idx = cc[0:m] * 1024
        badw = bad[0:m].copy()
        for j in range(1, 6):
            idx += cc[j:m + j] * _POW[j]
            badw |= bad[j:m + j]
        idx[badw] = INVALID
        return idx

    def _score_offsets(self, seq: str, idx: Optional[np.ndarray], offsets: Sequence[int]) -> Tuple[float, ...]:
        """Reference score of seq[k:] for each k in offsets (idx: _hexamer_index(seq) or None)."""
        n_full = len(seq)
        get = self.scores.get
        tail_stop = seq[n_full - 3:] in self.stop_codons
        out = []
        for k in offsets:
            n = n_full - k
            if n < 6:
                out.append(0.0)
                continue
            m = n - 3 if (n % 3 == 0 and tail_stop) else n  # positions scored: 0 .. m-1
            ms = min(m, 5)
            vals = np.empty(m + 1, dtype=np.float64)
            vals[0] = 0.0
            for i in range(ms):
                vals[i + 1] = get(f"{seq[k:k + i + 1]}-{i % 3}", 0.0)
            if m > 5:
                # positions i = 5 .. m-1 of seq[k:], i.e. 6-mers seq[k+i-5 : k+i+1] = idx[k+i-5]
                t = idx[k:k + m - 5]
                fr = (np.arange(5, m)) % 3
                vals[6:] = self.tbl[fr, t]
            out.append(float(np.add.accumulate(vals)[-1]))
        return tuple(out)

    # -- public --------------------------------------------------------------------
    def score(self, seq: str) -> float:
        if len(seq) < 6:
            return 0.0
        return self._score_offsets(seq, self._hexamer_index(seq), (0,))[0]

    def forward_frames(self, seq: str) -> Tuple[float, float, float]:
        """Scores of seq, seq[1:], seq[2:]."""
        idx = self._hexamer_index(seq) if len(seq) >= 6 else None
        return self._score_offsets(seq, idx, (0, 1, 2))

    def six_frames(self, seq: str, rev: str) -> Tuple[float, ...]:
        """Scores of seq, seq[1:], seq[2:], rev, rev[1:], rev[2:] (rev = reverse complement)."""
        return self.forward_frames(seq) + self.forward_frames(rev)


def compile_table(scores: Dict[str, float], stop_codons) -> Optional[MarkovTable]:
    """MarkovTable for a score dict, or None if any key is outside ^[ACGT]{1,6}-[012]$."""
    if not all(_KEY.match(k) for k in scores):
        return None
    return MarkovTable(scores, stop_codons)
