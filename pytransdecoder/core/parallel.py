"""
Multi-process scoring of candidate and decoy ORFs.

Both loops are CPU-bound pure Python, so they use worker processes (threads would be serialised by the
GIL). Each worker runs the same per-ORF arithmetic as the single-process code. Work is split into chunks
in input order, and results are put back together in that order, so the output does not depend on the
number of workers.
"""

from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from typing import Dict, FrozenSet, List, Sequence, Tuple

from pytransdecoder.core.sequence import reverse_complement

MARKOV_ORDER = 5

# per-worker state, set by _init_worker (inherited on fork, or passed on spawn)
_HEX: Dict[str, float] = {}
_STOPS: FrozenSet[str] = frozenset()
_FINDER_ARGS: dict = {}
_TABLE = None  # core/markov_fast.MarkovTable of _HEX, or None (reference loop)


def score_markov(sequence: str, markov_order: int, scores: Dict[str, float], stop_codons) -> float:
    """Sum of framed k-mer log-likelihoods; a terminal in-frame stop codon is not scored."""
    seq_length = len(sequence)
    if seq_length < markov_order + 1:
        return 0.0
    total_score = 0.0
    for i in range(seq_length):
        frame = i % 3
        markov_use = min(i, markov_order)
        kmer = sequence[i - markov_use:i + 1]
        if i == seq_length - 3 and frame == 0:
            if sequence[i:i + 3] in stop_codons:
                break
        total_score += scores.get(f"{kmer}-{frame}", 0.0)
    return total_score


def _init_worker(hex_scores, stop_codons, finder_args):
    global _HEX, _STOPS, _FINDER_ARGS, _TABLE
    from pytransdecoder.core.markov_fast import compile_table
    _HEX, _STOPS, _FINDER_ARGS = hex_scores, stop_codons, finder_args
    _TABLE = compile_table(hex_scores, stop_codons)


def _score_six_frames(sequence: str) -> Tuple[float, ...]:
    rev = reverse_complement(sequence)
    if _TABLE is not None:
        return _TABLE.six_frames(sequence, rev)
    return tuple(score_markov(s, MARKOV_ORDER, _HEX, _STOPS)
                 for s in (sequence, sequence[1:], sequence[2:], rev, rev[1:], rev[2:]))


def _six_frame_chunk(seqs: Sequence[str]) -> List[Tuple[float, ...]]:
    return [_score_six_frames(s) for s in seqs]


def _decoy_chunk(args):
    """Decoy ORFs of a chunk of (tid, seq) in the given order: (n_decoy, [(stratum, contrast), ...])."""
    from pytransdecoder.core import decoys
    from pytransdecoder.core.prng import Xoshiro128ss
    from pytransdecoder.core.orf_finder import ORFFinder
    items, n_shuf, seed = args
    finder = ORFFinder(**_FINDER_ARGS)
    n_decoy, out = 0, []
    for tid, seq in items:
        for k in range(n_shuf):
            shuf = decoys.dinucleotide_shuffle(seq, Xoshiro128ss(f"{seed}:{tid}:{k}"))
            for orf in finder.find_all_orfs(shuf, f"decoy{k}_{tid}"):
                n_decoy += 1
                cds = orf.sequence.upper()
                if _TABLE is not None:
                    s0, s1, s2 = (round(x, 2) for x in _TABLE.forward_frames(cds))
                else:
                    s0 = round(score_markov(cds, MARKOV_ORDER, _HEX, _STOPS), 2)
                    s1 = round(score_markov(cds[1:], MARKOV_ORDER, _HEX, _STOPS), 2)
                    s2 = round(score_markov(cds[2:], MARKOV_ORDER, _HEX, _STOPS), 2)
                c = s0 - max(s1, s2)
                if s0 > 0 and c > 0:
                    out.append((decoys.stratum(len(cds) // 3, len(seq), orf.orf_type == 'complete',
                                               'len_tx_type'), c))
    return n_decoy, out


def _chunks_by_size(items, sizes, n_chunks):
    """Consecutive chunks of roughly equal total size (order preserved)."""
    total = sum(sizes) or 1
    target = total / n_chunks
    chunks, cur, acc = [], [], 0
    for it, sz in zip(items, sizes):
        cur.append(it)
        acc += sz
        if acc >= target:
            chunks.append(cur)
            cur, acc = [], 0
    if cur:
        chunks.append(cur)
    return chunks


def _pool(cpu, hex_scores, stop_codons, finder_args=None):
    ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else None
    return ProcessPoolExecutor(max_workers=cpu, mp_context=ctx, initializer=_init_worker,
                               initargs=(hex_scores, stop_codons, finder_args or {}))


def score_six_frames(seqs: List[str], hex_scores, stop_codons, cpu: int) -> List[Tuple[float, ...]]:
    """Six-frame scores of each sequence, in input order."""
    if cpu <= 1 or len(seqs) < 2:
        _init_worker(hex_scores, stop_codons, {})
        return _six_frame_chunk(seqs)
    # many small chunks so slow chunks (long ORFs) do not leave workers idle
    chunks = _chunks_by_size(seqs, [len(s) for s in seqs], cpu * 16)
    with _pool(cpu, hex_scores, stop_codons) as ex:
        return [r for part in ex.map(_six_frame_chunk, chunks) for r in part]


def decoy_scores(transcripts: List[Tuple[str, str]], n_shuf: int, seed: int, hex_scores, stop_codons,
                 finder_args: dict, cpu: int):
    """(n_decoy, eligible decoy (stratum, contrast) list) over transcripts in the given order."""
    if cpu <= 1:
        _init_worker(hex_scores, stop_codons, finder_args)
        return _decoy_chunk((transcripts, n_shuf, seed))
    chunks = _chunks_by_size(transcripts, [len(s) for _, s in transcripts], cpu * 16)
    n_decoy, out = 0, []
    with _pool(cpu, hex_scores, stop_codons, finder_args) as ex:
        for n, part in ex.map(_decoy_chunk, [(c, n_shuf, seed) for c in chunks]):
            n_decoy += n
            out.extend(part)
    return n_decoy, out
