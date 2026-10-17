"""
Base frequency calculation

Computes nucleotide frequencies for background model.
Replaces compute_base_probs.pl
"""

from typing import Dict
from Bio import SeqIO


def compute_base_counts(fasta_file: str, top_strand_only: bool = False) -> Dict[str, int]:
    """
    Count nucleotides in a FASTA file, as compute_base_probs.pl does.

    Unless top_strand_only, each base also increments its complement so that
    the counts reflect both strands (A=T, C=G).
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    counts = {base: 0 for base in 'ACGT'}
    for record in SeqIO.parse(fasta_file, "fasta"):
        for base in str(record.seq).upper():
            if base in counts:
                counts[base] += 1
                if not top_strand_only:
                    counts[complement[base]] += 1
    return counts


def compute_base_frequencies(fasta_file: str, top_strand_only: bool = False) -> Dict[str, float]:
    """Nucleotide frequencies from a FASTA file (see compute_base_counts)."""
    counts = compute_base_counts(fasta_file, top_strand_only)
    total = sum(counts.values())
    if total == 0:
        return {base: 0.25 for base in 'ACGT'}
    return {base: counts[base] / total for base in 'ACGT'}


def write_base_frequencies(counts: Dict[str, int], output_file: str):
    """
    Write base_freqs.dat in TransDecoder's format: base, count, ratio (3 decimals).

    The rounded ratio is what both implementations read back as the background
    model, so it is written exactly as compute_base_probs.pl does.
    """
    total = sum(counts.values())
    with open(output_file, 'w') as f:
        for base in 'ACGT':
            ratio = counts[base] / total if total else 0.25
            f.write(f"{base}\t{counts[base]}\t{ratio:.3f}\n")


def read_base_frequencies(input_file: str) -> Dict[str, float]:
    """Read background base frequencies from base_freqs.dat (base, count, ratio)."""
    frequencies = {}
    with open(input_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 3:
                frequencies[parts[0]] = float(parts[2])
            elif len(parts) == 2:  # files written by PyTransDecoder <= 0.1.0
                frequencies[parts[0]] = float(parts[1])
    missing = set('ACGT') - set(frequencies)
    if missing:
        raise ValueError(f"{input_file}: missing base frequencies for {sorted(missing)}")
    return frequencies


def calculate_gc_content(frequencies: Dict[str, float]) -> float:
    """Calculate GC content from base frequencies"""
    return frequencies.get('G', 0) + frequencies.get('C', 0)
