"""
Base frequency calculation

Computes nucleotide frequencies for background model.
Replaces compute_base_probs.pl
"""

from typing import Dict
from collections import defaultdict
from Bio import SeqIO


def compute_base_frequencies(
    fasta_file: str,
    top_strand_only: bool = False
) -> Dict[str, float]:
    """
    Compute nucleotide base frequencies from FASTA file.
    
    Args:
        fasta_file: Path to FASTA file
        top_strand_only: If True, only count forward strand
    
    Returns:
        Dictionary mapping nucleotides to frequencies
    """
    counts = defaultdict(int)
    total = 0
    
    # Count bases
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).upper()
        
        for base in sequence:
            if base in 'ACGT':
                counts[base] += 1
                total += 1
                
                # If counting both strands, add reverse complement
                if not top_strand_only:
                    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                    counts[complement[base]] += 1
                    total += 1
    
    # Convert to frequencies
    if total == 0:
        return {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    
    frequencies = {base: counts[base] / total for base in 'ACGT'}
    
    return frequencies


def write_base_frequencies(frequencies: Dict[str, float], output_file: str):
    """
    Write base frequencies to file.
    
    Format matches TransDecoder's base_freqs.dat
    """
    with open(output_file, 'w') as f:
        for base in sorted(frequencies.keys()):
            f.write(f"{base}\t{frequencies[base]:.6f}\n")


def read_base_frequencies(input_file: str) -> Dict[str, float]:
    """Read base frequencies from file"""
    frequencies = {}
    
    with open(input_file) as f:
        for line in f:
            if line.strip():
                base, freq = line.strip().split('\t')
                frequencies[base] = float(freq)
    
    return frequencies


def calculate_gc_content(frequencies: Dict[str, float]) -> float:
    """Calculate GC content from base frequencies"""
    return frequencies.get('G', 0) + frequencies.get('C', 0)
