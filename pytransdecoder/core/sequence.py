"""
Sequence utilities
"""

from Bio.Seq import Seq


def reverse_complement(sequence: str) -> str:
    """
    Return reverse complement of a DNA/RNA sequence.
    
    Args:
        sequence: DNA or RNA sequence
        
    Returns:
        Reverse complement sequence
    """
    seq = Seq(sequence)
    return str(seq.reverse_complement())
