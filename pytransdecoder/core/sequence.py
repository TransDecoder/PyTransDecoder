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
    if sequence.isascii():
        return sequence.translate(_rc_table())[::-1]
    return str(Seq(sequence).reverse_complement())


_RC = None


def _rc_table():
    """Per-character complement of ASCII, taken from Biopython (whose reverse complement maps each
    character independently), so translate + reverse equals Seq.reverse_complement."""
    global _RC
    if _RC is None:
        _RC = str.maketrans({chr(o): str(Seq(chr(o)).reverse_complement()) for o in range(128)})
    return _RC
