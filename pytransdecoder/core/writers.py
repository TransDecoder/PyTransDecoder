"""
Output file writers for TransDecoder
"""

from typing import List, TextIO
from pathlib import Path
from .models import ORF


class GFF3Writer:
    """Write ORFs to GFF3 format"""
    
    def __init__(self, filepath: Path):
        self.filepath = filepath
        self.file: TextIO = None
    
    def __enter__(self):
        self.file = open(self.filepath, 'w')
        self._write_header()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file:
            self.file.close()
    
    def _write_header(self):
        """Write GFF3 header"""
        self.file.write("##gff-version 3\n")
    
    def write_orf(self, orf: ORF):
        """Write a single ORF"""
        self.file.write(orf.to_gff3() + "\n")


class FastaWriter:
    """Write sequences to FASTA format"""
    
    def __init__(self, filepath: Path):
        self.filepath = filepath
        self.file: TextIO = None
    
    def __enter__(self):
        self.file = open(self.filepath, 'w')
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file:
            self.file.close()
    
    def write_sequence(self, header: str, sequence: str, wrap: int = 60):
        """
        Write a sequence in FASTA format.
        
        Args:
            header: FASTA header (should start with >)
            sequence: Sequence string
            wrap: Line width for wrapping (0 = no wrap)
        """
        if not header.startswith('>'):
            header = '>' + header
        
        self.file.write(header + '\n')
        
        if wrap > 0:
            # Wrap sequence
            for i in range(0, len(sequence), wrap):
                self.file.write(sequence[i:i+wrap] + '\n')
        else:
            self.file.write(sequence + '\n')
