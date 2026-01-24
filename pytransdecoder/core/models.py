"""
Data models for PyTransDecoder

Minimal ORF representation - replaces the 5,588-line Gene_obj.pm
"""

from dataclasses import dataclass, field
from typing import Optional, Tuple, List


@dataclass
class ORF:
    """
    Open Reading Frame representation.
    
    This replaces Gene_obj.pm with only the essential functionality
    that TransDecoder actually uses.
    """
    # Identifiers
    transcript_id: str      # asmbl_id in Gene_obj
    gene_id: str           # TU_feat_name
    model_id: str          # Model_feat_name
    
    # Coordinates (1-based, as in original TransDecoder)
    start: int
    end: int
    strand: str            # '+' or '-'
    
    # Sequences
    sequence: str          # CDS nucleotide sequence
    protein: str           # Translated protein
    
    # Metadata
    orf_type: str         # complete, 5prime_partial, 3prime_partial, internal
    length: int           # Length in nucleotides
    frame: int = 0        # Reading frame (0, 1, 2)
    phase: int = 0        # CDS phase for GFF3 output
    description: str = "" # com_name in Gene_obj
    
    # For selection/ranking (populated in Predict phase)
    markov_scores: Optional[List[float]] = None
    blast_hits: List[str] = field(default_factory=list)
    pfam_hits: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        """Set defaults after initialization"""
        if not self.description:
            orient = '+' if self.strand == '+' else '-'
            self.description = f"ORF type:{self.orf_type} ({orient})"
    
    def get_cds_length(self) -> int:
        """Get CDS length in nucleotides (replaces Gene_obj::get_CDS_length)"""
        return self.length
    
    def get_model_span(self) -> Tuple[int, int]:
        """Get (min, max) coordinates (replaces Gene_obj::get_model_span)"""
        return (min(self.start, self.end), max(self.start, self.end))
    
    def to_gff3(self, source: str = "transdecoder") -> str:
        """
        Generate GFF3 representation (replaces Gene_obj::to_GFF3_format).
        
        Returns multi-line string with gene, mRNA, and CDS features.
        """
        lend, rend = self.get_model_span()
        
        # Attributes for each feature type
        gene_attrs = f"ID={self.gene_id};Name=ORF_{self.gene_id}"
        mrna_attrs = f"ID={self.model_id};Parent={self.gene_id};Name={self.model_id}"
        cds_attrs = f"ID=cds.{self.model_id};Parent={self.model_id}"
        
        lines = []
        
        # Gene feature
        lines.append("\t".join([
            self.transcript_id,
            source,
            "gene",
            str(lend),
            str(rend),
            ".",
            self.strand,
            ".",
            gene_attrs
        ]))
        
        # mRNA feature
        lines.append("\t".join([
            self.transcript_id,
            source,
            "mRNA",
            str(lend),
            str(rend),
            ".",
            self.strand,
            ".",
            mrna_attrs
        ]))
        
        # CDS feature (with phase)
        lines.append("\t".join([
            self.transcript_id,
            source,
            "CDS",
            str(lend),
            str(rend),
            ".",
            self.strand,
            str(self.phase),
            cds_attrs
        ]))
        
        return "\n".join(lines)
    
    def to_fasta_cds(self) -> str:
        """Generate FASTA format for CDS sequence"""
        header = f">{self.model_id} type:{self.orf_type} {self.transcript_id}:{self.start}-{self.end}({self.strand})"
        return f"{header}\n{self.sequence}"
    
    def to_fasta_protein(self, genetic_code: str = "universal") -> str:
        """Generate FASTA format for protein sequence"""
        header = f">{self.model_id} type:{self.orf_type} gc:{genetic_code} {self.transcript_id}:{self.start}-{self.end}({self.strand})"
        return f"{header}\n{self.protein}"
