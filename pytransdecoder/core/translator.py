"""
Nucleotide to protein translation using BioPython genetic codes.

Replaces Nuc_translator.pm functionality.
"""

from typing import List, Set
from Bio.Data import CodonTable
from Bio.Seq import Seq


class Translator:
    """
    Handles translation of nucleotide sequences to protein.
    
    Supports multiple genetic codes via BioPython's CodonTable.
    """
    
    # Mapping from TransDecoder genetic code names to BioPython table IDs
    GENETIC_CODE_MAP = {
        'universal': 1,
        'standard': 1,
        'vertebrate_mitochondrial': 2,
        'yeast_mitochondrial': 3,
        'mold_mitochondrial': 4,
        'invertebrate_mitochondrial': 5,
        'ciliate': 6,
        'dasycladacean': 6,
        'hexamita': 6,
        'echinoderm_mitochondrial': 9,
        'euplotid': 10,
        'bacterial': 11,
        'alternative_yeast': 12,
        'ascidian_mitochondrial': 13,
        'flatworm_mitochondrial': 14,
        'chlorophycean_mitochondrial': 16,
        'trematode_mitochondrial': 21,
        'scenedesmus_obliquus_mitochondrial': 22,
        'thraustochytrium_mitochondrial': 23,
        'pterobranchia_mitochondrial': 24,
        'sr1_gracilibacteria': 25,
        'pachysolen_tannophilus': 26,
        'candida': 12,  # CUG encodes Ser instead of Leu
        'acetabularia': 6,
        'tetrahymena': 6,
        'mesodinium': 6,
        'peritrich': 6,
    }
    
    def __init__(self, genetic_code: str = "universal"):
        """
        Initialize translator with specified genetic code.
        
        Args:
            genetic_code: Name of genetic code (case-insensitive)
        """
        code_name = genetic_code.lower().replace(' ', '_').replace('-', '_')
        
        if code_name not in self.GENETIC_CODE_MAP:
            available = ', '.join(sorted(self.GENETIC_CODE_MAP.keys()))
            raise ValueError(
                f"Unknown genetic code: {genetic_code}. "
                f"Available codes: {available}"
            )
        
        self.genetic_code_name = genetic_code
        self.table_id = self.GENETIC_CODE_MAP[code_name]
        self.codon_table = CodonTable.unambiguous_dna_by_id[self.table_id]
    
    def get_stop_codons(self) -> List[str]:
        """
        Get stop codons for current genetic code.
        
        Returns:
            List of stop codon sequences (uppercase)
        """
        return sorted(self.codon_table.stop_codons)
    
    def get_start_codons(self, allow_non_met: bool = False) -> List[str]:
        """
        Get start codons for current genetic code.
        
        Args:
            allow_non_met: If True, return all possible start codons.
                          If False, return only ATG (standard start).
        
        Returns:
            List of start codon sequences (uppercase)
        """
        if allow_non_met:
            # Use all start codons from the table
            return sorted(self.codon_table.start_codons)
        else:
            # Standard start codon only
            return ['ATG']
    
    def translate(self, sequence: str, frame: int = 0) -> str:
        """
        Translate nucleotide sequence to protein.
        
        Args:
            sequence: DNA/RNA sequence
            frame: Reading frame offset (0, 1, or 2)
        
        Returns:
            Protein sequence (stops are represented as *)
        """
        # Adjust for frame
        seq = sequence[frame:].upper()
        
        # Convert to Bio.Seq and translate
        bio_seq = Seq(seq)
        
        try:
            protein = str(bio_seq.translate(table=self.table_id, to_stop=False))
        except Exception as e:
            # Handle translation errors (e.g., incomplete codons)
            # Translate as much as possible
            complete_codons = len(seq) // 3
            truncated = seq[:complete_codons * 3]
            if truncated:
                protein = str(Seq(truncated).translate(table=self.table_id, to_stop=False))
            else:
                protein = ""
        
        return protein
    
    @classmethod
    def get_available_genetic_codes(cls) -> List[str]:
        """Get list of supported genetic code names"""
        return sorted(cls.GENETIC_CODE_MAP.keys())
