"""
ORF (Open Reading Frame) finder

Replicates the core algorithm from Longest_orf.pm
"""

from typing import List, Dict, Tuple
from .models import ORF
from .translator import Translator
from .sequence import reverse_complement


class ORFFinder:
    """
    Find all Open Reading Frames in nucleotide sequences.
    
    Replicates the functionality of Longest_orf.pm
    """
    
    def __init__(
        self,
        min_protein_length: int = 100,
        allow_5prime_partial: bool = True,
        allow_3prime_partial: bool = True,
        allow_non_met_starts: bool = False,
        forward_strand: bool = True,
        reverse_strand: bool = True,
        genetic_code: str = "universal",
        complete_orfs_only: bool = False
    ):
        """
        Initialize ORF finder.
        
        Args:
            min_protein_length: Minimum protein length in amino acids
            allow_5prime_partial: Allow ORFs without start codon
            allow_3prime_partial: Allow ORFs without stop codon
            allow_non_met_starts: Allow non-Met start codons
            forward_strand: Analyze forward strand
            reverse_strand: Analyze reverse strand
            genetic_code: Genetic code name
            complete_orfs_only: Only return complete ORFs (has start and stop)
        """
        self.min_protein_length = min_protein_length
        self.min_nt_length = min_protein_length * 3
        self.allow_5prime_partial = allow_5prime_partial
        self.allow_3prime_partial = allow_3prime_partial
        self.allow_non_met_starts = allow_non_met_starts
        self.forward_strand = forward_strand
        self.reverse_strand = reverse_strand
        self.genetic_code = genetic_code
        self.complete_orfs_only = complete_orfs_only
        
        self.translator = Translator(genetic_code)
        self.stop_codons = set(self.translator.get_stop_codons())
        self.start_codons = set(self.translator.get_start_codons(allow_non_met_starts))
    
    def find_all_orfs(self, sequence: str, seq_id: str) -> List[ORF]:
        """
        Find all ORFs in a sequence.
        
        Args:
            sequence: Nucleotide sequence
            seq_id: Sequence identifier
        
        Returns:
            List of ORF objects, sorted by length (longest first)
        """
        orfs = []
        original_length = len(sequence)
        
        # Forward strand
        if self.forward_strand:
            orfs.extend(self._find_orfs_strand(sequence, seq_id, '+', original_length))
        
        # Reverse strand
        if self.reverse_strand:
            rev_seq = reverse_complement(sequence)
            orfs.extend(self._find_orfs_strand(rev_seq, seq_id, '-', original_length))
        
        # Sort by length descending (longest first)
        orfs.sort(key=lambda x: x.length, reverse=True)
        
        return orfs
    
    def _find_orfs_strand(self, sequence: str, seq_id: str, strand: str, original_length: int) -> List[ORF]:
        """Find all ORFs on one strand"""
        sequence = sequence.upper()
        orfs = []
        
        # Find stops and starts for each frame
        for frame in range(3):
            # Find all stop codons in this frame
            stops, real_stops = self._find_stops_in_frame(sequence, frame)
            
            # Find all start codons in this frame (passing real_stops to avoid starts at stop positions)
            starts = self._find_starts_in_frame(sequence, frame, real_stops)
            
            # Build ORFs from start-stop pairs
            orfs.extend(
                self._build_orfs_from_coords(
                    starts, stops, real_stops, sequence, seq_id, strand, frame, original_length
                )
            )
        
        return orfs
    
    def _find_stops_in_frame(self, sequence: str, frame: int) -> Tuple[List[int], set]:
        """
        Find all stop codon positions in a specific frame.
        
        Returns:
            Tuple of (all_stops, real_stop_set) where:
                - all_stops: List of all stop positions (real + implicit)
                - real_stop_set: Set of positions that are actual stop codons
        """
        stops = []
        real_stops = set()
        
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon in self.stop_codons:
                stops.append(i)
                real_stops.add(i)
        
        # Add implicit stops at end for 3' partials (NOT real stop codons)
        # Perl adds seq_length, seq_length-1, seq_length-2 as implicit stops
        if self.allow_3prime_partial:
            seq_len = len(sequence)
            for offset in [0, -1, -2]:
                pos = seq_len + offset
                if pos >= 0 and pos not in real_stops:
                    stops.append(pos)
        
        return stops, real_stops
    
    def _find_starts_in_frame(self, sequence: str, frame: int, real_stops: set) -> List[int]:
        """
        Find all start codon positions in a specific frame.
        
        Args:
            sequence: DNA/RNA sequence
            frame: Reading frame (0, 1, or 2)
            real_stops: Set of positions that are actual stop codons (to avoid adding stops as starts)
            
        Returns:
            List of start positions
        """
        starts = []
        
        # Add implicit start at beginning for 5' partials
        # But don't add it if it's already a stop codon (Perl's logic)
        if self.allow_5prime_partial and frame not in real_stops:
            starts.append(frame)
        
        # Find explicit start codons
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon in self.start_codons:
                starts.append(i)
        
        return starts
    
    def _build_orfs_from_coords(
        self,
        starts: List[int],
        stops: List[int],
        real_stops: set,
        sequence: str,
        seq_id: str,
        strand: str,
        frame: int,
        original_length: int
    ) -> List[ORF]:
        """
        Build ORF objects from start and stop coordinate lists.
        
        Important: Replicates Perl's algorithm where:
        - Each stop codon is used only once
        - Each start codon is paired with only its first valid stop
        This prevents generating overlapping ORFs.
        
        Args:
            real_stops: Set of positions that are actual stop codons
                       (vs implicit end-of-sequence stops)
        """
        orfs = []
        last_used_stop = -1  # Track the last stop position used in this frame
        
        # Iterate through starts in order
        for start_pos in starts:
            # Skip starts that come before or at the last used stop
            # This prevents reusing the same stop codon
            if start_pos <= last_used_stop:
                continue
            
            # Find the first valid stop after this start (Perl's algorithm)
            valid_stop = None
            for stop_pos in stops:
                if (stop_pos > start_pos and 
                    (stop_pos - start_pos) % 3 == 0):
                    valid_stop = stop_pos
                    break  # Take first valid stop (replicates Perl's 'last')
            
            if valid_stop is None:
                continue
            
            # Mark this stop as used
            last_used_stop = valid_stop
            
            # Check if this is a real stop codon or implicit end-of-sequence
            is_real_stop = valid_stop in real_stops
            
            # Extract ORF sequence
            # For real stops: include the 3 bases of the stop codon
            # For implicit stops (3' partial): extract up to end of sequence
            if is_real_stop:
                # Real stop codon - include it (valid_stop points to start of stop codon)
                orf_end = valid_stop + 3
            else:
                # Implicit stop (3' partial) - extract to end of sequence
                orf_end = len(sequence)
            
            orf_seq = sequence[start_pos:orf_end]
            
            # Ensure CDS length is a multiple of 3 (trim partial codons)
            remainder = len(orf_seq) % 3
            if remainder != 0:
                orf_seq = orf_seq[:-remainder]
                orf_end = start_pos + len(orf_seq)
            
            # Length check
            if len(orf_seq) < self.min_nt_length:
                continue
            
            # Translate
            protein = self.translator.translate(orf_seq)
            
            # Verify no internal stops (sanity check from Perl code)
            stop_count = protein.count('*')
            if stop_count > 1:
                # Should not happen if algorithm is correct
                continue
            
            # Determine ORF type
            has_start = sequence[start_pos:start_pos+3] in self.start_codons
            has_stop = is_real_stop  # Only real stop codons count
            
            is_5prime_partial = not has_start
            is_3prime_partial = not has_stop
            
            # Determine ORF type string
            if not is_5prime_partial and not is_3prime_partial:
                orf_type = "complete"
            elif not is_5prime_partial:
                orf_type = "3prime_partial"
            elif not is_3prime_partial:
                orf_type = "5prime_partial"
            else:
                orf_type = "internal"
            
            # Skip if complete_orfs_only and not complete
            if self.complete_orfs_only and orf_type != "complete":
                continue
            
            # Calculate coordinates
            # For reverse strand, transform back to original sequence coordinates
            # Perl uses: revcomp_coord($pos, $seq_length) = $seq_length - $pos + 1
            # Perl: start_pos_adj = start_pos + 1, stop_pos_adj = stop_pos + 1 + 2
            # Perl ALWAYS adds +3 to stop position, even for implicit stops
            # This affects both coordinate reporting and length calculations
            if strand == '-':
                # Transform coordinates back to original sequence. For 3' partial
                # ORFs, use the trimmed sequence end rather than the implicit stop
                # position so the reported span matches the extracted CDS.
                start_adj = start_pos + 1  # 1-based position of start
                end_adj = orf_end

                # Transform to original sequence coordinates
                orf_start = original_length - start_adj + 1
                orf_end_coord = original_length - end_adj + 1
            else:
                # Forward strand: convert to 1-based. For 3' partial ORFs, use the
                # trimmed CDS endpoint rather than valid_stop + 3.
                orf_start = start_pos + 1  # 1-based position of start
                orf_end_coord = orf_end
            
            # Create ORF object
            orf = ORF(
                transcript_id=seq_id,
                gene_id="",  # Will be set later
                model_id="",  # Will be set later
                start=orf_start,
                end=orf_end_coord,
                strand=strand,
                sequence=orf_seq,
                protein=protein,
                orf_type=orf_type,
                length=len(orf_seq),
                frame=frame
            )
            
            orfs.append(orf)
        
        return orfs
