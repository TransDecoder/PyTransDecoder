# PyTransDecoder Architecture

## System Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                     PyTransDecoder System                        │
└─────────────────────────────────────────────────────────────────┘

                              │
                              ▼
                    ┌─────────────────┐
                    │  CLI Interface  │
                    │   (Click/arg)   │
                    └─────────────────┘
                              │
                ┌─────────────┴─────────────┐
                │                           │
                ▼                           ▼
    ┌───────────────────┐      ┌──────────────────┐
    │  longorfs.py      │      │   predict.py     │
    │  (Phase 1)        │      │   (Phase 2)      │
    └───────────────────┘      └──────────────────┘
                │                           │
                └─────────────┬─────────────┘
                              │
            ┌─────────────────┼─────────────────┐
            │                 │                 │
            ▼                 ▼                 ▼
    ┌──────────┐      ┌──────────┐     ┌──────────┐
    │   Core   │      │ Scoring  │     │Selection │
    │  Modules │      │  System  │     │  Logic   │
    └──────────┘      └──────────┘     └──────────┘
```

## Module Dependency Graph

```
pytransdecoder/
│
├── __main__.py  ──────┐
│                      │
├── longorfs.py  ──────┼───> core.fasta
│                      │     core.translator
│                      │     core.orf_finder
│                      │     core.gff_utils
│                      │     utils.cli
│                      │
├── predict.py   ──────┼───> core.fasta
│                      │     core.gff_utils
│                      │     scoring.markov
│                      │     scoring.pwm
│                      │     selection.selector
│                      │     selection.homology
│                      │     utils.cli
│                      │
├── core/              │
│   ├── fasta.py  ─────┼───> biopython
│   ├── translator.py ─┼───> biopython
│   ├── orf_finder.py ─┼───> translator
│   │                  │     sequence
│   ├── sequence.py  ──┼───> biopython
│   ├── gff_utils.py ──┼───> gene_obj
│   └── gene_obj.py ───┘
│
├── scoring/
│   ├── markov.py ─────────> numpy, pandas
│   ├── pwm.py ────────────> numpy
│   ├── hexamer.py ────────> numpy
│   └── base_freqs.py ─────> numpy
│
├── selection/
│   ├── selector.py ───────> pandas, overlap
│   ├── homology.py ───────> pandas
│   └── overlap.py ────────> numpy
│
└── utils/
    ├── cli.py ────────────> click
    ├── checkpoints.py
    └── logging.py ────────> logging
```

## Data Flow: Phase 1 (LongOrfs)

```
Input: transcripts.fasta
   │
   ├──> FastaReader
   │       │
   │       ├──> For each transcript:
   │       │       │
   │       │       ├──> ORFFinder
   │       │       │      │
   │       │       │      ├──> Translator.get_stop_codons()
   │       │       │      │
   │       │       │      ├──> Find all stops in 6 frames
   │       │       │      │
   │       │       │      ├──> Find all starts upstream of stops
   │       │       │      │
   │       │       │      ├──> Build ORF objects
   │       │       │      │
   │       │       │      └──> Filter by min_length
   │       │       │
   │       │       └──> Collect ORFs
   │       │
   │       └──> BaseFrequencyCalculator
   │              │
   │              └──> Compute nucleotide frequencies
   │
   ├──> Output Writers:
   │      │
   │      ├──> GFF3Writer ────> longest_orfs.gff3
   │      ├──> FastaWriter ───> longest_orfs.cds
   │      ├──> FastaWriter ───> longest_orfs.pep
   │      └──> FreqWriter ────> base_freqs.dat
   │
   └──> Checkpoint: longorfs.ok
```

## Data Flow: Phase 2 (Predict)

```
Input: 
   - transcripts.fasta
   - longest_orfs.{cds,pep,gff3}
   - base_freqs.dat
   - [optional] blast_hits.outfmt6
   - [optional] pfam_hits.domtblout

   │
   ├──> Training Set Selection
   │       │
   │       ├──> Get top N longest ORFs
   │       ├──> Remove redundant proteins
   │       └──> Select final training set
   │
   ├──> Model Training
   │       │
   │       ├──> MarkovModel.train(training_cds)
   │       │      │
   │       │      ├──> Extract hexamer frequencies
   │       │      ├──> Compute log-likelihoods
   │       │      └──> Save hexamer.scores
   │       │
   │       └──> PWM.train(training_starts)
   │              │
   │              └──> Save start_pwm.model
   │
   ├──> Scoring All ORFs
   │       │
   │       ├──> For each ORF:
   │       │      │
   │       │      └──> MarkovModel.score_all_frames()
   │       │             │
   │       │             └──> Output: longest_orfs.cds.scores
   │       │
   │       └──> If homology data provided:
   │              │
   │              ├──> Parse BLAST hits
   │              └──> Parse Pfam hits
   │
   ├──> ORF Selection
   │       │
   │       ├──> For each transcript:
   │       │      │
   │       │      ├──> Get all ORFs
   │       │      │
   │       │      ├──> Filter ORFs:
   │       │      │      │
   │       │      │      ├──> Retain if long (GC-dependent threshold)
   │       │      │      ├──> Retain if BLAST hit
   │       │      │      ├──> Retain if Pfam hit
   │       │      │      └──> Retain if positive frame 1 score
   │       │      │
   │       │      ├──> Rank ORFs:
   │       │      │      │
   │       │      │      ├──> By homology evidence count
   │       │      │      ├──> By frame 1 score
   │       │      │      └──> By length
   │       │      │
   │       │      ├──> Remove overlapping ORFs
   │       │      │
   │       │      └──> Select best (or single best)
   │       │
   │       └──> Collect selected ORFs
   │
   ├──> Start Codon Refinement (optional)
   │       │
   │       └──> For 5' partial ORFs:
   │              │
   │              ├──> Scan upstream region
   │              ├──> Score with PWM
   │              └──> Extend to best start
   │
   ├──> Output Writers:
   │      │
   │      ├──> GFF3Writer ────> transcripts.transdecoder.gff3
   │      ├──> FastaWriter ───> transcripts.transdecoder.cds
   │      ├──> FastaWriter ───> transcripts.transdecoder.pep
   │      └──> BEDWriter ─────> transcripts.transdecoder.bed
   │
   └──> Checkpoint: predict.ok
```

## Class Diagrams

### Core Classes

```python
# ORF Data Structure
@dataclass
class ORF:
    """Represents an Open Reading Frame"""
    sequence: str           # Nucleotide sequence
    protein: str           # Translated protein sequence
    start: int             # Start position (1-based)
    end: int               # End position (1-based)
    length: int            # Length in nucleotides
    strand: str            # '+' or '-'
    frame: int             # Reading frame (0, 1, 2)
    transcript_id: str     # Parent transcript ID
    orf_id: str           # Unique ORF identifier
    
    # Completeness flags
    is_5prime_partial: bool = False
    is_3prime_partial: bool = False
    is_complete: bool = False
    
    # Scores (populated in Phase 2)
    markov_scores: Optional[np.ndarray] = None  # 6-frame scores
    pwm_score: Optional[float] = None
    
    # Homology (populated in Phase 2)
    blast_hits: List[str] = field(default_factory=list)
    pfam_hits: List[str] = field(default_factory=list)
    
    def get_coords(self) -> Tuple[int, int]:
        """Return (min, max) coordinates"""
        return (min(self.start, self.end), max(self.start, self.end))
    
    def overlaps(self, other: 'ORF', max_pct: float = 10.0) -> bool:
        """Check if this ORF overlaps another by more than max_pct"""
        # Implementation
        pass
    
    def to_gff3_line(self) -> str:
        """Convert to GFF3 format line"""
        # Implementation
        pass
```

### ORF Finder

```python
class ORFFinder:
    """
    Finds all Open Reading Frames in nucleotide sequences.
    
    Implements the core algorithm from Longest_orf.pm
    """
    
    def __init__(self,
                 min_protein_length: int = 100,
                 allow_5prime_partial: bool = True,
                 allow_3prime_partial: bool = True,
                 allow_non_met_starts: bool = False,
                 forward_strand: bool = True,
                 reverse_strand: bool = True,
                 genetic_code: str = "Standard",
                 complete_orfs_only: bool = False):
        """Initialize ORF finder with parameters"""
        self.min_protein_length = min_protein_length
        self.allow_5prime_partial = allow_5prime_partial
        self.allow_3prime_partial = allow_3prime_partial
        self.allow_non_met_starts = allow_non_met_starts
        self.forward_strand = forward_strand
        self.reverse_strand = reverse_strand
        self.genetic_code = genetic_code
        self.complete_orfs_only = complete_orfs_only
        
        self.translator = Translator(genetic_code)
    
    def find_all_orfs(self, sequence: str, seq_id: str) -> List[ORF]:
        """
        Find all ORFs in a sequence.
        
        Args:
            sequence: Nucleotide sequence
            seq_id: Sequence identifier
            
        Returns:
            List of ORF objects, sorted by length (descending)
        """
        orfs = []
        
        # Forward strand
        if self.forward_strand:
            orfs.extend(self._find_orfs_one_strand(sequence, seq_id, '+'))
        
        # Reverse strand
        if self.reverse_strand:
            rev_seq = reverse_complement(sequence)
            orfs.extend(self._find_orfs_one_strand(rev_seq, seq_id, '-'))
        
        # Sort by length descending
        orfs.sort(key=lambda x: x.length, reverse=True)
        
        return orfs
    
    def _find_orfs_one_strand(self, sequence: str, 
                              seq_id: str, 
                              strand: str) -> List[ORF]:
        """Find ORFs on one strand"""
        # 1. Find all stop codons
        stops = self._find_stop_codons(sequence)
        
        # 2. Find all start codons
        starts = self._find_start_codons(sequence, stops)
        
        # 3. Build ORFs from start-stop pairs
        orfs = self._build_orfs(starts, stops, sequence, seq_id, strand)
        
        # 4. Filter by minimum length
        min_nt_length = self.min_protein_length * 3
        orfs = [orf for orf in orfs if orf.length >= min_nt_length]
        
        return orfs
    
    def _find_stop_codons(self, sequence: str) -> List[List[int]]:
        """
        Find all stop codons in all 3 frames.
        
        Returns:
            List of 3 lists, one per frame
        """
        stop_codons = self.translator.get_stop_codons()
        stops = [[], [], []]  # One list per frame
        
        for i in range(len(sequence) - 2):
            frame = i % 3
            codon = sequence[i:i+3].upper()
            if codon in stop_codons:
                stops[frame].append(i)
        
        # Handle 3' partials: add implicit stop at sequence end
        if self.allow_3prime_partial:
            seq_len = len(sequence)
            for frame in range(3):
                # Add stop at last position in frame
                last_pos = seq_len - ((seq_len - frame) % 3)
                stops[frame].append(last_pos)
        
        return stops
    
    def _find_start_codons(self, sequence: str, 
                           stops: List[List[int]]) -> List[List[int]]:
        """
        Find all start codons upstream of stops.
        
        Returns:
            List of 3 lists, one per frame
        """
        start_codons = self.translator.get_start_codons(
            allow_non_met=self.allow_non_met_starts
        )
        starts = [[], [], []]  # One list per frame
        
        # Find explicit start codons
        for i in range(len(sequence) - 2):
            frame = i % 3
            codon = sequence[i:i+3].upper()
            if codon in start_codons:
                starts[frame].append(i)
        
        # Handle 5' partials: add implicit start at sequence beginning
        if self.allow_5prime_partial:
            for frame in range(3):
                starts[frame].insert(0, frame)  # Start at position 0, 1, or 2
        
        return starts
    
    def _build_orfs(self, starts: List[List[int]], 
                    stops: List[List[int]], 
                    sequence: str,
                    seq_id: str,
                    strand: str) -> List[ORF]:
        """Build ORF objects from start-stop pairs"""
        orfs = []
        
        for frame in range(3):
            frame_starts = starts[frame]
            frame_stops = stops[frame]
            
            for stop_pos in frame_stops:
                # Find all starts upstream of this stop
                valid_starts = [s for s in frame_starts if s < stop_pos]
                
                for start_pos in valid_starts:
                    # Extract ORF sequence
                    orf_seq = sequence[start_pos:stop_pos+3]
                    
                    # Translate
                    protein = self.translator.translate(orf_seq)
                    
                    # Check completeness
                    is_5prime_partial = (start_pos < 3)
                    is_3prime_partial = not self._has_stop_codon(orf_seq)
                    is_complete = not (is_5prime_partial or is_3prime_partial)
                    
                    # Skip if complete_orfs_only and not complete
                    if self.complete_orfs_only and not is_complete:
                        continue
                    
                    # Create ORF object
                    orf = ORF(
                        sequence=orf_seq,
                        protein=protein,
                        start=start_pos + 1,  # Convert to 1-based
                        end=stop_pos + 3,
                        length=len(orf_seq),
                        strand=strand,
                        frame=frame,
                        transcript_id=seq_id,
                        orf_id=f"{seq_id}.p{len(orfs)+1}",
                        is_5prime_partial=is_5prime_partial,
                        is_3prime_partial=is_3prime_partial,
                        is_complete=is_complete
                    )
                    
                    orfs.append(orf)
        
        return orfs
    
    def _has_stop_codon(self, sequence: str) -> bool:
        """Check if sequence ends with a stop codon"""
        if len(sequence) < 3:
            return False
        last_codon = sequence[-3:].upper()
        return last_codon in self.translator.get_stop_codons()
```

### Markov Model

```python
class MarkovModel:
    """
    5th-order Markov model for scoring coding potential.
    
    Scores sequences based on hexamer (6-mer) frequencies
    learned from a training set of coding sequences.
    """
    
    def __init__(self, order: int = 5):
        """
        Args:
            order: Markov model order (default: 5 for hexamers)
        """
        self.order = order
        self.kmer_length = order + 1  # 6 for order 5
        self.scores = {}  # hexamer -> log-likelihood score
        self.trained = False
    
    def train(self, sequences: List[str], 
              base_frequencies: Dict[str, float]):
        """
        Train model on coding sequences.
        
        Args:
            sequences: List of coding sequences
            base_frequencies: Background nucleotide frequencies
        """
        # Count hexamer frequencies in training set
        coding_counts = self._count_kmers(sequences)
        coding_freqs = self._normalize_counts(coding_counts)
        
        # Calculate expected frequencies under background model
        background_freqs = self._calculate_background_freqs(
            base_frequencies, self.kmer_length
        )
        
        # Calculate log-likelihood ratios
        for kmer in coding_freqs:
            coding_prob = coding_freqs[kmer]
            background_prob = background_freqs.get(kmer, 1e-10)
            self.scores[kmer] = np.log(coding_prob / background_prob)
        
        self.trained = True
    
    def score_all_frames(self, sequence: str) -> np.ndarray:
        """
        Score sequence in all 6 frames.
        
        Args:
            sequence: Nucleotide sequence
            
        Returns:
            Array of 6 scores [+1, +2, +3, -1, -2, -3]
        """
        if not self.trained:
            raise ValueError("Model must be trained before scoring")
        
        scores = np.zeros(6)
        
        # Forward strand
        scores[0] = self.score(sequence, frame=0)
        scores[1] = self.score(sequence[1:], frame=0)
        scores[2] = self.score(sequence[2:], frame=0)
        
        # Reverse strand
        rev_seq = reverse_complement(sequence)
        scores[3] = self.score(rev_seq, frame=0)
        scores[4] = self.score(rev_seq[1:], frame=0)
        scores[5] = self.score(rev_seq[2:], frame=0)
        
        return scores
    
    def score(self, sequence: str, frame: int = 0) -> float:
        """
        Score a single sequence/frame.
        
        Args:
            sequence: Nucleotide sequence
            frame: Reading frame offset (0, 1, or 2)
            
        Returns:
            Log-likelihood score
        """
        sequence = sequence.upper()
        total_score = 0.0
        
        # Use variable-order Markov model (up to order 5)
        for i in range(len(sequence)):
            # Determine effective order based on position
            effective_order = min(i, self.order)
            kmer_start = max(0, i - effective_order)
            kmer = sequence[kmer_start:i+1]
            
            # Lookup score (default to 0 if not found)
            total_score += self.scores.get(kmer, 0.0)
        
        return total_score
    
    def _count_kmers(self, sequences: List[str]) -> Dict[str, int]:
        """Count k-mer occurrences in sequences"""
        counts = {}
        for seq in sequences:
            seq = seq.upper()
            for i in range(len(seq) - self.kmer_length + 1):
                kmer = seq[i:i+self.kmer_length]
                if 'N' not in kmer:  # Skip ambiguous bases
                    counts[kmer] = counts.get(kmer, 0) + 1
        return counts
    
    def _normalize_counts(self, counts: Dict[str, int]) -> Dict[str, float]:
        """Convert counts to frequencies"""
        total = sum(counts.values())
        return {k: v / total for k, v in counts.items()}
    
    def save(self, filepath: str):
        """Save model to file"""
        with open(filepath, 'w') as f:
            for kmer, score in sorted(self.scores.items()):
                f.write(f"{kmer}\t{score:.6f}\n")
    
    def load(self, filepath: str):
        """Load model from file"""
        self.scores = {}
        with open(filepath) as f:
            for line in f:
                kmer, score = line.strip().split('\t')
                self.scores[kmer] = float(score)
        self.trained = True
```

## State Management

```python
class WorkflowState:
    """
    Manages workflow state and checkpoints.
    
    Allows resuming interrupted runs.
    """
    
    def __init__(self, workdir: Path):
        self.workdir = workdir
        self.checkpoint_dir = workdir / "__checkpoints"
        self.checkpoint_dir.mkdir(exist_ok=True)
    
    def is_complete(self, step_name: str) -> bool:
        """Check if a step has been completed"""
        checkpoint = self.checkpoint_dir / f"{step_name}.ok"
        return checkpoint.exists()
    
    def mark_complete(self, step_name: str):
        """Mark a step as complete"""
        checkpoint = self.checkpoint_dir / f"{step_name}.ok"
        checkpoint.touch()
    
    def clear_checkpoint(self, step_name: str):
        """Remove a checkpoint"""
        checkpoint = self.checkpoint_dir / f"{step_name}.ok"
        if checkpoint.exists():
            checkpoint.unlink()
```

## Error Handling Strategy

```python
class TransDecoderError(Exception):
    """Base exception for TransDecoder errors"""
    pass

class InvalidSequenceError(TransDecoderError):
    """Raised when sequence contains invalid characters"""
    pass

class InvalidGeneticCodeError(TransDecoderError):
    """Raised when genetic code is not supported"""
    pass

class FileFormatError(TransDecoderError):
    """Raised when input file format is invalid"""
    pass

class NoORFsFoundError(TransDecoderError):
    """Raised when no ORFs are found"""
    pass
```

## Configuration Management

```python
@dataclass
class LongOrfsConfig:
    """Configuration for LongOrfs phase"""
    transcripts_file: Path
    min_protein_length: int = 100
    genetic_code: str = "Standard"
    strand_specific: bool = False
    output_dir: Path = Path.cwd()
    gene_trans_map: Optional[Path] = None
    complete_orfs_only: bool = False
    verbose: bool = False

@dataclass
class PredictConfig:
    """Configuration for Predict phase"""
    transcripts_file: Path
    output_dir: Path = Path.cwd()
    top_orfs_train: int = 500
    retain_long_orfs_mode: str = "dynamic"  # or "strict"
    retain_long_orfs_length: int = 1000000
    retain_pfam_hits: Optional[Path] = None
    retain_blastp_hits: Optional[Path] = None
    single_best_only: bool = False
    no_refine_starts: bool = False
    genetic_code: str = "Standard"
    verbose: bool = False
```

This architecture provides a clean, modular design that should be maintainable and testable. Each component has clear responsibilities and well-defined interfaces.
