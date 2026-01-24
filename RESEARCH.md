# TransDecoder to Python Port - Research Document

## Executive Summary

TransDecoder is a Perl-based tool for identifying candidate protein-coding regions within transcript sequences. This document analyzes the current Perl implementation and outlines the strategy for porting it to Python.

## Current Architecture Overview

> **📋 See [GENE_OBJ_ANALYSIS.md](GENE_OBJ_ANALYSIS.md) for detailed analysis of what code is actually needed from Gene_obj.pm (spoiler: only ~3% of the 5,588 lines!)**

### Two-Phase Workflow

TransDecoder operates in two distinct phases:

1. **TransDecoder.LongOrfs** - Initial ORF discovery
2. **TransDecoder.Predict** - ORF selection and refinement

### Key Components

#### 1. Main Scripts
- `TransDecoder.LongOrfs` (422 lines) - Extracts all potential ORFs from transcripts
- `TransDecoder.Predict` (502 lines) - Scores ORFs and selects best candidates

#### 2. Core Perl Libraries (PerlLib/)
- **Fasta_reader.pm** - FASTA file parsing
- **Nuc_translator.pm** - DNA/RNA to protein translation with multiple genetic codes
- **Longest_orf.pm** - ORF identification algorithm
- **Gene_obj.pm** (5588 lines) - Complex gene/exon/CDS data structure
- **PWM.pm** - Position Weight Matrix for start codon refinement
- **GFF3_utils2.pm** - GFF3 format handling
- **GTF_utils2.pm** - GTF format handling
- **Pipeliner.pm** - Command execution and checkpointing
- **DelimParser.pm** - Delimited file parsing
- **Overlap_piler.pm** - Handling overlapping features

#### 3. Utility Scripts (util/)
- `compute_base_probs.pl` - Calculate nucleotide frequencies
- `score_CDS_likelihood_all_6_frames.pl` - Markov model scoring
- `seq_n_baseprobs_to_loglikelihood_vals.pl` - Generate hexamer scores
- `select_best_ORFs_per_transcript.pl` - ORF selection logic
- `train_start_PWM.pl` - Position weight matrix training
- `start_codon_refinement.pl` - Start codon adjustment
- `exclude_similar_proteins.pl` - Redundancy removal
- `get_top_longest_fasta_entries.pl` - Extract longest sequences
- `cdna_alignment_orf_to_genome_orf.pl` - Map cDNA ORFs to genome
- `gff3_file_to_proteins.pl` - Extract proteins from GFF3

## Algorithm Details

### Phase 1: LongOrfs

**Input:** Transcripts FASTA file

**Process:**
1. Parse transcripts using Fasta_reader
2. For each transcript:
   - Identify all stop codons in all 6 frames
   - Identify all start codons (Met or alternative)
   - Pair starts with stops to form ORFs
   - Handle partial ORFs (5' and 3')
3. Filter by minimum protein length (default: 100 aa)
4. Compute base frequencies for later use
5. Output: CDS sequences, peptides, GFF3 annotations

**Key Features:**
- Strand-specific option (top strand only)
- Support for multiple genetic codes (Universal, Euplotes, Tetrahymena, Candida, etc.)
- Complete ORFs only option (require Met start and stop codon)
- Gene-to-transcript mapping support

### Phase 2: Predict

**Input:** 
- Transcripts FASTA
- LongOrfs output (CDS, peptides, GFF3)
- Optional: BLASTP results, Pfam domain hits

**Process:**
1. Select top N longest ORFs (default: 500) for training
2. Remove redundant sequences
3. Build 5th-order Markov model from training set
4. Score all ORFs in all 6 frames using hexamer frequencies
5. Apply retention rules:
   - Automatically retain long ORFs (dynamic threshold based on GC content)
   - Retain ORFs with BLAST hits
   - Retain ORFs with Pfam domain hits
   - Retain ORFs with positive coding score in frame 1
6. Select best ORF per transcript (optional: single best only)
7. Refine start codons using PWM (optional)
8. Output final predictions

**Scoring System:**
- **Hexamer (6-mer) Markov scoring** - Measures coding potential
- **GC-content-based thresholds** - Dynamic ORF length cutoffs (25-80% GC)
- **Homology prioritization** - BLAST/Pfam hits ranked highest
- **Length-based tiebreaking** - Longer ORFs preferred

### Key Algorithms

#### ORF Finding (Longest_orf.pm)
```
1. Find all stop codons in sequence
2. For each stop codon:
   - Look upstream for start codons (Met or alternatives)
   - Create ORF from start to stop
3. Handle partial ORFs:
   - 5' partial: no start codon (begins at sequence start)
   - 3' partial: no stop codon (ends at sequence end)
4. Sort ORFs by length descending
```

#### Markov Model Training
```
1. Extract hexamer frequencies from training CDS set
2. Calculate log-likelihood ratios vs background frequencies
3. Store scores for all 4^6 = 4096 hexamers
4. Score sequence: sum log-likelihood for overlapping hexamers
```

#### Start Codon Refinement (PWM)
```
1. Train Position Weight Matrix on known start sites
2. For 5' partial ORFs:
   - Scan upstream region
   - Score potential start codons using PWM
   - Select highest-scoring valid start
```

## Technical Requirements for Python Port

### 1. Core Dependencies

#### Bioinformatics Libraries
- **BioPython** - FASTA/GFF3 parsing, sequence manipulation, translation
  - `Bio.SeqIO` for FASTA reading
  - `Bio.Seq` for sequence operations and translation
  - `Bio.SeqUtils.CodonUsage` for codon statistics
  
- **gffutils** or **BCBio.GFF** - GFF3 file handling
  
#### Data Science Libraries
- **NumPy** - Numerical computations, array operations
- **Pandas** - Data structures for scores, hits, annotations
- **SciPy** - Statistical functions if needed

#### Utilities
- **Click** or **argparse** - Command-line interface
- **tqdm** - Progress bars
- **logging** - Better than print statements

### 2. File Format Support

Must handle:
- **FASTA** - Input transcripts, output peptides/CDS
- **GFF3** - Output ORF annotations
- **GTF** - Alternative annotation format
- **Tab-delimited** - Gene-transcript maps, scores, BLAST/Pfam results

### 3. Genetic Code Support

BioPython provides built-in genetic code tables:
- Standard (Universal)
- Vertebrate Mitochondrial
- Yeast Mitochondrial
- Mold/Protozoan Mitochondrial
- Invertebrate Mitochondrial
- Ciliate/Dasycladacean/Hexamita
- Echinoderm/Flatworm Mitochondrial
- Euplotid
- Bacterial/Plant Plastid
- Alternative Yeast
- Ascidian Mitochondrial
- Alternative Flatworm Mitochondrial
- Chlorophycean Mitochondrial
- Trematode Mitochondrial
- Scenedesmus obliquus Mitochondrial
- Thraustochytrium Mitochondrial
- Pterobranchia Mitochondrial
- SR1/Gracilibacteria
- And more...

## Proposed Python Architecture

### Directory Structure
```
PyTransDecoder/
├── pytransdecoder/
│   ├── __init__.py
│   ├── __main__.py              # Entry point
│   ├── longorfs.py              # Phase 1 implementation
│   ├── predict.py               # Phase 2 implementation
│   ├── core/
│   │   ├── __init__.py
│   │   ├── fasta.py            # FASTA I/O
│   │   ├── translator.py       # Translation engine
│   │   ├── orf_finder.py       # ORF discovery
│   │   ├── gene_obj.py         # Gene/ORF data structures
│   │   ├── gff_utils.py        # GFF3/GTF handling
│   │   └── sequence.py         # Sequence utilities
│   ├── scoring/
│   │   ├── __init__.py
│   │   ├── markov.py           # Markov model
│   │   ├── pwm.py              # Position weight matrix
│   │   ├── hexamer.py          # Hexamer scoring
│   │   └── base_freqs.py       # Nucleotide frequencies
│   ├── selection/
│   │   ├── __init__.py
│   │   ├── selector.py         # ORF selection logic
│   │   ├── homology.py         # BLAST/Pfam parsing
│   │   └── overlap.py          # Overlap detection
│   └── utils/
│       ├── __init__.py
│       ├── cli.py              # Command-line parsing
│       ├── checkpoints.py      # Progress tracking
│       └── logging.py          # Logging setup
├── scripts/                     # Utility scripts
│   ├── compute_base_probs.py
│   ├── score_cds_likelihood.py
│   ├── train_start_pwm.py
│   └── select_best_orfs.py
├── tests/
│   ├── test_orf_finder.py
│   ├── test_translator.py
│   ├── test_markov.py
│   └── test_integration.py
├── setup.py
├── pyproject.toml
├── README.md
└── requirements.txt
```

### Module Responsibilities

#### pytransdecoder.longorfs
- Main entry point for Phase 1
- CLI argument parsing
- Orchestrates ORF finding workflow
- Generates output files

#### pytransdecoder.predict
- Main entry point for Phase 2
- Orchestrates training and scoring
- ORF selection and filtering
- Final output generation

#### pytransdecoder.core.fasta
```python
class FastaReader:
    """Iterator over FASTA sequences."""
    def __init__(self, filepath: str)
    def __iter__(self) -> Iterator[SeqRecord]
    
class FastaWriter:
    """Write sequences to FASTA."""
    def write(self, seq_id: str, sequence: str, description: str = "")
```

#### pytransdecoder.core.translator
```python
class Translator:
    """Translate nucleotide sequences with configurable genetic codes."""
    def __init__(self, genetic_code: str = "Standard")
    def translate(self, sequence: str, frame: int = 0) -> str
    def get_stop_codons(self) -> List[str]
    def get_start_codons(self, allow_non_met: bool = False) -> List[str]
```

#### pytransdecoder.core.orf_finder
```python
@dataclass
class ORF:
    """Data structure for an Open Reading Frame."""
    sequence: str
    protein: str
    start: int
    end: int
    length: int
    strand: str
    frame: int
    is_5prime_partial: bool
    is_3prime_partial: bool
    is_complete: bool

class ORFFinder:
    """Find ORFs in nucleotide sequences."""
    def __init__(self, 
                 min_length: int = 100,
                 allow_5prime_partial: bool = True,
                 allow_3prime_partial: bool = True,
                 allow_non_met_starts: bool = False,
                 both_strands: bool = True,
                 genetic_code: str = "Standard")
    
    def find_orfs(self, sequence: str, seq_id: str) -> List[ORF]
    def _find_stop_codons(self, sequence: str) -> List[int]
    def _find_start_codons(self, sequence: str, stops: List[int]) -> List[int]
    def _build_orfs(self, starts: List[int], stops: List[int], 
                    sequence: str, strand: str) -> List[ORF]
```

#### pytransdecoder.scoring.markov
```python
class MarkovModel:
    """5th-order Markov model for coding potential."""
    def __init__(self, order: int = 5)
    def train(self, sequences: List[str])
    def score(self, sequence: str, frame: int = 0) -> float
    def score_all_frames(self, sequence: str) -> np.ndarray
    def save(self, filepath: str)
    def load(self, filepath: str)
```

#### pytransdecoder.scoring.pwm
```python
class PositionWeightMatrix:
    """PWM for start codon scoring."""
    def __init__(self, length: int = 20)
    def add_sequence(self, sequence: str)
    def build(self)
    def score(self, sequence: str) -> float
    def find_best_start(self, sequence: str, 
                        window: int = 50) -> Tuple[int, float]
```

#### pytransdecoder.selection.selector
```python
class ORFSelector:
    """Select best ORFs per transcript."""
    def __init__(self,
                 min_auto_accept_length: int = 1000000,
                 blast_hits: Dict = None,
                 pfam_hits: Dict = None,
                 single_best: bool = False,
                 max_overlap_pct: float = 10.0)
    
    def select_orfs(self, orfs: List[ORF], scores: pd.DataFrame) -> List[ORF]
    def _rank_orf(self, orf: ORF) -> Tuple[int, float, int]
    def _check_overlap(self, orf1: ORF, orf2: ORF) -> float
```

## Implementation Strategy

### Minimal Viable Port Philosophy

**Core Principle: Port only what's actively used**

Before implementing any module, verify it's actually needed:

1. **Grep for usage** - Search TransDecoder scripts for method/function calls
2. **Trace execution** - Follow the code path to see what's invoked
3. **Test-driven** - If a test passes without it, you don't need it
4. **Iterative** - Start minimal, add only when needed

**Example: Gene_obj.pm Analysis**

Gene_obj.pm is 5588 lines, but TransDecoder likely uses only:
- Basic ORF coordinate storage
- CDS sequence retrieval  
- GFF3 output formatting
- Maybe a few utility methods

**Estimated actual usage: ~5-10% of Gene_obj.pm code**

Instead of porting all 5588 lines, create a simple ORF dataclass with just:
- Coordinates (start, end, strand)
- Sequences (DNA, protein)
- Metadata (transcript_id, orf_id, completeness flags)
- to_gff3() method

This could be ~100-200 lines instead of 5588!

**Apply same analysis to all modules:**
- Fasta_reader.pm → Use BioPython SeqIO (already exists!)
- Nuc_translator.pm → Use BioPython translation (already exists!)
- Longest_orf.pm → Core algorithm, must port carefully
- PWM.pm → Needed for start codon refinement
- Others → Evaluate individually

### Phase 1: Core Functionality (Weeks 1-2)

1. **Set up project structure**
   - Initialize Python package
   - Configure testing framework (pytest)
   - Set up CI/CD if needed

2. **Implement core modules**
   - FASTA reader/writer using BioPython
   - Translator with genetic code support
   - Sequence utilities (reverse complement, etc.)

3. **Implement ORF finder**
   - Stop codon identification
   - Start codon identification
   - ORF assembly
   - Handle partial ORFs
   - Unit tests with known sequences

### Phase 2: LongOrfs Script (Week 3)

1. **Implement TransDecoder.LongOrfs**
   - CLI interface with Click
   - Main workflow orchestration
   - GFF3 output generation
   - Progress bars with tqdm

2. **Testing**
   - Run on test data from `sample_data/`
   - Compare outputs with Perl version
   - Validate GFF3 format

### Phase 3: Scoring System (Weeks 4-5)

1. **Base frequency calculation**
   - Compute nucleotide frequencies
   - Background model

2. **Markov model**
   - Hexamer frequency extraction
   - Log-likelihood calculation
   - 6-frame scoring

3. **PWM for start codon refinement**
   - Position-specific scoring
   - Start codon scanning

### Phase 4: Predict Script (Weeks 6-7)

1. **Implement TransDecoder.Predict**
   - CLI interface
   - Training set selection
   - Redundancy removal
   - Model training and scoring

2. **ORF selection logic**
   - Parse BLAST/Pfam results
   - Implement ranking system
   - Handle overlaps
   - Single best mode

3. **Start codon refinement**
   - PWM training
   - 5' ORF extension

### Phase 5: Testing & Validation (Week 8)

1. **Integration testing**
   - Run on all sample datasets
   - Compare outputs with Perl version
   - Performance benchmarking

2. **Edge cases**
   - Empty inputs
   - Very short sequences
   - Unusual genetic codes
   - Large datasets (memory/speed)

### Phase 6: Documentation & Polish (Week 9)

1. **User documentation**
   - README with examples
   - Installation instructions
   - Tutorial

2. **API documentation**
   - Docstrings for all public functions
   - Generate Sphinx docs

3. **Packaging**
   - PyPI-ready setup
   - Version management

## Key Challenges & Solutions

### Challenge 1: Gene_obj.pm Complexity
**Problem:** 5588-line Perl module with complex gene structure representation

**Solution:** 
- **Analyze actual usage first** - Determine which Gene_obj methods TransDecoder actually calls
- Start with minimal dataclass-based ORF structure
- Only port the specific Gene_obj functionality that's actively used
- Use Python's dataclasses and type hints for cleaner implementation
- **Key insight:** Much of Gene_obj.pm is likely unused by TransDecoder and can be skipped entirely

### Challenge 2: Performance
**Problem:** Processing large transcriptomes

**Solution:**
- Use NumPy for numerical operations
- Consider multiprocessing for embarrassingly parallel tasks
- Profile and optimize hot paths
- Consider Numba for JIT compilation if needed

### Challenge 3: Exact Output Matching
**Problem:** May need to match Perl output exactly for validation

**Solution:**
- Use same random seed (1234)
- Use same rounding/formatting
- Test incrementally with small datasets
- Document any intentional differences

### Challenge 4: Dependency on External Tools
**Problem:** Some Perl scripts call external binaries

**Solution:**
- Identify all external dependencies
- Consider whether to re-implement or call external tools
- Document required external tools clearly

## Testing Strategy

### Unit Tests
- Test each module independently
- Use pytest fixtures for common data
- Mock external dependencies

### Integration Tests
- Use data from `__testing/` directory
- Compare outputs with expected results
- Test full workflows end-to-end

### Validation Tests
- Run on `sample_data/` examples
- Compare with Perl TransDecoder outputs
- Check:
  - Number of ORFs found
  - ORF coordinates
  - Protein sequences
  - GFF3 format validity
  - Scores (within floating-point tolerance)

## Performance Considerations

### Memory
- Use iterators for FASTA parsing (don't load all into memory)
- Stream processing where possible
- Consider chunking for very large files

### Speed
- Critical path is ORF finding and scoring
- NumPy for vectorized operations
- Multiprocessing for independent transcripts
- Profile before optimizing

### Scalability
- Should handle millions of transcripts
- Progress indicators for long runs
- Checkpoint system to resume interrupted runs

## Compatibility Notes

### Input/Output Compatibility
- Must accept same command-line arguments (or provide migration guide)
- Must produce compatible GFF3 output
- Protein sequences should match exactly
- CDS sequences should match exactly

### Feature Parity
Must support:
- All genetic codes
- Strand-specific mode
- Complete ORFs only mode
- Gene-to-transcript mapping
- BLAST/Pfam integration
- Single best mode
- Start codon refinement
- Dynamic/strict long ORF retention modes

### Nice-to-Have Additions
- Better progress reporting
- JSON output option
- BED format output
- More flexible filtering options
- Better error messages
- Parallel processing built-in

## Dependencies Summary

### Required Python Packages
```
biopython>=1.81
numpy>=1.24
pandas>=2.0
click>=8.1
tqdm>=4.65
```

### Optional Packages
```
gffutils>=0.12  # Alternative GFF handling
pysam>=0.21     # If BAM/SAM support needed
pytest>=7.4     # Testing
pytest-cov      # Coverage
black           # Code formatting
mypy            # Type checking
```

### External Tools (Optional)
- HMMER (for Pfam searches)
- BLAST+ (for homology searches)

## Migration Guide for Users

### Command Equivalence
```bash
# Perl version
TransDecoder.LongOrfs -t transcripts.fasta -m 100

# Python version
pytransdecoder longorfs -t transcripts.fasta -m 100
# or
python -m pytransdecoder longorfs -t transcripts.fasta -m 100
```

### Output Compatibility
- GFF3 files should be compatible
- FASTA files should be identical
- Checkpoint directories may differ in structure

## Timeline Estimate

**With selective porting (only essential code):**

- **Week 1:** Code analysis + Core modules (FASTA/BioPython, minimal ORF model, ORF finder)
- **Week 2:** LongOrfs implementation and testing
- **Week 3:** Scoring system (Markov, PWM) 
- **Week 4:** Predict implementation (training, scoring)
- **Week 5:** ORF selection logic and homology parsing
- **Week 6:** Integration testing and validation
- **Week 7:** Documentation, polish, and packaging

**Total: ~7 weeks for focused port** (down from 9 weeks)

**Time savings by avoiding unnecessary code:**
- Gene_obj.pm: Skip ~5000 lines → Save ~1 week
- Use BioPython instead of porting FASTA/translator → Save ~3-4 days
- Skip unused utility scripts → Save ~2-3 days

## Success Criteria

1. ✅ All modules pass unit tests
2. ✅ Integration tests pass on sample data
3. ✅ Outputs match Perl version (within tolerance)
4. ✅ Performance comparable to Perl version
5. ✅ All genetic codes supported
6. ✅ Documentation complete
7. ✅ Easy to install and use

## Next Steps

1. Review this research document
2. Confirm requirements and scope
3. Set up Python project structure
4. Begin Phase 1 implementation
5. Establish testing framework early
6. Iterate with frequent validation against Perl version

## References

- TransDecoder Wiki: https://github.com/TransDecoder/TransDecoder/wiki
- BioPython Tutorial: http://biopython.org/DIST/docs/tutorial/Tutorial.html
- Genetic Codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
- GFF3 Specification: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
