# PyTransDecoder Implementation Guide

## Quick Start Implementation Checklist

### Phase 0: Project Setup & Code Analysis (Days 1-2)

**IMPORTANT: Analyze Before Porting**

Before writing any code, determine what's actually needed:

- [ ] Grep through TransDecoder scripts to find Gene_obj method calls
- [ ] Identify which PerlLib modules are actually imported and used
- [ ] Map out the minimal data structures needed for ORFs
- [ ] Document which utility scripts are essential vs. optional
- [ ] Create a "must-have" vs "nice-to-have" feature list

**Example Analysis:**
```bash
# Find all Gene_obj method calls
grep -r "Gene_obj" TransDecoder/TransDecoder.* TransDecoder/util/
grep -r "->" TransDecoder/TransDecoder.LongOrfs | grep -v "^#"

# Find which modules are actually used
grep "^use " TransDecoder/TransDecoder.LongOrfs
grep "^use " TransDecoder/TransDecoder.Predict
```

### Phase 0: Project Setup (Day 1)

- [ ] Initialize Python package structure
- [ ] Create pyproject.toml and setup.py
- [ ] Set up virtual environment
- [ ] Install dependencies (BioPython, NumPy, Pandas, Click)
- [ ] Configure pytest
- [ ] Create initial module files

### Phase 1: Core Modules (Week 1)

#### Day 1-2: FASTA and Sequence Utilities

**Files to create:**
- `pytransdecoder/core/__init__.py`
- `pytransdecoder/core/sequence.py`
- `pytransdecoder/core/fasta.py`
- `tests/test_fasta.py`
- `tests/test_sequence.py`

**Implementation priority:**
1. Reverse complement function
2. FastaReader class using BioPython
3. FastaWriter class
4. Basic tests with sample sequences

#### Day 3-4: Translator

**Files to create:**
- `pytransdecoder/core/translator.py`
- `tests/test_translator.py`

**Key functions:**
- Translation using BioPython's genetic code tables
- Stop codon identification for each genetic code
- Start codon identification (Met and alternatives)
- Frame-aware translation

**Test cases:**
- Test all supported genetic codes
- Test partial sequences
- Test sequences with N's
- Verify stop codon lists match Perl version

#### Day 5-7: ORF Finder

**Files to create:**
- `pytransdecoder/core/orf_finder.py`
- `pytransdecoder/core/models.py` (ORF dataclass)
- `tests/test_orf_finder.py`

**Critical algorithm:**
- Must exactly replicate Longest_orf.pm logic
- Test with TransDecoder test data

**Test strategy:**
- Create simple test sequences with known ORFs
- Compare against Perl TransDecoder output
- Edge cases: very short sequences, no ORFs, all stops

### Phase 2: LongOrfs Command (Week 2)

#### Day 8-10: Main LongOrfs Script

**Files to create:**
- `pytransdecoder/longorfs.py`
- `pytransdecoder/core/gff_utils.py`
- `tests/test_longorfs.py`

**Implementation steps:**
1. CLI argument parsing with Click
2. Main workflow orchestration
3. GFF3 output generation
4. Progress bars with tqdm
5. Base frequency calculation

**Integration test:**
```bash
python -m pytransdecoder longorfs \
  -t tests/data/Trinity.fasta \
  -m 100 \
  -O test_output/
```

Compare outputs:
- Number of ORFs
- ORF coordinates
- Protein sequences
- GFF3 format

#### Day 11-12: GFF3 Output

**Key requirements:**
- Must match Perl TransDecoder GFF3 format exactly
- Validate with GFF3 validators
- Test with downstream tools

### Phase 3: Scoring System (Weeks 3-4)

#### Day 13-15: Base Frequencies and Hexamer Scoring

**Files to create:**
- `pytransdecoder/scoring/__init__.py`
- `pytransdecoder/scoring/base_freqs.py`
- `pytransdecoder/scoring/hexamer.py`
- `tests/test_base_freqs.py`
- `tests/test_hexamer.py`

**Implementation:**
1. Nucleotide frequency calculation
2. K-mer counting
3. Background model calculation

#### Day 16-19: Markov Model

**Files to create:**
- `pytransdecoder/scoring/markov.py`
- `tests/test_markov.py`

**Critical features:**
- 5th-order Markov model
- Log-likelihood scoring
- 6-frame scoring
- Match Perl scores within tolerance

**Validation:**
- Use same training sequences as Perl version
- Compare hexamer scores file
- Compare per-sequence scores

#### Day 20-22: Position Weight Matrix

**Files to create:**
- `pytransdecoder/scoring/pwm.py`
- `tests/test_pwm.py`

**Features:**
- Train on start codon contexts
- Score potential start sites
- Find best alternative starts

### Phase 4: Selection and Predict (Weeks 5-6)

#### Day 23-26: Homology Data Parsing

**Files to create:**
- `pytransdecoder/selection/__init__.py`
- `pytransdecoder/selection/homology.py`
- `tests/test_homology.py`

**Parse formats:**
- BLAST outfmt6
- HMMER domtblout

#### Day 27-30: ORF Selection Logic

**Files to create:**
- `pytransdecoder/selection/selector.py`
- `pytransdecoder/selection/overlap.py`
- `tests/test_selector.py`

**Implementation:**
1. ORF ranking algorithm
2. Overlap detection
3. Single best mode
4. GC-content-based thresholds

#### Day 31-35: Predict Command

**Files to create:**
- `pytransdecoder/predict.py`
- `tests/test_predict.py`

**Integration test:**
```bash
python -m pytransdecoder predict \
  -t tests/data/Trinity.fasta \
  --retain_blastp_hits tests/data/blastp.outfmt6 \
  --retain_pfam_hits tests/data/pfam.domtblout
```

### Phase 5: Testing and Validation (Week 7)

#### Day 36-38: Integration Tests

**Test all sample datasets:**
- simple_transcriptome_target
- cufflinks_example
- pasa_example
- stringtie_example
- supertranscripts_example

#### Day 39-40: Performance Testing

**Benchmarks:**
- Memory usage profiling
- Speed comparison with Perl version
- Large dataset testing (millions of transcripts)

#### Day 41-42: Edge Case Testing

**Test scenarios:**
- Empty input files
- Very short sequences (< 100 bp)
- Sequences with all N's
- Very long transcripts (> 100kb)
- All genetic codes
- Gzipped inputs

### Phase 6: Documentation and Polish (Week 8-9)

#### Day 43-45: Documentation

- [ ] README with installation instructions
- [ ] Usage examples for all modes
- [ ] API documentation
- [ ] Migration guide from Perl version

#### Day 46-48: Utility Scripts

Port utility scripts from util/:
- compute_base_probs.py
- score_cds_likelihood.py
- select_best_orfs.py
- train_start_pwm.py

#### Day 49-50: Packaging

- [ ] Finalize setup.py/pyproject.toml
- [ ] Create entry points
- [ ] Test installation via pip
- [ ] Create Docker container (optional)

---

## Detailed Implementation Examples

### Example 1: Basic ORF Finding

```python
# pytransdecoder/core/orf_finder.py

from dataclasses import dataclass
from typing import List, Optional
from Bio.Seq import Seq
from .translator import Translator
from .sequence import reverse_complement

@dataclass
class ORF:
    """Open Reading Frame"""
    sequence: str
    protein: str
    start: int  # 1-based
    end: int    # 1-based
    length: int
    strand: str
    frame: int
    transcript_id: str
    orf_id: str
    is_5prime_partial: bool = False
    is_3prime_partial: bool = False
    is_complete: bool = False

class ORFFinder:
    """Find ORFs in nucleotide sequences"""
    
    def __init__(self, 
                 min_protein_length: int = 100,
                 allow_5prime_partial: bool = True,
                 allow_3prime_partial: bool = True,
                 genetic_code: str = "Standard"):
        self.min_protein_length = min_protein_length
        self.min_nt_length = min_protein_length * 3
        self.allow_5prime_partial = allow_5prime_partial
        self.allow_3prime_partial = allow_3prime_partial
        self.translator = Translator(genetic_code)
    
    def find_all_orfs(self, sequence: str, seq_id: str) -> List[ORF]:
        """
        Find all ORFs in both strands.
        
        Returns ORFs sorted by length (longest first).
        """
        orfs = []
        
        # Forward strand
        orfs.extend(self._find_orfs_strand(sequence, seq_id, '+'))
        
        # Reverse strand
        rev_seq = reverse_complement(sequence)
        orfs.extend(self._find_orfs_strand(rev_seq, seq_id, '-'))
        
        # Sort by length descending
        orfs.sort(key=lambda x: x.length, reverse=True)
        
        return orfs
    
    def _find_orfs_strand(self, sequence: str, 
                          seq_id: str, 
                          strand: str) -> List[ORF]:
        """Find all ORFs on one strand"""
        sequence = sequence.upper()
        stop_codons = set(self.translator.get_stop_codons())
        start_codons = set(self.translator.get_start_codons())
        
        orfs = []
        
        # Process each frame
        for frame in range(3):
            # Find stops in this frame
            stops = []
            for i in range(frame, len(sequence) - 2, 3):
                codon = sequence[i:i+3]
                if codon in stop_codons:
                    stops.append(i)
            
            # Add implicit stop at end for 3' partials
            if self.allow_3prime_partial:
                last_pos = frame + ((len(sequence) - frame) // 3) * 3
                if last_pos not in stops:
                    stops.append(last_pos)
            
            # For each stop, look for starts
            for stop_pos in stops:
                # Find starts upstream in same frame
                starts = []
                
                # Add implicit start at beginning for 5' partials
                if self.allow_5prime_partial:
                    starts.append(frame)
                
                # Find explicit starts
                for i in range(frame, stop_pos, 3):
                    codon = sequence[i:i+3]
                    if codon in start_codons:
                        starts.append(i)
                
                # Build ORF from each start to this stop
                for start_pos in starts:
                    orf_seq = sequence[start_pos:stop_pos+3]
                    
                    # Skip if too short
                    if len(orf_seq) < self.min_nt_length:
                        continue
                    
                    # Translate
                    protein = self.translator.translate(orf_seq)
                    
                    # Check completeness
                    has_start = sequence[start_pos:start_pos+3] in start_codons
                    has_stop = sequence[stop_pos:stop_pos+3] in stop_codons
                    
                    is_5prime_partial = not has_start
                    is_3prime_partial = not has_stop
                    is_complete = has_start and has_stop
                    
                    # Create ORF
                    orf_num = len(orfs) + 1
                    orf = ORF(
                        sequence=orf_seq,
                        protein=protein,
                        start=start_pos + 1,  # Convert to 1-based
                        end=stop_pos + 3,
                        length=len(orf_seq),
                        strand=strand,
                        frame=frame,
                        transcript_id=seq_id,
                        orf_id=f"{seq_id}.p{orf_num}",
                        is_5prime_partial=is_5prime_partial,
                        is_3prime_partial=is_3prime_partial,
                        is_complete=is_complete
                    )
                    
                    orfs.append(orf)
        
        return orfs

```

### Example 2: Command-Line Interface

```python
# pytransdecoder/cli.py

import click
from pathlib import Path
from . import longorfs, predict

@click.group()
@click.version_option(version='0.1.0')
def cli():
    """PyTransDecoder: Python port of TransDecoder"""
    pass

@cli.command()
@click.option('-t', '--transcripts', 
              required=True,
              type=click.Path(exists=True),
              help='Transcripts FASTA file')
@click.option('-m', '--min-protein-length',
              default=100,
              type=int,
              help='Minimum protein length (default: 100)')
@click.option('-G', '--genetic-code',
              default='Standard',
              help='Genetic code (default: Standard)')
@click.option('-S', '--strand-specific',
              is_flag=True,
              help='Analyze top strand only')
@click.option('-O', '--output-dir',
              type=click.Path(),
              default='.',
              help='Output directory')
@click.option('--complete-orfs-only',
              is_flag=True,
              help='Only output complete ORFs')
@click.option('-v', '--verbose',
              is_flag=True,
              help='Verbose output')
def longorfs_cmd(transcripts, min_protein_length, genetic_code,
                 strand_specific, output_dir, complete_orfs_only, verbose):
    """Extract long ORFs from transcripts (Phase 1)"""
    
    config = longorfs.LongOrfsConfig(
        transcripts_file=Path(transcripts),
        min_protein_length=min_protein_length,
        genetic_code=genetic_code,
        strand_specific=strand_specific,
        output_dir=Path(output_dir),
        complete_orfs_only=complete_orfs_only,
        verbose=verbose
    )
    
    longorfs.run(config)

@cli.command()
@click.option('-t', '--transcripts',
              required=True,
              type=click.Path(exists=True),
              help='Transcripts FASTA file')
@click.option('-O', '--output-dir',
              type=click.Path(),
              default='.',
              help='Output directory (from LongOrfs)')
@click.option('-T', '--top-orfs-train',
              default=500,
              type=int,
              help='Top ORFs for training (default: 500)')
@click.option('--retain-pfam-hits',
              type=click.Path(exists=True),
              help='Pfam domain hits file')
@click.option('--retain-blastp-hits',
              type=click.Path(exists=True),
              help='BLASTP hits file (outfmt 6)')
@click.option('--single-best-only',
              is_flag=True,
              help='Retain only single best ORF per transcript')
@click.option('--no-refine-starts',
              is_flag=True,
              help='Skip start codon refinement')
@click.option('-v', '--verbose',
              is_flag=True,
              help='Verbose output')
def predict_cmd(transcripts, output_dir, top_orfs_train,
                retain_pfam_hits, retain_blastp_hits,
                single_best_only, no_refine_starts, verbose):
    """Predict likely coding regions (Phase 2)"""
    
    config = predict.PredictConfig(
        transcripts_file=Path(transcripts),
        output_dir=Path(output_dir),
        top_orfs_train=top_orfs_train,
        retain_pfam_hits=Path(retain_pfam_hits) if retain_pfam_hits else None,
        retain_blastp_hits=Path(retain_blastp_hits) if retain_blastp_hits else None,
        single_best_only=single_best_only,
        no_refine_starts=no_refine_starts,
        verbose=verbose
    )
    
    predict.run(config)

if __name__ == '__main__':
    cli()
```

### Example 3: GFF3 Output

```python
# pytransdecoder/core/gff_utils.py

from typing import List, TextIO
from .models import ORF

class GFF3Writer:
    """Write ORFs to GFF3 format"""
    
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.file = None
    
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
        """Write a single ORF as a GFF3 gene/mRNA/CDS structure"""
        
        # Gene feature
        gene_line = self._format_feature(
            seqid=orf.transcript_id,
            source="transdecoder",
            feature_type="gene",
            start=orf.start,
            end=orf.end,
            score=".",
            strand=orf.strand,
            phase=".",
            attributes=f"ID={orf.orf_id}.gene;Name=ORF_{orf.orf_id}"
        )
        self.file.write(gene_line + "\n")
        
        # mRNA feature
        mrna_line = self._format_feature(
            seqid=orf.transcript_id,
            source="transdecoder",
            feature_type="mRNA",
            start=orf.start,
            end=orf.end,
            score=".",
            strand=orf.strand,
            phase=".",
            attributes=f"ID={orf.orf_id};Parent={orf.orf_id}.gene"
        )
        self.file.write(mrna_line + "\n")
        
        # CDS feature
        cds_line = self._format_feature(
            seqid=orf.transcript_id,
            source="transdecoder",
            feature_type="CDS",
            start=orf.start,
            end=orf.end,
            score=".",
            strand=orf.strand,
            phase="0",
            attributes=f"ID={orf.orf_id}.cds;Parent={orf.orf_id}"
        )
        self.file.write(cds_line + "\n")
    
    def _format_feature(self, seqid, source, feature_type, start, end,
                       score, strand, phase, attributes):
        """Format a GFF3 feature line"""
        return "\t".join([
            seqid,
            source,
            feature_type,
            str(start),
            str(end),
            score,
            strand,
            phase,
            attributes
        ])
```

### Example 4: Testing Strategy

```python
# tests/test_orf_finder.py

import pytest
from pytransdecoder.core.orf_finder import ORFFinder
from pytransdecoder.core.models import ORF

def test_simple_complete_orf():
    """Test finding a simple complete ORF"""
    # ATG...TAA with in-frame stop
    sequence = "ATGGCATAA"  # M A *
    
    finder = ORFFinder(min_protein_length=1)  # Allow short ORFs for testing
    orfs = finder.find_all_orfs(sequence, "test_seq")
    
    assert len(orfs) > 0
    assert orfs[0].is_complete
    assert orfs[0].protein == "MA"
    assert orfs[0].start == 1
    assert orfs[0].end == 9

def test_5prime_partial():
    """Test 5' partial ORF (no start codon)"""
    sequence = "GCATAA"  # A * (missing start)
    
    finder = ORFFinder(
        min_protein_length=1,
        allow_5prime_partial=True
    )
    orfs = finder.find_all_orfs(sequence, "test_seq")
    
    assert len(orfs) > 0
    assert orfs[0].is_5prime_partial
    assert not orfs[0].is_complete

def test_3prime_partial():
    """Test 3' partial ORF (no stop codon)"""
    sequence = "ATGGCA"  # M A (missing stop)
    
    finder = ORFFinder(
        min_protein_length=1,
        allow_3prime_partial=True
    )
    orfs = finder.find_all_orfs(sequence, "test_seq")
    
    assert len(orfs) > 0
    assert orfs[0].is_3prime_partial
    assert not orfs[0].is_complete

def test_multiple_frames():
    """Test ORFs in different frames"""
    # Create sequence with ORFs in different frames
    sequence = "ATGATGATGTAA"  # Frame 0: M M M *
    
    finder = ORFFinder(min_protein_length=1)
    orfs = finder.find_all_orfs(sequence, "test_seq")
    
    # Should find ORF in frame 0
    frame_0_orfs = [orf for orf in orfs if orf.frame == 0]
    assert len(frame_0_orfs) > 0

def test_both_strands():
    """Test finding ORFs on both strands"""
    # Forward: ATG...TAA
    # Reverse comp: TTA...CAT
    sequence = "ATGGCATAA"
    
    finder = ORFFinder(min_protein_length=1)
    orfs = finder.find_all_orfs(sequence, "test_seq")
    
    plus_orfs = [orf for orf in orfs if orf.strand == '+']
    minus_orfs = [orf for orf in orfs if orf.strand == '-']
    
    assert len(plus_orfs) > 0
    # May or may not have minus ORFs depending on sequence

def test_min_length_filter():
    """Test minimum length filtering"""
    # Short ORF: 9 nt = 3 aa
    sequence = "ATGGCATAA"
    
    # Should find with min_protein_length=1
    finder1 = ORFFinder(min_protein_length=1)
    orfs1 = finder1.find_all_orfs(sequence, "test_seq")
    assert len(orfs1) > 0
    
    # Should not find with min_protein_length=10
    finder2 = ORFFinder(min_protein_length=10)
    orfs2 = finder2.find_all_orfs(sequence, "test_seq")
    assert len(orfs2) == 0

def test_genetic_code():
    """Test different genetic codes"""
    # UAG is stop in Standard, but Gln in some codes
    sequence = "ATGTAGTAA"  # ATG UAG UAA
    
    # Standard code: UAG is stop
    finder_std = ORFFinder(min_protein_length=1, genetic_code="Standard")
    orfs_std = finder_std.find_all_orfs(sequence, "test_seq")
    
    # Find shortest ORF (should stop at UAG)
    shortest_std = min(orfs_std, key=lambda x: x.length)
    assert shortest_std.length == 6  # ATG TAG

def test_complete_orfs_only():
    """Test complete ORFs only mode"""
    # One complete, one partial
    sequence = "ATGGCATAA" + "GCA"  # Complete + partial
    
    finder = ORFFinder(
        min_protein_length=1,
        allow_5prime_partial=False,
        allow_3prime_partial=False
    )
    orfs = finder.find_all_orfs(sequence, "test_seq")
    
    # All ORFs should be complete
    for orf in orfs:
        assert orf.is_complete or len(orfs) == 0

@pytest.mark.integration
def test_against_perl_output():
    """Integration test: compare against Perl TransDecoder output"""
    # Load test sequence
    test_seq_file = "tests/data/test_transcript.fasta"
    expected_orfs_file = "tests/data/expected_orfs.txt"
    
    # Run ORF finder
    # ... compare results ...
    # This test requires actual test data
    pass
```

---

## Critical Implementation Notes

### 1. ORF Coordinate System

TransDecoder uses **1-based coordinates**. Ensure all position calculations account for this:

```python
# Python string indexing is 0-based
sequence[0:3]  # First codon in 0-based indexing

# But ORF.start and ORF.end should be 1-based
orf.start = 1  # First nucleotide
orf.end = 3    # Third nucleotide
```

### 2. Frame Calculation

Reading frames:
- Frame 0: positions 0, 3, 6, 9, ...
- Frame 1: positions 1, 4, 7, 10, ...
- Frame 2: positions 2, 5, 8, 11, ...

```python
frame = position % 3
```

### 3. Strand Orientation

For reverse strand:
1. Reverse complement the sequence
2. Find ORFs as normal
3. Coordinates are relative to reverse-complemented sequence
4. May need to convert back to forward strand coordinates for output

### 4. Stop Codon Handling

Stop codons at the END of an ORF are included in the nucleotide sequence but NOT in the protein sequence:

```python
orf_seq = "ATGGCATAA"  # ATG GCA TAA
protein = "MA"         # Stop (*) not included
```

### 5. Random Seed

For reproducibility, set random seed to 1234 (matching Perl version):

```python
import random
random.seed(1234)
```

### 6. Floating Point Comparison

When comparing scores with Perl output, use tolerance:

```python
import math
assert math.isclose(python_score, perl_score, rel_tol=1e-4)
```

---

## Development Workflow

### 1. Start with Tests

Write tests first for each module:
```bash
pytest tests/test_orf_finder.py -v
```

### 2. Incremental Implementation

Implement one method at a time, running tests after each:
```bash
pytest tests/test_orf_finder.py::test_simple_complete_orf -v
```

### 3. Compare with Perl

Regularly compare outputs with Perl TransDecoder:
```bash
# Run Perl version
cd TransDecoder
./TransDecoder.LongOrfs -t test.fasta

# Run Python version
cd PyTransDecoder
python -m pytransdecoder longorfs -t test.fasta

# Compare
diff TransDecoder/test.fasta.transdecoder_dir/longest_orfs.pep \
     PyTransDecoder/test.fasta.transdecoder_dir/longest_orfs.pep
```

### 4. Profile Performance

Use profiling tools to identify bottlenecks:
```bash
python -m cProfile -o profile.stats -m pytransdecoder longorfs -t large_test.fasta
python -m pstats profile.stats
```

### 5. Continuous Integration

Set up GitHub Actions or similar:
```yaml
# .github/workflows/test.yml
name: Test
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: actions/setup-python@v2
        with:
          python-version: '3.9'
      - run: pip install -e .[dev]
      - run: pytest tests/ -v --cov=pytransdecoder
```

This implementation guide provides concrete examples and a step-by-step approach to porting TransDecoder to Python.
