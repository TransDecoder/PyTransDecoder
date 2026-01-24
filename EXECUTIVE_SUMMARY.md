# TransDecoder Python Port - Executive Summary

## Overview

TransDecoder is a mature Perl-based bioinformatics tool for identifying protein-coding regions in transcript sequences. This document summarizes the research findings for porting it to Python.

## What TransDecoder Does

TransDecoder predicts likely coding sequences from RNA transcript assemblies using:
1. **ORF identification** - Finds all potential open reading frames
2. **Markov modeling** - Scores coding potential using hexamer frequencies
3. **Homology evidence** - Incorporates BLAST and Pfam domain matches
4. **Start codon refinement** - Optimizes 5' boundaries using position weight matrices

## Current Implementation Analysis

### Size and Complexity
- **Main scripts**: 2 Perl scripts (~900 lines total)
- **Core libraries**: 14 Perl modules (~7,800 lines total)
  - Largest: Gene_obj.pm (5,588 lines) - but much of this may not be needed
- **Utility scripts**: 23 helper scripts
- **Supported genetic codes**: 20+ variants
- **File formats**: FASTA, GFF3, GTF, tab-delimited

### Key Algorithms

#### Phase 1: LongOrfs (Initial Discovery)
1. Parse FASTA transcripts
2. Find all stop codons in 6 reading frames
3. Find start codons upstream of stops
4. Build ORFs from start-stop pairs
5. Handle partial ORFs (missing start/stop)
6. Filter by minimum length (default: 100 aa)
7. Output: CDS, peptides, GFF3, base frequencies

#### Phase 2: Predict (Selection)
1. Select top N longest ORFs for training (default: 500)
2. Remove redundant sequences
3. Train 5th-order Markov model on hexamers
4. Score all ORFs in all 6 frames
5. Apply retention criteria:
   - Long ORFs (GC-content-dependent threshold)
   - BLAST/Pfam homology hits
   - Positive coding score in correct frame
6. Rank and select best ORF(s) per transcript
7. Refine start codons using PWM
8. Output final predictions

## Python Port Strategy

### Recommended Tech Stack

**Core Dependencies:**
- **BioPython** (1.81+) - FASTA I/O, translation, genetic codes
- **NumPy** (1.24+) - Numerical operations, array handling
- **Pandas** (2.0+) - Data tables for scores and hits
- **Click** (8.1+) - Command-line interface

**Testing:**
- **pytest** - Unit and integration testing
- **pytest-cov** - Code coverage

### Proposed Architecture

```
pytransdecoder/
├── __init__.py
├── __main__.py              # CLI entry point
├── longorfs.py              # Phase 1 command
├── predict.py               # Phase 2 command
├── core/                    # Core functionality
│   ├── fasta.py            # FASTA I/O
│   ├── translator.py       # Translation with genetic codes
│   ├── orf_finder.py       # ORF discovery algorithm
│   ├── sequence.py         # Sequence utilities
│   ├── gff_utils.py        # GFF3/GTF handling
│   └── models.py           # Data structures (ORF class)
├── scoring/                 # Scoring system
│   ├── base_freqs.py       # Nucleotide frequencies
│   ├── markov.py           # 5th-order Markov model
│   ├── hexamer.py          # Hexamer scoring
│   └── pwm.py              # Position weight matrix
├── selection/               # ORF selection
│   ├── selector.py         # Selection logic
│   ├── homology.py         # BLAST/Pfam parsing
│   └── overlap.py          # Overlap detection
└── utils/                   # Utilities
    ├── cli.py              # CLI helpers
    ├── checkpoints.py      # Progress tracking
    └── logging.py          # Logging configuration
```

### Implementation Timeline

**Total estimated time: 7 weeks** (with selective porting)

- **Week 1**: Code analysis + Core modules (leverage BioPython, minimal ORF model)
- **Week 2**: LongOrfs command and testing
- **Week 3**: Scoring system (Markov, PWM)
- **Weeks 4-5**: Predict command and selection logic
- **Week 6**: Integration testing and validation
- **Week 7**: Documentation and packaging

**Time saved by selective porting:**
- Using BioPython instead of custom FASTA/translator: ~4 days
- Minimal ORF dataclass instead of full Gene_obj: ~7 days
- Skipping unused utility scripts: ~3 days

### Key Advantages of Python Port

1. **Better maintainability** - Modern language, cleaner syntax
2. **Rich ecosystem** - BioPython, NumPy, Pandas
3. **Easier installation** - pip install, no CPAN dependencies
4. **Better testing** - pytest framework
5. **Type hints** - Better code clarity and IDE support
6. **Parallel processing** - Built-in multiprocessing
7. **Container-friendly** - Easy Docker integration

## Critical Success Factors

### Must-Have Features (Parity with Perl)
✓ All genetic code tables (20+ variants)
✓ Both strands or strand-specific
✓ Partial ORF handling (5' and 3')
✓ Complete ORFs only mode
✓ Minimum length filtering
✓ GFF3/GTF output
✓ BLAST/Pfam integration
✓ Markov model scoring
✓ Start codon refinement
✓ Single best ORF mode
✓ GC-content-based thresholds

### Performance Requirements
- Handle millions of transcripts
- Memory-efficient streaming
- Comparable speed to Perl version
- Progress indicators for long runs

### Validation Strategy
1. **Unit tests** - Each module independently
2. **Integration tests** - Full workflows
3. **Comparison tests** - Match Perl outputs
4. **Performance tests** - Speed and memory benchmarks

## Key Challenges and Solutions

### Challenge 1: Gene_obj.pm Complexity (5,588 lines)
**Solution**: Analyze actual usage first - likely only need 5-10% of functionality. Implement simple ORF dataclass (~100-200 lines) with just what TransDecoder actually uses, not the full 5,588 lines

### Challenge 2: Exact Output Matching
**Solution**: Use same random seed (1234), same precision, incremental validation

### Challenge 3: Performance
**Solution**: NumPy for vectorization, multiprocessing for parallelization, profiling

### Challenge 4: Genetic Code Support
**Solution**: Leverage BioPython's built-in genetic code tables

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Outputs don't match Perl | High | Incremental testing, small datasets first |
| Performance issues | Medium | Profile early, optimize hot paths, consider Cython/Numba |
| Missing edge cases | Medium | Comprehensive test suite, run on diverse data |
| Scope creep | Low | Strict feature parity initially, enhancements later |

## Recommended Next Steps

1. **Approve scope** - Confirm feature requirements and timeline
2. **Set up project** - Initialize Python package structure
3. **Implement core** - Start with FASTA, translator, ORF finder
4. **Iterate with testing** - Build test suite alongside code
5. **Validate continuously** - Compare with Perl outputs frequently
6. **Document thoroughly** - Keep docs updated as you build

## Resource Requirements

### Development
- 1 developer, ~9 weeks full-time
- Access to test datasets (provided in sample_data/)
- Computing resources for testing on large datasets

### Dependencies
- Python 3.8+ 
- BioPython, NumPy, Pandas, Click
- pytest for testing
- Optional: HMMER, BLAST+ (for external searches)

## Expected Outcomes

### Deliverables
1. **pytransdecoder** Python package
2. Command-line tools (longorfs, predict)
3. Comprehensive test suite
4. Documentation (README, API docs, migration guide)
5. PyPI package (optional)
6. Docker container (optional)

### Success Metrics
- ✅ All test cases pass
- ✅ Outputs match Perl version (within tolerance)
- ✅ Performance within 2x of Perl version
- ✅ Easy installation (pip install)
- ✅ Clear documentation
- ✅ Used successfully by early adopters

## Long-Term Vision

### Phase 1: Feature Parity (Weeks 1-9)
Port all existing TransDecoder functionality to Python

### Phase 2: Enhancements (Future)
- Parallel processing built-in
- Better progress reporting
- JSON output format
- BED format output
- REST API for web services
- More flexible filtering options
- Better error messages
- GUI (optional)

### Phase 3: Advanced Features (Future)
- Machine learning-based scoring
- Integration with annotation pipelines
- Cloud-native deployment
- Real-time processing
- Support for more file formats

## Conclusion

Porting TransDecoder to Python is **feasible and valuable**. The Perl codebase is well-structured with clear algorithms that translate naturally to Python. BioPython provides excellent support for core bioinformatics operations. The main implementation effort involves:

1. **ORF finding algorithm** (most critical, ~2 weeks)
2. **Markov model scoring** (moderately complex, ~2 weeks)  
3. **Selection logic** (moderate complexity, ~2 weeks)
4. **Integration and testing** (essential, ~3 weeks)

The result will be a more maintainable, easier-to-install, and more extensible tool that preserves all the functionality of the original TransDecoder.

---

## References

- [TransDecoder Wiki](https://github.com/TransDecoder/TransDecoder/wiki)
- [BioPython Documentation](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
- [NCBI Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- [GFF3 Specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

---

**Prepared by**: GitHub Copilot  
**Date**: January 24, 2026  
**Version**: 1.0
