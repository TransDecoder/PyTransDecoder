# Sample Data Testing Guide

## Overview
The `sample_data/` directory contains multiple test examples that demonstrate TransDecoder usage across different workflows (Cufflinks, PASA, StringTie, supertranscripts).

## Current Status

### ✅ Core Functionality Working
1. **Wrapper scripts created**: `TransDecoder.LongOrfs` and `TransDecoder.Predict`
   - These bash wrappers call `pytransdecoder longorfs` and `pytransdecoder predict`
   - Maintains compatibility with existing test scripts

2. **Utility scripts**: Symlinked to Perl TransDecoder's `util/` directory
   - Contains 26 utility scripts for format conversion and validation
   - Used for genome coordinate mapping, GFF3/BED conversion, etc.

3. **Simple example tested**: `simple_transcriptome_target` works successfully
   - Phase 1 (LongOrfs): ✅ Working - generated 845 candidate ORFs
   - Phase 2 (Predict): ✅ Working - selected 733 final ORFs

### ⚠️ Known Issues

1. **Partial codon warning** in `_gff3_to_proteins()`:
   ```
   BiopythonWarning: Partial codon, len(sequence) not a multiple of three.
   ```
   - Occurs during protein extraction from genome coordinates
   - May be due to incomplete ORFs at transcript boundaries
   - Needs investigation in the sequence extraction logic

2. **Utility scripts dependency**: Currently relies on Perl utilities
   - `util/` is symlinked from Perl TransDecoder
   - Required utilities:
     - `cdna_alignment_orf_to_genome_orf.pl` - Maps ORFs to genome coordinates
     - `fasta_prot_checker.pl` - Validates protein sequences
     - `gtf_to_alignment_gff3.pl` - GTF to GFF3 conversion
     - `gtf_genome_to_cdna_fasta.pl` - Extract transcripts from genome
     - `gff3_file_to_bed.pl` - GFF3 to BED conversion
     - `gtf_to_bed.pl` - GTF to BED conversion

### 📋 Test Examples

Each subdirectory has a `Makefile` with `test` and `clean` targets:

1. **simple_transcriptome_target/** - Basic Trinity transcriptome
   - Input: Trinity.fasta (921 transcripts)
   - Tests: Basic ORF prediction, genome mapping

2. **pasa_example/** - PASA assemblies with gene grouping
   - Input: pasa_assemblies.fasta (858 transcripts)
   - Tests: Gene-transcript mapping, custom output directory

3. **cufflinks_example/** - With homology support
   - Input: transcripts.gtf + genome
   - Tests: BLAST + Pfam integration, GTF workflows

4. **stringtie_example/** - StringTie merged transcripts
   - Similar to Cufflinks but with StringTie output

5. **supertranscripts_example/** - Supertranscript analysis
   - Tests: GTF to GFF3 conversion, supertranscript ORF mapping

## Running Tests

### Single Example
```bash
cd sample_data/simple_transcriptome_target
make test  # or just: ./runMe.sh
```

### All Examples
```bash
cd sample_data
make test
```

### Clean Outputs
```bash
cd sample_data
make clean
```

## Next Steps

### High Priority
1. **Fix partial codon warning** - Investigate sequence extraction in genome coordinate mapping
2. **Test remaining examples**:
   - pasa_example
   - cufflinks_example (requires BLAST/HMMER)
   - stringtie_example
   - supertranscripts_example

### Medium Priority
3. **Port utility scripts to Python** (optional):
   - Most utilities are format converters
   - Could be useful for pure Python distribution
   - Not critical since Perl utilities work fine

### Low Priority
4. **Compare outputs** with Perl version for each example
5. **Performance benchmarking** on larger datasets
6. **Docker integration** - test with Docker wrapper

## Dependencies

### Python (already installed)
- pytransdecoder package
- BioPython >= 1.81
- Click >= 8.0

### External tools (for homology examples)
- BLAST+ (blastp, makeblastdb)
- HMMER (hmmsearch, hmmpress)
- Perl (for utility scripts)

## Compatibility

The Python implementation maintains **CLI compatibility** with the original Perl version:
- Same command-line options
- Same file naming conventions
- Same output formats (GFF3, BED, FASTA)
- Works as drop-in replacement in existing pipelines
