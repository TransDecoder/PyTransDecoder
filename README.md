# PyTransDecoder

Python port of [TransDecoder](https://github.com/TransDecoder/TransDecoder) - identify candidate coding regions within transcript sequences.

## Installation

```bash
cd PyTransDecoder
pip install -e .
```

## Quick Start

### Phase 1: Extract Long ORFs

```bash
pytransdecoder longorfs -t transcripts.fasta
```

This creates a directory `transcripts.fasta.transdecoder_dir/` with:
- `longest_orfs.pep` - Protein sequences
- `longest_orfs.cds` - CDS sequences  
- `longest_orfs.gff3` - ORF annotations
- `base_freqs.dat` - Nucleotide frequencies

### Options

```
-t, --transcripts PATH          Input transcripts FASTA file [required]
-m, --min-protein-length INT    Minimum protein length (default: 100 aa)
-G, --genetic-code TEXT         Genetic code (default: Standard)
-S, --strand-specific           Only analyze top strand
-O, --output-dir PATH           Output directory
--gene-trans-map PATH           Gene-to-transcript mapping file
--complete-orfs-only            Only output complete ORFs
-v, --verbose                   Verbose output
```

## Supported Genetic Codes

- universal/standard
- vertebrate_mitochondrial
- yeast_mitochondrial
- invertebrate_mitochondrial
- ciliate/tetrahymena/dasycladacean
- euplotid
- bacterial
- candida
- And 15+ more...

## Phase 2: Predict (Coming Soon)

Phase 2 (TransDecoder.Predict) will be implemented after Phase 1 is validated.

## Development Status

- ✅ Phase 1 (LongOrfs): Implemented, ready for testing
- ⏳ Phase 2 (Predict): Not yet implemented

## Testing Against Perl Version

To validate output matches the Perl version:

```bash
# Run Perl version
cd ../TransDecoder
./TransDecoder.LongOrfs -t test.fasta

# Run Python version
cd ../PyTransDecoder  
pytransdecoder longorfs -t test.fasta

# Compare outputs
diff ../TransDecoder/test.fasta.transdecoder_dir/longest_orfs.pep \
     test.fasta.transdecoder_dir/longest_orfs.pep
```

## Requirements

- Python 3.8+
- BioPython >= 1.81
- Click >= 8.0
- tqdm >= 4.65
