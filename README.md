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

## CLI Compatibility

PyTransDecoder's CLI accepts **both underscore and dash formats** for option names, ensuring compatibility with existing Perl TransDecoder workflows:

```bash
# Both of these work identically:
pytransdecoder predict -t transcripts.fasta --retain-pfam-hits pfam.domtblout
pytransdecoder predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout
```

This means you can use existing scripts and Makefiles without modification.

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

## Phase 2: Predict Likely Coding Regions

### Basic Usage

```bash
pytransdecoder predict -t transcripts.fasta
```

This creates final output files in the same directory as your transcripts:
- `transcripts.fasta.transdecoder.gff3` - Gene predictions
- `transcripts.fasta.transdecoder.bed` - BED format
- `transcripts.fasta.transdecoder.pep` - Protein sequences
- `transcripts.fasta.transdecoder.cds` - CDS sequences

### With Homology Support

```bash
# Run BLASTP against UniProt
blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep \
       -db uniprot_sprot.pep -max_target_seqs 1 \
       -outfmt 6 -evalue 1e-5 -num_threads 4 > blastp.outfmt6

# Run Pfam domain search
hmmscan --cpu 4 --domtblout pfam.domtblout Pfam-A.hmm \
        transcripts.fasta.transdecoder_dir/longest_orfs.pep

# Run predict with homology data
pytransdecoder predict -t transcripts.fasta \
    --retain-blastp-hits blastp.outfmt6 \
    --retain-pfam-hits pfam.domtblout
```

### Options

```
-t, --transcripts PATH            Input transcripts FASTA file [required]
-O, --output-dir PATH             Output directory
-T, --top-orfs-train INT          Training ORFs (default: 500)
--retain-long-orfs-mode           'dynamic' or 'strict' (default: dynamic)
--retain-pfam-hits PATH           Pfam domain hits (domtblout format)
--retain-blastp-hits PATH         BLAST hits (outfmt 6 format)
--single-best-only                Only best ORF per transcript
--no-refine-starts                Skip start codon refinement
-G, --genetic-code TEXT           Genetic code (default: Standard)
-v, --verbose                     Verbose output
```

### How It Works

1. **Training**: Selects top longest unique ORFs to train hexamer scoring model
2. **Scoring**: Scores all ORFs using Markov chain model (hexamer composition)
3. **Selection**: Chooses best ORFs based on:
   - Homology support (BLAST/Pfam hits)
   - Coding potential score
   - ORF length
   - GC content-based thresholds
4. **Output**: Generates final predictions in multiple formats

## Development Status

- ✅ Phase 1 (LongOrfs): Implemented and validated
- ✅ Phase 2 (Predict): Implemented and tested
- ⏳ Performance benchmarking on large datasets

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
- tqdm >= 4.65

**Note**: Ported by Claude.io Sonnet 4.5 under guidance by bhaas. Jan 24, 2026.

