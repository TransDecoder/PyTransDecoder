# TransDecoder Utility Scripts

This directory contains utility scripts for format conversion and validation.

It also contains the legacy `TransDecoder.LongOrfs` and `TransDecoder.Predict`
compatibility wrappers. The preferred user-facing entrypoint is `pyTransdecoder`,
but these wrappers remain available for phase-specific workflows and older scripts.

## Python Implementations

These scripts are pure Python with no Perl dependencies:

- **fasta_prot_checker.py** - Validates protein FASTA files
  - Checks for start codon (M)
  - Checks for stop codon (*)
  - Detects internal stop codons
  
- **gff3_file_to_bed.py** - Converts GFF3 to BED format
  - Simple format conversion
  - Preserves strand and score information

## Perl Scripts (Temporary)

These scripts still use Perl but will be ported to Python in future versions:

- **cdna_alignment_orf_to_genome_orf.pl** - Maps cDNA ORF coordinates to genome
- **gtf_to_alignment_gff3.pl** - Converts GTF to alignment GFF3 format
- **gtf_genome_to_cdna_fasta.pl** - Extracts transcript sequences from genome
- **gtf_to_bed.pl** - Converts GTF to BED format
- **gff3_gene_to_gtf_format.pl** - Converts GFF3 genes to GTF format

## Perl Library

The `PerlLib/` directory contains Perl modules required by the Perl scripts:
- Gene_obj.pm - Gene object model
- GTF_utils2.pm - GTF parsing utilities
- GFF3_utils2.pm - GFF3 parsing utilities
- Fasta_reader.pm - FASTA file parsing
- And others...

## Usage

All scripts can be called directly:

```bash
# Preferred full pipeline
../pyTransdecoder -t transcripts.fa

# Legacy phase-specific wrappers
./TransDecoder.LongOrfs -t transcripts.fa
./TransDecoder.Predict -t transcripts.fa

# Python scripts
./fasta_prot_checker.py proteins.fa
./gff3_file_to_bed.py annotations.gff3 > output.bed

# Perl scripts
./gtf_to_bed.pl transcripts.gtf > output.bed
./fasta_prot_checker.pl proteins.fa
```

## Future Work

The goal is to port all Perl utilities to Python to make PyTransDecoder
a completely self-contained Python package with no Perl dependencies.

Priority porting order:
1. gtf_to_bed.py - GTF to BED conversion
2. gtf_to_alignment_gff3.py - GTF to GFF3 conversion  
3. gtf_genome_to_cdna_fasta.py - Transcript extraction
4. cdna_alignment_orf_to_genome_orf.py - Coordinate mapping (most complex)
