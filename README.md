# PyTransDecoder

Python port of [TransDecoder](https://github.com/TransDecoder/TransDecoder) - identify candidate coding regions within transcript sequences.

## Installation

```bash
cd PyTransDecoder
pip install -e .
```

## Quick Start

Use `pyTransdecoder` as the primary interface when you want the standard end-to-end pipeline:

```bash
pyTransdecoder -t transcripts.fasta
```

This runs LongOrfs followed by Predict and writes the usual `*.transdecoder.{gff3,bed,pep,cds}` outputs.

You can also have the pipeline run homology support searches directly:

```bash
pyTransdecoder -t transcripts.fasta \
    --blast-search-pep uniprot_sprot.pep \
    --pfam-search-db Pfam-A.hmm
```

## Phase-Specific Commands

If you want to run the phases separately, use the subcommand CLI:

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
-m, --min-protein-length INT    Minimum protein length (default: 50 aa; ORFs under 100 codons
                                must also pass Predict's short-ORF frame-contrast test)
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

Legacy compatibility wrappers for `TransDecoder.LongOrfs` and `TransDecoder.Predict` now live under [`util/`](./util) instead of the repository root.

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

ORFs with a hit are kept even when their coding signal alone would not keep them. Hits to
transposable-element proteins (or to host proteins with a TE-derived segment) can therefore keep
short ORFs in TE-containing transcripts. These are kept on purpose; to find and remove them by
best-hit protein or by repeat overlap, see
[docs/homology_te_screen.md](docs/homology_te_screen.md).

For the full pipeline entrypoint, you can skip the separate Pfam step and let
`pyTransdecoder` prepare and search the HMM database for you:

```bash
pyTransdecoder -t transcripts.fasta --pfam-search-db Pfam-A.hmm
```

### Options

```
-t, --transcripts PATH            Input transcripts FASTA file [required]
-O, --output-dir PATH             Output directory
-T, --top-orfs-train INT          Training ORFs (default: 500)
--retain-long-orfs-mode           'dynamic' or 'strict' (default: dynamic)
--retain-pfam-hits PATH           Pfam domain hits (domtblout format)
--retain-blastp-hits PATH         BLAST hits (outfmt 6 format)
--retain-multiple-orfs            Report all non-overlapping ORFs passing selection
                                  (default: the best ORF per transcript)
--no-refine-starts                Skip start codon refinement
--no-homology-start-assist        Do not let BLASTP/Pfam hits decide 5'-partial vs trimmed starts
                                  (default: they do whenever hits are supplied)
--start-contrast-alpha FLOAT      Start refinement: UTR score per base of upstream segment (default: 0.4)
--start-kozak-weight FLOAT        Start refinement: bases of coding signal per bit of Kozak score (default: 17)
--max-qvalue FLOAT                Keep ORFs (without homology) with target-decoy q-value <= this
                                  (default: 0.05)
--precise                         Same as --max-qvalue 0.01: fewer calls on noncoding sequence
--sensitive                       Same as --max-qvalue 0.2: more short and partial ORFs
--decoy-shuffles INT              Shuffled copies per transcript (default: 0 = automatic)
--decoy-seed INT                  Seed of the decoy shuffles (default: 42)
--preset NAME|DIR                 Models and thresholds from a preset instead of the input
                                  (recommended below ~1,000 transcripts or for mostly partial sets)
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
   - GC content-based thresholds (ORFs long enough to be unlikely in random sequence of the
     transcriptome's GC content are kept regardless of score; applied only when the genetic
     code's stop codons are TAA/TAG/TGA, the set the thresholds were simulated with)
   - A false discovery rate against decoys: each input transcript is shuffled keeping its
     dinucleotide counts (so length and base composition are kept), decoy ORFs are found in the
     shuffles with the LongOrfs settings and genetic code, and scored with the same model. Within
     strata of ORF length, transcript length and ORF type (complete or partial), an ORF's q-value
     is the lowest estimated FDR (decoys / targets at or above its frame contrast) of any cutoff
     that accepts it. Without homology support an ORF needs in-frame score > 0, a higher score
     than the shifted frames, and q <= `--max-qvalue` (0.05; `--precise` 0.01, `--sensitive` 0.2).
     Short ORFs, which are far more frequent by chance, therefore need more coding evidence, with
     no fixed length cutoff. The decoys are reproducible (`--decoy-seed`) and identical in
     TransDecoder (Perl).
   - Translated noncanonical ORFs (Ribo-seq ORFs) mostly lack codon-usage signal, so finding
     them needs translation evidence.
   - **Small inputs**: below ~1,000 transcripts the coding model and the decoy statistics learned
     from the input are unreliable. `--preset NAME` uses a coding model, start model and q-value
     thresholds learned from a reference transcriptome instead (bundled: see `--help`;
     `pytransdecoder build-preset` makes one from a Predict run on a related species).
4. **Start refinement**: an ORF that runs to the 5' end of its transcript (no start codon) is
   either a truncated (5'-partial) ORF or one that extends into the 5' UTR. For the segment
   before its first in-frame ATG, two scores are compared, both in bases' worth of typical coding
   signal (the training ORFs' median per-base frame contrast):
   - **coding score**: the segment's frame contrast (coding score in frame minus the better
     shifted frame) divided by that typical per-base value;
   - **UTR score**: `--start-contrast-alpha` (0.4) x segment length + `--start-kozak-weight` (17)
     x the ATG's Kozak score in bits (log2-odds of -6..-1, +4..+5 against in-frame ATGs inside
     the training ORFs).
   If the UTR score is higher, the ORF is trimmed to start at the ATG; otherwise it stays 5'
   partial. Each bit of Kozak score lets the ATG overrule 17 bases of coding signal. See
   [docs/start_codon_refinement.md](docs/start_codon_refinement.md) for a step-by-step
   walkthrough with examples and advice on choosing alpha. The Kozak model is written to the
   work directory as `start_refinement.kozak.tsv` (bits), `.kozak.counts.tsv` and a sequence
   logo, `.kozak.logo.svg`.
   **With homology hits** (`--retain-blastp-hits` and/or `--retain-pfam-hits`), the hits of the
   ORF's longest-ORF peptide decide first (PyTransDecoder only; turn off with
   `--no-homology-start-assist`):
   - BLASTP, trim: the best hit's subject N-terminus (subject start <= 3) aligns within 5 aa of the
     first ATG;
   - BLASTP, keep 5' partial: the alignment starts >= 5 aa upstream of the ATG and the subject
     continues further N-terminal (subject start > 1);
   - Pfam, keep 5' partial: a domain envelope starts upstream of the ATG.
   Where no rule fires, or a trim and a keep call conflict, the coding-vs-UTR rule above decides.
   Each rule is used only if its hits file was given. Per-ORF evidence and decisions are written
   to `start_refinement.homology.tsv` in the work directory, and ORFs trimmed by homology carry
   `start_homology=R2` in the GFF3. See [docs/start_codon_refinement.md](docs/start_codon_refinement.md).
5. **Output**: Generates final predictions in multiple formats

## Annotating a Genome from Transcript Structures

Transcripts reconstructed from RNA-seq (StringTie, LRAA, Cufflinks, PASA, or Trinity aligned to
the genome) can be annotated with coding regions on the genome: ORFs predicted on the transcript
sequences are propagated through each transcript's exon structure, giving gene / mRNA / exon /
CDS / five_prime_UTR / three_prime_UTR models in GFF3.

From a genome and a GTF (transcript sequences are extracted with gffread):

```bash
pyTransdecoder --genome genome.fa --gtf stringtie.gtf -O out
# -> out/stringtie.cDNA.fasta.transdecoder.genome.gff3 (and .genome.bed)
```

From transcript sequences and their structures (GTF, alignment GFF3 from GMAP/PASA with
`Target=`, or minimap2 SAM/BAM spliced alignments):

```bash
pyTransdecoder -t Trinity.fasta --structures Trinity.minimap2.bam -O out
# or, after LongOrfs + Predict:
pytransdecoder genome-map --orfs Trinity.fasta.transdecoder.gff3 --structures Trinity.minimap2.bam \
    -t Trinity.fasta -o Trinity.fasta.transdecoder.genome.gff3
```

Notes:
- ORFs on the opposite strand of a spliced transcript are not propagated; single-exon
  transcripts take the ORF's orientation.
- Alignment blocks carry transcript coordinates (GFF3 `Target`, SAM CIGAR), so transcripts with
  indels or clipped ends map base-accurately; an ORF whose CDS is not within the aligned part of
  its transcript is not propagated. Only primary SAM/BAM alignments are used.
- GTF transcripts with strand `.` are taken as `+`, as gffread extracts them; use
  `--unstranded-as minus` for sequences extracted with TransDecoder's
  `util/gtf_genome_to_cdna_fasta.pl`.
- Isoforms are grouped by the GTF `gene_id` (or the TransDecoder gene id, e.g. with
  `--gene-trans-map`), with the chromosome and strand appended.

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
