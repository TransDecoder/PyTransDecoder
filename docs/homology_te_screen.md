# Screening homology-retained ORFs for transposable-element origin

When BLASTP or Pfam hits are supplied (`--retain-blastp-hits`, `--retain-pfam-hits`), an ORF with a
hit is kept even when its coding signal alone would not keep it. Many transcripts contain
transposable-element (TE) sequence (Alu, LINE-1, endogenous retroviruses), and TE-encoded proteins,
or host proteins with a TE-derived segment, are in Swiss-Prot. Homology to them can therefore keep
short ORFs in TE-containing transcripts, mostly in UTRs and noncoding RNAs.

PyTransDecoder keeps these ORFs on purpose: TE-derived ORFs can be translated and are of interest, for
example, as a source of tumour antigens. If you do not want them, the steps below find them. The
example is human (GENCODE/UCSC `hg38` chromosome names) with transcripts mapped to the genome; adapt
the names to your data.

Inputs: the homology run (`transcripts.fa`, its `transcripts.fa.transdecoder_dir/`,
`transcripts.fa.transdecoder.pep` and, for step 4, `transcripts.fa.transdecoder.genome.gff3` from
`pytransdecoder genome-map`), its `blastp.outfmt6`, and the Swiss-Prot FASTA used for the search.
The commands need `awk`, `sort`, `comm`, `join` and, for step 4, [bedtools](https://bedtools.readthedocs.io).

## 1. ORFs kept only with homology evidence

Run Predict once more without the hits, on a copy of the LongOrfs output, and compare the ORF IDs
(LongOrfs IDs, the same in both runs):

```bash
mkdir abinitio && cd abinitio
ln -s ../transcripts.fa .
cp -r ../transcripts.fa.transdecoder_dir .
pytransdecoder predict -t transcripts.fa          # same options as the homology run, minus the hits
cd ..
ids() { grep '^>' "$1" | cut -c2- | cut -d' ' -f1 | sort -u; }
ids transcripts.fa.transdecoder.pep          > homology.ids
ids abinitio/transcripts.fa.transdecoder.pep > abinitio.ids
comm -23 homology.ids abinitio.ids > homology_only.ids
```

## 2. Best hit and its protein name

BLAST and DIAMOND list each query's hits best first, so the first line per ORF is its best hit.

```bash
awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' blastp.outfmt6 | sort -k1,1 > best_hit.tsv
# Swiss-Prot names: ">sp|Q6B4Z3|UTY_PANTR Histone demethylase UTY OS=Pan troglodytes ..."
zcat uniprot_sprot.fasta.gz | awk '/^>/ {id = substr($1, 2); $1 = ""; sub(/^ /, ""); print id"\t"$0}' |
  sort -k1,1 > sprot_names.tsv
join -t $'\t' homology_only.ids best_hit.tsv | sort -t $'\t' -k2,2 |
  join -t $'\t' -1 2 -2 1 - sprot_names.tsv |
  awk -F'\t' -v OFS='\t' '{print $2, $1, $3}' | sort -k1,1 > homology_only.best_hit.tsv
# columns: ORF id, best-hit subject, subject name and species
```

## 3. Screen by best-hit protein

```bash
grep -iE 'LINE-1|transpos|reverse transcriptase|endogenous retrovirus|retrovirus-related' \
  homology_only.best_hit.tsv | cut -f1 > te_by_name.ids
```

The keyword list is an example; check the most frequent best hits
(`cut -f2,3 homology_only.best_hit.tsv | sort | uniq -c | sort -k1,1nr | head -30`) and extend it.
A name screen misses host proteins with a TE-derived segment. In a human lymphoblastoid
transcriptome (GM12878, Trinity) searched against Swiss-Prot without human entries, the most frequent
best hit of homology-only ORFs was the C-terminus of chimpanzee UTY (`UTY_PANTR`), which is encoded
by an Alu element; its name gives no hint. The repeat-overlap screen (step 4) flagged 97% of those
that mapped to the genome.

## 4. Screen by repeat overlap (genome-mapped predictions)

Flag ORFs whose CDS lies mostly (here at least half its bases) in annotated transposable elements,
for example the UCSC RepeatMasker track:

```bash
# TE intervals: RepeatMasker classes SINE, LINE, LTR, DNA, RC, Retroposon
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
zcat rmsk.txt.gz | awk -F'\t' -v OFS='\t' '$12 ~ /^(SINE|LINE|LTR|DNA|RC|Retroposon)/ {print $6, $7, $8}' |
  sort -k1,1 -k2,2n | bedtools merge > te.bed

# CDS segments of the predictions, named by ORF id
awk -F'\t' -v OFS='\t' '$3 == "CDS" { match($9, /Parent=[^;]+/); print $1, $4 - 1, $5, substr($9, RSTART + 7, RLENGTH - 7) }' \
  transcripts.fa.transdecoder.genome.gff3 | sort -k1,1 -k2,2n > cds.bed

# fraction of each ORF's CDS bases in TEs
bedtools coverage -a cds.bed -b te.bed |
  awk -F'\t' -v OFS='\t' '{cov[$4] += $6; len[$4] += $7} END {for (o in len) print o, cov[o] / len[o]}' |
  sort -k1,1 > cds_te_fraction.tsv
join -t $'\t' homology_only.ids cds_te_fraction.tsv | awk -F'\t' '$2 >= 0.5 {print $1}' > te_by_overlap.ids
```

bedtools may warn about an "inconsistent naming convention" for unplaced contigs; this is harmless.
ORFs on transcripts that could not be mapped to the genome have no CDS segments and are not
screened by this step.

Without a genome, run RepeatMasker (or `dfamscan.pl` with Dfam) on the predicted CDSs
(`transcripts.fa.transdecoder.cds`) and use the masked fraction of each CDS instead.

## 5. Remove the flagged ORFs

```bash
sort -u te_by_name.ids te_by_overlap.ids > te_flagged.ids
# GFF3 / BED (transcript or genome coordinates): drop every line that names a flagged ORF id
# (-w: "X.p1" matches "X.p1", "cds.X.p1", "X.p1.exon1", but not "X.p10")
grep -v -w -F -f te_flagged.ids transcripts.fa.transdecoder.gff3 > filtered.gff3
grep -v -w -F -f te_flagged.ids transcripts.fa.transdecoder.genome.gff3 > filtered.genome.gff3
# proteins and CDSs (the ORF id is the first word of each header)
seqkit grep -v -f te_flagged.ids transcripts.fa.transdecoder.pep > filtered.pep
```

Whether to drop these ORFs depends on the question: for a reference annotation, removing them gives
fewer unannotated short ORFs; for studies of TE-derived translation, they are the candidates.
