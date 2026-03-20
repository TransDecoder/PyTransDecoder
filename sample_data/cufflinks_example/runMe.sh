#!/bin/bash

set -ev

export PERL_HASH_SEED=0

if [ ! -e test.genome.fasta ]; then
    gunzip -c test.genome.fasta.gz > test.genome.fasta
fi


if [ ! -e transcripts.gtf ]; then
    gunzip -c transcripts.gtf.gz > transcripts.gtf
fi

if [ ! -e mini_Pfam-A.hmm ]; then
    gunzip -c mini_Pfam-A.hmm.gz > mini_Pfam-A.hmm
fi

if [ ! -e mini_sprot.db.pep ]; then
    gunzip -c mini_sprot.db.pep.gz > mini_sprot.db.pep
fi

cmd=""
## Run the primary PyTransdecoder workflow.
if [ "$1" == "" ]; then
    cmd="../../pyTransdecoder --genome test.genome.fasta --gtf transcripts.gtf"
else
    cmd="../../pyTransdecoder --genome test.genome.fasta --gtf transcripts.gtf \
        --blast-search-pep mini_sprot.db.pep \
        --pfam-search-db mini_Pfam-A.hmm \
        -v"
fi

eval $cmd

## make bed files for viewing with GenomeView

# covert cufflinks gtf to bed
../../util/gtf_to_bed.pl transcripts.gtf > transcripts.bed

# convert the genome-based gene-gff3 file to bed
../../util/gff3_file_to_bed.pl transcripts.cDNA.fasta.transdecoder.genome.gff3 > transcripts.cDNA.fasta.transdecoder.genome.bed


# ensure no fatal problems w/ pep file
../../util/fasta_prot_checker.pl transcripts.cDNA.fasta.transdecoder.pep

# Done!  Coding region genome annotations provided as: transcripts.cDNA.fasta.transdecoder.genome.*


exit 0
