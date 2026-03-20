#!/bin/bash -ve

export PERL_HASH_SEED=0


if [ ! -f "chrX.fa" ]; then
    gunzip -c chrX.fa.gz > chrX.fa
fi


../../pyTransdecoder --genome chrX.fa --gtf stringtie_merged.gtf -S $ARGS


exit 0
