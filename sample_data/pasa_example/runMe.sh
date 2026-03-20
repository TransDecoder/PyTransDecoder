#!/bin/bash -ve

if [ ! -e genome.fasta ]; then
    gunzip -c genome.fasta.gz > genome.fasta
fi

if [ ! -e pasa_assemblies.fasta ]; then
    gunzip -c pasa_assemblies.fasta.gz > pasa_assemblies.fasta
fi

if [ ! -e pasa_assemblies.GTF ]; then
    gunzip -c pasa_assemblies.GTF.gz > pasa_assemblies.GTF
fi

if [ ! -e pasa_assemblies_described.txt ]; then
    gunzip -c pasa_assemblies_described.txt.gz > pasa_assemblies_described.txt
fi


# get the gene-to-transcript relationships
cut -f2,3 pasa_assemblies_described.txt > pasa.gene_trans_map.txt


../../pyTransdecoder --gtf pasa_assemblies.GTF --genome genome.fasta --gene-trans-map pasa.gene_trans_map.txt -S


../../util/fasta_prot_checker.pl pasa_assemblies.cDNA.fasta.transdecoder.pep


echo "Done.  See pasa_assemblies.fasta.transdecoder.\*"


exit 0
