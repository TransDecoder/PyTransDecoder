#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (cleanme.pl 
                        runMe.sh
                        test.genome.fasta
                        test.genome.fasta.fai
                        test.genome.fasta.gz
                        test.tophat.sam.gz
                        transcripts.gtf
                        transcripts.gtf.gz
                        Makefile
                        mini_Pfam-A.hmm
                        mini_Pfam-A.hmm.gz
                        mini_sprot.db.pep
                        mini_sprot.db.pep.gz
                                         );


my %keep = map { + $_ => 1 } @files_to_keep;


foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}

`rm -rf ./transcripts.fasta.transdecoder_dir/`;
`rm -rf ./transcripts.fasta.transdecoder_dir.__checkpoints`;
`rm -rf ./transcripts.fasta.transdecoder_dir.__checkpoints_longorfs/`;
`rm -rf ./transcripts.cDNA.fasta.transdecoder_dir/`;
`rm -rf ./transcripts.cDNA.fasta.transdecoder_dir.__checkpoints`;
`rm -rf ./transcripts.cDNA.fasta.transdecoder_dir.__checkpoints_longorfs/`;

foreach my $generated (
    qw(
        blastp.outfmt6
        mini_Pfam-A.hmm.h3f
        mini_Pfam-A.hmm.h3i
        mini_Pfam-A.hmm.h3m
        mini_Pfam-A.hmm.h3p
        mini_sprot.db.pep.pdb
        mini_sprot.db.pep.phr
        mini_sprot.db.pep.pin
        mini_sprot.db.pep.pot
        mini_sprot.db.pep.psq
        mini_sprot.db.pep.ptf
        mini_sprot.db.pep.pto
        pfam.domtblout
        transcripts.bed
        transcripts.cDNA.fasta
        transcripts.cDNA.fasta.transdecoder.bed
        transcripts.cDNA.fasta.transdecoder.cds
        transcripts.cDNA.fasta.transdecoder.genome.bed
        transcripts.cDNA.fasta.transdecoder.genome.gff3
        transcripts.cDNA.fasta.transdecoder.gff3
        transcripts.cDNA.fasta.transdecoder.pep
        transcripts.fasta
        transcripts.fasta.transdecoder.bed
        transcripts.fasta.transdecoder.cds
        transcripts.fasta.transdecoder.gff3
        transcripts.fasta.transdecoder.pep
        transcripts.gff3
    )
) {
    unlink($generated) if -e $generated;
}

exit(0);
