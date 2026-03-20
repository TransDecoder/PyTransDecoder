#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (cleanme.pl 
                        runMe.sh
stringtie_merged.gff3
stringtie_merged.gtf
chrX.fa
chrX.fa.fai
chrX.fa.gz
Makefile

                                         );


my %keep = map { + $_ => 1 } @files_to_keep;


foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}

`rm -rf ./stringtie_merged.transcripts.fasta.transdecoder_dir`;
`rm -rf ./stringtie_merged.transcripts.fasta.transdecoder_dir.__checkpoints`;
`rm -rf ./stringtie_merged.transcripts.fasta.transdecoder_dir.__checkpoints_longorfs`;
`rm -rf ./stringtie_merged.cDNA.fasta.transdecoder_dir`;
`rm -rf ./stringtie_merged.cDNA.fasta.transdecoder_dir.__checkpoints`;
`rm -rf ./stringtie_merged.cDNA.fasta.transdecoder_dir.__checkpoints_longorfs`;

foreach my $generated (
    qw(
        stringtie_merged.cDNA.fasta
        stringtie_merged.cDNA.fasta.transdecoder.bed
        stringtie_merged.cDNA.fasta.transdecoder.cds
        stringtie_merged.cDNA.fasta.transdecoder.genome.gff3
        stringtie_merged.cDNA.fasta.transdecoder.gff3
        stringtie_merged.cDNA.fasta.transdecoder.pep
    )
) {
    unlink($generated) if -e $generated;
}

exit(0);
