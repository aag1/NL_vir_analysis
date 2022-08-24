#!/usr/bin/perl
use strict;
use warnings;



### DTR detection is based on https://github.com/simroux/VirSorter/blob/master/Scripts/Step_1_contigs_cleaning_and_gene_prediction.pl



### INPUT PARAMETERS
my $fastaF = $ARGV[0];
my $outF = $ARGV[1];
my $min_overlap = $ARGV[2];



### READING INPUT FASTA
open FILE1, '<', $fastaF or die "Failed to open $fastaF: $!\n";

my %seq_base;
my $id_seq = "";

while (my $line = <FILE1>) {

	chomp($line);

	if ($line =~ /^>(\S*)/) { $id_seq = $1 }
	else { $seq_base{$id_seq} .= $line }

}

close FILE1 or die "Failed to close $fastaF: $!\n";



### DETECTION OF SEQUENCES WITH TERMINAL REPEATS
open FILE2, '>', $outF or die "Failed to open $outF: $!\n";
print FILE2 "seq_id\tterm_rep_type\tterm_rep_len\n";

for my $k (keys %seq_base) {

    my $s = $seq_base{$k};

    my $prefix = substr($s, 0, $min_overlap);


    # DTR
    if ($s =~ /.+($prefix.*?)$/) {

        my $suffix = $1;

        my $l = length($suffix);

        my $test = substr($s, 0, $l);

        if ($suffix eq $test) {

            print FILE2 "$k\tDTR\t$l\n";

        }
    }


    # ITR
    my $r = reverse $s;
    $r =~ tr/ATGCatgc/TACGtacg/;

    if ($prefix eq substr($r, 0, $min_overlap)) {

        my @S = split '', $s;
        my @R = split '', $r;

        my $i = 0;

        while ($i < @S) {

            if ($S[$i] eq $R[$i]) { $i++ } else { last }

        }

        print FILE2 "$k\tITR\t$i\n";

    }

}

close FILE2 or die "Failed to close $outF: $!\n";

# Note that for $s = 'AAbcAAb....AAbcAAb' and $min_overlap = 2, suffix AAb will be detected.
