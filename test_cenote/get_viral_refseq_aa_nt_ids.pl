use strict;
use warnings;
use Bio::SeqIO;
use Bio::Annotation::Collection;
use Data::Dumper;


my $seqio = Bio::SeqIO -> new(-file => "/data/umcg-tifn/DATABASES/viral_refseq_209/viral_refseq_209_protein.gpff", -format => "genbank" );


print "genome_id\tprotein_id\tcoded_by\tprotein_desc\n";

while (my $entry = $seqio -> next_seq) {

    my $protein_id = $entry -> accession_number . "." . $entry -> version;


    my $protein_desc = $entry -> desc;


    my $genome_id = "";

    my @dblinks = $entry -> annotation -> get_Annotations("dblink");

    foreach my $x (@dblinks) {

        if ($x->{"database"} eq "REFSEQ") { $genome_id = $x->{"primary_id"} . "." . $x->{"version"} }

    }


    my $coded_by = "";

    for my $f ($entry -> all_SeqFeatures) {

        if ($f -> primary_tag eq "CDS") {

            my @arr = $f -> each_tag_value("coded_by");

            $coded_by = join ";", @arr;

        }

    }


    print $genome_id . "\t" . $protein_id . "\t" . $coded_by . "\t" . $protein_desc . "\n";

}
