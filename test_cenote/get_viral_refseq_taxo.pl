use strict;
use warnings;
use Bio::SeqIO;
use Bio::Annotation::Collection;


my $seqio = Bio::SeqIO -> new(-file => "/data/umcg-tifn/DATABASES/viral_refseq_209/viral_refseq_209_genomic.gbff", -format => "genbank" );


print "genome_id\tgenome_desc\tgenome_length\thost\tlab_host\tisolation_source\ttaxonomy\n";

while (my $entry = $seqio -> next_seq) {

    my $genome_id = $entry -> accession_number . "." . $entry -> version;
    
    my $genome_desc = $entry -> desc;

    my $genome_length = $entry -> length;

    my $host = "";
    my $lab_host = "";
    my $isolation_source = "";

    for my $feat_object ($entry -> get_SeqFeatures) {

            my $prt = $feat_object -> primary_tag;
            if ($prt ne "source") { next }

            if ($feat_object -> has_tag("host"))             { ($host)             = $feat_object -> get_tag_values("host") }
            if ($feat_object -> has_tag("lab_host"))         { ($lab_host)         = $feat_object -> get_tag_values("lab_host") }
            if ($feat_object -> has_tag("isolation_source")) { ($isolation_source) = $feat_object -> get_tag_values("isolation_source") }

    }

    my $taxonomy = join ";", reverse($entry -> species -> classification);

    print $genome_id . "\t" . $genome_desc . "\t" . $genome_length . "\t" . $host . "\t" . $lab_host . "\t" . $isolation_source . "\t" . $taxonomy . "\n";

}
