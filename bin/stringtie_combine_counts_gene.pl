#!/usr/bin/perl

=head1

stringtie_combine_rawcounts.pl -

=head1 SYNOPSIS

    stringtie_combine_rawcounts.pl -l [sample names comma sep list] -r [reference gene file] -g [outgenefile] -f [out trait file] -t [out pheno file]

=head1 COMMAND-LINE OPTIONS
  ARGUMENTS
 -l sample names comma separated (required)
 -r reference gene file (required)
 -g path to genes outfile (required)
 -f path to traits outfile (required)
 -t out pheno file (required)

=head1 DESCRIPTION

=head1 AUTHOR

 Nicolas Morales (nm529@cornell.edu)

=cut

use strict;

use Getopt::Std;
use Data::Dumper;
use Carp qw /croak/ ;
use Pod::Usage;

our ($opt_l, $opt_r, $opt_g, $opt_f, $opt_t);

getopts('l:r:g:f:t:');

if (!$opt_l || !$opt_r || !$opt_g || !$opt_f || !$opt_t) {
    die "Must provide options -l (sample name list) -r (ref gene file) -g (genes outfile) -f (traits outfile) -t (out pheno file)\n";
}

my @samples = split ',', $opt_l;

my @result;
my %unique_genes;
my %data_hash;

my $accession = ['39671', 'TMEB419'];

my %ref_hash;
open(my $fh, '<', $opt_r)
    or die "Could not open file '$opt_r' $!";

while (my $row = <$fh>) {
    my @columns = split ' ', $row;
    my $ref_transcript = $columns[4];
    my @ref_val = split "=", $ref_transcript;
    $ref_hash{$ref_val[1]}++;
}
print STDERR Dumper \%ref_hash;

my $accession = 1;

foreach my $sample (@samples){
    my $filename = $sample.".txt";
    open(my $fh, '<', $filename)
	or die "Could not open file '$filename' $!";
    while (my $row = <$fh>) {
	my @data = split '\t',$row;
	#print STDERR Dumper \@data;
	my $chr = $data[0];
	my $start = $data[1];
	my $end = $data[2];
	my $gene_name = $data[3];
	my $transcript_name = $data[4];
	my $count = $data[5];
        #Remove semicolons from end of strings and quotes
        my $chop = chop($gene_name);
        my $chop = chop($transcript_name);
        $chop = chop($count);
        $chop = chop($count);
        $gene_name = substr($gene_name, 1, length($gene_name)-2);
        $transcript_name = substr($transcript_name, 1, length($transcript_name)-2);
        $count = substr($count, 1, length($count)-2);

        if (exists($ref_hash{$transcript_name})) {
            #my $full_name = $gene_name."_".$chr."_".$start."_".$end;
            $unique_genes{$transcript_name}++;
            $data_hash{$transcript_name}->{$sample} = $count;
        }
    }
}

print STDERR "Transcripts: ".scalar(keys %unique_genes)."\n";

my %unique_full_trait;
open(my $fh, '>', $opt_g)
    or die "Could not open file '$opt_g' $!";

foreach my $g (keys %unique_genes){
    my $accession_string = sprintf("%07d",$accession);
    my $gene_full_name = $g."|ST2:".$accession_string;
    $unique_full_trait{$g} = $gene_full_name;
    #if ($unique_genes{$g} != 56 && $unique_genes{$g} != 112 && $unique_genes{$g} != 168){
    #    die "ERROR: ".$g."  ".$unique_genes{$g}."\n";
    #}
    print $fh $gene_full_name."\n";
    $accession++;
}

my %trait_names;
my %final_hash;
my @traits;

foreach my $gene (keys %unique_full_trait){
    foreach my $sample (@samples){
        if ($data_hash{$gene}->{$sample}){
	    #print STDERR $sample."\n";
            my $trait_name;
            my $gene_trait = $unique_full_trait{$gene};
	    if (index($sample, 'FR1') != -1) {
                $trait_name = $gene_trait."||cass fibrous root|CASSTISS:0000011||cass end of day|CASSTIME:0000003||cass day 30|CASSTIME:0000103||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
            }
	    if (index($sample, 'FR3') != -1) {
		$trait_name = $gene_trait."||cass fibrous root|CASSTISS:0000011||cass end of day|CASSTIME:0000003||cass day 42|CASSTIME:0000115||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
	    }
	    if (index($sample, 'FR7') != -1) {
		$trait_name = $gene_trait."||cass fibrous root|CASSTISS:0000011||cass end of day|CASSTIME:0000003||cass day 60|CASSTIME:0000133||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
	    }
	    if (index($sample, 'SR1') != -1) {
		$trait_name = $gene_trait."||cass pre-storage root|CASSTISS:0000012||cass end of day|CASSTIME:0000003||cass day 30|CASSTIME:0000103||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
	    }
	    if (index($sample, 'SR3') != -1) {
		$trait_name = $gene_trait."||cass pre-storage root|CASSTISS:0000012||cass end of day|CASSTIME:0000003||cass day 42|CASSTIME:0000115||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
	    }
	    if (index($sample, 'SR7') != -1) {
		$trait_name = $gene_trait."||cass storage root|CASSTISS:0000013||cass end of day|CASSTIME:0000003||cass day 60|CASSTIME:0000133||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
	    }
	    if ((index($sample, 'SIL1') != -1)||(index($sample, 'SILb') != -1)) {
		$trait_name = $gene_trait."||cass sink leaf|CASSTISS:0000004||cass end of day|CASSTIME:0000003||cass day 30|CASSTIME:0000103||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
	    }
	    if (index($sample, 'SIL3') != -1) {
		$trait_name = $gene_trait."||cass sink leaf|CASSTISS:0000004||cass end of day|CASSTIME:0000003||cass day 42|CASSTIME:0000115||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
	    }
	    if (index($sample, 'SIL7') != -1) {
		$trait_name = $gene_trait."||cass sink leaf|CASSTISS:0000004||cass end of day|CASSTIME:0000003||cass day 60|CASSTIME:0000133||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
	    }
	    if (index($sample, 'SOL1') != -1) {
		$trait_name = $gene_trait."||cass source leaf|CASSTISS:0000005||cass end of day|CASSTIME:0000003||cass day 30|CASSTIME:0000103||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
	    }
	    if (index($sample, 'SOL3') != -1) {
		$trait_name = $gene_trait."||cass source leaf|CASSTISS:0000005||cass end of day|CASSTIME:0000003||cass day 42|CASSTIME:0000115||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
	    }
	    if (index($sample, 'SOL7') != -1) {
		$trait_name = $gene_trait."||cass source leaf|CASSTISS:0000005||cass end of day|CASSTIME:0000003||cass day 60|CASSTIME:0000133||TPM|CASSUNIT:0000005||IITA|CASSINST:0000006";
	    }
            
            if (!$trait_name){
                die "no trait name constructed\n";
            }
            $trait_names{$trait_name}++;
            $final_hash{$sample}->{$trait_name} = $data_hash{$gene}->{$sample};

        }
    }
}

open(my $fh, '>', $opt_f)
    or die "Could not open file '$opt_f' $!";

foreach my $t (keys %trait_names){
    print $fh $t."\n";
    push @traits, $t;
}

open(my $fh, '>', $opt_t)
    or die "Could not open file '$opt_t' $!";
print $fh '"studyYear","studyDbId","studyName","studyDesign","locationDbId","locationName","germplasmDbId","germplasmName","germplasmSynonyms","observationLevel","observationUnitDbId","observationUnitName","replicate","blockNumber","plotNumber",';
print $fh join ',', map { qq!"$_"! } @traits;
print $fh "\n";

foreach my $sample (@samples){
    print $fh '"2016","1623","RNA_Andreas","CRD","3","Ibadan","39671","TMEB419","","plot","X","'.$sample.'","1","1","1",';
    my $count = 1;
    foreach my $t (keys %trait_names){
        print $fh '"'.$final_hash{$sample}->{$t}.'"';
        if ($count < scalar(keys(%trait_names))){
            print $fh ",";
        }
        $count++;
    }
    print $fh "\n";
}

#foreach my $g (keys %unique_genes){
#    print $fh $g."\t";
#my $count = 1;
#foreach my $s (@samples){
#print $fh $data_hash{$g}->{$s};
#if ($count < scalar(@samples)){
#print $fh "\t";
#}
#$count++;
#}
#print $fh "\n";
#}

print STDERR "Complete!\n";

