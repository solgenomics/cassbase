#!/usr/bin/perl

=head1

transform_wolfgang_3_proteomics_files_to_pheno_file.pl 

=head1 SYNOPSIS

    transform_wolfgang_3_proteomics_files_to_pheno_file.pl  -i [infile1] -j [infile2] -k [infile3] -o [outfile]
    
    perl bin/transform_wolfgang_3_proteomics_files_to_pheno_file.pl -i /home/vagrant/Downloads/proteomics1.csv -j /home/vagrant/Downloads/proteomics2.csv -k /home/vagrant/Downloads/proteomics3.csv -o /home/vagrant/cxgn/cassbase/bin/pheno.csv

=head1 COMMAND-LINE OPTIONS
  ARGUMENTS
 -i proteomics data from experiment 1 CASS 2018
 -j proteomics data from experiment 2 CASS 2018
 -k proteomics data from experiment 3 CASS 2018
 -o output lucy.tsv file that can be indexed using TEA script

=head1 DESCRIPTION


=head1 AUTHOR

 Nicolas Morales (nm529@cornell.edu)

=cut

use strict;

use Getopt::Std;
use Data::Dumper;
use Carp qw /croak/ ;
use Pod::Usage;
use Text::CSV;
use JSON;

our ($opt_i, $opt_o);

getopts('i:o:');

if (!$opt_i || !$opt_o) {
    die "Must provide options -i (input file) -o (out file)\n";
}

my %genotype_id_hash = (
    'EXP1' => 'I050128',
    'EXP3' => 'I050128',
    'EXP2' => {
        '2' => 'Z010116',
        '6' => 'SIMONYE'
    }
);

my %tissue_type_hash = (
    'EXP1' => ['Phloem', 'Cambium', 'Xylem'],
    'EXP2' => ['Xylem'],
    'EXP3' => ['Phloem', 'Outer Xylem', 'Central Xylem']
);

my $csv = Text::CSV->new({ sep_char => ',' });

open(my $fh1, '<', $opt_i)
    or die "Could not open file '$opt_i' $!";

my $header_row = <$fh1>;
my @header_columns;
if ($csv->parse($header_row)) {
    @header_columns = $csv->fields();
} else {
    die "Could Not Parse Line: $header_row\n";
}
my $col_max = scalar(@header_columns)-1;

my @outfile_header_row = ("studyYear","programDbId","programName","programDescription","studyDbId","studyName","studyDescription","studyDesign","plotWidth","plotLength","fieldSize","fieldTrialIsPlannedToBeGenotyped","fieldTrialIsPlannedToCross","plantingDate","harvestDate","locationDbId","locationName","germplasmDbId","germplasmName","germplasmSynonyms","observationLevel","observationUnitDbId","observationUnitName","replicate","blockNumber","plotNumber","rowNumber","colNumber","entryType","plantNumber","plantedSeedlotStockDbId","plantedSeedlotStockUniquename","plantedSeedlotCurrentCount","plantedSeedlotCurrentWeightGram","plantedSeedlotBoxName","plantedSeedlotTransactionCount","plantedSeedlotTransactionWeight","plantedSeedlotTransactionDescription","availableGermplasmSeedlotUniquenames");

my %data_hash;
my %seen_composed_traits;

while ( my $row = <$fh1> ){
    my @columns;
    if ($csv->parse($row)) {
        @columns = $csv->fields();
    } else {
        die "Could not parse line $row\n";
    }

    my $gene_name = $columns[0];
    my $gene_annotation = $columns[1];

    foreach my $i (2 .. $col_max) {
        my $sample = $header_columns[$i];
        my $genotype_id = substr($sample, 0, 3);
        my $tissue_type = substr($sample, 3, 1);
        my $rep = substr($sample, 4, 1);
        my $value = $columns[$i];
        
        my $accession_name = $genotype_id_hash{$genotype_id};
        if (!$accession_name) {
            die "EXITED: $genotype_id NOT FOUND\n";
        }

        my $tissue = $tissue_type_hash{$tissue_type};
        if (!$tissue) {
            die "EXITED: $tissue_type NOT FOUND\n";
        }

        my $composed_trait = $gene_name."|X:0000000||".$tissue."|X:0000000||NA|X:0000000||NA|X:0000000||DESeq2NormaledCounts|X:0000000||ER|X:0000000|XX:0000000||||".$gene_annotation;
        $seen_composed_traits{$composed_trait}++;
        $data_hash{$accession_name}->{$rep}->{$composed_trait} = $value;
    }
}

my @traits = keys %seen_composed_traits;
push @outfile_header_row, @traits;

open(my $fh_o, ">", $opt_o) || die("\nERROR:\n");
    print STDERR $opt_o."\n";

    print $fh_o join (',', @outfile_header_row), "\n";

    while (my ($accession, $rep_hash) = each %data_hash) {
        while (my ($rep, $composed_trait_hash) = each %$rep_hash) {
            print $fh_o '"2018","NA","CASS","CASS","NA","CASS_2018_DEglobal","RNAseq 12 genotypes","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",'.$accession.',"NA","tissue_sample","NA",'.$accession."_".$rep.',"NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"';
            
            foreach (@traits) {
                my $value = $composed_trait_hash->{$_};
                print $fh_o ',"'.$value.'"';
            }
            print $fh_o "\n";
        }
    }

print STDERR "Script complete!\n";
