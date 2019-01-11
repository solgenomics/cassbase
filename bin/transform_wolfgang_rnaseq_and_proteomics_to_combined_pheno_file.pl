#!/usr/bin/perl

=head1

transform_wolfgang_rnaseq_and_proteomics_to_combined_pheno_file.pl 

=head1 SYNOPSIS

    transform_wolfgang_rnaseq_and_proteomics_to_combined_pheno_file.pl  -i [infile] -o [outfile]
    
    perl bin/transform_wolfgang_rnaseq_and_proteomics_to_combined_pheno_file.pl -r /home/vagrant/Downloads/rnaseq.csv -i /home/vagrant/Downloads/proteomics1.csv -j /home/vagrant/Downloads/proteomics2.csv -k /home/vagrant/Downloads/proteomics3.csv -o /home/vagrant/cxgn/cassbase/bin/pheno.csv

=head1 COMMAND-LINE OPTIONS
  ARGUMENTS
 -r rnaseq file from wolfgang for CASS 2018 rnaseq DESeq2NormaledCounts
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

our ($opt_r, $opt_i, $opt_j, $opt_k, $opt_o);

getopts('r:i:j:k:o:');

if (!$opt_r || !$opt_i || !$opt_j || !$opt_k || !$opt_o) {
    die "Must provide options -r (RNAseq input file) -i (proteomics input file1) -j (proteomics input file2) -k (proteomics input file3) -o (out file)\n";
}

my %genotype_id_hash = (
    'G02' => 'Z010116',
    'G04' => 'I082418',
    'G06' => 'SIMONYE',
    'G07' => 'KALESO',
    'G11' => 'I011231',
    'G12' => 'I071378',
    'G14' => 'I090485',
    'G15' => 'I990304',
    'G16' => 'I051599',
    'G18' => 'I011663',
    'G20' => 'I090576',
    'G22' => 'I050128',
    'EXP1' => 'I050128',
    'EXP3' => 'I050128',
    'EXP2' => {
        '2' => 'Z010116',
        '6' => 'SIMONYE'
    }
);

my %tissue_type_hash = (
    'p' => 'Phloem',
    'c' => 'Cambium',
    'x' => 'Xylem',
    'EXP1' => ['Phloem', 'Cambium', 'Xylem'],
    'EXP2' => ['Xylem'],
    'EXP3' => ['Phloem', 'Outer Xylem', 'Central Xylem']
);

my $csv = Text::CSV->new({ sep_char => ',' });

open(my $fh0, '<', $opt_r)
    or die "Could not open file '$opt_r' $!";

my $header_row0 = <$fh0>;
my @header_columns0;
if ($csv->parse($header_row0)) {
    @header_columns0 = $csv->fields();
} else {
    die "Could Not Parse Line: $header_row0\n";
}
my $col_max0 = scalar(@header_columns0)-1;

my @outfile_header_row = ("studyYear","programDbId","programName","programDescription","studyDbId","studyName","studyDescription","studyDesign","plotWidth","plotLength","fieldSize","fieldTrialIsPlannedToBeGenotyped","fieldTrialIsPlannedToCross","plantingDate","harvestDate","locationDbId","locationName","germplasmDbId","germplasmName","germplasmSynonyms","observationLevel","observationUnitDbId","observationUnitName","replicate","blockNumber","plotNumber","rowNumber","colNumber","entryType","plantNumber","plantedSeedlotStockDbId","plantedSeedlotStockUniquename","plantedSeedlotCurrentCount","plantedSeedlotCurrentWeightGram","plantedSeedlotBoxName","plantedSeedlotTransactionCount","plantedSeedlotTransactionWeight","plantedSeedlotTransactionDescription","availableGermplasmSeedlotUniquenames");

my %data_hash;
my %seen_composed_traits;

while ( my $row = <$fh0> ){
    my @columns;
    if ($csv->parse($row)) {
        @columns = $csv->fields();
    } else {
        die "Could not parse line $row\n";
    }

    my $gene_name = $columns[0];
    my $gene_annotation = $columns[1];

    foreach my $i (2 .. $col_max0) {
        my $sample = $header_columns0[$i];
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

open(my $fh1, '<', $opt_i)
    or die "Could not open file '$opt_i' $!";

print STDERR $opt_i."\n";

my $header_row1 = <$fh1>;
my @header_columns1;
if ($csv->parse($header_row1)) {
    @header_columns1 = $csv->fields();
} else {
    die "Could Not Parse Line: $header_row1\n";
}
my $col_max1 = scalar(@header_columns1)-1;

while ( my $row = <$fh1> ){
    my @columns;
    if ($csv->parse($row)) {
        @columns = $csv->fields();
    } else {
        die "Could not parse line $row\n";
    }

    my $gene_name = $columns[0];
    my $gene_annotation = $columns[1];
    my $accession_name = $genotype_id_hash{'EXP1'};

    foreach my $i (4 .. $col_max1) {
        my $sample = $header_columns1[$i];
        my ($tissue_type, $rep) = split '\|', $sample;
        my $value = $columns[$i];
        
        #print STDERR Dumper [$accession_name, $tissue_type, $rep, $gene_name, $value];
        my $composed_trait = $gene_name."|X:0000000||".$tissue_type."|X:0000000||NA|X:0000000||NA|X:0000000||Log2Intensity|X:0000000||ER|X:0000000|XX:0000000||||".$gene_annotation;
        $seen_composed_traits{$composed_trait}++;
        $data_hash{$accession_name}->{$rep}->{$composed_trait} = $value;
    }
}

open(my $fh2, '<', $opt_j)
    or die "Could not open file '$opt_j' $!";

print STDERR $opt_j."\n";

my $header_row2 = <$fh2>;
my @header_columns2;
if ($csv->parse($header_row2)) {
    @header_columns2 = $csv->fields();
} else {
    die "Could Not Parse Line: $header_row2\n";
}
my $col_max2 = scalar(@header_columns2)-1;

while ( my $row = <$fh2> ){
    my @columns;
    if ($csv->parse($row)) {
        @columns = $csv->fields();
    } else {
        die "Could not parse line $row\n";
    }

    my $gene_name = $columns[0];
    my $gene_annotation = $columns[1];
    my $tissue_type = 'Xylem';

    foreach my $i (4 .. $col_max2) {
        my $sample = $header_columns2[$i];
        my ($accession_name, $rep) = split '\|', $sample;
        my $value = $columns[$i];
        
        #print STDERR Dumper [$accession_name, $tissue_type, $rep, $gene_name, $value];
        my $composed_trait = $gene_name."|X:0000000||".$tissue_type."|X:0000000||NA|X:0000000||NA|X:0000000||Log2Intensity|X:0000000||ER|X:0000000|XX:0000000||||".$gene_annotation;
        $seen_composed_traits{$composed_trait}++;
        $data_hash{$accession_name}->{$rep}->{$composed_trait} = $value;
    }
}

open(my $fh3, '<', $opt_k)
    or die "Could not open file '$opt_k' $!";

print STDERR $opt_k."\n";

my $header_row3 = <$fh3>;
my @header_columns3;
if ($csv->parse($header_row3)) {
    @header_columns3 = $csv->fields();
} else {
    die "Could Not Parse Line: $header_row3\n";
}
my $col_max3 = scalar(@header_columns3)-1;

while ( my $row = <$fh3> ){
    my @columns;
    if ($csv->parse($row)) {
        @columns = $csv->fields();
    } else {
        die "Could not parse line $row\n";
    }

    my $gene_name = $columns[0];
    my $gene_annotation = $columns[1];
    my $accession_name = $genotype_id_hash{'EXP3'};

    foreach my $i (4 .. $col_max3) {
        my $sample = $header_columns3[$i];
        my ($tissue_type, $rep) = split '\|', $sample;
        my $value = $columns[$i];

        #print STDERR Dumper [$accession_name, $tissue_type, $rep, $gene_name, $value];
        my $composed_trait = $gene_name."|X:0000000||".$tissue_type."|X:0000000||NA|X:0000000||NA|X:0000000||Log2Intensity|X:0000000||ER|X:0000000|XX:0000000||||".$gene_annotation;
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
            print $fh_o '"2018","NA","CASS","CASS","NA","CASS_2018_DEglobal_Proteomics_Combined","RNAseq 12 genotypes combined with Proteomics Experiments. Experiment 1 PCX, Experiment 2 Carotenoids, Experiment 3 PX1X2","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",'.$accession.',"NA","tissue_sample","NA",'.$accession."_".$rep.',"NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"';
            
            foreach (@traits) {
                my $value = $composed_trait_hash->{$_};
                print $fh_o ',"'.$value.'"';
            }
            print $fh_o "\n";
        }
    }

print STDERR "Script complete!\n";
