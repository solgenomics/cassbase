#!/usr/bin/perl

=head1

transform_phenotypes_to_expression_atlas_lucy_index.pl 

=head1 SYNOPSIS

    transform_phenotypes_to_expression_atlas_lucy_index.pl  -i [infile]

=head1 COMMAND-LINE OPTIONS
  ARGUMENTS
 -i phenotype file downloaded directly from website

=head1 DESCRIPTION


=head1 AUTHOR

 Nicolas Morales (nm529@cornell.edu)

=cut

use strict;

use Getopt::Std;
use Data::Dumper;
use Carp qw /croak/ ;
use Pod::Usage;
use Spreadsheet::ParseExcel;
use Bio::Chado::Schema;
use Statistics::Basic qw(:all);

our ($opt_i, $opt_o, $opt_c);

getopts('i:o:c:');

if (!$opt_i || !$opt_o || !$opt_c) {
    pod2usage(-verbose => 2, -message => "Must provide options -i (input file) -o (lucy out file) -c (corr out file)\n");
}

my $parser   = Spreadsheet::ParseExcel->new();
my $excel_obj = $parser->parse($opt_i);

my $worksheet = ( $excel_obj->worksheets() )[0]; #support only one worksheet
my ( $row_min, $row_max ) = $worksheet->row_range();
my ( $col_min, $col_max ) = $worksheet->col_range();

my @data_out;

my @traits;
for my $col ( 14 .. $col_max) {
    my $multiterm_trait = $worksheet->get_cell(3,$col)->value();
    my @component_terms = split /\|\|/, $multiterm_trait;
    my $chebi_term = @component_terms[0];
    my $tissue_term = @component_terms[1];
    my $collection_term = @component_terms[2];
    my $age_term = @component_terms[3];
    $tissue_term =~ s/cass //;
    $collection_term =~ s/cass //;
    $age_term =~ s/cass //;
    push @traits, [$chebi_term, $tissue_term, $collection_term, $age_term];
}

my %intermed;
my %step2s;
for my $row ( 4 .. $row_max ) {

	my $accession_name = $worksheet->get_cell($row,7)->value();
    #print STDERR $accession_name."\n";
    if ($accession_name eq 'IITA-TMS-IBA011412'){
    
    for( my $i=0; $i<scalar(@traits); $i++) {
        my $trait_col = $i + 14;
        #print STDERR "$row $trait_col\n";
        my $value = '';
        if ($worksheet->get_cell($row,$trait_col)) {
            $value = $worksheet->get_cell($row,$trait_col)->value();
        }
        
        my $chebi_term = @traits[$i]->[0];
        my $tissue_term = @traits[$i]->[1];
        my $collection_term = @traits[$i]->[2];
        my $age_term = @traits[$i]->[3];
        
        my $temp_key = "$chebi_term, $accession_name, $tissue_term, $collection_term, $age_term";
        my $step2 = "$accession_name $tissue_term, $collection_term, $age_term";
        if (exists($intermed{$temp_key})) {
            my $values = $intermed{$temp_key}->[3];
            push @$values, $value;
            $intermed{$temp_key}->[3] = $values;
        } else {
            if (index($tissue_term, 'leaf') != -1) {
                $intermed{$temp_key} = [$chebi_term, "$accession_name Leaf, $collection_term, $age_term", $step2, [$value]];
            }
            if (index($tissue_term, 'root') != -1) {
                $intermed{$temp_key} = [$chebi_term, "$accession_name Root, $collection_term, $age_term", $step2, [$value]];
            }
            if (index($tissue_term, 'stem') != -1) {
                $intermed{$temp_key} = [$chebi_term, "$accession_name Stem, $collection_term, $age_term", $step2, [$value]];
            }
        }
        $step2s{$step2} = 1;
    }

    }

}
#print STDERR Dumper \%intermed;
#print STDERR Dumper keys %intermed;

my @step2s_sorted;
foreach (sort keys %step2s) {
    push @step2s_sorted, $_;
}

my %corr_out;
foreach (keys %intermed) {
    my $chebi_term = $intermed{$_}->[0];
    my $step1 = $intermed{$_}->[1];
    my $step2 = $intermed{$_}->[2];
    my $values = $intermed{$_}->[3];
    my $average = mean($values);
    my $stddev = stddev($values);
    
    push @data_out, [$chebi_term, $step1, $step2, $average->query, $stddev->query, $values];

    $corr_out{$chebi_term}->{$step2} = $average->query;
}
#print STDERR Dumper \@data_out;
#print STDERR Dumper \%corr_out;

my $fh;
open($fh, ">", $opt_o);
    print STDERR $opt_o."\n";
    foreach (@data_out) {
        my $values = $_->[5];
        print $fh "$_->[0]\t$_->[1]\t$_->[2]\t$_->[3]\t$_->[4]\t".join(',', @$values),"\n";
    }
close $fh;

open($fh, ">", $opt_c);
    print STDERR $opt_c."\n";
    print $fh "Metabolites\t", join("\t", @step2s_sorted), "\n";
    foreach my $chebi (sort keys %corr_out) {
        print $fh "$chebi\t";
        my $vals = $corr_out{$chebi};
        my $step = 1;
        foreach my $step2 (@step2s_sorted) {
            my $c = $vals->{$step2} ? $vals->{$step2} : '';
            print $fh "$c";
            if ($step < scalar(@step2s_sorted)) {
                print $fh "\t";
            }
            $step++;
        }
        print $fh "\n";
    }
close $fh;

print STDERR "Script Complete.\n";
