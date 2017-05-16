#!/usr/bin/perl

=head1

construct_expression_file.pl 

=head1 SYNOPSIS

	awk '$3 == "mRNA" {print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' mesculenta_305_v6.1.gene_exons.gff3 > mesculenta_305_v6.1_mRNA.txt

    construct_expression_file.pl  -i [infile created from stringtie_combine_counts_gene.pl combined_pheno.csv] -g [infile created from awk command above] -a [annotation infile] -o [outfile_prefix]
	
	perl construct_expression_file.pl -i /data2/Mueller/Nicolas/andreas_2017/stringtie_tpm/combined_pheno.csv -g /data2/Mueller/Nicolas/andreas_2017/mesculenta_305_v6.1_mRNA.txt -a /data2/Mueller/Nicolas/andreas_2017/Mesculenta_305_v6.1.annotation_info.txt -o expression_file

    perl construct_expression_file.pl -i ~/test_combined.csv -g /data2/Mueller/Nicolas/andreas_2017/mesculenta_305_v6.1_mRNA.txt -a /data2/Mueller/Nicolas/andreas_2017/Mesculenta_305_v6.1.annotation_info.txt -o expression_file

=head1 COMMAND-LINE OPTIONS
  ARGUMENTS
 -i input expression value file
 -g input gff file
 -a input annotation file
 -o output expression file

=head1 DESCRIPTION


=head1 AUTHOR

 Nicolas Morales (nm529@cornell.edu)

=cut

use strict;

use Getopt::Std;
use Data::Dumper;
use Carp qw /croak/ ;
use Pod::Usage;
use Statistics::Basic qw(:all);
use Statistics::R;
use Text::CSV;

our ($opt_i, $opt_g, $opt_a, $opt_o);

getopts('i:g:a:o:');

if (!$opt_i || !$opt_g || !$opt_a || !$opt_o){
	die "Must give -i, -g, -a, -o\n";
}

my %pacid_hash;
open(my $fh, '<', $opt_g)
    or die "Could not open file '$opt_g' $!";

while (my $row = <$fh>) {
	my @columns = split '\t', $row;
	my @atrr_columns = split ';', $columns[5];
	my @atrr_col0 = split '=', $atrr_columns[0];
	my @atrr_col1 = split '=', $atrr_columns[1];
	my @atrr_col2 = split '=', $atrr_columns[2];
	my $id = $atrr_col0[1];
	my $name = $atrr_col1[1];
	my $pacid = $atrr_col2[1];
	$pacid_hash{$pacid} = {'chr'=>$columns[0], 'type'=>$columns[1], 'start'=>$columns[2], 'end'=>$columns[3], 'strand'=>$columns[4], 'transid'=>$id, 'transname'=>$name};
}
#print STDERR Dumper \%pacid_hash;

my %annotation_hash;
open(my $fh, '<', $opt_a)
    or die "Could not open file '$opt_a' $!";

my $header = <$fh>;
while (my $row = <$fh>) {
	my @columns = split '\t', $row;
	#pacId	locusName	transcriptName	peptideName	Pfam	Panther	KOG	KEGG/ec	KO	GO	Best-hit-arabi-name	arabi-symbol	arabi-defline
	my $pacid_col = $columns[0];#
	my @pacid_terms = split ',',$pacid_col;
	if (scalar(@pacid_terms)>1){
		die "More than one pacid term on a single line\n";
	}
	my $gene_name_col = $columns[1];#
	my @gene_name_terms = split ',',$gene_name_col;
	if (scalar(@gene_name_terms)>1){
		die "More than one gene_name term on a single line\n";
	}
	my $transname_col = $columns[2];#
	my @transname_terms = split ',',$transname_col;
	if (scalar(@transname_terms)>1){
		die "More than one transname term on a single line\n";
	}
	my $pepname = $columns[3];#
	my $pfam = $columns[4];#
	my $panther = $columns[5];#
	my $kog = $columns[6];#
	my $keggec = $columns[7];#
	my $ko = $columns[8];#
	my $go = $columns[9];#
	my $arabi_name = $columns[10];#
	my $arabi_sym = $columns[11];#
	my $arabi_def = $columns[12];#
	chomp($arabi_def);

	my $pacid = $pacid_terms[0];
	my $transname = $transname_terms[0];
	my $gene_name = $gene_name_terms[0];

	my $chr = $pacid_hash{$pacid}->{chr};#
	my $type = $pacid_hash{$pacid}->{type};
	my $start = $pacid_hash{$pacid}->{start};#
	my $end = $pacid_hash{$pacid}->{end};#
	my $strand = $pacid_hash{$pacid}->{strand};
	my $transname_gff = $pacid_hash{$pacid}->{transname};#

	if ($transname_gff ne $transname){
		die "names not matching\n";
	}
	if (!$transname || !$start || !$end){
		die "missing info\n";
	}

	$annotation_hash{$transname.".v6.1"} = {chr=>$chr, type=> $type, start=>$start, end=>$end, strand=>$strand, pac_id=>$pacid, gene_name=>$gene_name, pep_name=>$pepname, pfam=>$pfam, panther=>$panther, kog=>$kog, keggec=>$keggec, ko=>$ko, go=>$go, arabi_name=>$arabi_name, arabi_sym=>$arabi_sym, arabi_def=>$arabi_def};
}
#print STDERR Dumper \%annotation_hash;

my $csv = Text::CSV->new({ sep_char => ',' });

open(my $fh, '<', $opt_i)
    or die "Could not open file '$opt_i' $!";

my $header_row = <$fh>;
my @columns;
if ($csv->parse($header_row)) {
    @columns = $csv->fields();
} else {
    die "Could not parse header.\n";
}
my $num_cols = scalar(@columns);

my @traits;
my $num_col_before = 15;
for my $col_num ($num_col_before .. $num_cols-1) {
    my $trait = $columns[$col_num];
    push @traits, $trait;
}
#print STDERR Dumper \@traits;

my %gene_sample_hash;
while ( my $row = <$fh> ){
    my @columns;
    if ($csv->parse($row)) {
        @columns = $csv->fields();
    } else {
        die "Could not parse row $row.\n";
    }
    my $i = 0;
    my $sample_name = $columns[11];
    my $accession = $columns[7];
    for my $col_num ($num_col_before .. $num_cols-1) {
        my $trait = $traits[$i];
	#print STDERR Dumper $trait;
        my @trait_parts = split /\|\|/, $trait;
	#print STDERR Dumper \@trait_parts;
        my @transname_parts = split /\|/, $trait_parts[0];
        my $transname = $transname_parts[0];
        my $sample = $sample_name.'||'.$accession.'||'.$trait_parts[1].'||'.$trait_parts[2].'||'.$trait_parts[3].'||'.$trait_parts[4].'||'.$trait_parts[5];
        $i++;
        my $val = $columns[$col_num];
        $gene_sample_hash{$transname}->{$sample} = $val;
    }
}
#print STDERR Dumper \%gene_sample_hash;

my @out_array;
my @out_header = ('TranscriptName', 'Chromosome', 'Start', 'End', 'Type', 'Strand', 'PacID', 'GeneName', 'PepName', 'PFam', 'Panther', 'KOG', 'KEGGEC', 'KO', 'GO', 'ArabiName', 'ArabiSym', 'ArabiDef');

my %unique_samples;
foreach my $transname (keys %gene_sample_hash) {
    my $sh = $gene_sample_hash{$transname};
    foreach my $s (sort keys %$sh) {
        $unique_samples{$s}++;
    }
}
print STDERR Dumper \%unique_samples;
foreach (sort keys %unique_samples){
    push @out_header, $_;
}
push @out_array, \@out_header;

foreach my $transname (sort keys %gene_sample_hash) {
    #print STDERR Dumper $transname;
    #print STDERR Dumper $annotation_hash{$transname};
    my @outline = ($transname, $annotation_hash{$transname}->{chr}, $annotation_hash{$transname}->{start}, $annotation_hash{$transname}->{end}, $annotation_hash{$transname}->{type}, $annotation_hash{$transname}->{strand}, $annotation_hash{$transname}->{pac_id}, $annotation_hash{$transname}->{gene_name}, $annotation_hash{$transname}->{pep_name}, $annotation_hash{$transname}->{pfam}, $annotation_hash{$transname}->{panther}, $annotation_hash{$transname}->{kog}, $annotation_hash{$transname}->{keggec}, $annotation_hash{$transname}->{ko}, $annotation_hash{$transname}->{go}, $annotation_hash{$transname}->{arabi_name}, $annotation_hash{$transname}->{arabi_sym}, $annotation_hash{$transname}->{arabi_def});
    my $sh = $gene_sample_hash{$transname};
    foreach my $s (sort keys %$sh) {
        push @outline, $sh->{$s};
    }
    push @out_array, \@outline;
}

open(my $fh, '>', $opt_o)
    or die "Could not open file '$opt_o' $!";

foreach (@out_array){
    print $fh join "\t", @$_;
    print $fh "\n";
}

print STDERR "Script complete!\n";

