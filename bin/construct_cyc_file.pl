#!/usr/bin/perl

=head1

construct_cyc_file.pl 

=head1 SYNOPSIS

	awk '$3 == "mRNA" {print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' mesculenta_305_v6.1.gene_exons.gff3 > mesculenta_305_v6.1_mRNA.txt

    construct_cyc_file.pl  -i [infile created from gff above] -a [annotation infile] -o [outfile_prefix]
	
	perl construct_cyc_file.pl -i ~/Downloads/Mesculenta/v6.1/annotation/Mesculenta_305_v6.1_mRNA.txt -a ~/Downloads/Mesculenta/v6.1/annotation/Mesculenta_305_v6.1.annotation_info.txt -o cyc_file

=head1 COMMAND-LINE OPTIONS
  ARGUMENTS
 -i input gff file
 -a input annotation file
 -o output cyc file prefix. will create a file for each chr

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

our ($opt_i, $opt_a, $opt_o);

getopts('i:a:o:');

if (!$opt_i || !$opt_a || !$opt_o){
	die "Must give -i, -a, -o\n";
}

my %pacid_hash;
open(my $fh, '<', $opt_i)
    or die "Could not open file '$opt_i' $!";

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

my @out_array;
my %chrs;
open(my $fh, '<', $opt_a)
    or die "Could not open file '$opt_a' $!";

my $header = <$fh>;
while (my $row = <$fh>) {
	my @columns = split '\t', $row;
	#pacId	locusName	transcriptName	peptideName	Pfam	Panther	KOG	KEGG/ec	KO	GO	Best-hit-arabi-name	arabi-symbol	arabi-defline
	my $pacid = $columns[0];#
	my $gene_name = $columns[1];#
	my $transname = $columns[2];#
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

	my @glines;

	#ID,NAME,SYNONYM,STARTBASE,ENDBASE,GENE-COMMENT,FUNCTION,PRODUCT-TYPE,EC,FUNCTION-COMMENT,DBLINK,GO,INTRON
	my $id_line = "ID\t$transname\n";
	push @glines, $id_line;
	my $name_line = "NAME\t$transname\n";
	push @glines, $name_line;
	my $syn_line = $pepname ? "SYNONYM\t$pepname\n" : '';
	if ($syn_line){ push @glines, $syn_line; }
	my $start_line = "STARTBASE\t$start\n";
	push @glines, $start_line;
	my $end_line = "ENDBASE\t$end\n";
	push @glines, $end_line;
	my $function_line = $arabi_def ? "FUNCTION\t$arabi_def\n" : '';
	if ($function_line){ push @glines, $function_line; }
	my $product_type_line = "PRODUCT-TYPE\tP\n";
	push @glines, $product_type_line;
	my $ec_line = $keggec ? "EC\t$keggec\n" : '';
	if ($ec_line){ push @glines, $ec_line; }
	my $go_line = $go ? "GO\t$go\n" : '';
	if ($go_line){ push @glines, $go_line; }
	my $gcomment_line1 = $gene_name ? "GENE-COMMENT\tGENENAME:$gene_name\n" : '';
	if ($gcomment_line1){ push @glines, $gcomment_line1; }
	my $gcomment_line2 = $pfam ? "GENE-COMMENT\tPFAM:$pfam\n" : '';
	if ($gcomment_line2){ push @glines, $gcomment_line2; }
	my $gcomment_line3 = $panther ? "GENE-COMMENT\tPANTHER:$panther\n" : '';
	if ($gcomment_line3){ push @glines, $gcomment_line3; }
	my $gcomment_line4 = $kog ? "GENE-COMMENT\tKOG:$kog\n" : '';
	if ($gcomment_line4){ push @glines, $gcomment_line4; }
	my $gcomment_line5 = $ko ? "GENE-COMMENT\tKO:$ko\n" : '';
	if ($gcomment_line5){ push @glines, $gcomment_line5; }
	my $fcomment_line1 = $arabi_name ? "FUNCTION-COMMENT\tARABINAME:$arabi_name\n" : '';
	if ($fcomment_line1){ push @glines, $fcomment_line1; }
	my $fcomment_line2 = $arabi_sym ? "FUNCTION-COMMENT\tARABISYM:$arabi_sym\n" : '';
	if ($fcomment_line2){ push @glines, $fcomment_line2; }
	push @glines, "//\n";
	push @out_array, { $chr => \@glines };
	$chrs{$chr}++;
}
#print STDERR Dumper \%pacid_hash;

foreach my $chr (keys %chrs){
	my $file_name = $opt_o."_".$chr.".pf";
	print STDERR $file_name."\n";

	open(my $fh, '>', $file_name)
    	or die "Could not open file $!";

	foreach my $chrline (@out_array){
		if (exists($chrline->{$chr})){
			my $lines = $chrline->{$chr};
			foreach (@$lines){
				print $fh $_;
			}
		}
	}

}

print STDERR "Script complete!\n";

