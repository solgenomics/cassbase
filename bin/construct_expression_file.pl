#!/usr/bin/perl

=head1

construct_expression_file.pl 

=head1 SYNOPSIS

	awk '$3 == "mRNA" {print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' mesculenta_305_v6.1.gene_exons.gff3 > mesculenta_305_v6.1_mRNA.txt

    construct_expression_file.pl  -i [infile created from stringtie_combine_counts_gene.pl combined_pheno.csv] -a [annotation infile] -o [outfile_prefix]
	
	perl construct_expression_file.pl -i /data2/Mueller/Nicolas/andreas_2017/stringtie_tpm/combined_pheno.csv -a /data2/Mueller/Nicolas/andreas_2017/Mesculenta_305_v6.1.annotation_info.txt -o expression_file

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
	my @pepname_terms = split ',',$pepname;
	my $pfam = $columns[4];#
	my @pfam_terms = split ',',$pfam;
	my $panther = $columns[5];#
	my @panther_terms = split ',',$panther;
	my $kog = $columns[6];#
	my @kog_terms = split ',',$kog;
	my $keggec = $columns[7];#
	my @keggec_terms = split ',',$keggec;
	my $ko = $columns[8];#
	my @ko_terms = split ',',$ko;
	my $go = $columns[9];#
	my @go_terms = split ',', $go;
	my $arabi_name = $columns[10];#
	my @arabi_name_terms = split ',', $arabi_name;
	my $arabi_sym = $columns[11];#
	my @arabi_sym_terms = split ',', $arabi_sym;
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

	my @glines;

	#ID,NAME,SYNONYM,STARTBASE,ENDBASE,GENE-COMMENT,FUNCTION,PRODUCT-TYPE,EC,FUNCTION-COMMENT,DBLINK,GO,INTRON
	my $id_line = "ID\t$transname\n";
	push @glines, $id_line;
	my $name_line = "NAME\t$transname\n";
	push @glines, $name_line;
	foreach (@pepname_terms){
		push @glines, "SYNONYM\t$_\n"; 
	}
	my $start_line = "STARTBASE\t$start\n";
	push @glines, $start_line;
	my $end_line = "ENDBASE\t$end\n";
	push @glines, $end_line;
	my $function_line = $arabi_def ? "FUNCTION\t$arabi_def\n" : '';
	if ($function_line){ push @glines, $function_line; }
	my $product_type_line = "PRODUCT-TYPE\tP\n";
	push @glines, $product_type_line;
	foreach (@keggec_terms){
		push @glines, "EC\t$_\n"; 
	}
	foreach (@go_terms){
		push @glines, "GO\t$_\n"; 
	}
	my $gcomment_line1 = $gene_name ? "GENE-COMMENT\tGENENAME:$gene_name\n" : '';
	if ($gcomment_line1){ push @glines, $gcomment_line1; }
	foreach (@pfam_terms){
		push @glines, "GENE-COMMENT\tPFAM:$_\n"; 
	}
	foreach (@panther_terms){
		push @glines, "GENE-COMMENT\tPANTHER:$_\n";
	}
	foreach (@kog_terms){
		push @glines, "GENE-COMMENT\tKOG:$_\n";
	}
	foreach (@ko_terms){
		push @glines, "GENE-COMMENT\tKO:$_\n";
	}
	foreach (@arabi_name_terms){
		push @glines, "FUNCTION-COMMENT\tARABINAME:$_\n";
	}
	foreach (@arabi_sym_terms){
		push @glines, "FUNCTION-COMMENT\tARABISYM:$_\n";
	}
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

