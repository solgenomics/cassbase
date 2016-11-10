#!/usr/bin/perl

=head1

transform_phenotypes_to_expression_atlas_lucy_index.pl 

=head1 SYNOPSIS

    transform_phenotypes_to_expression_atlas_lucy_index.pl  -i [infile]
    
    perl bin/transform_phenotypes_to_expression_atlas_lucy_index.pl -i /home/vagrant/Downloads/cass_phenotype.xls -o /home/vagrant/cxgn/cassbase/bin/lucy.tsv -c /home/vagrant/cxgn/cassbase/bin/pre_corr.tsv -f /home/vagrant/cxgn/cassbase/bin/corr.tsv -p /home/vagrant/cxgn/cassbase/bin/project.txt -d /home/vagrant/cxgn/cassbase/bin/desc.tsv

=head1 COMMAND-LINE OPTIONS
  ARGUMENTS
 -i phenotype file downloaded directly from website
 use absolute paths for output files

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
use Statistics::R;

our ($opt_i, $opt_p, $opt_o, $opt_c, $opt_f, $opt_d);

getopts('i:p:o:c:f:d:');

if (!$opt_i || !$opt_p || !$opt_o || !$opt_c || !$opt_f || !$opt_d) {
    pod2usage(-verbose => 2, -message => "Must provide options -i (input file) -p (project file out) -o (lucy out file) -c (corr pre-3col out file) -f (corr out file) -d (metabolite description oufile )\n");
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
    $tissue_term =~ s/cass //g;
    $collection_term =~ s/cass //g;
    $age_term =~ s/cass //g;
    my ($tiss, $tiss_ont) = split /\|/, $tissue_term; #/#
    my ($metabolite, $chebi_ont) = split /\|/, $chebi_term; #/#
    $tiss =~ s/ /_/g;
    $tiss =~ s/\s/_/g;
    $metabolite =~ s/ /_/g;
    $metabolite =~ s/\s/_/g;
    $metabolite =~ s/\(//g;
    $metabolite =~ s/\)//g;
    push @traits, [$metabolite, $tiss, $collection_term, $age_term];
}

my %intermed;
my %corr_steps;
my %project_info;
my %accession_info_hash;
for my $row ( 4 .. $row_max ) {

    my $accession_name = $worksheet->get_cell($row,7)->value();
    my $project_name = $worksheet->get_cell($row,2)->value();
    my $project_design = $worksheet->get_cell($row,3)->value();
    my $project_location = $worksheet->get_cell($row,5)->value();
    my $project_year = $worksheet->get_cell($row,0)->value();
    
    $project_info{$project_name} = {
        design => $project_design,
        location => $project_location,
        year => $project_year
    };

    #print STDERR $accession_name."\n";
    #if ($accession_name eq 'IITA-TMS-IBA011412'){
    
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

        my $stage;
        if (index($tissue_term, 'leaf') != -1) {
            $stage = 'leaf';
        } elsif (index($tissue_term, 'root') != -1) {
            $stage = 'root';
        } elsif (index($tissue_term, 'stem') != -1) {
            $stage = 'stem';
        }
        if (!$stage) {
            die "Tissue not leaf, root, or stem";
        }

        $accession_info_hash{$project_name}->{$accession_name}->{$stage}->{$tissue_term} = 1;

        my $temp_key = "$chebi_term, $accession_name, $tissue_term, $collection_term, $age_term";
        my $step2 = "$accession_name";
        my $corr_step = "$accession_name, $tissue_term";
        if (exists($intermed{$temp_key})) {
            my $values = $intermed{$temp_key}->[3];
            push @$values, $value;
            $intermed{$temp_key}->[3] = $values;
        } else {
            $intermed{$temp_key} = [$chebi_term, $tissue_term, $step2, [$value], $corr_step];
        }
        $corr_steps{$corr_step} = 1;
    }

    #}

}
#print STDERR Dumper \%intermed;
#print STDERR Dumper keys %intermed;
#print STDERR Dumper \%project_info;
#print STDERR Dumper \%accession_info_hash;

foreach my $project_name (keys %accession_info_hash) {
    $project_info{$project_name}->{accessions} = $accession_info_hash{$project_name};
}
#print STDERR Dumper \%project_info;

my @corr_steps_sorted;
foreach (sort keys %corr_steps) {
    push @corr_steps_sorted, $_;
}

my %corr_out;
my %unique_metabolites;
foreach (sort keys %intermed) {
    my $chebi_term = $intermed{$_}->[0];
    my $step1 = $intermed{$_}->[1];
    my $step2 = $intermed{$_}->[2];
    my $values = $intermed{$_}->[3];
    my $corr_step = $intermed{$_}->[4];

    my @non_empty_values = grep($_ ne "", @$values);
    my $average = mean(\@non_empty_values);
    my $stddev = stddev(\@non_empty_values);

    my $display_average = sprintf("%.2f", $average->query);
    my $display_stddev = sprintf("%.2f", $stddev->query);
    my @non_empty_values_formatted;
    foreach (@non_empty_values) {
        push @non_empty_values_formatted, sprintf("%.2f", $_);
    }

    push @data_out, [$chebi_term, $step1, $step2, $display_average, $display_stddev, \@non_empty_values_formatted];

    $corr_out{$chebi_term}->{$corr_step} = $display_average;
    $unique_metabolites{$chebi_term}++;
}

#print STDERR Dumper \@data_out;
#print STDERR Dumper \%corr_out;

open(my $fh, ">", $opt_d);
    print STDERR $opt_d."\n";
    foreach (keys %unique_metabolites) {
        print $fh "1\t$_\t$_\n";
    }
close $fh;

open(my $fh, ">", $opt_o);
    print STDERR $opt_o."\n";
    foreach (@data_out) {
        my $values = $_->[5];
        my $metabolite = $_->[0];
        print $fh "$metabolite\t$_->[2]\t$_->[1]\t$_->[3]\t$_->[4]\t".join(',', @$values),"\n";
    }
close $fh;

open(my $fh, ">", $opt_c) || die("\nERROR:\n");
    print STDERR $opt_c."\n";
    print $fh "Metabolites\t", join("\t", @corr_steps_sorted), "\n";
    foreach my $chebi (sort keys %corr_out) {
        print $fh "$chebi\t";
        my $vals = $corr_out{$chebi};
        my $step = 1;
        foreach my $corr_step (@corr_steps_sorted) {
            my $c = $vals->{$corr_step};
            if (!$c) { $c = 0;}
            print $fh "$c";
            if ($step < scalar(@corr_steps_sorted)) {
                print $fh "\t";
            }
            $step++;
        }
        print $fh "\n";
    }
close $fh;

my $R = Statistics::R->new();
my $out1 = $R->run(
    qq`data<-read.delim(exp <- "$opt_c", header=TRUE)`,
    qq`rownames(data)<-data[,1]`,
    qq`data<-data[,-1]`,
    qq`minexp <- 0`,
    qq`CorrMat<-as.data.frame(cor(t(data[rowSums(data)>minexp,]), method="spearman"))`,
    qq`write.table(CorrMat, file="$opt_f", sep="\t")`,
);

open (my $file_fh, "<", "$opt_f") || die ("\nERROR: the file $opt_f could not be found\n");
    my $header = <$file_fh>;
    chomp($header);
    my @gene_header = split("\t",$header);
    unshift(@gene_header,"");

    my %pairs;
    my @final_out;
    while (my $line = <$file_fh>) {
        chomp($line);
        my @line = split("\t",$line);

        for (my $n = 1; $n < @gene_header; $n++) {
            $line[0] =~ s/\"//g;
            $gene_header[$n] =~ s/\"//g;

            my $hash_key = $line[0]."_".$gene_header[$n];
            my $hash_key2 = $gene_header[$n]."_".$line[0];

            if ($line[$n] >= 0.10 && $line[0] ne $gene_header[$n]) {
                if (!$pairs{$hash_key}) {
                    push @final_out, "$line[0]\t$gene_header[$n]\t".sprintf("%.2f",$line[$n]);
                    $pairs{$hash_key} = 1;
                    $pairs{$hash_key2} = 1;
                }
            }
        }
    }
close $file_fh;

open (my $file_fh, ">", "$opt_f") || die ("\nERROR:\n");
    print STDERR $opt_f."\n";
    foreach (@final_out) {
        print $file_fh $_;
        print $file_fh "\n";
    }
close $file_fh;

open (my $file_fh, ">", "$opt_p") || die ("\nERROR:\n");
    print STDERR $opt_p."\n";
    print $file_fh "#organism\norganism_species: Manihot esculenta\norganism_variety: \norganism_description: Manihot esculenta\n# organism - end\n\n";
    foreach my $project_name (keys %project_info) {
        my $project_year = $project_info{$project_name}->{year};
        my $project_design = $project_info{$project_name}->{design};
        my $project_location = $project_info{$project_name}->{location};

        print $file_fh "#project\nproject_name: $project_name\nproject_contact: \nproject_description: Accessions used in Trial '$project_name', Location: '$project_location', Breeding Program: 'CASS', Project Design: '$project_design'\nexpr_unit: ug/gDW\nindex_dir_name: cass_index\n# project - end\n\n";

        my $accession_hash = $project_info{$project_name}->{accessions};
        foreach my $accession (keys %$accession_hash) {

            print $file_fh "# figure --- All info needed for a cluster of images (usually includes a stage and all its tissues). Copy this block as many times as you need (including as many tissue layer blocks as you need).\nfigure_name: $accession\ncube_stage_name: $accession\nconditions:\n# write figure metadata\n\n";

            my $stage_hash = $accession_hash->{$accession};
            foreach my $stage (keys %$stage_hash) {
                print $file_fh "#stage layer\nlayer_name: $accession\nlayer_description:\nlayer_type: stage\nbg_color:\nlayer_image: $accession.png\nimage_width: 250\nimage_height: 500\ncube_ordinal: 10\nimg_ordinal: 10\norgan: $stage\n# layer - end\n\n";

                my $tissue_hash = $stage_hash->{$stage};
                foreach my $tissue (keys %$tissue_hash) {
                    print $file_fh "#tissue layer\nlayer_name: $tissue\nlayer_description: $tissue\nlayer_type: tissue\nbg_color:\nlayer_image: $tissue.png\nimage_width: 250\nimage_height: 500\ncube_ordinal: 100\nimg_ordinal: 100\norgan: $stage\n# layer - end\n\n";
                }
            }

            print $file_fh "# figure - end\n\n";
        }
    }
close $file_fh;

print STDERR "Script Complete.\n";
