#!/usr/bin/perl

=head1

transform_phenotypes_to_expression_atlas_lucy_index.pl 

=head1 SYNOPSIS

    transform_phenotypes_to_expression_atlas_lucy_index.pl  -i [infile]
    
    perl bin/transform_phenotypes_to_expression_atlas_lucy_index.pl -i /home/vagrant/Downloads/cass_phenotype.csv -o /home/vagrant/cxgn/cassbase/bin/lucy.tsv -c /home/vagrant/cxgn/cassbase/bin/pre_corr.tsv -f /home/vagrant/cxgn/cassbase/bin/corr.tsv -p /home/vagrant/cxgn/cassbase/bin/project.txt -d /home/vagrant/cxgn/cassbase/bin/desc.tsv -v 1 -n project_name

=head1 COMMAND-LINE OPTIONS
  ARGUMENTS
 -i phenotype csv file downloaded directly from website
 -o output lucy.tsv file that can be indexed using TEA script
 -c output pre-correlation file. This is the file fed into the R script below.
 -f output correlation.tsv file that can be indexed using TEA script
 -p output project.txt file that can be loaded into database using TEA script
 -d output description file that can be indexed by TEA script
 -v version number for how to group variables (e.g. grouping traits with tissues). currently 1 to 5
 -n project name

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
use Statistics::Basic qw(:all);
use Statistics::R;
use Text::CSV;
use JSON;
use LWP::Simple;

our ($opt_i, $opt_p, $opt_o, $opt_c, $opt_f, $opt_d, $opt_v, $opt_n, $opt_t);

getopts('i:p:o:c:f:d:v:n:t:');

if (!$opt_i || !$opt_p || !$opt_o || !$opt_c || !$opt_f || !$opt_d || !$opt_v || !$opt_n || !$opt_t) {
    die "Must provide options -i (input file) -p (project file out) -o (lucy out file) -c (corr pre-3col out file) -f (corr out file) -d (metabolite description oufile -v (script_version) -n (project name) -t (temp dir)\n";
}

open(my $fh, '<', $opt_i)
    or die "Could not open file '$opt_i' $!";

my $brapi_json_response = <$fh>;
my $brapi_response = decode_json $brapi_json_response;
print STDERR Dumper $brapi_response;
my $remote_file_path = $brapi_response->{metadata}->{datafiles}->[0];
print STDERR $remote_file_path."\n";
my $phenotype_download = $opt_t."/phenotype_download.csv";
my $status_code = mirror($remote_file_path, $phenotype_download);
print STDERR $status_code."\n";

my $csv = Text::CSV->new({ sep_char => ',' });

open(my $fh, '<', $phenotype_download)
    or die "Could not open file '$phenotype_download' $!";

my $trait_row = <$fh>;
my @columns;
if ($csv->parse($trait_row)) {
    @columns = $csv->fields();
} else {
    die "Could Not Parse Line: $trait_row\n";
}
my $col_max = scalar(@columns)-1;

#my $parser   = Spreadsheet::ParseExcel->new();
#my $excel_obj = $parser->parse($opt_i);

#my $worksheet = ( $excel_obj->worksheets() )[0]; #support only one worksheet
#my ( $row_min, $row_max ) = $worksheet->row_range();
#my ( $col_min, $col_max ) = $worksheet->col_range();

my @data_out;

my @traits;
for my $col ( 15 .. $col_max) {
    my $multiterm_trait = $columns[$col];
    my @component_terms = split /\|\|/, $multiterm_trait; #/#
    if ($opt_v == 1){
        my $chebi_term;
        my $tissue_term;
        my $collection_term;
        my $age_term;
        my $unit_term;
        my $final_term;
        if (scalar(@component_terms) == 6){
            $chebi_term = $component_terms[0];
            $tissue_term = $component_terms[1];
            $collection_term = $component_terms[2];
            $age_term = $component_terms[3];
            $unit_term = $component_terms[4];
            $final_term = $component_terms[5];
        }
        if (scalar(@component_terms) == 3){
            $chebi_term = $component_terms[0];
            $age_term = $component_terms[1];
            $final_term = $component_terms[2];
        }
        my $i=rindex($final_term, /\|/); #/#
        my $institute_term=substr($final_term,0,$i);

        $tissue_term =~ s/cass //g;
        $collection_term =~ s/cass //g;
        $age_term =~ s/cass //g;
        my ($tiss, $tiss_ont) = split /\|/, $tissue_term; #/#
        my ($metabolite, $chebi_ont) = split /\|/, $chebi_term; #/#
        my ($collection, $collection_ont) = split /\|/, $collection_term; #/#
        my ($age, $age_ont) = split /\|/, $age_term; #/#
        my ($unit, $unit_ont) = split /\|/, $unit_term; #/#
        my ($institute, $inst_ont) = split /\|/, $institute_term; #/#
        $tiss =~ s/ /_/g;
        $tiss =~ s/\s/_/g;
        $tiss =~ s/\(//g;
        $tiss =~ s/\)//g;
        $metabolite =~ s/ /_/g;
        $metabolite =~ s/\s/_/g;
        $metabolite =~ s/\(//g;
        $metabolite =~ s/\)//g;
        $collection =~ s/ /_/g;
        $collection =~ s/\s/_/g;
        $collection =~ s/\(//g;
        $collection =~ s/\)//g;
        $age =~ s/ /_/g;
        $age =~ s/\s/_/g;
        $age =~ s/\(//g;
        $age =~ s/\)//g;
        $institute =~ s/ /_/g;
        $institute =~ s/\s/_/g;
        $institute =~ s/\(//g;
        $institute =~ s/\)//g;
        $unit =~ s/ /_/g;
        $unit =~ s/\s/_/g;
        $unit =~ s/\(//g;
        $unit =~ s/\)//g;
        push @traits, [$metabolite."_".$institute, $tiss, $collection, $age, $unit];
    } elsif ($opt_v == 2){
        if (scalar(@component_terms) == 6){
            my $chebi_term = $component_terms[0];
            my $tissue_term = $component_terms[1];
            my $collection_term = $component_terms[2];
            my $age_term = $component_terms[3];
            my $unit_term = $component_terms[4];
            my $final_term = $component_terms[5];
            my $i=rindex($final_term, /\|/); #/#
            my $institute_term=substr($final_term,0,$i);

            $tissue_term =~ s/cass //g;
            $collection_term =~ s/cass //g;
            $age_term =~ s/cass //g;
            my ($tiss, $tiss_ont) = split /\|/, $tissue_term; #/#
            my ($metabolite, $chebi_ont) = split /\|/, $chebi_term; #/#
            my ($collection, $collection_ont) = split /\|/, $collection_term; #/#
            my ($age, $age_ont) = split /\|/, $age_term; #/#
            my ($unit, $unit_ont) = split /\|/, $unit_term; #/#
            my ($institute, $inst_ont) = split /\|/, $institute_term; #/#
            $tiss =~ s/ /_/g;
            $tiss =~ s/\s/_/g;
            $tiss =~ s/\(//g;
            $tiss =~ s/\)//g;
            $metabolite =~ s/ /_/g;
            $metabolite =~ s/\s/_/g;
            $metabolite =~ s/\(//g;
            $metabolite =~ s/\)//g;
            $collection =~ s/ /_/g;
            $collection =~ s/\s/_/g;
            $collection =~ s/\(//g;
            $collection =~ s/\)//g;
            $age =~ s/ /_/g;
            $age =~ s/\s/_/g;
            $age =~ s/\(//g;
            $age =~ s/\)//g;
            $institute =~ s/ /_/g;
            $institute =~ s/\s/_/g;
            $institute =~ s/\(//g;
            $institute =~ s/\)//g;
            $unit =~ s/ /_/g;
            $unit =~ s/\s/_/g;
            $unit =~ s/\(//g;
            $unit =~ s/\)//g;
            push @traits, [$metabolite."_".$tiss."_".$institute, $tiss, $collection, $age, $unit];
        } elsif (scalar(@component_terms) == 3){
            my $agronomic_term = $component_terms[0];
            my $age_term = $component_terms[1];
            my $final_term = $component_terms[2];
            my $i=rindex($final_term, /\|/); #/#
            my $institute_term=substr($final_term,0,$i);

            $age_term =~ s/cass //g;
            my ($agro, $agro_ont) = split /\|/, $agronomic_term; #/#
            my ($age, $age_ont) = split /\|/, $age_term; #/#
            my ($institute, $inst_ont) = split /\|/, $institute_term; #/#
            $agro =~ s/ /_/g;
            $agro =~ s/\s/_/g;
            $agro =~ s/\(//g;
            $agro =~ s/\)//g;
            $age =~ s/ /_/g;
            $age =~ s/\s/_/g;
            $age =~ s/\(//g;
            $age =~ s/\)//g;
            $institute =~ s/ /_/g;
            $institute =~ s/\s/_/g;
            $institute =~ s/\(//g;
            $institute =~ s/\)//g;
            push @traits, [$agro."_".$institute, '', '', $age, ''];
        }
    } elsif ($opt_v == 3 || $opt_v == 4 || $opt_v == 5){
        if (scalar(@component_terms) == 6){
            my $chebi_term = $component_terms[0];
            my $tissue_term = $component_terms[1];
            my $collection_term = $component_terms[2];
            my $age_term = $component_terms[3];
            my $unit_term = $component_terms[4];
            my $final_term = $component_terms[5];
            my $i=rindex($final_term, /\|/); #/#
            my $institute_term=substr($final_term,0,$i);

            $tissue_term =~ s/cass //g;
            $collection_term =~ s/cass //g;
            $age_term =~ s/cass //g;
            my ($tiss, $tiss_ont) = split /\|/, $tissue_term; #/#
            my ($metabolite, $chebi_ont) = split /\|/, $chebi_term; #/#
            my ($collection, $collection_ont) = split /\|/, $collection_term; #/#
            my ($age, $age_ont) = split /\|/, $age_term; #/#
            my ($unit, $unit_ont) = split /\|/, $unit_term; #/#
            my ($institute, $inst_ont) = split /\|/, $institute_term; #/#
            $tiss =~ s/ /_/g;
            $tiss =~ s/\s/_/g;
            $tiss =~ s/\(//g;
            $tiss =~ s/\)//g;
            $metabolite =~ s/ /_/g;
            $metabolite =~ s/\s/_/g;
            $metabolite =~ s/\(//g;
            $metabolite =~ s/\)//g;
            $collection =~ s/ /_/g;
            $collection =~ s/\s/_/g;
            $collection =~ s/\(//g;
            $collection =~ s/\)//g;
            $age =~ s/ /_/g;
            $age =~ s/\s/_/g;
            $age =~ s/\(//g;
            $age =~ s/\)//g;
            $institute =~ s/ /_/g;
            $institute =~ s/\s/_/g;
            $institute =~ s/\(//g;
            $institute =~ s/\)//g;
            $unit =~ s/ /_/g;
            $unit =~ s/\s/_/g;
            $unit =~ s/\(//g;
            $unit =~ s/\)//g;
            push @traits, [$metabolite."_".$tiss."_".$institute, $tiss, $collection, $age, $unit];
        } elsif (scalar(@component_terms) == 3){
            my $agronomic_term = $component_terms[0];
            my $age_term = $component_terms[1];
            my $final_term = $component_terms[2];
            my $i=rindex($final_term, /\|/); #/#
            my $institute_term=substr($final_term,0,$i);

            $age_term =~ s/cass //g;
            my ($agro, $agro_ont) = split /\|/, $agronomic_term; #/#
            my ($age, $age_ont) = split /\|/, $age_term; #/#
            my ($institute, $inst_ont) = split /\|/, $institute_term; #/#
            $agro =~ s/ /_/g;
            $agro =~ s/\s/_/g;
            $agro =~ s/\(//g;
            $agro =~ s/\)//g;
            $age =~ s/ /_/g;
            $age =~ s/\s/_/g;
            $age =~ s/\(//g;
            $age =~ s/\)//g;
            $institute =~ s/ /_/g;
            $institute =~ s/\s/_/g;
            $institute =~ s/\(//g;
            $institute =~ s/\)//g;
            push @traits, [$agro."_".$institute, '', '', $age, ''];
        }
    }
}

#print STDERR Dumper \@traits;
#print STDERR scalar(@traits)."\n";

my %intermed;
my %corr_steps;
my %project_info;
my %accession_info_hash;
my $project_name = $opt_n;

while ( my $row = <$fh> ){
    my @columns;
    if ($csv->parse($row)) {
        @columns = $csv->fields();
    } else {
        die "Could not parse line $row\n";
    }
    
    my $accession_name = $columns[7];
    #my $project_name = $columns[2];
    my $project_design = $columns[3];
    my $project_location = $columns[5];
    my $project_year = $columns[0];
    
    $project_info{$project_name} = {
        design => $project_design,
        location => $project_location,
        year => $project_year
    };

    #print STDERR $accession_name."\n";

    for( my $i=0; $i<scalar(@traits); $i++) {
        my $trait_col = $i + 15;
        #print STDERR "$row $trait_col\n";
        my $value = '';
        if ($columns[$trait_col]) {
            $value = $columns[$trait_col];
        }
        my $trait_term = $traits[$i]->[0];
        my $tissue_term = $traits[$i]->[1];
        my $collection_term = $traits[$i]->[2];
        my $age_term = $traits[$i]->[3];
        my $unit_term = $traits[$i]->[4];

        #Normalize mg/g and ng/g to ug/g
        if ($value && $unit_term eq 'mg/gDW'){
            $value = $value * 1000;
        }
        if ($value && $unit_term eq 'ng/gDW'){
            $value = $value / 1000;
        }

        my $stage;
        my $temp_key;
        my $step2;
        my $corr_step;
        if ($opt_v == 1){
            if (index($tissue_term, 'leaf') != -1) {
                $stage = 'leaf';
            } elsif (index($tissue_term, 'root') != -1) {
                $stage = 'root';
            } elsif (index($tissue_term, 'stem') != -1) {
                $stage = 'stem';
            }
            if (!$stage) {
                print STDERR Dumper $traits[$i];
                print STDERR "Tissue not leaf, root, or stem\n";
                next;
            }

            $temp_key = "$trait_term, $accession_name, $tissue_term, $collection_term, $age_term";
            $step2 = "$accession_name-$collection_term-$age_term";
            $corr_step = "$accession_name, $tissue_term, $collection_term, $age_term";
            $accession_info_hash{$project_name}->{$step2}->{$stage}->{$tissue_term} = 1;

            if (exists($intermed{$temp_key})) {
                my $values = $intermed{$temp_key}->[3];
                push @$values, $value;
                $intermed{$temp_key}->[3] = $values;
            } else {
                $intermed{$temp_key} = [$trait_term, $tissue_term, $step2, [$value], $corr_step];
            }
            $corr_steps{$corr_step} = 1;
        } elsif ($opt_v == 2){
            $temp_key = "$trait_term, $accession_name, $age_term";
            $step2 = "$accession_name-$age_term";
            $corr_step = "$accession_name, $age_term";
            $accession_info_hash{$project_name}->{$step2}->{'plant'}->{'plant'} = 1;

            if (exists($intermed{$temp_key})) {
                my $values = $intermed{$temp_key}->[3];
                push @$values, $value;
                $intermed{$temp_key}->[3] = $values;
            } else {
                $intermed{$temp_key} = [$trait_term, 'plant', $step2, [$value], $corr_step];
            }
            $corr_steps{$corr_step} = 1;
        } elsif ($opt_v == 3){
            if ($age_term ne 'week_16'){
                next;
            }
            $temp_key = "$trait_term, $accession_name";
            $step2 = "$accession_name";
            $corr_step = "$accession_name";
            $accession_info_hash{$project_name}->{$step2}->{'plant'}->{'plant'} = 1;

            if (exists($intermed{$temp_key})) {
                my $values = $intermed{$temp_key}->[3];
                push @$values, $value;
                $intermed{$temp_key}->[3] = $values;
            } else {
                $intermed{$temp_key} = [$trait_term, 'plant', $step2, [$value], $corr_step];
            }
            $corr_steps{$corr_step} = 1;
        } elsif ($opt_v == 4){
            if ($age_term eq 'week_16' && ($trait_term eq 'fresh_root_weight_per_plant' || $trait_term eq 'fresh_shoot_weight_per_plant' || $trait_term eq 'fresh_total_weight_per_plant' || $trait_term eq 'harvest_index')){
                next;
            }
            $temp_key = "$trait_term, $accession_name";
            $step2 = "$accession_name";
            $corr_step = "$accession_name";
            $accession_info_hash{$project_name}->{$step2}->{'plant'}->{'plant'} = 1;

            if (exists($intermed{$temp_key})) {
                my $values = $intermed{$temp_key}->[3];
                push @$values, $value;
                $intermed{$temp_key}->[3] = $values;
            } else {
                $intermed{$temp_key} = [$trait_term, 'plant', $step2, [$value], $corr_step];
            }
            $corr_steps{$corr_step} = 1;
        } elsif ($opt_v == 5){
            if (($age_term eq 'week_16' && ($trait_term eq 'fresh_root_weight_per_plant' || $trait_term eq 'fresh_shoot_weight_per_plant' || $trait_term eq 'fresh_total_weight_per_plant' || $trait_term eq 'harvest_index')) || $accession_name eq 'TMEB419'){
                next;
            }
            $temp_key = "$trait_term, $accession_name";
            $step2 = "$accession_name";
            $corr_step = "$accession_name";
            $accession_info_hash{$project_name}->{$step2}->{'plant'}->{'plant'} = 1;

            if (exists($intermed{$temp_key})) {
                my $values = $intermed{$temp_key}->[3];
                push @$values, $value;
                $intermed{$temp_key}->[3] = $values;
            } else {
                $intermed{$temp_key} = [$trait_term, 'plant', $step2, [$value], $corr_step];
            }
            $corr_steps{$corr_step} = 1;
        }
    }

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
my %unique_traits;
foreach (sort keys %intermed) {
    my $trait_term = $intermed{$_}->[0];
    my $step1 = $intermed{$_}->[1];
    my $step2 = $intermed{$_}->[2];
    my $values = $intermed{$_}->[3];
    my $corr_step = $intermed{$_}->[4];

    my @non_empty_values = grep($_ ne "", @$values);
    #print STDERR Dumper \@non_empty_values;
    my $average = mean(\@non_empty_values);
    my $stddev = stddev(\@non_empty_values);

    my $display_average = sprintf("%.2f", $average->query);
    my $display_stddev = sprintf("%.2f", $stddev->query);
    my @non_empty_values_formatted;
    foreach (@non_empty_values) {
        push @non_empty_values_formatted, sprintf("%.2f", $_);
    }

    push @data_out, [$trait_term, $step1, $step2, $display_average, $display_stddev, \@non_empty_values_formatted];

    $corr_out{$trait_term}->{$corr_step} = $display_average;
    $unique_traits{$trait_term}++;
}

#print STDERR Dumper \@data_out;
#print STDERR Dumper \%corr_out;

open(my $fh, ">", $opt_d);
    print STDERR $opt_d."\n";
    foreach (keys %unique_traits) {
        print $fh "1\t$_\t$_\n";
    }
close $fh;

open(my $fh, ">", $opt_o);
    print STDERR $opt_o."\n";
    foreach (@data_out) {
        my $values = $_->[5];
        my $trait_term = $_->[0];
        print $fh "$trait_term\t$_->[2]\t$_->[1]\t$_->[3]\t$_->[4]\t".join(',', @$values),"\n";
    }
close $fh;

open(my $fh, ">", $opt_c) || die("\nERROR:\n");
    print STDERR $opt_c."\n";
    print $fh "Metabolites\t", join("\t", @corr_steps_sorted), "\n";
    foreach my $trait_term (sort keys %corr_out) {
        print $fh "$trait_term\t";
        my $vals = $corr_out{$trait_term};
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
        my $project_name_index_dir = $project_name;
        $project_name_index_dir =~ s/ //g;
        $project_name_index_dir =~ s/\s//g;
        my $project_year = $project_info{$project_name}->{year};
        my $project_design = $project_info{$project_name}->{design};
        my $project_location = $project_info{$project_name}->{location};

        print $file_fh "#project\nproject_name: $project_name\nproject_contact: \nproject_description: Accessions used in Trial '$project_name', Location: '$project_location', Breeding Program: 'CASS', Project Design: '$project_design'\nexpr_unit: ug/gDW\nindex_dir_name: cass_index_$project_name_index_dir\n# project - end\n\n";

        my $accession_hash = $project_info{$project_name}->{accessions};
        foreach my $accession (keys %$accession_hash) {

            print $file_fh "# figure --- All info needed for a cluster of images (usually includes a stage and all its tissues). Copy this block as many times as you need (including as many tissue layer blocks as you need).\nfigure_name: $accession\ncube_stage_name: $accession\nconditions:\n# write figure metadata\n\n";

            #if ($opt_v == 1){
                my $stage_hash = $accession_hash->{$accession};
                foreach my $stage (keys %$stage_hash) {
                    print $file_fh "#stage layer\nlayer_name: $accession\nlayer_description:\nlayer_type: stage\nbg_color:\nlayer_image: plant_background.png\nimage_width: 250\nimage_height: 500\ncube_ordinal: 10\nimg_ordinal: 10\norgan: $stage\n# layer - end\n\n";

                    my $tissue_hash = $stage_hash->{$stage};
                    foreach my $tissue (keys %$tissue_hash) {
                        print $file_fh "#tissue layer\nlayer_name: $tissue\nlayer_description: $tissue\nlayer_type: tissue\nbg_color:\nlayer_image: $tissue.png\nimage_width: 250\nimage_height: 500\ncube_ordinal: 100\nimg_ordinal: 100\norgan: $stage\n# layer - end\n\n";
                    }
                }
            #}

            print $file_fh "# figure - end\n\n";
        }
    }
close $file_fh;

print STDERR "Script Complete.\n";
