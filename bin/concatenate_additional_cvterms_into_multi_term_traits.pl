#!/usr/bin/env perl

=head1
concatenate_additional_cvterms_into_multi_term_traits.pl

=head1 SYNOPSIS

    this script is very specific to cassbase ontologies, but could possibly be generalized through more sophisticated recursion.

    perl concatenate_additional_cvterms_into_multi_term_traits.pl -H localhost -D fixture2 -l 'linamarin|CHEBI:16441[OR]sn-glycerol 3-phosphate|CHEBI:15978[OR]trans-4-coumaric acid|CHEBI:32374,cass_tissues|CASSTISS:0000000,cass time of day|CASSTIME:0000001,cass number of weeks|CASSTIME:0000005,cass_units|CASSUNIT:0000000' -d CASSFT -c postcomposed_terms

=head1 COMMAND-LINE OPTIONS

 -H  host name
 -D  database name
 -l  comma separated list of parent cvterms. the first term can take [OR] separated specific cvterms. all children of the parent term (is_a relationship) will be concatenated together and saved.
 -d  the db name that the new cvterms will be stored under. Will increment the accession
 -c  the cv name that the new cvterms will be stored under.

=head2 DESCRIPTION


=head2 AUTHOR

Nicolas Morales (nm529@cornell.edu)

April 2014

=head2

=cut

use strict;
use warnings;

use lib 'lib';
use Getopt::Std;
use Bio::Chado::Schema;
use CXGN::DB::InsertDBH;
use CXGN::DB::Connection;
use Data::Dumper;
use SGN::Model::Cvterm;
use Try::Tiny;

our ($opt_H, $opt_D, $opt_l, $opt_d, $opt_c);
getopts('H:D:l:d:c:');



if (!$opt_D || !$opt_H || !$opt_l || !$opt_d || !$opt_c ) {
  die("Exiting: options missing\nRequires: -D -H -l -d -c");
}

my $dbh = CXGN::DB::InsertDBH
  ->new({
	 dbname => $opt_D,
	 dbhost => $opt_H,
	 dbargs => {AutoCommit => 1,
		    RaiseError => 1},
	});

my $schema = Bio::Chado::Schema->connect(  sub { $dbh->get_actual_dbh() } );

my $db = $schema->resultset("General::Db")->find({name=>$opt_d});
my $db_id = $db->db_id();
my $cv = $schema->resultset("Cv::Cv")->find({name=>$opt_c});
my $cv_id = $cv->cv_id();

my $accession_q = "SELECT accession from dbxref WHERE db_id=? ORDER BY accession::int DESC LIMIT 1;";
my $h = $dbh->prepare($accession_q);
$h->execute($db_id);
my ($last_accession) = $h->fetchrow_array();
my $accession = $last_accession + 1;
print STDERR $accession."\n";

my @parent_trait_names = split /,/, $opt_l;

my $first_element = splice @parent_trait_names, 0, 1;
my @first_terms = split /\[OR\]/, $first_element;
foreach my $i (@first_terms) {
    my @children_array;

    push (@children_array, [$i]);

    foreach my $j (@parent_trait_names) {
        my $children;
        if ($j eq 'cass_tissues|CASSTISS:0000000') {
            my @sub_nodes = ('cass leaf|CASSTISS:0000001', 'cass stem|CASSTISS:0000002', 'cass root|CASSTISS:0000003');
            foreach my $t (@sub_nodes) {
                my $sub_children = get_children($schema, $t);
                push @$children, @$sub_children;
            }
        } else {
            $children = get_children($schema, $j);
        }
        push (@children_array, $children);
    }

    print Dumper \@children_array;

    my $count = 0;
    my @concatenated_terms;
    my $first_term = $children_array[0];
    foreach my $a (@$first_term) {
        my $a_concat_term = $a;
        my $second_term = $children_array[1];
        foreach my $b (@$second_term) {
            my $b_concat_term = $a_concat_term.'||'.$b;
            my $third_term = $children_array[2];
            foreach my $c (@$third_term) {
                my $c_concat_term = $b_concat_term.'||'.$c;
                my $fourth_term = $children_array[3];
                foreach my $d (@$fourth_term) {
                    my $d_concat_term = $c_concat_term.'||'.$d;
                    my $fifth_term = $children_array[4];
                    foreach my $e (@$fifth_term) {
                        my $e_concat_term = $d_concat_term.'||'.$e;
                        push @concatenated_terms, $e_concat_term;
                    }
                }
            }
        }
    }



    #print Dumper \@concatenated_terms;
    print scalar(@{$children_array[0]}) * scalar(@{$children_array[1]}) * scalar(@{$children_array[2]}) * scalar(@{$children_array[3]}) * scalar(@{$children_array[4]})."\n";
    print scalar(@concatenated_terms)."\n";

    foreach (@concatenated_terms) {
        my $accession_string = sprintf("%07d",$accession);
        my $dbxref = $schema->resultset("General::Dbxref")->create({db_id=>$db_id, accession=>$accession_string});
        my $dbxref_id = $dbxref->dbxref_id();
        my $cvterm = $schema->resultset("Cv::Cvterm")->create({cv_id=>$cv_id, name=>$_, dbxref_id=>$dbxref_id});
        $accession++;
        $count++;
    }

    print STDERR "Added $count new terms.\n";
}

print STDERR "Complete.\n";


sub get_children {
    my $schema = shift;
    my $term = shift;
    print $term."\n";
    my $parent_node_cvterm_id = SGN::Model::Cvterm->get_cvterm_row_from_trait_name($schema, $term)->cvterm_id();
    my $rel_cvterm_id = SGN::Model::Cvterm->get_cvterm_row($schema, 'is_a', 'relationship')->cvterm_id();

    my $children_ref = $schema->resultset("Cv::CvtermRelationship")->search({type_id => $rel_cvterm_id, object_id => $parent_node_cvterm_id})->search_related('subject');
    my @children;
    while (my $child = $children_ref->next() ) {
        my $dbxref_info = $child->search_related('dbxref');
        my $accession = $dbxref_info->first()->accession();
        my $db_info = $dbxref_info->search_related('db');
        my $db_name = $db_info->first()->name();
        push @children, $child->name."|".$db_name.":".$accession;
    }
    return \@children;
}
