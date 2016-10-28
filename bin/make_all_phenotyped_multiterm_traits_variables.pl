#!/usr/bin/env perl

=head1
make_all_phenotyped_multiterm_traits_variables.pl

=head1 SYNOPSIS

    perl make_all_phenotyped_multiterm_traits_variables.pl -H localhost -D fixture2 -c postcomposed_terms

=head1 COMMAND-LINE OPTIONS

 -H  host name
 -D  database name
 -c  the cv name that the multiterm traits are under.

=head2 DESCRIPTION


=head2 AUTHOR

Nicolas Morales (nm529@cornell.edu)

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

our ($opt_H, $opt_D, $opt_c);
getopts('H:D:c:');



if (!$opt_D || !$opt_H || !$opt_c ) {
  die("Exiting: options missing\nRequires: -D -H -c ");
}

my $dbh = CXGN::DB::InsertDBH
  ->new({
	 dbname => $opt_D,
	 dbhost => $opt_H,
	 dbargs => {AutoCommit => 1,
		    RaiseError => 1},
	});

my $chado_schema = Bio::Chado::Schema->connect(  sub { $dbh->get_actual_dbh() } );

my $coderef = sub {

    my $is_a_term_id = SGN::Model::Cvterm->get_cvterm_row($chado_schema, 'is_a', 'relationship')->cvterm_id();

    my $variable_term_id;
    try {
        $variable_term_id = SGN::Model::Cvterm->get_cvterm_row($chado_schema, 'VARIABLE_OF', $opt_c)->cvterm_id();
    } catch {
        my $variable_term = $chado_schema->resultset("Cv::Cvterm")->create_with({
            name => 'VARIABLE_OF',
            cv   => $opt_c,
        });
        $variable_term_id = $variable_term->cvterm_id();
    };

    my $is_a_term_id = SGN::Model::Cvterm->get_cvterm_row($chado_schema, 'is_a', 'relationship')->cvterm_id();

    my $root_term_id;
    try {
        $root_term_id = SGN::Model::Cvterm->get_cvterm_row($chado_schema, $opt_c, $opt_c)->cvterm_id();
    } catch {
        my $variable_term = $chado_schema->resultset("Cv::Cvterm")->create_with({
            name => $opt_c,
            cv   => $opt_c,
        });
        $root_term_id = $variable_term->cvterm_id();
    };

    my $q = "SELECT cvterm.cvterm_id FROM phenotype JOIN cvterm on (phenotype.cvalue_id=cvterm.cvterm_id) JOIN cv on (cvterm.cv_id=cv.cv_id) WHERE cv.name=?;";
    my $h = $chado_schema->storage->dbh()->prepare($q);
    $h->execute($opt_c);
    my %cvterm_ids;
    while (my ($cvterm) = $h->fetchrow_array()) {
        $cvterm_ids{$cvterm} = 1;
    }


    my $count_q = "SELECT count(*) FROM cvterm_relationship WHERE subject_id=? and (type_id=$is_a_term_id or type_id=$variable_term_id);";
    my $count_h = $chado_schema->storage->dbh()->prepare($count_q);

    foreach (keys %cvterm_ids) {

        my $count_stmt = $count_h->execute($_);
        my $count = $count_h->fetchrow;
        print STDERR "$_ ::: $count\n";
        if (!$count || $count == 0){
            my $create_q = "INSERT INTO cvterm_relationship (type_id, subject_id, object_id) VALUES ($variable_term_id, ?, ?);";
            my $create_h = $chado_schema->storage->dbh()->prepare($create_q);
            $create_h->execute($_, $root_term_id);
        } elsif ($count == 1) {
            my $update_q = "UPDATE cvterm_relationship SET type_id=$variable_term_id WHERE subject_id=? and (type_id=$is_a_term_id or type_id=$variable_term_id);";
            my $update_h = $chado_schema->storage->dbh()->prepare($update_q);
            $update_h->execute($_);
        } else {
            die "More than one is_a relationship already saved...\n";
        }

    }
};

try {
    $chado_schema->txn_do($coderef);
} catch {
    die "Patch failed! Transaction exited." . $_ .  "\n" ;
};

print "You're done!\n";
