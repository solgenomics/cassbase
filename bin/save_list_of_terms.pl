#!/usr/bin/env perl

=head1
save_list_of_terms.pl

=head1 SYNOPSIS

    perl save_list_of_terms.pl -H localhost -D fixture2 -i input_list.txt -c 'stringtie_transcripts'

=head1 COMMAND-LINE OPTIONS

 -H  host name
 -D  database name
 -i  input list file. should be a flat list of traits, each on their own line. E.g. Manes.05G094400.1.v6.1|ST:0302847
 -c  cv name to save these terms under

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

our ($opt_H, $opt_D, $opt_i, $opt_c);
getopts('H:D:i:c:');



if (!$opt_D || !$opt_H || !$opt_i || !$opt_c) {
  die("Exiting: options missing\nRequires: -D -H -i -c");
}

my $dbh = CXGN::DB::InsertDBH
  ->new({
	 dbname => $opt_D,
	 dbhost => $opt_H,
	 dbargs => {AutoCommit => 1,
		    RaiseError => 1},
	});

my $schema = Bio::Chado::Schema->connect(  sub { $dbh->get_actual_dbh() } );

my $cv = $schema->resultset("Cv::Cv")->find_or_create({name=>$opt_c});
my $cv_id = $cv->cv_id();

open(my $fh, '<', $opt_i)
    or die "Could not open file '$opt_i' $!";

while (my $row = <$fh>) {
    my ($term, $db_acc) = split '\|', $row;
    my ($db, $acc) = split '\:', $db_acc;

    my $db = $schema->resultset("General::Db")->find_or_create({name=>$db});
    my $db_id = $db->db_id();
    my $dbxref = $schema->resultset("General::Dbxref")->create({db_id=>$db_id, accession=>$acc});
    my $dbxref_id = $dbxref->dbxref_id();
    my $cvterm = $schema->resultset("Cv::Cvterm")->create({cv_id=>$cv_id, name=>$term, dbxref_id=>$dbxref_id});
}


print STDERR "Complete.\n";
