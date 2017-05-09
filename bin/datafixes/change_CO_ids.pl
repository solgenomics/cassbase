#!/usr/bin/perl

=head1

change_CO_ids.pl - CO terms for "Total fresh weight|CO:0000308", "total fresh weight by plot|CO:0000309", "total fresh weight by plant|CO:0000310" were changed to "Total fresh weight|CO:0000324", "total fresh weight by plot|CO:0000325", "total fresh weight by plant|CO:0000326". The postcomposed_terms are updated by this script.

=head1 SYNOPSIS

    change_CO_ids.pl -H [dbhost] -D [dbname]

=head1 COMMAND-LINE OPTIONS
  ARGUMENTS
 -H host name (required) e.g. "localhost"
 -D database name (required) e.g. "cxgn_cassava"

=head1 DESCRIPTION

=head1 AUTHOR

 Nicolas Morales (nm529@cornell.edu)

=cut

use strict;

use Getopt::Std;
use Data::Dumper;
use Carp qw /croak/ ;
use Bio::Chado::Schema;
use CXGN::DB::InsertDBH;

our ($opt_H, $opt_D);

getopts('H:D:');

if (!$opt_H || !$opt_D) {
    die "Must provide options -H (hostname), -D (database name) \n";
}

my $dbhost = $opt_H;
my $dbname = $opt_D;

my $dbh = CXGN::DB::InsertDBH->new({ 
	dbhost=>$dbhost,
	dbname=>$dbname,
	dbargs => {AutoCommit => 1, RaiseError => 1}
});

my $schema= Bio::Chado::Schema->connect(  sub { $dbh->get_actual_dbh() } );
$dbh->do('SET search_path TO public,sgn');

my $q = "select cvterm_id, name from cvterm where (name like '%CO:0000308%' OR name like '%CO:0000309%' OR name like '%CO:0000310%');";
my $h = $dbh->prepare($q);
$h->execute();

my $q2 = "update cvterm set name=? where cvterm_id=?;";

my $changed_count = 0;
while (my ($cvterm_id, $cvterm) = $h->fetchrow_array()) {
    $cvterm =~ s/CO:0000308/CO:0000324/;
    $cvterm =~ s/CO:0000309/CO:0000325/;
    $cvterm =~ s/CO:0000310/CO:0000326/;
	print STDERR $cvterm."\n";
	print STDERR "ID:$cvterm_id\n";
    my $h2 = $dbh->prepare($q2);
    $h2->execute($cvterm, $cvterm_id);
    $changed_count++;
}
print STDERR "Changed: $changed_count\n";
print STDERR "Script Complete.\n";
