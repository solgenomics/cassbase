#!/usr/bin/perl

=head1

change_week24_to_week16.pl - phenotypes were wrongly loaded using week 24. this script swaps those phenotypes to week 16.

=head1 SYNOPSIS

    change_week24_to_week16.pl -H [dbhost] -D [dbname] -l 12,123,42

=head1 COMMAND-LINE OPTIONS
  ARGUMENTS
 -H host name (required) e.g. "localhost"
 -D database name (required) e.g. "cxgn_cassava"
 -l comma separated list of trial ids

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

our ($opt_H, $opt_D, $opt_l);

getopts('H:D:l:');

if (!$opt_H || !$opt_D || !$opt_l) {
    die "Must provide options -H (hostname), -D (database name), -l (list of trial ids) \n";
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

my @trial_ids = split ',', $opt_l;

my $q = "select phenotype_id, cvterm.name from phenotype join cvterm on (phenotype.cvalue_id=cvterm.cvterm_id) join nd_experiment_phenotype using(phenotype_id) join nd_experiment_project using(nd_experiment_id) where project_id=? and cvterm.name like '%week 24%';";

my $q2 = "select cvterm_id from cvterm where name =?;";

my $q3 = "update phenotype set cvalue_id=? where phenotype_id=?;";

my $changed_count = 0;
foreach my $trial_id (@trial_ids) {
    print STDERR $trial_id."\n";
    my $h = $dbh->prepare($q);
    $h->execute($trial_id);
    while (my ($phenotype_id, $cvterm) = $h->fetchrow_array()) {
        $cvterm =~ s/week 24/week 16/;
        $cvterm =~ s/0000029/0000021/;
        my $h2 = $dbh->prepare($q2);
        $h2->execute($cvterm);
        print STDERR $cvterm."\n";
        while (my ($cvterm_id) = $h2->fetchrow_array()) {
            print STDERR "ID:$cvterm_id\n";
            my $h3 = $dbh->prepare($q3);
            $h3->execute($cvterm_id, $phenotype_id);
            $changed_count++;
        }
    }
}
print STDERR "Changed: $changed_count\n";
print STDERR "Script Complete.\n";
