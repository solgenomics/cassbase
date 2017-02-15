#!/usr/bin/env perl

=head1
append_institute_to_trait.pl

=head1 SYNOPSIS

    perl append_institute_to_trait.pl -H localhost -D fixture2

=head1 COMMAND-LINE OPTIONS

 -H  host name
 -D  database name

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

our ($opt_H, $opt_D);
getopts('H:D:');

if (!$opt_D || !$opt_H) {
  die("Exiting: options missing\nRequires: -D -H\n");
}

my $dbh = CXGN::DB::InsertDBH
  ->new({
	 dbname => $opt_D,
	 dbhost => $opt_H,
	 dbargs => {AutoCommit => 1,
		    RaiseError => 1},
	});

my $schema = Bio::Chado::Schema->connect(  sub { $dbh->get_actual_dbh() } );

my $cv = $schema->resultset("Cv::Cv")->find({name=>'postcomposed_terms'});
my $cv_id = $cv->cv_id();

my $db = $schema->resultset("General::Db")->find({name=>'CASSFT'});
my $db_id = $db->db_id();

my %all_cvterms;
my $q = "SELECT cvterm.cvterm_id, cvterm.name from cvterm where cv_id=$cv_id;";
my $h = $dbh->prepare($q);
$h->execute();
while( my ($cvterm_id, $cvterm) = $h->fetchrow_array() ){
    if (index($cvterm, 'Manes') == -1) {
        $all_cvterms{$cvterm} = $cvterm_id;
    }
}

foreach (keys %all_cvterms){
    my $cvterm_id = $all_cvterms{$_};
    my $q = "UPDATE cvterm set name=? where cvterm_id=?;";
    my $h = $dbh->prepare($q);
    my $new_name = $_."||ER|CASSINST:0000001";
    $h->execute($new_name, $cvterm_id);
}

my @all_cvterms;
foreach (keys %all_cvterms){
    push @all_cvterms, $_;
}

my @additional_institutes = (
    'KL|CASSINST:0000002',
    'GO_AF|CASSINST:0000003',
    'GO_MS|CASSINST:0000004',
    'GO_AS|CASSINST:0000005',
    'IITA|CASSINST:0000006',
    'BTI|CASSINST:0000007'
);

$q = "SELECT accession FROM dbxref WHERE db_id=$db_id ORDER BY dbxref_id DESC LIMIT 1;";
$h = $dbh->prepare($q);
$h->execute();
my ($last_accession) = $h->fetchrow_array();
my $accession = $last_accession + 1;

my @new_terms;
foreach my $a (@all_cvterms){
    foreach my $b (@additional_institutes){
        push @new_terms, $a."||".$b;
    }
}

foreach (@new_terms){
    my $dbxref = $schema->resultset("General::Dbxref")->create({db_id=>$db_id, accession=>$accession});
    my $dbxref_id = $dbxref->dbxref_id();
    my $cvterm = $schema->resultset("Cv::Cvterm")->create({cv_id=>$cv_id, name=>$_, dbxref_id=>$dbxref_id});
    $accession++;
}

print STDERR "Complete.\n";

