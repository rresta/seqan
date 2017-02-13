#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use constant USAGE =><<END;
use Bio::SearchIO;
use Bio::SeqIO;
use CMSequence;
use CMResults;

SYNOPSIS:

Input is a stockholm format sequence alignment with a structure annotation and
a database. Builds a CM and uses this to search the database for putative
homologues.

DESCRIPTION:


OPTIONS:

EXAMPLES:


DEPENDENCIES:



AUTHOR:


COPYRIGHT:

This program is free software. You may copy, modify, and redistribute
it under the same terms as Perl itself.

END

my $help = 0;
my $db = "";
my $q = "";

GetOptions(
           "help" => \$help,
	   "h"    => \$help,
	   "q=s"    => \$q,
	   "db=s"   => \$db
          ) or die USAGE;

$help and die USAGE;

if (!$q){
    die "Missing stockholm file!";
}

if (!$db){
    die "Missing database!";
}

my @dbn = split /\//, $db;
my $dbn = pop(@dbn);

#my $cmd = "/home/pgardner/bin/rsearch -n 1000 -E 10 -m /home/pgardner/INST/rsearch-1.1/matrices $q $db >$q\.$dbn";
my $cmd = "/home/pgardner/bin/rsearch -n 1000 -E 100000 -m /home/pgardner/INST/rsearch-1.1/matrices/RIBOSUM85-60.mat $q $db >$q\.$dbn";

print "Running:\n$cmd\n\n";
system $cmd;


print "\nSee:\n$q\.$dbn\n\nHit:\n";
system "grep -c \^\"\>\" $q\.$dbn";
print "sequences.\n";
