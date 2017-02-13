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

Input is an epn format sequence alignment with a structure annotation and a
database. Searches the database for putative homologues.

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
my $pcw = 0.1;
GetOptions(
           "help" => \$help,
	   "h"    => \$help,
	   "q=s"    => \$q,
	   "db=s"   => \$db,
	   "pcw=s" => \$pcw
          ) or die USAGE;

$help and die USAGE;

if (!$q){
    die "Missing epn file!";
}

if (!$db){
    die "Missing database!";
}


my @dbn = split /\//, $db;
my $dbn = pop(@dbn);

open my $fin, "$q" or die "$q: $!\n";

my $sc = 0;
my $begin = "";
my $end = "";
while ( my $l1 = <$fin> ) {
    
    chomp($l1);
    #print "l1 = $l1\n";
    my @l1 = split /\t/, $l1;
    if (substr($l1,0,1) =~ '>'){
	$sc++;
    }
    elsif ($sc == 1) {
	my $l = length($l1);
	my $b = substr($l1,0,1);
	my $e = substr($l1,$l-1,1);
	#print " b = $b\n e = $e\n";
	$begin = $begin . $b;
	$end = $end . $e;
	
    }

}

my $outfile = "$q\.$dbn\_pcw". $pcw . ".out";
my $command = "/home/pgardner/bin/erpin $q $db   \-$begin,$end \-nomask \-fwd \-sum /home/pgardner/INST/erpin4.2.5.server/sum/SUM.dat -pcw $pcw >$outfile";

print "Searching sequences and structure in $q against $dbn with:
$command
\n";

system "$command";
system "tail -n5 $outfile";

