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

Input is a CM and a database. RaveNnA uses the CM to exactly filter the
database.

NB. Use full path names...

DESCRIPTION:


OPTIONS:

EXAMPLES:


DEPENDENCIES:

Making the database:
cp \$db data/
ls \$db > data/\$db\.list
release/cmzasha --no-stderr --gc-count data/\$db\.list
To "data/ravenna.config.tab" add (tab delimited):
db \$db data/\$db\.list data/EmblIdAndOrganismCompact_Barrick3.tab data/default.heurcreationspec

#

AUTHOR:


COPYRIGHT:

This program is free software. You may copy, modify, and redistribute
it under the same terms as Perl itself.

END

my $help = 0;
my $db = "";
my $q = "";
my $heur = 0;
my $heurl = 0; #local+heuristic
GetOptions(
           "help" => \$help,
	   "h"    => \$heur,
	   "hl"   => \$heurl,
	   "q=s"    => \$q,
	   "db=s"   => \$db
          ) or die USAGE;

$help and die USAGE;

if (!$q){
    die "Missing cm file!";
}

if (!$db){
    die "Missing database!";
}

my @qn = split /\//, $q;
my $qn1 = pop(@qn);
my $path = join("\/",@qn);

#Ugly hack as RaveNnA won't run if a filename begins with tRNA!!!!:
my $req = $path . "/not_" . $qn1;

my $qn2 = pop(@qn);

my @dbn = split /\//, $db;
my $dbn = pop(@dbn);

chdir("/home/pgardner/INST/ravenna-0.2f");

system "cp $q $req";

#system "perl ravenna.pl -scoreThreshold 0 -database $dbn -cmFileName $q -workDir data/t 200 -addToBaseName $db\_$qn2";

#my $cmd = "perl ravenna.pl -onlyForwardStrand -global -scoreThreshold 0.00001 -database $dbn -cmFileName $req -workDir $path 200";
my $cmd = "";
if($heur){
    $cmd = "perl src/ravenna.pl -filter ML -onlyForwardStrand -global -scoreThreshold 0.00001 -database $dbn -cmFileName $req -workDir $path 200";
}
elsif($heurl){
    $cmd = "perl src/ravenna.pl -filter ML -onlyForwardStrand -local -scoreThreshold 0.00001 -database $dbn -cmFileName $req -workDir $path 200";
}
else {
    $cmd = "perl src/ravenna.pl -onlyForwardStrand -scoreThreshold 0.00001 -database $dbn -cmFileName $req -workDir $path 200";
}
print "Running:\n$cmd\n\n";
system $cmd;


