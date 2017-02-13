#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use constant USAGE =><<END;

SYNOPSIS:

    run_asr.pl [OPTIONS] [INFILE]
    
DESCRIPTION:

INFILE is an erpin format RNA sequence alignment annotated with a secondary
structure.  Calls stomb, mb and rnanc and formats the output to stockholm
files which can be viewed using RALEE.

OPTIONS:

    --runmb 
           Run MrBayes. Only use if you require updated tree files. By default
           the program handles this sensibly (i.e. only runs if required files
           are absent).
    --n [0..inf]
           Number of sequences to sample from internal nodes and simulated
           branch tips. Default is 20.
    --l [0..inf]
           Lambda parameter, sets the branch-lengths for the positive predictive
           sequences (PSR). Default is 0.5.
    --h [0/1]
           When h=1 the harmonic-mean of all branch-lengths is used for the 
           PSR sequences. Lambda then becomes a scaling factor. [yes=1,no=0] 
           Default is 0.
    --i [0/1]
           Include original data in the output. [yes=1,no=0] Default is 1.
    --c [0/1]
           Condition on canonical pair frequencies. [yes=1,no=0] Default is 1.
    --debug
           Print debug info.

EXAMPLES:

run_asr.pl --n 5 --b 1 --l 0.2 tRNA00.erpin  
Samples 5 sequences from each internal node and each imaginary tip. Imaginary branch lengths are 0.2. 

DEPENDENCIES:

MrBayes v3.1.1 
sreformat (ships with Sean Eddy's squid package)

AUTHOR:

Jon P. Bollback
Paul P. Gardner

COPYRIGHT:

This program is free software. You may copy, modify, and redistribute
it under the same terms as Perl itself.

END

my $help = 0;
my $n = 20;
my $l = 0.5;
my $runmb = 0;
my $i = 1;
my $c = 1;
my $h = 0;
my $debug = 0;


GetOptions(
           "help" => \$help,
           "h"    => \$help,
	   "runmb" => \$runmb,
	   "n=s"    => \$n,
           "h=s"    => \$h,
	   "l=s"    => \$l,
	   "i=s"    => \$i,
	   "c=s"    => \$c,
	   "debug"  => \$debug
          ) or die USAGE;

$help and die USAGE;

@ARGV = ('-') unless @ARGV;
my $q = shift @ARGV;

if (!$q){
    die "Missing erpin file!";
}

if ($runmb || !(-e $q . ".stomb.con" && -e $q . ".stomb.stat")){
    print "\nRunning MrBayes....this may take a while...\n\n";
    system "stomb $q $q.stomb";
    system "mb $q.stomb";
    system "rm $q.stomb.t $q.stomb.trprobs $q.stomb.mcmc $q.stomb.parts $q.stomb.p ";
}

print "\nRunning rnanc....\nCommand: rnanc $q $q.stomb.con $q.stomb.stat $q.rnanc $l $n $i $c $debug\n
Input files: $q $q.stomb.con $q.stomb.stat
Output file: $q.rnanc
PSR branch-lengths (lambda): $l
Number of sequences to sample per-node/tip: $n
Include original data in output [yes=1,no=0]: $i
Condition on canonical pair frequencies [yes=1,no=0]: $c
Print debugging information [yes=1,no=0]: $debug
\n";



system "rnanc $q $q.stomb.con $q.stomb.stat $q.rnanc $h $l $n $i $c $debug";
system "cat $q.rnanc";

system "perl -p -n -e 's/\>SS_cons/\>\#\=GCXSS_cons/g;' $q.rnanc >$q.bak";
system "sreformat --pfam stockholm $q.bak | tr \"X()-\"  \" <>.\" > $q.rnanc.stk";
system "grep -v ^anc $q.rnanc.stk | tr \"()-\" \"<>.\" > $q.rnanc.seqNpsr.stk";
system "grep -v ^psr $q.rnanc.stk | tr \"()-\" \"<>.\" > $q.rnanc.seqNasr.stk";
system "rm $q.bak";

print "
For stockholm format output see:
$q.rnanc.stk
$q.rnanc.seqNpsr.stk
$q.rnanc.seqNasr.stk\n\n";

