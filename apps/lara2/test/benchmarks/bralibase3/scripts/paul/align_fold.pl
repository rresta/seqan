#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use constant USAGE =><<END;
#use Bio::SearchIO;
#use Bio::SeqIO;
use CMSequence;
use CMResults;

SYNOPSIS:

Take as input a fasta filename containing unaligned fasta files - align the
sequences using proalign and then fold using RNAalifold (both must be
installed and in your path -> ~/bin/ProAlign_0.5a1.jar). 

Returns Stockholm formatted output for Infernal:
filename.infernal.stk

Returns Stockholm formatted output for RSearch (refold here and/or use MFE 
structures):
filename.num.rsearch.stk

Returns erpin format:
filename.erpin.epn

DESCRIPTION:


OPTIONS:

EXAMPLES:


DEPENDENCIES:

use Bio::SearchIO;
use Bio::SeqIO;

These must be in your path (~/bin/):
ProAlign_0.5a1.jar
RNAalifold
sreformat
parent2epn.pl (from erpin4.2.5.server/scripts) 

AUTHOR:

Paul P. Gardner

COPYRIGHT:

This program is free software. You may copy, modify, and redistribute
it under the same terms as Perl itself.

END

my $help = 0;
my $q = "";

GetOptions(
           "help" => \$help,
	   "h"    => \$help,
	   "q=s"    => \$q
          ) or die USAGE;

$help and die USAGE;

if (!$q){
    die "Missing query sequence(s)!";
}

print "\#"x70 . "\n"."Aligning sequences in $q with proalign\n";


system "java -jar ~/bin/ProAlign_0.5a1.jar -nogui  -seqfile=$q  -newtree  -outfile=$q.pir -outformat=pir";

open my $fin, "$q.pir" or die "$q.pir: $!\n";

open my $fileout, ">".$q.".proalign" or die "$q.proalign: $!\n";

#Converting from fasta-like "pir" format to "fasta"
while ( my $l1 = <$fin> ) {

    chomp($l1);
    my @l1 = split /\t/, $l1;
        if (substr($l1,0,1) =~ '>'){
	$l1 =~ s/DL;//g;
    }
    else {
	$l1 =~ s/\*//g;
    }
    if ($l1) {
	printf $fileout "$l1\n";
	#print "$l1\n";
    }
}

close $fin;
close $fileout;
system "rm $q\.pir $q\.pir.min";

system "sreformat --pfam stockholm  $q.proalign";

system "sreformat --pfam stockholm $q.proalign \>$q.proalign.stk";

print "\n\nFolding $q\.proalign alignment with RNAalifold\n";
system "sreformat clustal $q.proalign | RNAalifold \>$q.proalign.aliout";

open $fin, "$q.proalign.aliout" or die "$q.proalign.aliout: $!\n";
my $cnt = 0;
my $ss = "";
my $rs_ss = "";
my $ess = "";

print "Structure:\n";
while ( my $l1 = <$fin> ) {

    print "$l1";
    
    if ($cnt){
	chomp($l1);
	my @l1 = split (/ /,$l1);

	$_ = $l1[0];
	tr/\(\)/\<\>/;
	$ss = $_;

	$_ = $l1[0];
	tr/\(\)/\>\</;
	$rs_ss = $_;

	$_ = $l1[0];
	tr/\./\-/;
	$ess = $_;
	
    }
    $cnt++;
}



#Add Alifold prediction to stockholm format file
open $fin, "$q.proalign.stk" or die "$q.proalign.stk: $!\n";
open $fileout, ">".$q.".proalign.struct.stk.tmp" or die "$q.proalign.struct.stk: $!\n";
open my $erpfileout, ">".$q.".erpin.fa" or die "$q.erpin.fa: $!\n";
printf $erpfileout "\>SS_cons\n$ess\n";

my $ssl = "#=GC SS_cons              ";
my $scount = 0;
my $pad = "";
while ( my $l1 = <$fin> ) {
    
        if (substr($l1,0,1) =~ /\#/ || substr($l1,0,1) =~ /\n/){
	    printf $fileout "$l1";
	}
        elsif (substr($l1,0,2) !~ /\/\//){
	    
	    my @l1 = split /\s+/, $l1;
	    my $diff = length($ssl)-length($l1[0]);
	    
	    if ($diff>=0){
		
		$pad = " " x $diff;
		printf $fileout $l1[0] .$pad . $l1[1] . "\n";

	    }
	    else {
		#printf STDERR "WARNING: long (>12 chars) sequence names aren't supported!\n\n";

		$pad = " ";
		printf $fileout $l1[0] .$pad . $l1[1] . "\n";
	    }
	    
	    #RSearch output
	    open my $rfileout, ">"."$q\.$scount\.rsearch.stk" or die "$q\.$scount\.rsearch.stk: $!\n";
	    my $rs_str = "";
	    my $rs_ss2 = "";
	    ($rs_str,$rs_ss2) = degap($l1[1],$rs_ss);
	    printf $rfileout "# STOCKHOLM 1.0

$l1[0]$pad$rs_str\n$ssl$rs_ss2\n//\n";
	    
		$scount++;
	    close $rfileout;
	    
	    #Erpin output
	    printf $erpfileout "\>$l1[0]\n$l1[1]\n";

	}
	else {
	    printf $fileout "$ssl$ss\n//\n";
	}
}
close $fileout;
close $erpfileout;

system "sreformat --pfam stockholm $q.proalign.struct.stk.tmp \>$q.infernal.stk; rm $q.proalign.struct.stk.tmp";

system "parent2epn.pl $q\.erpin.fa \>$q.erpin.epn";

print "

Created files:
$q\.proalign
$q.proalign.stk
$q\.\*\.rsearch.stk
$q.infernal.stk
$q.erpin.epn\n";

print "\#"x70 . "\n";

sub degap {
    
    my $str = shift;
    
    $str =~ s/\_\.\~/\-/g;
    my $ss = shift;
    
    chomp($str);
    chomp($ss);

    
    my @bp_tbl = make_pair_table( $ss );
    
    my @str = split //, $str;
    my @ss = split //, $ss;
    
    my $new_str = "";
    my $new_ss = "";
    my $i=0;

    
    if ($str =~ /\-/){
	while ( $i <= $bp_tbl[0] ){
	    
	    if ($str[$i] && $ss[$i]){
		if ($str[$i] =~ '\-'){
		    
		    if ($ss[$i] =~ '\>') {
			my $idx = $bp_tbl[$i+1];
			$ss[$idx-1] = '.';
			$ss[$i] = '.';
			
		    }
		    elsif ($ss[$i] =~ '\<') {
			my $idx = $bp_tbl[$i+1];
			$ss[$idx-1] = '.';
			$ss[$i] = '.';
		    }
		}
		
	    }	
	    
	$i++;
	    
	}
	
	$i = 0;
	while ( $i <= $bp_tbl[0] ){
	    
	    if ($str[$i] && $ss[$i]){
		if ($str[$i] =~ '\-'){
		    $str[$i] = "";
		    $ss[$i] = "";
		}
		
	    }	
	    
	    $i++;
	    
	}
	
	$new_str = join("",@str);
	$new_ss = join("",@ss);
    }
    else {
	$new_str = $str;
	$new_ss  = $ss;
    }
    
    printf STDOUT "\n$str\n$ss\nTO:\n$new_str\n$new_ss\n\n";
    return $new_str, $new_ss;
}

sub make_pair_table {
    my $str = shift;
    my $i;
    my @bps = ();
    my @bps_square = ();
    my @bps_curly = ();
    my @bps_angle = ();
    my @pair_table;
    my $prime5;
    my $prime3;
    my $count = 0;
    
    my $len = length($str);
    
    #print "make_pair_table: len = $len\n";

    $pair_table[0] = $len;

    for (my $j = 0; $j < $len; $j++ ) {
	$pair_table[$j+1] = 0;
	if ( substr($str,$j,1) =~ '\(' ){
	    push(@bps, $j);
	    ++$count;
	}	
	elsif ( substr($str,$j,1) =~ '\)' ){
	    $prime5 = pop(@bps);
	    $prime3 = $j;
	    $pair_table[$prime3+1] = $prime5+1;
	    $pair_table[$prime5+1] = $prime3+1;
	    --$count;
	    #print "base-pair: $pair_table[$i][$prime3] . $pair_table[$i][$prime5]\n"
	}
	elsif ( substr($str,$j,1) =~ '\[' ){
	    push(@bps_square, $j);
	    ++$count;
	}	
	elsif ( substr($str,$j,1) =~ '\]' ){
	    $prime5 = pop(@bps_square);
	    $prime3 = $j;
	    $pair_table[$prime3+1] = $prime5+1;
	    $pair_table[$prime5+1] = $prime3+1;
	    --$count;
	}
	elsif ( substr($str,$j,1) =~ '\{' ){
	    push(@bps_curly, $j);
	    ++$count;
	}	
	elsif ( substr($str,$j,1) =~ '\}' ){
	    $prime5 = pop(@bps_curly);
	    $prime3 = $j;
	    $pair_table[$prime3+1] = $prime5+1;
	    $pair_table[$prime5+1] = $prime3+1;
	    --$count;
	}
	elsif ( substr($str,$j,1) =~ '\>' ){
	    push(@bps_angle, $j);
	    ++$count;
	}	
	elsif ( substr($str,$j,1) =~ '\<' ){
	    $prime5 = pop(@bps_angle);
	    $prime3 = $j;
	    $pair_table[$prime3+1] = $prime5+1;
	    $pair_table[$prime5+1] = $prime3+1;
	    --$count;
	}
    }
    
#({[<.>]})
    
#print "pair table:\n";
#for ($i = 1; $i < $len+1; $i++ ) {
#    print "$pair_table[$i] ";
#}
    
#print "\n";
#
#print "pair_table = $#pair_table\n";    

if ($count){
    print "Unbalanced brackets in:
$str\n";
    die;
}    


    return @pair_table;
    
}
