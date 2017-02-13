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

Runs Infernal ver 0.7

Input is a stockholm format sequence alignment with a structure annotation and
a database. Builds a CM and uses this to search the database for putative
homologues.

DESCRIPTION:

OPTIONS:

EXAMPLES:

run_infernal7.pl -q query.stk -db database.fasta
run_infernal7.pl -local -q query.stk -db database.fasta


DEPENDENCIES:



AUTHOR:


COPYRIGHT:

This program is free software. You may copy, modify, and redistribute
it under the same terms as Perl itself.

END

my $help = 0;
my $db = "";
my $q = "";
my $local = 0;

GetOptions(
           "help" => \$help,
	   "h"    => \$help,
	   "local" => \$local,
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

print "Building CM using $q\n";
#system "touch $q.7.cm; rm $q.7.cm";
system "/home/pgardner/INST/infernal-0.7/src/cmbuild -F $q.7.cm $q >/dev/null";

print "Performing homology search of $q versus database $db\n";

my $outfile1 = "";
my $outfile2 = "";

if ($local){
    $outfile1 = "$q.$dbn.cmsearch7-local";
    $outfile2 = "$q.$dbn.infernal7-local";
    system "/home/pgardner/INST/infernal-0.7/src/cmsearch --toponly --local $q.7.cm $db >$outfile1";
}
else {
    $outfile1 = "$q.$dbn.cmsearch7";
    $outfile2 = "$q.$dbn.infernal7";
    system "/home/pgardner/INST/infernal-0.7/src/cmsearch --toponly $q.7.cm $db >$outfile1";
}


open( RES, "$outfile1" ) or die "can't read $outfile1";
    my $res = CMResults->new();
    $res -> parse_infernal( \*RES );
    $res = $res -> remove_overlaps();

open my $fileout, ">$outfile2" or die "$outfile2 $!\n";


my $count = 0;
    foreach my $unit ( sort { $b->bits <=> $a->bits } $res->eachUnit() ) {
	    my $outstring = sprintf( "%10s%8d%8d%8d%8d%10s", 
				  $unit->seqname, 
				  $unit->start_seq, 
				  $unit->end_seq, 
				  $unit->start_mod, 
				  $unit->end_mod, 
				  $unit->bits );

	    printf $fileout $outstring."\n";
	    
	    $count++;
	}

printf $fileout "Hit: $count sequences\n";

print "Hit: $count sequences\n";
print "Output:\n\t$outfile1\n\t$outfile2\n";

######################################################################

######### 

package CMResults;

sub new {
    my $ref = shift;
    my $class = ref($ref) || $ref;
    my $self = {
        'units' => [],
        'seq'   => {},
        'name'  => undef };
    bless $self, $class;
    return $self;
}

sub name {
    my $self  = shift;
    my $value = shift;
    $self->{'name'} = $value if( defined $value );
    return $self->{'name'};
}

sub addUnit {
    my $self = shift;
    my $unit = shift;
    my $name = $unit->seqname();

    if( !exists $self->{'seq'}->{$name} ) {
        warn "Adding a domain of $name but with no CMSequence. Will be kept in domain array but not added to a CMSequence";
    } else {
        $self->{'seq'}->{$name}->addUnit($unit);
    }
    push( @{$self->{'units'}}, $unit );
}

sub eachUnit {
    my $self = shift;
    return @{$self->{'units'}};
}

sub addSequence {
    my $self = shift;
    my $seq  = shift;
    my $name = $seq->name();
    if( exists $self->{'seq'}->{$name} ) {
        warn "You already have $name in CMResults. Replacing by a new entry!";
    }
    $self->{'seq'}->{$name} = $seq;
}

sub eachSequence {
    my $self = shift;
    my (@array,$name);
    foreach $name ( keys %{$self->{'seq'}} ) {
        push( @array, $self->{'seq'}->{$name} );
    }
    return @array;
}

sub getSequence {
    my $self = shift;
    my $name = shift;
    return $self->{'seq'}->{$name};
}

sub parse_infernal {
    my $self = shift;
    my $file = shift;

    my( $id, $start, $end, $ready, $modst, $moden );
    my $unit;  # this should always be the last added Unit

    while( <$file> ) {
        chomp;
        if( /^sequence:\s+(\S+)\s*/ ) {
            if( $1 =~ /^(\S+)\/(\d+)-(\d+)/ ) {
                ( $id, $start, $end ) = ( $1, $2, $3 );
            }
            elsif( ($id) = $1 =~ /^(\S+)/ ) {
                $start = 1;
            }
            else { 
                die "Don't recognise cmsearch output line [$_]";
            }
            unless( $self -> getSequence( $id ) ) {
                my $seq = CMSequence->new();
                $seq    -> name( $id );
                $self   -> addSequence( $seq );
            }
        }
        elsif( /^\s+$/ ) {
            $ready = 1;
        }
        elsif( /^hit\s+\d+\s*:\s+(.*)\s+bits/ ) {
            my $rest = $1;
            my( $st, $en, $bits );
            if( $rest =~ /(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/ ) {
                ( $st, $en, $modst, $moden, $bits ) = ( $1, $2, $3, $4, $5 );
            }
            elsif( $rest =~ /(\d+)\s+(\d+)\s+(\S+)/ ) {
                ( $st, $en, $bits ) = ( $1, $2, $3 );
            }
            else {
                warn "Don't recognise cmsearch output line [$_]";
            }

            $ready = 1;
        
            $st += $start - 1;
            $en += $start - 1;

            $unit = CMUnit->new();
            $unit -> seqname( $id );
            $unit -> modname( " " );
            $unit -> modacc( " " );
            $unit -> start_seq( $st );
            $unit -> end_seq( $en );
            $unit -> start_mod( $modst ) if $modst;
            $unit -> end_mod( $moden ) if $moden;
            $unit -> bits( $bits );
            $unit -> evalue( " " );

            $self -> addUnit( $unit );
        }
        elsif( /^\s+(\d+)\s+.*\s+(\d+)\s*$/ and $ready ) {
            # unit is already in results object, but this should still
            # get to where it needs to be
            $ready = 0;
            $unit -> start_mod( $1 ) unless $unit -> start_mod();
            $unit -> end_mod( $2 );
        }
    }
    return $self;
}    

sub remove_overlaps {
    my $self = shift;
    my $new = CMResults->new();
    foreach my $seq ( $self -> eachSequence() ) {
        my $newseq = CMSequence->new();
        $newseq -> name( $seq -> name() );
        $new -> addSequence( $newseq );

      UNIT:
	foreach my $unit1 ( sort { $b->bits <=> $a->bits } $seq -> eachUnit() ) {
	    foreach my $unit2 ( $newseq -> eachUnit() ) {
		if( ( $unit1->start_seq >= $unit2->start_seq and $unit1->start_seq <= $unit2->end_seq ) or
		    ( $unit1->end_seq   >= $unit2->start_seq and $unit1->end_seq   <= $unit2->end_seq ) or
		    ( $unit1->start_seq <= $unit2->start_seq and $unit1->end_seq   >= $unit2->end_seq ) ) {
		    next UNIT;
		}
	    }
	    $new -> addUnit( $unit1 );
	}
    }
    return $new;
}

sub filter_on_cutoff {
    my $self = shift;
    my $thr  = shift;
    my ($new,$seq,$unit,@array,@narray);

    if( !defined $thr ) {
        carp("CMResults: filter on cutoff needs an argument");
    }

    $new = CMResults->new();
    foreach $seq ( $self->eachSequence()) {
        my $newseq = CMSequence->new();
        $newseq->name($seq->name);
        $new->addSequence($newseq);
	foreach $unit ( $seq->eachUnit() ) {
	    if( $unit->bits() < $thr ) {
		next;
	    }
	    $new->addUnit($unit);
	}
    }
    return $new;
}

############## 

package CMSequence;

sub new {
    my $ref = shift;
    my $class = ref($ref) || $ref;
    my $self = {
        'name'   => undef,
        'units'  => [] };
    bless $self, $class;
    return $self;
}

sub name {
    my $self = shift;
    my $name = shift;

    if( defined $name ) {
        $self->{'name'} = $name;
    }
    return $self->{'name'};
}

sub addUnit {
    my $self = shift;
    my $unit = shift;
    push(@{$self->{'units'}},$unit); 
}

sub eachUnit {
    my $self = shift;
    return @{$self->{'units'}};
}

##############

package CMUnit;

sub new {
    my $ref = shift;
    my $class = ref($ref) || $ref;
    my $self = {
        seqname    => undef,
        start_seq  => undef,
	end_seq    => undef,
        modname    => undef,
        modacc     => undef,
        start_mod  => undef,
	end_mod    => undef,
        bits       => undef,
        evalue     => undef,
	};

    bless $self, $class;
    return $self;
}

sub seqname {
    my $self  = shift;
    my $value = shift;
    $self->{'seqname'} = $value if( defined $value );
    return $self->{'seqname'};
}

sub modname {
    my $self  = shift;
    my $value = shift;
    $self->{'modname'} = $value if( defined $value );
    return $self->{'modname'};
}

sub modacc {
    my $self  = shift;
    my $value = shift;
    $self->{'modacc'} = $value if( defined $value );
    return $self->{'modacc'};
}

sub bits {
    my $self  = shift;
    my $value = shift;
    $self->{'bits'} = $value if( defined $value );
    return $self->{'bits'};
}

sub evalue {
    my $self  = shift;
    my $value = shift;
    $self->{'evalue'} = $value if( defined $value );
    return $self->{'evalue'};
}

sub start_seq {
    my $self = shift;
    my $value = shift;
    $self->{'start_seq'} = $value if( defined $value );
    return $self->{'start_seq'};
}

sub end_seq {
    my $self = shift;
    my $value = shift;
    $self->{'end_seq'} = $value if( defined $value );
    return $self->{'end_seq'};
}

sub start_mod {
    my $self = shift;
    my $value = shift;
    $self->{'start_mod'} = $value if( defined $value );
    return $self->{'start_mod'};
}

sub end_mod {
    my $self = shift;
    my $value = shift;
    $self->{'end_mod'} = $value if( defined $value );
    return $self->{'end_mod'};
}

##############

