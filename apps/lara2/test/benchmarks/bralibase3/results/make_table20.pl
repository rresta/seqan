#!/usr/bin/perl

use strict;
use warnings;


my @seq5_shortnames = (
"NCBIblast_W11_top_24",
"NCBIblast_W7_top_22",
"NCBIblast_W7_541010_top_26",
"WUblast_W11_top_72",
"WUblast_W7_top_103",
"WUblast_W3_top_108",
"WUblast_W7_541010_top_86",
"WUblast_W7_pupy_top_92",
"fasta_z11_top_89",
"fasta_541010_z11_top_ktup6_85",
"fasta_U_z11_152",
"fasta_z11_R9540_1010_6",
"fasta_z11_sdf_1010_top_392",
"paralign_b1000_top_11",
"paralign_541010_b1000_top_91",
"ssearch_z11_top_97",
"ssearch_541010_z11_top_91",
"ssearch_U_z11_159",
"HMMer1_hmmfs_-6.14",
"HMMer1_hmmls_-1.16",
"HMMer2_hmmfs_-3",
"HMMer2_hmmls_-11.7",
"SAM35_sw2_-19.46",
"SAM35_sw0_-18.61",
"SAM35-HMMer2_-17.6",		       
"erp",
"infernal7_opt",
"infernal7_loc_opt",
 "infernal",
"ravennaML",
"ravennaML-local",
"rsearchD100s10n1000",
"RSmatch12_15"
);	

my @seq5_shortnamesWU = (
			 "NCBIblast_W7_541010_top_26",
			 "WUblast_W7_541010_top_86",
			 "fasta_541010_z11_top_ktup6_85",
			 "paralign_541010_b1000_top_91",
			 "ssearch_541010_z11_top_91"
			 );

my @seq5_shortnamesRNA = (
			  "WUblast_W11_top_72",
			  "WUblast_W7_pupy_top_92",
			  "fasta_z11_top_89",
			  "fasta_U_z11_152",
			  "fasta_z11_R9540_1010_6",
			  "fasta_z11_sdf_1010_top_392",
			  "ssearch_z11_top_97",
			  "ssearch_U_z11_159"
			  );

my %seq5_tps_files = ( 
"fasta_541010_z11_top_85" => "20/fasta_541010_z11_top_85",
"fasta_541010_z11_top_ktup6_85" => "20/fasta_541010_z11_top_ktup6_85",
"fasta_U_z11_152" => "20/fasta_U_z11_152",
"fasta_z11_R9540_1010_6" => "20/fasta_z11_R9540_1010_top_6",
"fasta_z11_sdf_1010_top_392" => "20/fasta_z11_sdf_1010_top_392",
"fasta_z11_top_89" => "20/fasta_z11_top_89",
"fasta_z11_top_ktup2_157" => "20/fasta_z11_top_ktup2_157",
"fasta_z11_top_ktup3_157" => "20/fasta_z11_top_ktup3_157",
"fasta_z11_top_ktup4_159" => "20/fasta_z11_top_ktup4_159",
"fasta_z11_top_ktup5_152" => "20/fasta_z11_top_ktup5_152",
"fasta_z11_top_ktup6_152" => "20/fasta_z11_top_ktup6_152",
"HMMer1_hmmfs_-6.14" => "20/HMMer1_hmmfs_-6.14",
"HMMer1_hmmls_-1.16" => "20/HMMer1_hmmls_-1.16",
"HMMer2_hmmfs_-3" => "20/HMMer2_hmmfs_-3",
"HMMer2_hmmls_-11.7" => "20/HMMer2_hmmls_-11.7",
"NCBIblast_W11_top_24" => "20/NCBIblast_W11_top_24",
"NCBIblast_W7_541010_top_26" => "20/NCBIblast_W7_541010_top_26",
"NCBIblast_W7_top_22" => "20/NCBIblast_W7_top_22",
"paralign_541010_b1000_top_91" => "20/paralign_541010_b1000_top_91",
"paralign_b1000_top_11" => "20/paralign_b1000_top_11",
"RSmatch12_15" => "20/RSmatch12_15",
"SAM35-HMMer2_-17.6"  => "20/SAM35-HMMer2_-17.6",
"SAM34-HMMer2_-18.2" => "20/SAM34-HMMer2_-18.2",
"SAM34_sw0_-18.43" => "20/SAM34_sw0_-18.43",
"SAM34_sw2_-18.53" => "20/SAM34_sw2_-18.53",
"SAM35_sw0_-18.61" => "20/SAM35_sw0_0306_-18.05",
"SAM35_sw2_-19.46" => "20/SAM35_sw2_0306_-6.28",
"ssearch_541010_z11_top_91" => "20/ssearch_541010_z11_top_91",
"ssearch_U_z11_159" => "20/ssearch_U_z11_159",
"ssearch_z11_top_97" => "20/ssearch_z11_top_97",
"WUblast_W11_top_72" => "20/WUblast_W11_top_72",
"WUblast_W3_top_108" => "20/WUblast_W3_top_108",
"WUblast_W7_541010_top_86" => "20/WUblast_W7_541010_top_86",
"WUblast_W7_ncbi_top_11" => "20/WUblast_W7_ncbi_top_11",
"WUblast_W7_pupy_top_92" => "20/WUblast_W7_pupy_top_92",
"WUblast_W7_top_103" => "20/WUblast_W7_top_103",
"erp" => "20/erpin_PCW0.1_20",
 "infernal" => "20/infernal",
"infernal7_opt" => "20/infernal7_7.47",
"infernal7_loc_opt" => "20/infernal7-local_11.74",
"ravenna" => "20/ravenna2f",
"ravennaML" => "20/ravenna2f-ML",
"ravennaML-local" => "20/ravenna2f-ML-local",
"rsearchD100s10n1000" => "20/rsearchD100s10n1000"
);

my %seq5_longnames = (
"fasta_541010_z11_top_85" => "\\fasta (65\\%%)",
"fasta_541010_z11_top_ktup6_85" => "\\fasta (ktup6,65\\%%)",
"fasta_U_z11_152" => "\\fasta (U)",
"fasta_z11_R9540_1010_6" => "\\fasta (RIBOSUM)",
"fasta_z11_sdf_1010_top_392" => "\\fasta (FOLDALIGN)",
"fasta_z11_top_89" => "\\fasta ",
"fasta_z11_top_ktup2_157" => "\\fasta (ktup2)",
"fasta_z11_top_ktup3_157" => "\\fasta (ktup3)",
"fasta_z11_top_ktup4_159" => "\\fasta (ktup4)",
"fasta_z11_top_ktup5_152" => "\\fasta (ktup5)",
"fasta_z11_top_ktup6_152" => "\\fasta (ktup6)",
"HMMer1_hmmfs_-6.14" => "\\hmmer (1.8.4,local)",
"HMMer1_hmmls_-1.16" => "\\hmmer (1.8.4,global)",
"HMMer2_hmmfs_-3" => "\\hmmer (2.3.2,local)",
"HMMer2_hmmls_-11.7" => "\\hmmer (2.3.2,global)",
"NCBIblast_W11_top_24" => "\\ncbiblast",
"NCBIblast_W7_541010_top_26" => "\\ncbiblast (W7,65\\%%)",
"NCBIblast_W7_top_22" => "\\ncbiblast (W7)",
"paralign_541010_b1000_top_91" => "\\paralign (65\\%%)",
"paralign_b1000_top_11" => "\\paralign",
"RSmatch12_15" => "\\rsmatch",
"SAM35-HMMer2_-17.6"  => "\\minitab[l]{\\sam (model) + \\\\ \\hmmer (2.3.2,search)}",
"SAM34-HMMer2_-18.2" => "\\minitab[l]{\\sam (3.4,model) + \\\\ \\hmmer (2.3.2,search)}",
"SAM34_sw0_-18.43" => "\\sam (3.4,global)",
"SAM34_sw2_-18.53" => "\\sam (3.4,local)",
"SAM35_sw0_-18.61" => "\\sam (global)",
"SAM35_sw2_-19.46" => "\\sam (local)",
"ssearch_541010_z11_top_91" => "\\ssearch (65\\%%)",
"ssearch_U_z11_159" => "\\ssearch (U)",
"ssearch_z11_top_97" => "\\ssearch ",
"WUblast_W11_top_72" => "\\wublast",
"WUblast_W3_top_108" => "\\wublast (W3)",
"WUblast_W7_541010_top_86" => "\\wublast (W7,65\\%%)",
"WUblast_W7_ncbi_top_11" => "\\wublast (W7,99\\%%)",
"WUblast_W7_pupy_top_92" => "\\wublast (W7,PUPY)",
"WUblast_W7_top_103" => "\\wublast (W7)",
"erp" => "\\erpin",
"infernal" => "\\infernal (0.55)",
"infernal7_opt" => "\\infernal (0.7)",
"infernal7_loc_opt" => "\\infernal (0.7,local)",
"ravenna" => "\\ravenna",
"ravennaML" => "\\ravenna (ML)",
"ravennaML-local" => "\\ravenna (ML,local)",
"rsearchD100s10n1000"  => "\\rsearch"
);


 my %seq5_toprint = %seq5_tps_files;


my %dp_size = qw(
	       rRNA	602
	       tRNA	1114
	       U5	235
	       );

my %seq5_tps_vals = ();
my %seq5_sens_vals = ();
my %seq5_fps_vals = ();
my %seq5_spec_vals = ();
my %seq5_mcc_vals = ();


#while (my $file = shift(@seq5_tps_files)) {

foreach my $filekey ( keys %seq5_tps_files ) {
    
    my $file = $seq5_tps_files{$filekey};

    my $tpfile = $file . ".tps";
    my $fpfile = $file . ".fps";

#    printf STDERR "Opening $tpfile & $fpfile\n";
    open my $tpfileptr, "$tpfile" or die "Couldn't open $tpfile\n";
    open my $fpfileptr, "$fpfile" or die "Couldn't open $fpfile\n";
    
    while ( (my $tpl = <$tpfileptr>) && (my $fpl = <$fpfileptr>) ){
	
	chomp($tpl);
	chomp($fpl);
	
	my @tpline = split /\s+/, $tpl;

	if (!isnumeric($tpline[3])){
	    print "ERROR: $tpfile has erroneous entry:\n$tpl\n";
	}
      	my $sens = $tpline[3];
	my $tp = $sens * $dp_size{$tpline[0]};
	my $fn = $dp_size{$tpline[0]} - $tp;
	
	push @{ $seq5_tps_vals{$filekey} }, $tp;
	push @{ $seq5_sens_vals{$filekey} }, $sens;

	my @fpline = split /\s+/, $fpl;
      	my $spec = (1 - $fpline[3]);
	my $fp = $fpline[3] * $dp_size{$fpline[0]} * 10;
	my $tn = $dp_size{$fpline[0]} * 10 - $fp;
		
	push @{ $seq5_fps_vals{$filekey} }, $fp;
	push @{ $seq5_spec_vals{$filekey} }, $spec;
	
	my $numerator = ($tp + $fp)*($tp + $fn)*($tn + $fp)*($tn + $fn);
	my $mcc = 0;
	if ($numerator){
	    $mcc = ($tp * $tn - $fp * $fn)/sqrt($numerator );
	}
	push @{ $seq5_mcc_vals{$filekey} }, $mcc;
	
	
	
	#printf STDOUT "$tpl\t$spec\t%0.0f\t%0.0f\t%0.0f\t%0.0f\t%0.2f\n", $tp, $tn, $fp, $fn, $mcc;
	
	if ($tpline[2] ne $fpline[2] ){
	    die "$tpline[2] is not equal to $fpline[2]\n in $tpfile & $fpfile\n";
	}
    }

}

my %seq5_mean_sens = ();
my %seq5_mean_spec = ();
my %seq5_mean_mcc = ();
my %seq5_var_sens = ();
my %seq5_var_spec = ();
my %seq5_var_mcc = ();

my %seq5_median_mcc = ();
my %seq5_uq_mcc = ();
my %seq5_lq_mcc = ();

my %seq5_median_sens = ();
my %seq5_uq_sens = ();
my %seq5_lq_sens = ();

my %seq5_median_spec = ();
my %seq5_uq_spec = ();
my %seq5_lq_spec = ();


foreach my $family ( keys %seq5_mcc_vals ) {
    #print "$family: ";
    my @mcc_sorted = ();
    my @sens_sorted = ();
    my @spec_sorted = ();
    my $mean_mcc = 0;
    my $mean_sens = 0;
    my $mean_spec = 0;
    my $var_mcc = 0;
    my $var_sens = 0;
    my $var_spec = 0;
    my $cnt = 0;
    foreach my $i ( 0 .. $#{ $seq5_mcc_vals{$family} } ) {
        $mean_mcc += $seq5_mcc_vals{$family}[$i];
	$mean_sens += $seq5_sens_vals{$family}[$i];
	$mean_spec += $seq5_spec_vals{$family}[$i];
        $var_mcc += $seq5_mcc_vals{$family}[$i]*$seq5_mcc_vals{$family}[$i];
	$var_sens += $seq5_sens_vals{$family}[$i]*$seq5_sens_vals{$family}[$i];
	$var_spec += $seq5_spec_vals{$family}[$i]*$seq5_spec_vals{$family}[$i];
	$cnt++;
	
	push(@mcc_sorted, $seq5_mcc_vals{$family}[$i]);
	push(@sens_sorted, $seq5_sens_vals{$family}[$i]);
	push(@spec_sorted, $seq5_spec_vals{$family}[$i]);

    }
    #$mean_mcc = $mean_mcc/($cnt);
    $seq5_mean_mcc{$family} = $mean_mcc/($cnt);
    $seq5_mean_sens{$family} = $mean_sens/($cnt);
    $seq5_mean_spec{$family} = $mean_spec/($cnt);

    $seq5_var_mcc{$family}  = $var_mcc/($cnt)  - ($mean_mcc/($cnt)) * ($mean_mcc/($cnt));
    $seq5_var_sens{$family} = $var_sens/($cnt) - ($mean_sens/($cnt)) * ($mean_sens/($cnt));
    $seq5_var_spec{$family} = $var_spec/($cnt) - ($mean_spec/($cnt)) * ($mean_spec/($cnt));
    
    
    #sample size & quartiles are even, 90, 180, 270:
    @mcc_sorted = sort{$a <=> $b} @mcc_sorted;
    my $mid = sprintf("%.0f", 0.5 * ($#{ $seq5_mcc_vals{$family} }+1) );
    $seq5_median_mcc{$family} = ($mcc_sorted[$mid-1]+$mcc_sorted[$mid])/2;
    my $uq = sprintf("%.0f", 0.75 * ($#{ $seq5_mcc_vals{$family} }+1) );
    $seq5_uq_mcc{$family} = ($mcc_sorted[$uq-1]+$mcc_sorted[$uq])/2;
    my $lq = sprintf("%.0f", 0.25 * ($#{ $seq5_mcc_vals{$family} }+1) );
    $seq5_lq_mcc{$family} = ($mcc_sorted[$lq-1]+$mcc_sorted[$lq])/2;
   
    @sens_sorted = sort{$a <=> $b} @sens_sorted;
    $seq5_median_sens{$family} = ($sens_sorted[$mid-1]+$sens_sorted[$mid])/2;
    $seq5_uq_sens{$family} = ($sens_sorted[$uq-1]+$sens_sorted[$uq])/2;
    $seq5_lq_sens{$family} = ($sens_sorted[$lq-1]+$sens_sorted[$lq])/2;

    @spec_sorted = sort{$a <=> $b} @spec_sorted;
    $seq5_median_spec{$family} = ($spec_sorted[$mid-1]+$spec_sorted[$mid])/2;
    $seq5_uq_spec{$family} = ($spec_sorted[$uq-1]+$spec_sorted[$uq])/2;
    $seq5_lq_spec{$family} = ($spec_sorted[$lq-1]+$spec_sorted[$lq])/2;

#    printf STDOUT "mcc: %0.3f\t%0.3f\t%0.3f\t$family\n", $seq5_lq_mcc{$family}, $seq5_median_mcc{$family},$seq5_uq_mcc{$family};
#    printf STDOUT "sens:%0.3f\t%0.3f\t%0.3f\t$family\n", $seq5_lq_sens{$family}, $seq5_median_sens{$family},$seq5_uq_sens{$family};
#    printf STDOUT "spec:%0.3f\t%0.3f\t%0.3f\t$family\n", $seq5_lq_spec{$family}, $seq5_median_spec{$family},$seq5_uq_spec{$family};
}

#print "\n";

my %tmp_mcc_hash = ();
my %seq5_mcc_ranks = ();
my %seq5_mcc_ranks_save = ();
my $sample_size = $#{ $seq5_mcc_vals{"erp"} };
my $old_family = "";

for (my $i=0; $i<=$sample_size; $i++){
  
    foreach my $family ( keys %seq5_mcc_vals ) {
	$tmp_mcc_hash{$family} = $seq5_mcc_vals{$family}[$i];
    }
    my @sorted_by_mcc = sort by_mcc @seq5_shortnames;
    
    my $cnt = 1;
    my $old_cnt = 1;
    foreach my $family (@sorted_by_mcc) {
	
#	if((defined $seq5_mcc_vals{$family}[$i]) && (defined $seq5_mcc_vals{$old_family}[$i])){
	if ($seq5_mcc_vals{$family}[$i] == $seq5_mcc_vals{$old_family}[$i]){
	    $seq5_mcc_ranks{$family} += $old_cnt;
	    push @{ $seq5_mcc_ranks_save{$family} }, $old_cnt;
#	    	print "$family\t$old_cnt\n";

	}
	else {
	    if (defined($seq5_mcc_ranks{$family})) {
		$seq5_mcc_ranks{$family} += $cnt;
		$old_cnt = $cnt;
		push @{ $seq5_mcc_ranks_save{$family} }, $cnt;
#		print "$family\t$cnt\n";

	    }
	    else {
		$seq5_mcc_ranks{$family} = $cnt;
		$old_cnt = $cnt;
		push @{ $seq5_mcc_ranks_save{$family} }, $cnt;
#		print "$family\t$cnt\n";

	    }
	}
 #   }
# 	if ($seq5_mcc_vals{$family}[$i] == $seq5_mcc_vals{$old_family}[$i]){
# 	    print "MCC $old_family = [$seq5_mcc_vals{$old_family}[$i]]\n";
# 	    print "MCC $family = [$seq5_mcc_vals{$family}[$i]]\n";
# 	}
	$old_family = $family;
	
	$cnt++;

    }
    
#    foreach my $family (@seq5_shortnames) {
#	print "$family\t$seq5_mcc_vals{$family}[$i]\t$seq5_mcc_ranks{$family} \n";
#    }
    
}


#     foreach my $family (@seq5_shortnames) {
# 	print "$family\t$seq5_mcc_ranks{$family} \n";
#     }

#printf STDOUT "Name\tSensitivity\tSpecificity\tMCC\tMean MCC Rank\n\n";

my $phmm = 1;
my $smeth = 1;
my $ssmeth = 1;


open my $out, ">table20.tex" or die "table20.tex: $!\n";




	printf $out "\\begin{table}[hbt]
\\centering
{\\footnotesize
\\begin{tabular}{||l|c|c|c|c|c||}
\\hline
\\hline
{\\bf Program} & {\\bf Sensitivity}   & {\\bf Specificity}   & {\\bf MCC}   & \\multicolumn{2}\{|c||}{\\bf Ave. MCC Rank} \\\\
                 &                      &                      &              & Median & Mean \\\\
\\hline

";

foreach my $family ( @seq5_shortnames ) {
    
    #printf STDOUT "$seq5_longnames{$family} \& %0.2f \$\\pm\$ %0.2f \& %0.4f \$\\pm\$ %0.2f \& %0.2f \$\\pm\$ %0.2f \& %0.2f\\\\ \n", $seq5_mean_sens{$family}, sqrt($seq5_var_sens{$family}), $seq5_mean_spec{$family}, sqrt($seq5_var_spec{$family}), $seq5_mean_mcc{$family}, sqrt($seq5_var_mcc{$family}), $seq5_mcc_ranks{$family}/$sample_size;

    if (($family =~ /NCBI/ || $family =~ /fasta/ ) && $ssmeth){
	printf $out "\\hline
\\multicolumn\{6\}\{\|\|c\|\|\}\{\\bf Sequence based methods\}\\\\
\\hline
";
	$ssmeth = 0;
    }
    elsif (($family =~ /hmmer/ || $family =~ /HMMer/) && $phmm){
	printf $out "\\hline
\\multicolumn\{6\}\{\|\|c\|\|\}\{\\bf Profile HMM methods\}\\\\
\\hline
";
	$phmm = 0;
    }
    elsif ($family =~ /erp/ && $smeth){
	printf $out "\\hline
\\multicolumn\{6\}\{\|\|c\|\|\}\{\\bf Structure enhanced methods\}\\\\
\\hline
";
	$smeth = 0;
    }
    @{ $seq5_mcc_ranks_save{$family} } = sort{$a <=> $b} @{ $seq5_mcc_ranks_save{$family} };

my $mid = sprintf("%.0f", 0.5 * $sample_size );
my $med = (${ $seq5_mcc_ranks_save{$family} }[$mid-1] + ${ $seq5_mcc_ranks_save{$family} }[$mid])/2;

#    if ($seq5_mcc_ranks{$family}/$sample_size <= 10){
    if ($med <= 10){
printf $out "$seq5_longnames{$family} \& (%0.2f,\\textbf\{%0.2f,\}%0.2f) \& (%0.3f,\\textbf\{%0.3f\},%0.3f) \& (%0.2f,\\textbf\{%0.2f\},%0.2f) \& \\textbf\{%0.1f\} \& %0.2f \\\\ \n", 
$seq5_lq_sens{$family}, $seq5_median_sens{$family},$seq5_uq_sens{$family}, 
$seq5_lq_spec{$family}, $seq5_median_spec{$family},$seq5_uq_spec{$family}, 
$seq5_lq_mcc{$family}, $seq5_median_mcc{$family},$seq5_uq_mcc{$family}, 
$med,
$seq5_mcc_ranks{$family}/$sample_size;
}
    else {
printf $out "$seq5_longnames{$family} \& (%0.2f,\\textbf\{%0.2f,\}%0.2f) \& (%0.3f,\\textbf\{%0.3f\},%0.3f) \& (%0.2f,\\textbf\{%0.2f\},%0.2f) \& %0.1f \& %0.2f\\\\ \n", 
$seq5_lq_sens{$family}, $seq5_median_sens{$family},$seq5_uq_sens{$family}, 
$seq5_lq_spec{$family}, $seq5_median_spec{$family},$seq5_uq_spec{$family}, 
$seq5_lq_mcc{$family}, $seq5_median_mcc{$family},$seq5_uq_mcc{$family}, 
$med,
$seq5_mcc_ranks{$family}/$sample_size;
	
    }

}

	printf $out "\\hline
\\hline
\\end{tabular}
\\caption[]{ Tabulated results of searches with {\\bf twenty} sequences. From
left to right column one contains program names and settings, column two sensitivity,
column three specificity, column four Matthew's correlation coefficient (MCC)
and column five contains both a median and mean ranking determined
by the MCC (see Box 2 for definitions). Sensitivity, specificity and MCC are summarised as 3-tuples
displaying the lower quartiles, medians and upper quartiles. Median MCC
rankings below ten (one being the best) are indicated with bold font.
See Methods section and Supplementary Tables~\\ref{table:homol:alg}\\&\\ref{table:homol:alg2} 
for algorithm settings.
}\\label{table:20seqs}
}
\\end{table}
";

close($out);

system "cat table20.tex";

foreach my $family ( keys %seq5_toprint ) {
my $outfile = $seq5_toprint{$family} . ".mcc";
open my $outf, ">$outfile" or die "$outfile: $!\n";
my $outfiler = $seq5_toprint{$family} . ".ranks";
open my $outrank, ">$outfiler" or die "$outfiler: $!\n";

#printf $outf "$family\n";
foreach my $imcc (@{ $seq5_mcc_vals{$family} }) {
printf $outf "$imcc\n";
}
foreach my $irank (@{ $seq5_mcc_ranks_save{$family} }) {
printf $outrank "$irank\n";
}

}


##############################
open $out, ">tableWU20.tex" or die "tableWU20.tex: $!\n";

$ssmeth = 1;
foreach my $family ( @seq5_shortnamesWU ) {
    
    if ($family =~ /NCBI/ && $ssmeth){
	printf $out "\\hline
\\multicolumn\{6\}\{\|\|c\|\|\}\{\\bf Twenty input sequences \}\\\\
\\hline
";
	$ssmeth = 0;
    }

my $mid = sprintf("%.0f", 0.5 * $sample_size );
my $med = (${ $seq5_mcc_ranks_save{$family} }[$mid-1] + ${ $seq5_mcc_ranks_save{$family} }[$mid])/2;

    if ($seq5_mcc_ranks{$family}/$sample_size <= 10){
printf $out "$seq5_longnames{$family} \& (%0.2f,\\textbf\{%0.2f,\}%0.2f) \& (%0.3f,\\textbf\{%0.3f\},%0.3f) \& (%0.2f,\\textbf\{%0.2f\},%0.2f) \& \\textbf\{%0.1f\} \& %0.2f \\\\ \n", 
$seq5_lq_sens{$family}, $seq5_median_sens{$family},$seq5_uq_sens{$family}, 
$seq5_lq_spec{$family}, $seq5_median_spec{$family},$seq5_uq_spec{$family}, 
$seq5_lq_mcc{$family}, $seq5_median_mcc{$family},$seq5_uq_mcc{$family}, 
$med,
$seq5_mcc_ranks{$family}/$sample_size;
}
    else {
printf $out "$seq5_longnames{$family} \& (%0.2f,\\textbf\{%0.2f,\}%0.2f) \& (%0.3f,\\textbf\{%0.3f\},%0.3f) \& (%0.2f,\\textbf\{%0.2f\},%0.2f) \& %0.2f \& %0.2f\\\\ \n", 
$seq5_lq_sens{$family}, $seq5_median_sens{$family},$seq5_uq_sens{$family}, 
$seq5_lq_spec{$family}, $seq5_median_spec{$family},$seq5_uq_spec{$family}, 
$seq5_lq_mcc{$family}, $seq5_median_mcc{$family},$seq5_uq_mcc{$family}, 
$med,
$seq5_mcc_ranks{$family}/$sample_size;
	
    }
}

	printf $out "\\hline
\\hline
\\end{tabular}
\\caption[]{ 
Tabulated results of the sequence based methods with identical scoring
parameters. From left to right column one contains program names, column
two sensitivity, column three specificity, column four Matthew's correlation
coefficient (MCC) and column five a median and mean ranking determined by the
MCC. Sensitivity, specificity and MCC are summarised as 3-tuples displaying
the lower quartiles, medians and upper quartiles.
}\\label{table:WU}
}
\\end{table}
";

system "cat tableWU5.tex tableWU20.tex > tableWU.tex";
#############################

open $out, ">tableRNA20.tex" or die "tableRNA20.tex: $!\n";

$ssmeth = 1;
foreach my $family ( @seq5_shortnamesRNA ) {
    

    if ($family =~ /WUblast/ && $ssmeth){
	printf $out "\\hline
\\multicolumn\{6\}\{\|\|c\|\|\}\{\\bf Twenty input sequences \}\\\\
\\hline
";
	$ssmeth = 0;
    }

my $mid = sprintf("%.0f", 0.5 * $sample_size );
my $med = (${ $seq5_mcc_ranks_save{$family} }[$mid-1] + ${ $seq5_mcc_ranks_save{$family} }[$mid])/2;

    if ($seq5_mcc_ranks{$family}/$sample_size <= 10){
printf $out "$seq5_longnames{$family} \& (%0.2f,\\textbf\{%0.2f,\}%0.2f) \& (%0.3f,\\textbf\{%0.3f\},%0.3f) \& (%0.2f,\\textbf\{%0.2f\},%0.2f) \& \\textbf\{%0.1f\} \& %0.2f\\\\ \n", 
$seq5_lq_sens{$family}, $seq5_median_sens{$family},$seq5_uq_sens{$family}, 
$seq5_lq_spec{$family}, $seq5_median_spec{$family},$seq5_uq_spec{$family}, 
$seq5_lq_mcc{$family}, $seq5_median_mcc{$family},$seq5_uq_mcc{$family}, 
$med,
$seq5_mcc_ranks{$family}/$sample_size;
}
    else {
printf $out "$seq5_longnames{$family} \& (%0.2f,\\textbf\{%0.2f,\}%0.2f) \& (%0.3f,\\textbf\{%0.3f\},%0.3f) \& (%0.2f,\\textbf\{%0.2f\},%0.2f) \& %0.1f \& %0.2f\\\\ \n", 
$seq5_lq_sens{$family}, $seq5_median_sens{$family},$seq5_uq_sens{$family}, 
$seq5_lq_spec{$family}, $seq5_median_spec{$family},$seq5_uq_spec{$family}, 
$seq5_lq_mcc{$family}, $seq5_median_mcc{$family},$seq5_uq_mcc{$family}, 
$med,
$seq5_mcc_ranks{$family}/$sample_size;
	
    }
}

	printf $out "\\hline
\\hline
\\end{tabular}
\\caption[]{ 
Tabulated results comparing sequence based methods with and without
RNA specific scoring parameters.  From left to right column one contains
program names and settings, column two sensitivity, column three specificity, column four
Matthew's correlation coefficient (MCC) and column five a median and mean 
ranking determined by the MCC. Sensitivity, specificity and MCC are summarised
as 3-tuples displaying the lower quartiles, medians and upper quartiles.
}\\label{table:RNA}
}
\\end{table}
";


close($out);

system "cat tableRNA5.tex tableRNA20.tex > tableRNA.tex";

#############################

sub by_mcc {
    $tmp_mcc_hash{$b} <=> $tmp_mcc_hash{$a};
}


sub isnumeric {
    my $num = shift;
    
    if ($num  !~  /^[0-9|.|,]*$/) {
    		return 0;
	    }
    else {
	return 1;	
    }	
}
