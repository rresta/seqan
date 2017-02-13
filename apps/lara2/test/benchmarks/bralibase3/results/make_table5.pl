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
"fasta34_top_89.rnanc.seq10asr",
"fasta_541010_z11_top_ktup6_85",
"fasta_U_z11_152",
"fasta_z11_R9540_1010_6",
"fasta_z11_sdf_1010_top_392",
"paralign_b1000_top_11",
"paralign_541010_b1000_top_91",
"ssearch_z11_top_97",
"ssearch34_top_97.rnanc.seq10asr",
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
"erpPSR",
		       "infernal7_opt",
		       "infernal7_loc_opt",
 		       "infernal",
"ravennaML",
"ravennaML-local",
"rsearch",
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
"fasta_541010_z11_top_85" => "5/fasta_541010_z11_top_85",
"fasta_541010_z11_top_ktup6_85" => "5/fasta_541010_z11_top_ktup6_85",
"fasta_U_z11_152" => "5/fasta_U_z11_152",
"fasta_z11_R9540_1010_6" => "5/fasta_z11_R9540_1010_top_6",
"fasta_z11_sdf_1010_top_392" => "5/fasta_z11_sdf_1010_top_392",
"fasta_z11_top_89" => "5/fasta_z11_top_89",
"fasta34_top_89.rnanc.seq10asr" => "5/fasta34_top_89.rnanc.seq10asr",
"fasta_z11_top_ktup2_157" => "5/fasta_z11_top_ktup2_157",
"fasta_z11_top_ktup3_157" => "5/fasta_z11_top_ktup3_157",
"fasta_z11_top_ktup4_159" => "5/fasta_z11_top_ktup4_159",
"fasta_z11_top_ktup5_152" => "5/fasta_z11_top_ktup5_152",
"fasta_z11_top_ktup6_152" => "5/fasta_z11_top_ktup6_152",
"HMMer1_hmmfs_-6.14" => "5/HMMer1_hmmfs_-6.14",
"HMMer1_hmmls_-1.16" => "5/HMMer1_hmmls_-1.16",
"HMMer2_hmmfs_-3" => "5/HMMer2_hmmfs_-3",
"HMMer2_hmmls_-11.7" => "5/HMMer2_hmmls_-11.7",
"NCBIblast_W11_top_24" => "5/NCBIblast_W11_top_24",
"NCBIblast_W7_541010_top_26" => "5/NCBIblast_W7_541010_top_26",
"NCBIblast_W7_top_22" => "5/NCBIblast_W7_top_22",
"paralign_541010_b1000_top_91" => "5/paralign_541010_b1000_top_91",
"paralign_b1000_top_11" => "5/paralign_b1000_top_11",
"RSmatch12_15" => "5/RSmatch12_15",
"SAM35-HMMer2_-17.6"  => "5/SAM35-HMMer2_-17.6",
"SAM34-HMMer2_-18.2" => "5/SAM34-HMMer2_-18.2",
"SAM34_sw0_-18.43" => "5/SAM34_sw0_-18.43",
"SAM34_sw2_-18.53" => "5/SAM34_sw2_-18.53",
"SAM35_sw0_-18.61" => "5/SAM35_sw0_0306_-18.05",
"SAM35_sw2_-19.46" => "5/SAM35_sw2_0306_-6.28",
"ssearch_U_z11_159" => "5/ssearch_U_z11_159",
"ssearch_541010_z11_top_91" => "5/ssearch_541010_z11_top_91",
"ssearch_z11_top_97" => "5/ssearch_z11_top_97",
"ssearch34_top_97.rnanc.seq10asr" => "5/ssearch34_top_97.rnanc.seq10asr",
"WUblast_W11_top_72" => "5/WUblast_W11_top_72",
"WUblast_W3_top_108" => "5/WUblast_W3_top_108",
"WUblast_W7_541010_top_86" => "5/WUblast_W7_541010_top_86",
"WUblast_W7_ncbi_top_11" => "5/WUblast_W7_ncbi_top_11",
"WUblast_W7_pupy_top_92" => "5/WUblast_W7_pupy_top_92",
"WUblast_W7_top_103" => "5/WUblast_W7_top_103",
		       "erp" => "5/erpin_PCW0.1",
		       "erpPSR" => "5/erpin_PSR",
 		       "infernal" => "5/infernal",
		       "infernal7_opt" => "5/infernal7_7.47",
		       "infernal7_loc_opt" => "5/infernal7-local_11.74",
"ravennaML" => "5/ravenna2f-ML",
"ravennaML-local" => "5/ravenna2f-ML-local",
"rsearch" => "5/rsearch",
		       );

my %seq5_longnames = (
"fasta_541010_z11_top_85" => "\\fasta (65\\%%)",
"fasta_541010_z11_top_ktup6_85" => "\\fasta (ktup6,65\\%%)",
"fasta_U_z11_152" => "\\fasta (U)",
"fasta_z11_R9540_1010_6" => "\\fasta (RIBOSUM)",
"fasta_z11_sdf_1010_top_392" => "\\fasta (FOLDALIGN)",
"fasta_z11_top_89" => "\\fasta",
"fasta34_top_89.rnanc.seq10asr" => "\\fasta [ASR]",
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
"SAM34_sw0_-18.43" => "\\sam (3.4,global)",
"SAM34_sw2_-18.53" => "\\sam (3.4,local)",
"SAM35_sw0_-18.61" => "\\sam (global)",
"SAM35_sw2_-19.46" => "\\sam (local)",
"SAM35-HMMer2_-17.6"  => "\\minitab[l]{\\sam (model) + \\\\ \\hmmer (2.3.2,search)}",
"SAM34-HMMer2_-18.2" => "\\minitab[l]{\\sam (3.4,model) + \\\\ \\hmmer (2.3.2,search)}",
"ssearch_U_z11_159" => "\\ssearch (U)",
"ssearch_541010_z11_top_91" => "\\ssearch (65\\%%)",
"ssearch_z11_top_97" => "\\ssearch",
"ssearch34_top_97.rnanc.seq10asr" => "\\ssearch [ASR]",
"WUblast_W11_top_72" => "\\wublast",
"WUblast_W3_top_108" => "\\wublast (W3)",
"WUblast_W7_541010_top_86" => "\\wublast (W7,65\\%%)",
"WUblast_W7_ncbi_top_11" => "\\wublast (W7,99\\%%)",
"WUblast_W7_pupy_top_92" => "\\wublast (W7,PUPY)",
"WUblast_W7_top_103" => "\\wublast (W7)",
		      "erp" => "\\erpin",
		      "erpPSR" => "\\erpin [PSR]",
 		      "infernal" => "\\infernal (0.55)",
		       "infernal7_opt" => "\\infernal (0.7)",
		      "infernal7_loc_opt" => "\\infernal (0.7,local)",
"ravennaML" => "\\ravenna (ML)",
"ravennaML-local" => "\\ravenna (ML,local)",
"rsearch" => "\\rsearch",
"RSmatch12_15" => "\\rsmatch"
		      );


my %seq5_toprint = %seq5_tps_files;



my %dp_size = qw(
	       rRNA	602
	       tRNA	1114
	       U5	235
	       );

#my $name_idx = 0;
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

# 	push @{ $seq5_tps_vals{$seq5_shortnames[$name_idx]} }, $tp;
# 	push @{ $seq5_sens_vals{$seq5_shortnames[$name_idx]} }, $sens;

	my @fpline = split /\s+/, $fpl;
	if ($fpline[3]<0 || $fpline[3]>1){
	    print "ERROR: $tpfile has erroneous entry:\n$tpl\n";
	}
      	my $spec = (1 - $fpline[3]);
	my $fp = $fpline[3] * $dp_size{$fpline[0]} * 10;
	my $tn = $dp_size{$fpline[0]} * 10 - $fp;
	
	
	push @{ $seq5_fps_vals{$filekey} }, $fp;
	push @{ $seq5_spec_vals{$filekey} }, $spec;
# 	push @{ $seq5_fps_vals{$seq5_shortnames[$name_idx]} }, $fp;
# 	push @{ $seq5_spec_vals{$seq5_shortnames[$name_idx]} }, $spec;
	
	my $numerator = ($tp + $fp)*($tp + $fn)*($tn + $fp)*($tn + $fn);
	my $mcc = 0;
	if ($numerator){
	    $mcc = ($tp * $tn - $fp * $fn)/sqrt($numerator);
	}


	push @{ $seq5_mcc_vals{$filekey} }, $mcc;
#	push @{ $seq5_mcc_vals{$seq5_shortnames[$name_idx]} }, $mcc;
	
	
	#print "$filekey\t$sens\t$spec\t$mcc\n";
	
	#printf STDOUT "$tpl\t$spec\t%0.0f\t%0.0f\t%0.0f\t%0.0f\t%0.2f\t$filekey\n", $tp, $tn, $fp, $fn, $mcc;
	
	if ($tpline[2] ne $fpline[2] ){
	    die "$tpline[2] is not equal to $fpline[2]\n in $tpfile & $fpfile\n";
	}
    }
#    $name_idx++;
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
    
    
    @mcc_sorted = sort{$a <=> $b} @mcc_sorted;
    
    #sample size (583) is odd (mid = 291):
    my $mid = sprintf("%.0f", 0.5 * ($#{ $seq5_mcc_vals{$family} } + 1) );
    $seq5_median_mcc{$family} = $mcc_sorted[$mid];
    
    #quartiles (1) are even:
    my $uq = sprintf("%.0f", 0.75 * ($#{ $seq5_mcc_vals{$family} } + 1) );
    $seq5_uq_mcc{$family} = ($mcc_sorted[$uq]+$mcc_sorted[$uq+1])/2;
    my $lq = sprintf("%.0f", 0.25 * ($#{ $seq5_mcc_vals{$family} } + 1) );
    $seq5_lq_mcc{$family} = ($mcc_sorted[$lq-1]+$mcc_sorted[$lq])/2;
    
    
    
    @sens_sorted = sort{$a <=> $b} @sens_sorted;
    $seq5_median_sens{$family} = $sens_sorted[$mid];
    $seq5_uq_sens{$family} = ($sens_sorted[$uq]+$sens_sorted[$uq+1])/2;
    $seq5_lq_sens{$family} = ($sens_sorted[$lq-1]+$sens_sorted[$lq])/2;

    @spec_sorted = sort{$a <=> $b} @spec_sorted;
    $seq5_median_spec{$family} = $spec_sorted[$mid];
    $seq5_uq_spec{$family} = ($spec_sorted[$uq]+$spec_sorted[$uq+1])/2;
    $seq5_lq_spec{$family} = ($spec_sorted[$lq-1]+$spec_sorted[$lq])/2;

#     printf STDOUT "mcc: %0.3f\t%0.3f\t%0.3f\t$family\n", $seq5_lq_mcc{$family}, $seq5_median_mcc{$family},$seq5_uq_mcc{$family};
#     printf STDOUT "sens:%0.3f\t%0.3f\t%0.3f\t$family\n", $seq5_lq_sens{$family}, $seq5_median_sens{$family},$seq5_uq_sens{$family};
#     printf STDOUT "spec:%0.3f\t%0.3f\t%0.3f\t$family\n", $seq5_lq_spec{$family}, $seq5_median_spec{$family},$seq5_uq_spec{$family};
}

#print "\n";

my %tmp_mcc_hash = ();
my %seq5_mcc_ranks = ();
my %seq5_mcc_ranks_save = ();
my $sample_size = $#{ $seq5_mcc_vals{"erp"} };
my $old_family = "";

# foreach my $family ( keys %seq5_mcc_vals ) {
#     system("rm $seq5_toprint{$family}" . ".mcc");
#     system("rm $seq5_toprint{$family}" . ".ranks");
# }

for (my $i=0; $i<=$sample_size; $i++){
  
    foreach my $family ( keys %seq5_mcc_vals ) {
	$tmp_mcc_hash{$family} = $seq5_mcc_vals{$family}[$i];
    }
    my @sorted_by_mcc = sort by_mcc @seq5_shortnames;
    
    my $cnt = 1;
    my $old_cnt = 1;
    my $test = 1;
#    print "ROUND $i \n";
    foreach my $family (@sorted_by_mcc) {
	
       	if ($seq5_mcc_vals{$family}[$i] == $seq5_mcc_vals{$old_family}[$i]){
	    #print "$family\t$old_family\n";
	    $seq5_mcc_ranks{$family} += $old_cnt;
	    push @{ $seq5_mcc_ranks_save{$family} }, $old_cnt;
	    $test = $old_cnt;
	}
	else {
	    if ($seq5_mcc_ranks{$family}) {
		$seq5_mcc_ranks{$family} += $cnt;
		$old_cnt = $cnt;
		push @{ $seq5_mcc_ranks_save{$family} }, $cnt;
		$test = $cnt;
	    }
	    else {
		$seq5_mcc_ranks{$family} = $cnt;
		$old_cnt = $cnt;
		push @{ $seq5_mcc_ranks_save{$family} }, $cnt;
		$test = $cnt;
	    }

	}
	
	$old_family = $family;
	$cnt++;
	
	#CHECK IN R:
	#T<-0.953191; F<-0.0012766

	#cnt <- 602; #rRNA
	#cnt <- 1114; #tRNA
	#cnt <- 235; #U5
	#tp <- T*cnt; fn <- cnt-tp; fp <- F*cnt*10; tn <- cnt*10 - fp; mcc <- (tp*tn - fp*fn)/sqrt( (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) ); mcc
	#print "$family: MCC = [$seq5_mcc_vals{$family}[$i]] rank = $test sens = $seq5_sens_vals{$family}[$i] spec = $seq5_spec_vals{$family}[$i]\n";

my $outfile = $seq5_toprint{$family} . ".mcc";
open my $outf, ">>$outfile" or die "$outfile: $!\n";
my $outfiler = $seq5_toprint{$family} . ".ranks";
open my $outrank, ">>$outfiler" or die "$outfiler: $!\n";
printf $outf "$seq5_mcc_vals{$family}[$i]\n";
printf $outrank "$test\n";
	
    }
    
#     foreach my $family (@seq5_shortnames) {
# 	print "$family\t$seq5_mcc_vals{$family}[$i]\t$seq5_mcc_ranks{$family} \n";
#     }
    
}

#printf STDOUT "Name\tSensitivity\tSpecificity\tMCC\tMean MCC Rank\n\n";

my $phmm = 1;
my $smeth = 1;
my $ssmeth = 1;


open my $out, ">table5.tex" or die "table5.tex: $!\n";

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

    if ($family =~ /NCBI/ && $ssmeth){
	printf $out "\\hline
\\multicolumn\{6\}\{\|\|c\|\|\}\{\\bf Sequence based methods\}\\\\
\\hline
";
	$ssmeth = 0;
    }
    elsif ($family =~ /HMMer/ && $phmm){
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
my $med = ${ $seq5_mcc_ranks_save{$family} }[$mid];
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
printf $out "$seq5_longnames{$family} \& (%0.2f,\\textbf\{%0.2f,\}%0.2f) \& (%0.3f,\\textbf\{%0.3f\},%0.3f) \& (%0.2f,\\textbf\{%0.2f\},%0.2f) \& %0.1f \& %0.2f \\\\ \n", 
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
\\caption[]{ Tabulated results of searches with {\\bf five} sequences. From
left to right column one contains program names and settings, column two sensitivity,
column three specificity, column four Matthew's correlation coefficient (MCC)
and column five contains both a median and mean ranking determined
by the MCC (see Box 2 for definitions). Sensitivity, specificity and MCC are 
summarised as 3-tuples displaying the lower quartiles, medians and upper quartiles. 
Median MCC rankings below ten (one being the best) are indicated with bold font. 
See Methods section and Supplementary Tables~\\ref{table:homol:alg}\\&\\ref{table:homol:alg2} for algorithm settings.
}\\label{table:5seqs}
}
\\end{table}
";
close($out);
system "cat table5.tex";

# foreach my $family ( keys %seq5_toprint ) {
# my $outfile = $seq5_toprint{$family} . ".mcc";
# open my $outf, ">$outfile" or die "$outfile: $!\n";
# my $outfiler = $seq5_toprint{$family} . ".ranks";
# open my $outrank, ">$outfiler" or die "$outfiler: $!\n";
# #printf $outf "$family\n";
# foreach my $imcc (@{ $seq5_mcc_vals{$family} }) {
# printf $outf "$imcc\n";
# }
# foreach my $irank (@{ $seq5_mcc_ranks_save{$family} }) {
# printf $outrank "$irank\n";
# }
# }
##############################
open $out, ">tableWU5.tex" or die "tableWU5.tex: $!\n";

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

$ssmeth = 1;
foreach my $family ( @seq5_shortnamesWU ) {
    

    if ($family =~ /NCBI/ && $ssmeth){
	printf $out "\\hline
\\multicolumn\{6\}\{\|\|c\|\|\}\{\\bf Five input sequences \}\\\\
\\hline
";
	$ssmeth = 0;
    }

    @{ $seq5_mcc_ranks_save{$family} } = sort{$a <=> $b} @{ $seq5_mcc_ranks_save{$family} };

my $mid = sprintf("%.0f", 0.5 * $sample_size );
my $med = ${ $seq5_mcc_ranks_save{$family} }[$mid];

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

close($out);

system "cat tableWU5.tex tableWU20.tex > tableWU.tex";
#############################

open $out, ">tableRNA5.tex" or die "tableRNA5.tex: $!\n";

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

$ssmeth = 1;
foreach my $family ( @seq5_shortnamesRNA ) {
    

    if ($family =~ /WUblast/ && $ssmeth){
	printf $out "\\hline
\\multicolumn\{6\}\{\|\|c\|\|\}\{\\bf Five input sequences \}\\\\
\\hline
";
	$ssmeth = 0;
    }

    
my $mid = sprintf("%.0f", 0.5 * $sample_size );
my $med = ${ $seq5_mcc_ranks_save{$family} }[$mid];

    if ($seq5_mcc_ranks{$family}/$sample_size <= 10){
printf $out "$seq5_longnames{$family} \& (%0.2f,\\textbf\{%0.2f,\}%0.2f) \& (%0.3f,\\textbf\{%0.3f\},%0.3f) \& (%0.2f,\\textbf\{%0.2f\},%0.2f) \& \\textbf\{%0.1f\} \& %0.2f\\\\ \n", 
$seq5_lq_sens{$family}, $seq5_median_sens{$family},$seq5_uq_sens{$family}, 
$seq5_lq_spec{$family}, $seq5_median_spec{$family},$seq5_uq_spec{$family}, 
$seq5_lq_mcc{$family}, $seq5_median_mcc{$family},$seq5_uq_mcc{$family}, 
#$seq5_mcc_ranks{$family}/$sample_size;
$med,
$seq5_mcc_ranks{$family};
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
