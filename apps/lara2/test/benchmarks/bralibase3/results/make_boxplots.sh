#!/bin/sh
R CMD BATCH --nosave make_boxplots.R

#Fix eps borders:
cat boxplot5_nocol.eps | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > boxplot5_2.eps
mv boxplot5_2.eps boxplot5_nocol.eps
cat boxplot20_nocol.eps | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > boxplot20_2.eps
mv boxplot20_2.eps boxplot20_nocol.eps

cat boxplot5_ranks_nocol.eps | perl -ne 's/0 0 504 503/0 \-20 544 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 543.83 503.47 cl/g; s/499.51 258.93 \(MCC\) .5 0 90 t/515.51 258.93 \(MCC\) .5 0 90 t
0 0 1 rgb
515.51 178.93 \(Sensitivity\) .5 0 90 t 
0 1 0 rgb
515.51 338.93 \(Specificity\) .5 0 90 t/g; print $_;' > boxplot5_2.eps
mv boxplot5_2.eps boxplot5_ranks_nocol.eps

cat boxplot20_ranks_nocol.eps | perl -ne 's/0 0 504 503/0 \-20 544 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 543.83 503.47 cl/g; s/499.51 258.93 \(MCC\) .5 0 90 t/515.51 258.93 \(MCC\) .5 0 90 t
0 0 1 rgb
515.51 178.93 \(Sensitivity\) .5 0 90 t 
0 1 0 rgb
515.51 338.93 \(Specificity\) .5 0 90 t/g; print $_;' > boxplot5_2.eps
mv boxplot5_2.eps boxplot20_ranks_nocol.eps

cat WU-scoring_ranks_nocol5.eps  | perl -ne 's/0 0 504 503/0 \-20 544 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 543.83 503.47 cl/g; s/499.51 258.93 \(MCC\) .5 0 90 t/515.51 258.93 \(MCC\) .5 0 90 t
0 0 1 rgb
515.51 178.93 \(Sensitivity\) .5 0 90 t 
0 1 0 rgb
515.51 338.93 \(Specificity\) .5 0 90 t/g; print $_;' > WU-scoring_2.eps
mv WU-scoring_2.eps WU-scoring_ranks_nocol5.eps

cat WU-scoring_ranks_nocol20.eps  | perl -ne 's/0 0 504 503/0 \-20 544 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 543.83 503.47 cl/g; s/499.51 258.93 \(MCC\) .5 0 90 t/515.51 258.93 \(MCC\) .5 0 90 t
0 0 1 rgb
515.51 178.93 \(Sensitivity\) .5 0 90 t 
0 1 0 rgb
515.51 338.93 \(Specificity\) .5 0 90 t/g; print $_;' > WU-scoring_2.eps
mv WU-scoring_2.eps WU-scoring_ranks_nocol20.eps

cat RNA-scoring_ranks_nocol.eps  | perl -ne 's/0 0 504 503/0 \-65 544 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-65.00 543.83 503.47 cl/g; s/499.51 258.93 \(MCC\) .5 0 90 t/515.51 258.93 \(MCC\) .5 0 90 t
0 0 1 rgb
515.51 178.93 \(Sensitivity\) .5 0 90 t 
0 1 0 rgb
515.51 338.93 \(Specificity\) .5 0 90 t/g; print $_;' > RNA-scoring_2.eps
mv RNA-scoring_2.eps RNA-scoring_ranks_nocol.eps

cat ASR_results_nocol.eps | perl -ne 's/0 0 504 503/0 \-65 544 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-65.00 543.83 503.47 cl/g; s/499.51 258.93 \(MCC\) .5 0 90 t/515.51 258.93 \(MCC\) .5 0 90 t
0 0 1 rgb
515.51 178.93 \(Sensitivity\) .5 0 90 t 
0 1 0 rgb
515.51 338.93 \(Specificity\) .5 0 90 t/g; print $_;' > boxplot5_2.eps
mv boxplot5_2.eps ASR_results_nocol.eps

###################################

cat boxplot5_sens_nocol.eps | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > boxplot5_sens_2.eps
mv boxplot5_sens_2.eps boxplot5_sens_nocol.eps

cat boxplot20_sens_nocol.eps | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > boxplot20_sens_2.eps
mv boxplot20_sens_2.eps boxplot20_sens_nocol.eps

cat boxplot5_spec_nocol.eps | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > boxplot5_spec_2.eps
mv boxplot5_spec_2.eps boxplot5_spec_nocol.eps
cat boxplot20_spec_nocol.eps | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > boxplot20_spec_2.eps
mv boxplot20_spec_2.eps boxplot20_spec_nocol.eps

########
cat WU-scoring_ranks_nocol.eps  | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > WU-scoring_2.eps
mv WU-scoring_2.eps WU-scoring_ranks_nocol.eps

cat WU-scoring_ranks_nocol20.eps  | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > WU-scoring_2.eps
mv WU-scoring_2.eps WU-scoring_ranks_nocol20.eps

##5
cat WU-scoring_nocol.eps  | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > WU-scoring_2.eps
mv WU-scoring_2.eps WU-scoring_nocol.eps

cat WU-scoring_sens_nocol.eps  | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > WU-scoring_2.eps
mv WU-scoring_2.eps WU-scoring_sens_nocol.eps

cat WU-scoring_spec_nocol.eps  | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > WU-scoring_2.eps
mv WU-scoring_2.eps WU-scoring_spec_nocol.eps

##20
cat WU-scoring20_nocol.eps  | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > WU-scoring_2.eps
mv WU-scoring_2.eps WU-scoring20_nocol.eps

cat WU-scoring20_sens_nocol.eps  | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > WU-scoring_2.eps
mv WU-scoring_2.eps WU-scoring20_sens_nocol.eps

cat WU-scoring20_spec_nocol.eps  | perl -ne 's/0 0 504 503/0 \-20 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-20.00 503.83 503.47 cl/g; print $_;' > WU-scoring_2.eps
mv WU-scoring_2.eps WU-scoring20_spec_nocol.eps


#######
cat RNA-scoring_nocol.eps  | perl -ne 's/0 0 504 503/0 \-65 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-65.00 503.83 503.47 cl/g; print $_;' > RNA-scoring_2.eps
mv RNA-scoring_2.eps RNA-scoring_nocol.eps

cat RNA-scoring_sens_nocol.eps  | perl -ne 's/0 0 504 503/0 \-65 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-65.00 503.83 503.47 cl/g; print $_;' > RNA-scoring_2.eps
mv RNA-scoring_2.eps RNA-scoring_sens_nocol.eps

cat RNA-scoring_spec_nocol.eps  | perl -ne 's/0 0 504 503/0 \-65 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-65.00 503.83 503.47 cl/g; print $_;' > RNA-scoring_2.eps
mv RNA-scoring_2.eps RNA-scoring_spec_nocol.eps

cat ASR_results_nocol.eps | perl -ne 's/0 0 504 503/0 \-30 504 503/g; s/0.00 0.00 503.83 503.47 cl/0.00 \-30.00 503.83 503.47 cl/g; print $_;' > boxplot5_2.eps
mv boxplot5_2.eps ASR_results_nocol.eps


