#!/bin/sh
# getgff.sh searchmethod cutoff opt identifier
# the identifier is used for the output

#method = $1
#cutoff = $2
#options = $3
#identifier = $4

for RNA in tRNA rRNA U5; do
  ls -1 ${RNA}/ |grep id | while read id; do 
    echo $RNA;
    echo $id;
    ls -1 ${RNA}/${id}/${RNA}*.fa | while read x; do
      getgff.py $1 $x Hsapiens.fasta $2 $3
    done > ${RNA}/${id}/${RNA}_$4.$1.gff
    ls -1 ${RNA}/${id}/${RNA}*.fa | while read x; do
	getgff.py $1 $x Hsapiens_reverse.fasta $2 $3
    done > ${RNA}/${id}/${RNA}_$4.$1_rev.gff
  done
done

