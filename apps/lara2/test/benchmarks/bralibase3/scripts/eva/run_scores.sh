#!/bin/bash
# $1 method
# $2 options
# $3 identifier

for RNA in tRNA rRNA U5; do
    for id in id4060 id5070 id6080 id7085 id8095; do
	scores.py $1 ${RNA}/${id}/${RNA}00.fa ${RNA}/${RNA}_pn.db $2
    done
done > $1_$3_pn.scores

cat $1_$3_pn.scores | grep _ >  $1_$3_pn.allsh
cat $1_$3_pn.scores | grep -v _ >  $1_$3_pn.all