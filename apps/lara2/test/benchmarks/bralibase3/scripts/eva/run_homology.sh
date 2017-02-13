#!/bin/sh
# run_homology.sh searchmethod cutoff opt identifier
# the identifier is used for the output

#method = $1
#cutoff = $2
#opt = $3
#ident = $4

for RNA in tRNA rRNA U5; do
  ls -1 ${RNA}/ |grep id | while read id; do 
    echo $RNA;
    echo $id;
    ls -1 ${RNA}/${id}/${RNA}*.fa | grep -v erpin | grep -v rnanc | while read x; do
      printf '%s\t' $x;
      homology.py $1 $x ${RNA}/${RNA} $2 $3
    done > ${RNA}/${id}/${RNA}_$4.$1.sens
    ls -1 ${RNA}/${id}/${RNA}*.fa |grep -v erpin | grep -v rnanc | while read x; do
	printf '%s\t' $x;
	homology.py $1 $x ${RNA}/sh${RNA} $2 $3
    done > ${RNA}/${id}/sh${RNA}_$4.$1.sens
  done
done


ls -1 */id*/*_$4.$1.sens |grep -v sh | xargs -n 40 cat | tr '/' '\t'| tr ' ' '\t' |sed -e 's/\.fa//g' | sort -f > $1_$4.tps
ls -1 */id*/*_$4.$1.sens |grep  sh | xargs -n 40 cat | tr '/' '\t'| tr ' ' '\t' |sed -e 's/\.fa//g' | sort -f > $1_$4.fps
