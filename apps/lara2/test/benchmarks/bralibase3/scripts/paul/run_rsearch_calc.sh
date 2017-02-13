#!/bin/sh
echo -ne "$1\t$2\t";
cat $1.*.rsearch.stk.$2 | grep ^">seq" | egrep -v '(>seq.>)|(>seq..>)|(>seq>)|(>seq...>)' | sort | uniq | wc -l 
#Had to add some nasty fixes to work around strange intermittent 
#sequence output errors. 
#See U5/id4060/U500.fa.3.rsearch.stk.U5 and U5/id4060/U5.fa.1.rsearch.stk.U5