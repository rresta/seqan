# laragu
The new version of Lara tool for the structural-alignment of RNA sequences 

To test the creation of the lamb string that will store the aligned interections that will be updated at each iteration
-i ./apps/lara2/demo/input.fa -ir ./apps/lara2/demo/input.fa -t 2 -v 3 -a -lsm ./apps/lara2/ribosum_matrices/RIBOSUM45-30_N.txt -lgo -2 -lge -1

To test the first local alignment
-i ./apps/lara2/demo/input.fa -ir ./apps/lara2/demo/input.fa -t 2 -v 3 -a

To run with the acquisition of ribosum matrices
-i ./apps/lara2/demo/input.fa -ir ./apps/lara2/demo/input.fa -t 2 -v 1 -lsm ./apps/lara2/ribosum_matrices/RIBOSUM45-30_N.txt

To run with a single fasta file run the following command
-i ./apps/lara2/demo/input.fa -t 2 -v 2

To run with two fasta files run the following command
-i ./apps/lara2/demo/input.fa -ir ./apps/lara2/demo/input.fa -t 2 -v 2

To test the dbn files run the following command
-i ./apps/lara2/demo/input.fa -ir ./apps/lara2/demo/rnafold.dbn -t 2 -v 2
