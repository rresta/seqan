# laragu
The new version of Lara tool for the structural-alignment of RNA sequences 

To test the first local alignment
-i ./demo/input.fa -ir ./demo/input.fa -t 2 -v 3 -a

To run with the acquisition of ribosum matrices
-i ./demo/input.fa -ir ./demo/input.fa -t 2 -v 1 -lsm ./ribosum_matrices/RIBOSUM45-30_N.txt

To run with a single fasta file run the following command
-i ./demo/input.fa -t 2 -v 2

To run with two fasta files run the following command
-i ./demo/input.fa -ir ./demo/input.fa -t 2 -v 2

To test the dbn files run the following command
-i ./demo/input.fa -ir ./demo/rnafold.dbn -t 2 -v 2
