#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "precomp_definitions.h"
#include "errorexit.h"
#include "matrices.h"
#include "structures.h"
#include "neighborjoining.h"

DataSet		*Data;
TreeNode		*Nodes, *root;
TreeNode		**downPass, **upPass;
TreeLabels		*Names;

double CalculateKimura80Distances(double p, double q)
{
	return (log(1.0/(1.0-2.0*(p-q))) + log(1.0/(1.0-(2.0*q))));
}

void CalculateNeighborJoiningTree(int nTaxa, int nSites)
{
	int i, j, n;
	
	// Loop over pairwise comparisons and get distances
	for(i = 0; i < nTaxa-1; i++)
		{
		for(j = 0; j < nTaxa; j++)
			{
			// First get number of transitions and transversions for comparison
			
		
			}
		}

	// Make tree
	
	
}

void CalculateProportionTiTv(int t1, int t2, int nSites, double *p, double *q)
{


}

void GetNJTree(char *dataName, char *treeName)
{
	int *structure, nTaxa, nSites, nNodes;
	
	// Read data and consensus secondary structure
	ReadDataFile(dataName, &nTaxa, &nSites, 0, &structure, Names, Data);
	nNodes = (2 * nTaxa) - 1;

	// Allocate memory for structures and pointer arrays
	Nodes = AllocateMemoryForNodeStruct(nNodes);
	downPass = AllocateMemoryForPointerArray(nNodes);
	upPass = AllocateMemoryForPointerArray(nNodes);

	// Get NJ Tree
	CalculateNeighborJoiningTree(int nTaxa, int nSites)


}


