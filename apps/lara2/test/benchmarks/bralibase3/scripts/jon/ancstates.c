#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "precomp_definitions.h"
#include "globals.h"
#include "errorexit.h"
#include "complex.h"
#include "gamma.h"
#include "linalg.h"
#include "matrices.h"
#include "structures.h"
#include "ancstates.h"

// Global Variables
int			 show_debug_info=0, isStructure=NO;
char			 numbers[11] = "0123456789";
double		 singlet_q_matrix[NSTATES][NSTATES], doublet_q_matrix[NDBLSTATES][NDBLSTATES];

void AllocateMemory(int x, int y, int num_nodes, TreeNode *nodes)
{
	int i;
	TreeNode *p;
	
	for(i = 0; i < num_nodes; i++)
		{
		p = nodes + i;
		p->frac_like = AllocateMemoryForFractionalLikelihoods(x);
		p->posterior = AllocateMemoryForPosteriorProbs(x);
		}
}

extern void	 CalcAncestralStates(int is_harmonic, double lambda, int debug_status, int condition_pairs, int output_original, int nsamples, char *dataName, char *treeName, char *outName, char *paramName)
{
	int			 nTaxa, nSites, numStems, numLoops, nNodes, length, *structure, *dmatrix, *anc_seqs;
	double	         *parameters, *log_P;
	long int		seed;
	TreeNode		*Nodes, *root;
	TreeNode		**downPass, **upPass;
	TreeLabels		*Names;
	TreeStrings	*Trees;
	SecondStruct	*SS;
	
	show_debug_info = debug_status;
	
	// Set seed value from clock
	seed = GetSeedFromClock();

	// Get number of taxa and sites
	length = GetNumberOfTaxaSitesAndAllocate(dataName, &nTaxa, &nSites);
	
	// Allocate some things based on nTaxa and nSites
	structure = (int *)malloc(sizeof(int) * (nSites));
	if(!structure) { printf("ReadDataFile Error: Problem allocating memory for structure.\n"); exit(1); }
	Names = (TreeLabels *)malloc(sizeof(TreeLabels) * (nTaxa));
	if(!Names) { printf("ReadDataFile Error: Problem allocating memory for Names.\n"); exit(1); }
	SS = (SecondStruct *)malloc(sizeof(SecondStruct) * (nSites));
	if(!SS) { printf("ReadDataFile Error: Problem allocating memory for Data.\n"); exit(1); }
	dmatrix = (int *)malloc(sizeof(int) * (nSites) * (nTaxa));
	if(!dmatrix) { printf("ReadDataFile Error: Problem allocating memory for dmatrix.\n"); exit(1); }

	// Get data file contents -- data and structure
	ReadDataFile(dataName, nTaxa, nSites, 0, dmatrix, structure, Names, SS, length);
	nNodes = (2 * nTaxa) -1;
	
	// Intialize stems and loops in Data
	InitializeStemsAndLoops(nTaxa, nSites, &numStems, &numLoops, structure, condition_pairs, dmatrix, SS);

	if(show_debug_info)
		printf("Number doublet sites = %d\nNumber of singleton sites = %d\n", numStems, numLoops);

	// Allocate memory for structures and pointer arrays
	Trees = AllocateMemoryForTreeStringStructure(1);
	Nodes = AllocateMemoryForNodeStruct(nNodes);
	downPass = AllocateMemoryForPointerArray(nNodes);
	upPass = AllocateMemoryForPointerArray(nNodes);
	parameters = (double *)malloc(sizeof(double) * 32);
	if(!parameters)
		PrintErrorAndExitSafely("CalcAncestralStates Error: Failure to allocate parameters.");

	anc_seqs = (int *)malloc(sizeof(int) * (nTaxa-2) * nsamples * nSites * 2);
	if(!anc_seqs)
		PrintErrorAndExitSafely("CalcAncestralStates Error: Failure to allocate anc_seqs.");
	log_P = (double *)malloc(sizeof(double) * (nTaxa-2) * nsamples * 2);
	if(!log_P)
		PrintErrorAndExitSafely("CalcAncestralStates Error: Failure to allocate log_P.");

	// Get trees and parameters from file (parameters should include both doublet model and singlet model)
	ReadTreeFile(treeName, Trees);
	ReadParameterFile(paramName, parameters);	
	
	// Allocate memory for fractional likelihoods and lists at each node
	AllocateMemory(20 * (numStems + numLoops), 0, nNodes, Nodes);

	// Set structure to defaults, i.e., nulls and undefined
	SetTreeStructureToDefaults(Nodes, upPass, downPass, root, nNodes, (numStems + numLoops));
	
	// Read tree into structure
	SetUpTree(Nodes, upPass, downPass, root, Names, Trees[0].t_string, nTaxa);	
	
	//GetAverageLengthToTips(upPass, downPass, counter);
	GetHarmonicMeanOfBranchLengths(is_harmonic, lambda, downPass, (2*nTaxa)-1, nTaxa);
	
	if(show_debug_info)
		ShowTreeStructure(nNodes, Nodes, downPass, upPass);
		
	// Assign Q-matrix for HKY85 model (Singlet model)
	SetQMatrix(parameters);
	
	// Assign Q-matrix for doublet model
	SetDoubletQMatrix(parameters);
	
	// Assign transition probabilities to each node for singelton sites
	SetSingeltonTransitionProbabilities(upPass, nNodes);

	// Assign transition probabilities to each node for doublet sites
	SetDoubletTransitionProbabilities(upPass, nNodes);

	// Get fractional likelihoods across rate categories
	GetFractionalLikelihoods(downPass, SS, structure, dmatrix, nSites, (numStems + numLoops), nNodes);
	
	// Calculate Empirical Bayes ancestral state reconstructions for each node, and predictive lineages
	GetAncestralStatesFromPosterior(upPass, structure, parameters, nSites, (numStems + numLoops), nNodes);

	// Sample ancestral nodes from marginal posterior distribution to generate ancestral sequences
	SampleAncestralSequences(downPass, SS, structure, log_P, &seed, anc_seqs, nsamples, nSites, nNodes, condition_pairs);

	//SampleHiddenAncestralSequences(downPass, root, SS, log_P, parameters, structure, &seed, anc_seqs, nsamples, nSites, nTaxa, nNodes);
	SimulateHiddenAncestralSequences(downPass, root, SS, log_P, parameters, structure, &seed, anc_seqs, nsamples, nSites, nTaxa, nNodes, lambda, condition_pairs);

	// Stochastically re-introduce gaps based on the frequency observed in the original data matrix columns
	ReIntroduceGapsInAncestralMatrix((nTaxa-2) * nsamples * 2, nSites, &seed, SS, structure, anc_seqs);

	if(show_debug_info)
		ShowAncestralSequences(log_P, structure, anc_seqs, (nTaxa-2) * nsamples * 2, nSites);

	//WriteResultsToFile_Old(outName, log_P, structure, anc_seqs, output_original, Names, dmatrix, nSites, nTaxa, nsamples);
	WriteResultsToFile_New(Nodes, downPass, outName, log_P, structure, anc_seqs, output_original, Names, dmatrix, nSites, nTaxa, nsamples, lambda, is_harmonic);

	//printf("Preparing to quit and exit the program...\n");
	free(structure);
	free(dmatrix);
	free(parameters);
	free(anc_seqs);
	free(log_P);
	free(upPass);
	free(downPass);
	FreeMemory(nNodes, Nodes);
	free(Nodes);
	free(Names);
	free(SS);
	//printf("\tMemory deallocated successfully.\n");
	//printf("Program exited successfully...\n");
}

void CalcCijk(double *c_ijk, int n, double **u, double **v)
{
	int 	i, j, k;
	double 	*pc;

		/* precalculate values for faster matrix mult in GTRChangeMatrix and GTRDerivatives */
		pc = c_ijk;
		for (i=0; i<n; i++)
			for (j=0; j<n; j++)
				for (k=0; k<n; k++)
				 	*pc++ = u[i][k] * v[k][j];	/* (note: pc = &c[i][j][k]) */
}

void CalcPij(double *c_ijk, int n, double *eigenValues, double v, double r, double **p)
{
	int		i, j, k;
	
	for (i=0; i<n; i++)
		{
		for (j=0; j<n; j++)
			{
			p[i][j] = 0.0;
			for (k=0; k<n; k++)
				{
				p[i][j] += c_ijk[i*n*n+j*n+k]*exp(eigenValues[k]*v*r);
				}
			if (p[i][j] < 0.0)
				p[i][j] = 0.0;
			}
		}
}

int CompareChar(char c)
{
	int i;

	for (i = 0; i < 10; i++)
		{
		if (c == numbers[i])
			return YES;
		}
	return NO;
}

int ComplexChangeMatrix(double *eigenValues, double *eigvalsImag, complex **Ceigvecs, complex **CinverseEigvecs, int n, double **p, double t, double r)
{
	int			i, j, k, rc;
	complex		**Cwork;

	rc = NO_ERROR;
	Cwork = pscmatrix(n);
	for(i=0; i<n; i++)
		{
		for(j=0; j<n; j++)
			{
			Cwork[i][j].re = 0.0;
			Cwork[i][j].im = 0.0;
			for(k = 0; k < n; k++)
				{
				Cwork[i][j] = Cadd(Cwork[i][j], Cmul(Ceigvecs[i][k], Cmul(CinverseEigvecs[k][j], Cexp(RCmul(t*r, Complex(eigenValues[k], eigvalsImag[k]))))));
				}
			}
		}
	for(i=0; i<n; i++)
		{
		for(j=0; j<n; j++)
			{
			p[i][j] = Cwork[i][j].re;
			if(fabs(Cwork[i][j].im) > 0.000000000001)
				{
	    		rc = ERROR;
				}
			}
		}
	free_pscmatrix (Cwork);
	return (rc);
}

void FreeMemory(int num_nodes, TreeNode *nodes)
{
	int i;
	TreeNode *p;
	
	for(i = 0; i < num_nodes; i++)
		{
		p = nodes + i;
		free(p->frac_like);
		}
}

void GetAncestralStatesFromPosterior(TreeNode **upPass, int *structure, double *parameters, int nSites, int compNSites, int nNodes)
{
	int i, k, n, m, counter;
	double sum, sumL, sumR, sumA;
	TreeNode *p, *l, *r, *a;

	for(i = 0; i < nNodes; i++)
		{
		p = upPass[i];
		if(p->left != NULL && p->right != NULL && p->anc != NULL)
			{
			// Set left, right, and anc pointers
			l = p->left;
			r = p->right;
			if(p->anc->anc == NULL)	// At nodes left and/or right of arbitrary root node
				{
				if(p == p->anc->right)
					a = p->anc->left;
				else
					a = p->anc->right;
				}
			else
				a = p->anc;
			
			// Get ancestral states from posterior for internal nodes
			counter = 0;			
			for(k = 0; k < nSites; k++)
				{
				if(structure[k] == 0)
					{
					sum = 0.0;
					for(m = 0; m < NSTATES; m++)
						{
						sumL = sumR = sumA = 0.0;					
						for(n = 0; n < NSTATES; n++)
							{
							sumL += l->t_probs[m][n] *l->frac_like[(counter*20)+16+n];
							sumR += r->t_probs[m][n] * r->frac_like[(counter*20)+16+n];
							sumA += a->t_probs[m][n] * a->frac_like[(counter*20)+16+n];
							}
						sum += p->posterior[(counter*20)+16+m] = (sumL * sumR * sumA) * parameters[12 + m];
						}
					// Normalize posteriors
					for(n = 0; n < NSTATES; n++)
						p->posterior[(counter*20)+16+n] /= sum;
					counter++;
					}
				else if (structure[k] == 1)
					{
					sum = 0.0;
					for(m = 0; m < NDBLSTATES; m++)
						{
						sumL = sumR = sumA = 0.0;					
						for(n = 0; n < NDBLSTATES; n++)
							{
							sumL += l->t_probs_rna[m][n] * l->frac_like[(counter*20) + n];
							sumR += r->t_probs_rna[m][n] * r->frac_like[(counter*20) + n];
							sumA += a->t_probs_rna[m][n] * a->frac_like[(counter*20) + n];
							}
						sum += p->posterior[(counter*20)+m] = (sumL * sumR * sumA) * parameters[12 + m];
						}
					// Normalize posteriors
					for(n = 0; n < NDBLSTATES; n++)
						p->posterior[(counter*20)+n] /= sum;
					counter++;
					}
				}
			}
		}
		
	// If debug output is requested, show fractional likelihoods
	if(show_debug_info)
		{
		printf("\nPosterior probabilities --\n");
		for(i = 0; i < nNodes; i++)
			{
			p = upPass[i];
			if(p->left != NULL && p->right != NULL && p->anc != NULL)
				{
				printf("\nNode %d:\n", p->index);
				counter = 0;
				for(k = 0; k < nSites; k++)
					{
					if(structure[k] == 0 || structure[k] == 1)
						{
						printf("\tSite %d (%d)\n\t\t", counter, structure[k]);
						for(n = 0; n < 20; n++)
							{
							printf("%lf ", p->posterior[(counter * 20) + n]);
							if(n == 3 || n == 7 || n == 11 || n == 15)
								printf("\n\t\t");
							}
						printf("\n");
						counter++;
						}
					}
				}
			}
		}

}

void GetAverageLengthToTips(TreeNode **upPass, TreeNode **downPass, int nNodes)
{
	int i, j, n;
	double length;
	TreeNode *p, *q;

	for(i = 0; i < nNodes; i++)
		{
		p = upPass[i];
		if(p->left != NULL && p->right != NULL)
			{
			n=0;
			for(j = 0; j < nNodes-1; j++)
				{
				q = downPass[j];
				if(q->left == NULL && q->right == NULL)
					{
					length=0.0;
					if(GetSumLengthToTipsFromNode(q, p, &length))
						{
						//printf("p=%d q=%d (%lf)\n", p->index, q->index, length);
						p->avg_length_to_tips += length;
						n++;
						}
					}
				}
			p->avg_length_to_tips /= (double)n;
			}
		}
}

int GetSumLengthToTipsFromNode(TreeNode *p, TreeNode *target, double *length)
{
	double t_len=0.0;
	while (p != target)
		{
		t_len += p->current_length;
		if(p->anc == NULL)
			{
			return 0;
			}
		p = p->anc;
		}
	(*length) = t_len;
	return 1;
}

void GetDoubletMembers(int *pos1, int *pos2, int doublet)
{
	if(doublet == 0)
		{
		*(pos1) = 0;
		*(pos2) = 0;
		}
	else if(doublet == 1)
		{
		*(pos1) = 0;
		*(pos2) = 1;
		}
	else if(doublet == 2)
		{
		*(pos1) = 0;
		*(pos2) = 2;
		}
	else if(doublet == 3)
		{
		*(pos1) = 0;
		*(pos2) = 3;
		}
	else if(doublet == 4)
		{
		*(pos1) = 1;
		*(pos2) = 0;
		}
	else if(doublet == 5)
		{
		*(pos1) = 1;
		*(pos2) = 1;
		}
	else if(doublet == 6)
		{
		*(pos1) = 1;
		*(pos2) = 2;
		}
	else if(doublet == 7)
		{
		*(pos1) = 1;
		*(pos2) = 3;
		}
	else if(doublet == 8)
		{
		*(pos1) = 2;
		*(pos2) = 0;
		}
	else if(doublet == 9)
		{
		*(pos1) = 2;
		*(pos2) = 1;
		}
	else if(doublet == 10)
		{
		*(pos1) = 2;
		*(pos2) = 2;
		}
	else if(doublet == 11)
		{
		*(pos1) = 2;
		*(pos2) = 3;
		}
	else if(doublet == 12)
		{
		*(pos1) = 3;
		*(pos2) = 0;
		}
	else if(doublet == 13)
		{
		*(pos1) = 3;
		*(pos2) = 1;
		}
	else if(doublet == 14)
		{
		*(pos1) = 3;
		*(pos2) = 2;
		}
	else if(doublet == 15)
		{
		*(pos1) = 3;
		*(pos2) = 3;
		}
}

int GetEigens(int n, double **q, double *eigenValues, double *eigvalsImag, double **eigvecs, double **inverseEigvecs, complex **Ceigvecs, complex **CinverseEigvecs)
{
	int			i, j, rc, *iwork, isComplex;
	double		**mwork, *dwork;
	complex 	**Cwork, *Ccol;

	/* allocate memory */
	dwork           = (double *)malloc((size_t) (n * sizeof(double)));
	iwork           = (int *)malloc((size_t) (n * sizeof(int)));

	/* calculate eigenvalues and eigenvectors */
	isComplex = NO;
	rc = EigenRealGeneral (n, q, eigenValues, eigvalsImag, eigvecs, iwork, dwork);
	if (rc != NO_ERROR)
		{
		if (rc == RC_COMPLEX_EVAL)
			{
			isComplex = YES;
			}
		}

	if (isComplex == NO)
		{
		mwork = psdmatrix (n);
		copy_psdmatrix (eigvecs, mwork, n);
		InvertMatrix (mwork, n, dwork, iwork, inverseEigvecs);
		free_psdmatrix (mwork);
		}
	else
		{
		for(i=0; i<n; i++)
			{
			if (eigvalsImag[i] == 0.0)
				{ 
				for(j=0; j<n; j++)
					{
					Ceigvecs[j][i].re = eigvecs[j][i];
					Ceigvecs[j][i].im = 0.0;
					}
				}
			else if (eigvalsImag[i] > 0.0)
				{ 
				for (j=0; j<n; j++)
					{
					Ceigvecs[j][i].re = eigvecs[j][i];
					Ceigvecs[j][i].im = eigvecs[j][i + 1];
					}
				}
			else if (eigvalsImag[i] < 0.0)
				{ 
				for (j=0; j<n; j++)
					{
					Ceigvecs[j][i].re =  eigvecs[j][i-1];
					Ceigvecs[j][i].im = -eigvecs[j][i];
					}
				}
			}
		Ccol = (complex *)malloc((size_t) (n * sizeof(complex)));
		Cwork = pscmatrix (n);
		copy_pscmatrix (Ceigvecs, Cwork, n);
		ComplexInvertMatrix (Cwork, n, dwork, iwork, CinverseEigvecs, Ccol);
		free (Ccol);
		free_pscmatrix (Cwork);
		}

	free (dwork);
	free (iwork);
	return (isComplex);
}

void GetFractionalLikelihoods(TreeNode **downPass, SecondStruct *SS, int *structure, int *data, int nSites, int compNSites, int nNodes)
{
	int i, j, n, m, counter=0;
	double sumL, sumR, x, y;
	TreeNode *p;
	SecondStruct *str;

	for(i = 0; i < nNodes; i++)
		{
		p = downPass[i];
		counter = 0;
		for(j = 0; j < nSites; j++)
			{
			if(structure[j] == 0 || structure[j] == 1)
				{
				if(p->left == NULL && p->right == NULL && p->anc != NULL)
					{
					str = SS + j;
					// Get nucleotide state for current site and initialize tips
					if(structure[j] == 0)
						IntializeSingeltonFracLikelihoodsAtTips(p, counter*20, data[(p->d_matrix_pos*nSites)+j]);
					else
						IntializeDoubletFracLikelihoodsAtTips(p, counter*20, data[(p->d_matrix_pos*nSites)+j],  data[(p->d_matrix_pos*nSites)+str->partner->site]);
					}
				else
					{
					// Calculate conditional likelihoods depending on whether they are singeltons or doublets
					if(structure[j] == 0)
						{
						for(n = 0; n < NSTATES; n++)
							{
							sumL = sumR = 0.0;
							for(m = 0; m < NSTATES; m++)
								{
								x = p->left->t_probs[n][m] * p->left->frac_like[(counter*20)+16+m];
								sumL += x;
								y = p->right->t_probs[n][m] * p->right->frac_like[(counter*20)+16+m];
								sumR += y;
								}
							p->frac_like[(counter*20)+16+n] = sumL * sumR;
							}
						}
					else
						{
						for(n = 0; n < NDBLSTATES; n++)
							{
							sumL = sumR = 0.0;
							for(m = 0; m < NDBLSTATES; m++)
								{
								x = p->left->t_probs_rna[n][m] * p->left->frac_like[(counter*20)+m];
								sumL += x;
								y = p->right->t_probs_rna[n][m] * p->right->frac_like[(counter*20)+m];
								sumR += y;
								}
							p->frac_like[(counter*20)+n] = sumL * sumR;
							}
						}
					}
				counter++;
				}
			}
		}

	// If debug output is requested, show fractional likelihoods
	if(show_debug_info)
		{
		printf("\nFractional likelihoods --\n");
		for(i = 0; i < nNodes; i++)
			{
			p = downPass[i];
			printf("\nNode %d:\n", p->index);
			for(j = 0; j < compNSites; j++)
				{
				printf("\tSite %d\n\t\t", j);
				for(n = 0; n < 20; n++)
					{
					printf("%lf ", p->frac_like[(j * 20) + n]);
					if(n == 3 || n == 7 || n == 11 || n == 15)
						printf("\n\t\t");
					}
				printf("\n");
				}
			}
		}
}

void GetHarmonicMeanOfBranchLengths(int is_harmonic, double lambda, TreeNode **downPass, int nNodes, int nTaxa)
{
	int i, n;
	double root_interval=0.0, sum_branches=0.0, harm_mean=0.0;
	TreeNode *p, *rt;

	if(is_harmonic == 1)
		{
		rt = downPass[nNodes-1];

		// Get number of branches in an unrooted tree
		n = (2*nTaxa)-3;

		// Get the sum of the denominator, sum(1/x)
		for(i = 0; i < nNodes-1; i++)
			{
			p = downPass[i];
			if(p->anc == rt)			
				root_interval += p->current_length;
			else
				sum_branches += 1.0/p->current_length;
			}

		// Get inverse of root_interval branch and add to sum_branches
		sum_branches += 1.0/root_interval;
		
		// Get reciprocal of sum_branches
		sum_branches = 1.0/sum_branches;

		// Calculate harmonic mean
		harm_mean = sum_branches * (double)n;
		}
	else if(is_harmonic == 0)	
		{
		harm_mean = 1.0;
		}
	else
		{
		PrintErrorAndExitSafely("GetHarmonicMeanOfBranchLengths Error: Parameter 'is_harmonic' has an incorrect value.");
		}
	
	// Assign harmonic mean to nodes
	for(i = 0; i < nNodes-1; i++)
		{
		p = downPass[i];
		p->avg_length_to_tips = harm_mean;
		//printf("%d: %lf\n", p->index, p->avg_length_to_tips);
		}

	//exit(1);
}

int GetNumberOfTaxaSitesAndAllocate(char *fileName, int *nTaxa, int *nSites)
{
	int n, nt=0, ns=0, len;
	char string[1000], x;
	FILE *in;

	//printf("%s\n", fileName);

	// Open the file in which data and structure is to be read from
	in = fopen(fileName,"r");
	if(!in)
		PrintErrorAndExitSafely("GetNumberOfTaxaSitesAndAllocate Error: Unable to open file.");

	n = 0;
	while(!feof(in))
		{
		x=fgetc(in);
		n++;
		}
	rewind(in);

	while(fgets(string, n, in))
		{
		len =strlen(string);
		if((char)string[len-1] == '\r' || (char)string[len-1] == '\n')
			string[len-1] = '\0';
		if((char)string[0] == '>')
			nt++;
		else
			{
			if(ns == 0)	
				{
				ns = strlen(string);
				}
			else
				{
				if(ns != strlen(string))
					PrintErrorAndExitSafely("GetNumberOfTaxaSitesAndAllocate Error: Sequences appear to be different lengths.");
				}
			}
		}

	*(nTaxa) = nt-1;
	*(nSites) = ns;

	fclose(in);

	return n;
}

long int GetSeedFromClock(void)
{
	time_t		 curTime;
	long int		seed;

	curTime = time(NULL);
	seed  = (long int)curTime;
	if (seed < 0)
		seed = -seed;
		
	return seed;
}

void IntializeDoubletFracLikelihoodsAtTips(TreeNode *p, int offset, int pos1, int pos2)
{
	int i, which;
	double pos[NDBLSTATES];

	for(i = 0; i < NDBLSTATES; i++)
		pos[i] = 0.0;
	
	// If both positions are defined, i.e., not ambiguous (N, ?, -, etc.) then proceed, otherwise the code needs to deal with ambiguity
	if(pos1 < 4 && pos2 < 4)
		{
		which = WhichDoublet(pos1, pos2);
		pos[which] = 1.0;
		}
	else
		{
		if(pos1 < 4 && pos2 > 4)		// First position is known, second is uncertain
			{
			which = WhichDoublet(pos1, 0);
			for(i = 0; i < NSTATES; i++)
				pos[which+i] = 1.0;
			}
		else if (pos1 > 4 && pos2 < 4)	// First position is uncertain, second position is known
			{
			which = WhichDoublet(0, pos2);
			for(i = 0; i < NSTATES; i++)
				pos[which + (i * 4)] = 1.0;
			}
		else						// Both positions are uncertain
			{
			for(i = 0; i < NDBLSTATES; i++)
				pos[i] = 1.0;
			}
		}

	for(i = 0; i < NDBLSTATES; i++)
		p->frac_like[offset + i] = pos[i];
}

void IntializeSingeltonFracLikelihoodsAtTips(TreeNode *p, int offset, int state)
{
	double pos_a=0.0, pos_c=0.0, pos_g=0.0, pos_t=0.0;

	if(state == 0)
		pos_a = 1.0;
	else if(state == 1)
		pos_c = 1.0;
	else if(state == 2)
		pos_g = 1.0;
	else if(state == 3)
		pos_t = 1.0;
	else if(state == 20)
		pos_a = pos_g = 1.0;
	else if(state == 21)
		pos_c = pos_g = 1.0;
	else if(state == 10)
		pos_a = pos_c = 1.0;
	else if(state == 32)
		pos_g = pos_t = 1.0;
	else if(state == 30)
		pos_a = pos_t = 1.0;
	else if(state == 31)
		pos_c = pos_t = 1.0;
	else if(state == 310)
		pos_a = pos_c = pos_t = 1.0;
	else if(state == 321)
		pos_c = pos_g = pos_t = 1.0;
	else if(state == 320)
		pos_a = pos_g = pos_t = 1.0;
	else if(state == 210)
		pos_a = pos_c = pos_t = 1.0;
	else if(state == 3210)
		pos_a = pos_c = pos_g = pos_t = 1.0;
	else if(state == 3211)
		pos_a = pos_c = pos_g = pos_t = 1.0;
	else if(state == 3212)
		pos_a = pos_c = pos_g = pos_t = 1.0;
	
	p->frac_like[offset + 16 + 0] = pos_a;
	p->frac_like[offset + 16 + 1] = pos_c;
	p->frac_like[offset + 16 + 2] = pos_g;
	p->frac_like[offset + 16 + 3] = pos_t;
}

void InitializeStemsAndLoops(int nTaxa, int nSites, int *numStems, int *numLoops, int *structure, int condition_pair, int *dmatrix, SecondStruct *SS)
{
	int i, j, ns=0, nl=0, tracker, num_noncon=0, num_gaps=0;
	SecondStruct *s1, *s2;
	
	for(i = 0; i < nSites; i++)
		{
		s1 = SS + i;
		if(i == 0)
			{
			s1->prev = NULL;
			s1->next = SS + i + 1;
			}
		else if(i == (nSites-1))
			{
			s1->prev = SS + i - 1;
			s1->next = NULL;
			}
		else
			{
			s1->prev = SS + i - 1;
			s1->next = SS + i + 1;
			}
		s1->site = i;
		s1->paired = NO;
		s1->freq_canonical = 0.0;
		s1->freq_gaps = 0.0;
		s1->partner = NULL;
		}

	for(i = 0; i < nSites; i++)
		{
		s1 = SS + i;
		if(structure[i] == 0)
			{
			num_gaps = 0;
			for(j = 0; j < nTaxa; j++)
				{
				if(dmatrix[(j*nSites)+s1->site] == 3211)
					num_gaps++;
				}
			s1->freq_gaps = (double)num_gaps/(double)nTaxa;
			nl++;
			}
		else if(structure[i] == 1)
			{
			s1->paired = YES;
			tracker = 0;
			isStructure = YES;
			// Find partner, j, of site, i
			for(j = i+1; j < nSites; j++)
				{
				s2 = SS + j;
				tracker += structure[j];
				if(tracker == (-1))
					{
					s2->paired = YES;
					s1->partner = s2;
					s2->partner = s1;
					j =  nSites;
					}
				}
			if(condition_pair == YES)
				{
				// Dermine frequency of canonical dinculeotides
				num_noncon = 0;
				num_gaps = 0;
				for(j = 0; j < nTaxa; j++)
					{
					if(IsCanonicalDoublet(dmatrix[(j*nSites)+s1->site], dmatrix[(j*nSites)+s2->site]) == NO)
						num_noncon++;
					if(dmatrix[(j*nSites)+s1->site] == 3211)
						num_gaps++;
					}
				s1->freq_canonical = s2->freq_canonical = (1.0 - ((double)num_noncon/(double)nTaxa));
				s1->freq_gaps = (double)num_gaps/(double)nTaxa;
				}
			ns++;
			}
		else if(structure[i] == -1)
			{
			num_gaps = 0;
			for(j = 0; j < nTaxa; j++)
				{
				if(dmatrix[(j*nSites)+s1->site] == 3211)
					num_gaps++;
				}
			s1->freq_gaps = (double)num_gaps/(double)nTaxa;
			}
		}
	
	*(numStems) = ns;
	*(numLoops) = nl;
	
	if(show_debug_info)
		ShowSecondaryStructure(nSites, SS);
}

int IsCanonicalDoublet(int one, int two)
{
	if((one + two) == 3 || (one+two) == 5)
		return YES;

	return NO;
}

int IsDoubletChangeCanonicalChange(int which)
{
	// 3, 6, 9, 11, 12, 14
	if(which == 3)
		return YES;
	else if(which == 6)
		return YES;
	else if(which == 9)
		return YES;
	else if(which == 11)
		return YES;
	else if(which == 12)
		return YES;
	else if(which == 14)
		return YES;
		
	return NO;
}

int NucCharToInt(char nuc)
{
    if(nuc == 'A')
        return 0;
    else if(nuc == 'C')
        return 1;
    else if(nuc == 'G')
        return 2;
    else if(nuc == 'T')
        return 3;
    else if(nuc == 'U')
        return 3;
    else if(nuc == 'R')
        return 20;
    else if(nuc == 'Y')
        return 31;
    else if(nuc == 'S')
        return 21;
    else if(nuc == 'W')
        return 30;
    else if(nuc == 'K')
        return 32;
    else if(nuc == 'M')
        return 10;
    else if(nuc == 'B')
        return 321;
    else if(nuc == 'D')
        return 320;
    else if(nuc == 'H')
        return 310;
    else if(nuc == 'V')
        return 210;
    else if(nuc == 'N')
        return 3210;
    else if(nuc == '-')
        return 3211;
    else if(nuc == '?')
        return 3212;
        
    return 3210;
}

char NucIntToChar(int nuc)
{
	if(nuc == 0)
		return 'A';
	else if(nuc == 1)
		return 'C';
	else if(nuc == 2)
		return 'G';
	else if(nuc == 3)
		return 'U';
	else if(nuc == 20)
		return 'R';
	else if(nuc == 31)
		return 'Y';
	else if(nuc == 21)
		return 'S';
	else if(nuc == 30)
		return 'W';
	else if(nuc == 32)
		return 'K';
	else if(nuc == 10)
		return 'M';
	else if(nuc == 321)
		return 'B';
	else if(nuc == 320)
		return 'D';
	else if(nuc == 310)
		return 'H';
	else if(nuc == 210)
		return 'V';
	else if(nuc == 3210)
		return 'N';
	else if(nuc == 3211)
		return '-';
	else if(nuc == 3212)
		return '?';
        
    return 'N';
}

void ReadDataFile(char *fileName, int nTaxa, int nSites, int *nNodes, int *dmatrix, int *structure, TreeLabels *Names, SecondStruct *SS, int n)
{
	int i, j, nt, ns, getStructure=NO, getTaxa=NO, getSites=NO, isFirstSeq=0;
	char *string, x;
	FILE *in;
	TreeLabels *nm;

	// Open the file in which data and structure is to be read from
	in = fopen(fileName,"r");
	if(!in)
		PrintErrorAndExitSafely("ReadDataFile Error: Unable to open file.");

	ns = nSites;
	nt = 0;
	
	string = (char *)malloc(sizeof(char) * (n+10));
	if(!string) { printf("ReadDataFile Error: Problem allocating memory for string.\n"); exit(1); }
	
	if(show_debug_info)		
		printf("Data file contents (filename = %s [nTaxa=%d, nSites=%d):\n", fileName, nTaxa, nSites);
	while(fgets(string, n, in))
		{
		if(show_debug_info)		
			printf("%s", string);
		if(isFirstSeq == 0)
			{
			getStructure = YES;
			}
		else if(getStructure)
			{
			i = 0;
			while(i < ns)
				{
				if(string[i] == '(')
					structure[i] = 1;
				else if(string[i] == ')')
					structure[i] = -1;
				else if(string[i] == '-')
					structure[i] = 0;
				i++;
				}
			getStructure = NO;
			getTaxa = YES;
			}
		else if(getTaxa)
			{
			nm = Names + nt;
			nm->index = nt;
			if(nt == 0)
				nm->prev = NULL;
			else
				nm->prev = Names + nt - 1;
			if(nt == nTaxa-1)
				nm->next = NULL;
			else
				nm->next = Names + nt + 1;
			i = 0;
			while(string[i] == ' ')
				i++;
			if(string[i] == '>')
				i++;
			j=0;
			while(string[i] != '\n')
				{
				nm->nm_label[j++] = string[i++];
				}
			nm->nm_label[j] = '\0';
			nt++;
			getTaxa=NO;
			getSites=YES;
			}
		else if(getSites)
			{
			i=0;
			while(string[i] == ' ')
				i++;
			j = 0;
			while(j < ns)
				{
				x = toupper(string[i++]);
				dmatrix[((nt-1) * ns) + (j++)] = NucCharToInt(x);
				}
			getSites=NO;
			getTaxa=YES;
			}
		isFirstSeq++;
		}
		
	free(string);
	fclose(in);
}

void ReadParameterFile(char *fileName, double *parameters)
{
	int i, n, slen, start, finish, insideCommentBlock, GetValues;
	double y;
	char *string, *token, x, delimiters[6];
	FILE *in;
	
	delimiters[0] = '\n';
	delimiters[1] = '\t';
	delimiters[2] = '\r';
	delimiters[3] =  '=';
	delimiters[4] =  ' ';
	delimiters[5] =  '\0';
	
	insideCommentBlock = NO;
	GetValues = NO;

	for(i = 0; i < 32; i++)
		parameters[i] = 1.0;

	in = fopen(fileName,"r");
	if(!in)
		{
		printf("ReadTreeFile Error: Unable to open file.\n");
		exit(1);
		}

	n = 0;
	while(!feof(in))
		{
		x=fgetc(in);
		n++;
		}
	rewind(in);
	
	string = (char *)malloc(sizeof(char) * (n+10));
	if(!string) { printf("ReadTreeFile Error: Problem allocating memory for string.\n"); exit(1); }

	token = (char *)malloc(sizeof(char) * (n+10));

	n = 0;
	while(!feof(in))
		{
		x = tolower(fgetc(in));
		if(insideCommentBlock == NO)
			{
			if(x == LSQRBRACKET) // Entering comment block
				{
				insideCommentBlock = YES;
				}
			else
				{
				string[n++] = x;
				}
			}
		else if(x == RSQRBRACKET)
			{
			insideCommentBlock = NO;
			}
		}
	string[--n] = '\0';
	slen = strlen(string);

	// Parse string
	if(isStructure)
		{
		n = 0;
		start = finish = 0;
		StringToken(&start, &finish, slen, string, delimiters, token);
		while(token[0] != '\0')
			{
			//printf("token: %s\n", token);
			if(!strcmp(token, "r(a<->c){1}"))
				{
				GetValues = YES;
				}
			else if (GetValues)
				{
				//printf("token: %s\n", token);
				sscanf(token, "%lf", &y);
				parameters[n++] = y;
				for(i = 0; i < 5; i++)
					StringToken(&start, &finish, slen, string, delimiters, token);
				}
			StringToken(&start, &finish, slen, string, delimiters, token);
			if(!strcmp(token, "end;") || !strcmp(token, "end"))
				break;
			}
		}
	else
		{
		n = 21;
		start = finish = 0;
		StringToken(&start, &finish, slen, string, delimiters, token);
		while(token[0] != '\0')
			{
			//printf("token: %s\n", token);
			if(!strcmp(token, "r(a<->c)"))
				{
				GetValues = YES;
				}
			else if (GetValues)
				{
				//printf("token: %s\n", token);
				sscanf(token, "%lf", &y);
				parameters[n++] = y;
				for(i = 0; i < 5; i++)
					StringToken(&start, &finish, slen, string, delimiters, token);
				}
			StringToken(&start, &finish, slen, string, delimiters, token);
			if(!strcmp(token, "end;") || !strcmp(token, "end"))
				break;
			}
		}
		
	free(token);
	free(string);
	fclose(in);
}

void ReadTreeFile(char *fileName, TreeStrings *TStrings)
{
	int n, slen, start, finish, insideCommentBlock, GetTree, isFirstPass;
	char *string, *token, x, delimiters[6];
	TreeStrings *tr;
	FILE *in;
	
	delimiters[0] = '\n';
	delimiters[1] = '\t';
	delimiters[2] = '\r';
	delimiters[3] =  '=';
	delimiters[4] =  ' ';
	delimiters[5] =  '\0';

	insideCommentBlock = NO;
	GetTree = NO;
	isFirstPass = YES;
	
	in = fopen(fileName,"r");
	if(!in)
		{
		printf("ReadTreeFile Error: Unable to open file.\n");
		exit(1);
		}

	n = 0;
	while(!feof(in))
		{
		x=fgetc(in);
		n++;
		}
	rewind(in);
	
	string = (char *)malloc(sizeof(char) * (n+10));
	if(!string) { printf("ReadTreeFile Error: Problem allocating memory for string.\n"); exit(1); }

	token = (char *)malloc(sizeof(char) * (n+10));

	n = 0;
	while(!feof(in))
		{
		x = tolower(fgetc(in));
		if(insideCommentBlock == NO)
			{
			if(x == LSQRBRACKET) // Entering comment block
				{
				insideCommentBlock = YES;
				}
			else
				{
				string[n++] = x;
				}
			}
		else if(x == RSQRBRACKET)
			{
			insideCommentBlock = NO;
			}
		}
	string[--n] = '\0';
	slen = strlen(string);
	
	// Parse string
	n = 0;
	start = finish = 0;
	StringToken(&start, &finish, slen, string, delimiters, token);
	while(token[0] != '\0')
		{
		if(!strcmp(token, "tree"))
			{
			if(!isFirstPass)
				{
				StringToken(&start, &finish, slen, string, delimiters, token);				
				GetTree = YES;
				}
			isFirstPass = NO;
			}
		else if(GetTree)
			{
			n = strlen(token);
			if(token[n-1] == ';')
				{
				token[n-1] = '\0';
				start++;
				}
			tr = TStrings + 0;
			tr->t_string = (char*)malloc(sizeof(char) * n+1);
			if(!tr->t_string)
				PrintErrorAndExitSafely("ReadTreeFile Error: Failure allocating memory for tr->t_string.");
			strcpy(tr->t_string, token);
			GetTree = NO;
			}
		StringToken(&start, &finish, slen, string, delimiters, token);
		if(!strcmp(token, "end;") || !strcmp(token, "end"))
			break;
		}

	free(token);
	free(string);
	fclose(in);
}

void ReIntroduceGapsInAncestralMatrix(int nsamples, int nSites, long int *seed, SecondStruct *SS, int *structure, int *anc_seqs)
{
	int i, j;
	double ran;
	SecondStruct *str;
	
	// Loop over structure and sample gaps accordingly
	for(i = 0; i < nSites; i++)
		{
		str = SS + i;
		if(str->freq_gaps > 0.0)
			{
			for(j = 0; j < nsamples; j++)
				{
				ran = RandomNumber(seed);
				if(ran < str->freq_gaps)
					{
					anc_seqs[(j*nSites) + i] = 3211;
					}
				}
			}
		}
}

/*-------| RandomNumber |------------------------------------------------
|   This pseudorandom number generator is described in:
|   Park, S. K. and K. W. Miller.  1988.  Random number generators: good
|      ones are hard to find.  Communications of the ACM, 31(10):1192-1201.
*/
static double RandomNumber (long int *seed)
{
	long int	lo, hi, test;

	hi = (*seed) / 127773;
	lo = (*seed) % 127773;
	test = 16807 * lo - 2836 * hi;
	if (test > 0)
		*seed = test;
	else
		*seed = test + 2147483647;
	return (double)(*seed) / (double)2147483647;
}

void ReIndexNodesForTree(TreeNode *p, int *i)
{
	if (p != NULL)
		{
		p->index = (*i)++;
		ReIndexNodesForTree(p->left, i);
		ReIndexNodesForTree(p->right, i);
		}
}

void SampleAncestralSequences(TreeNode **downPass, SecondStruct *SS, int *structure, double *log_P, long int *seed, int *anc_seqs, int nsamples, int nSites, int nNodes, int condition_pairs)
{
	int i, j, n, m, k, counter=0, matrix_counter=0, which, pos1, pos2;
	double ran, sum, prob=0.0;
	TreeNode *p;
	SecondStruct *str;
	
	// Pass down tree visiting internal nodes and sampling ancestral reconstructions, ### ((2*nTaxa-2) * nsamples) * nSites
	for(i = 0; i < nNodes; i++)
		{
		p = downPass[i];
		if(p->left != NULL && p->right != NULL)
			{
			for(m = 0; m < nsamples; m++)
				{
				counter = 0;
				prob = 0.0;
				for(j = 0; j < nSites; j++)
					{
					str = SS + j;
					if(structure[j] == 0)
						{
						sum = 0.0;
						ran = RandomNumber(seed);
						for(n = 0; n < NSTATES; n++)
							{
							sum += p->posterior[(counter*20)+16+n];
							if(ran < sum) // State has been selected
								{
								which = n;
								prob += log(p->posterior[(counter*20)+16+n]);
								n = NSTATES;
								}
							}
						anc_seqs[(matrix_counter * nSites)+j] = which;
						counter++;
						}
					else if(structure[j] == 1)
						{
						// If user has requested that the reconstructions be conditional on the frequency of unpaired doublets in the real data at paired sites
						if(condition_pairs == YES && str->freq_canonical == 1.0)
							{
							sum = 0.0;
							for(n = 0; n < NSTATES; n++)
								{
								for(k = 0; k < NSTATES; k++)
									{
									if(IsCanonicalDoublet(n, k) == NO)
										p->posterior[(counter*20)+(n*NSTATES)+k] = 0.0;
									sum += p->posterior[(counter*20)+(n*NSTATES)+k];
									}
								}
							for(n = 0; n < NDBLSTATES; n++)
								p->posterior[(counter*20)+n] /= sum;
							}
						else if(condition_pairs == YES && str->freq_canonical < 1.0)
							{
							sum = 0.0;
							for(n = 0; n < NSTATES; n++)
								{
								for(k = 0; k < NSTATES; k++)
									{
									if(IsCanonicalDoublet(n, k) == YES)
										sum += p->posterior[(counter*20)+(n*NSTATES)+k] *= str->freq_canonical;
									else
										sum += p->posterior[(counter*20)+(n*NSTATES)+k] *= (1.0 - str->freq_canonical);
									}
								}
							for(n = 0; n < NDBLSTATES; n++)
								p->posterior[(counter*20)+n] /= sum;
							}
						sum = 0.0;
						ran = RandomNumber(seed);
						for(n = 0; n < NDBLSTATES; n++)
							{
							sum += p->posterior[(counter*20)+n];
							if(ran < sum) // State has been selected
								{
								which = n;
								prob = log(p->posterior[(counter*20)+n]);
								n = NDBLSTATES;
								}
							}
						// Need to break doublet into singeltons and put them in the correct place in the matrix
						GetDoubletMembers(&pos1, &pos2, which);
						anc_seqs[(matrix_counter * nSites)+j] = pos1;
						anc_seqs[(matrix_counter * nSites)+(str->partner->site)] = pos2;
						counter++;
						}
					}
				log_P[matrix_counter] = prob;
				matrix_counter++;
				}
			}
		}
}

void SampleHiddenAncestralSequences(TreeNode **downPass, TreeNode *root, SecondStruct *SS, double *log_P,  double *parameters, int *structure, long int *seed, int *anc_seqs, int nsamples, int nSites, int nTaxa, int nNodes, int condition_pairs)
{
	int i, j, n, m, k, counter=0, matrix_counter=0, which, pos1, pos2;
	double ran, sum, sumL, sumA, *fracs, *posts, prob;
	double  tprobs[4][4];
	double tprobs_rna[16][16];
	TreeNode *p;
	SecondStruct *str;
	
	matrix_counter = nsamples*(nTaxa-2);
	
	// This is completely LAME, not sure why the pointer root is NULL
	for(i = 0; i < nNodes; i++)
		{
		p = downPass[i];
		if(p->anc == NULL)
			root = p;
		}
	
	// Allocate some memory
	fracs = (double *)malloc(sizeof(double) * 20);
	if(!fracs)
		PrintErrorAndExitSafely("SampleHiddenAncestralSequences Error: Failure to allocate fracs.");
	posts = (double *)malloc(sizeof(double) * 20);
	if(!posts)
		PrintErrorAndExitSafely("SampleHiddenAncestralSequences Error: Failure to allocate posts.");
		
	// Pass down tree visiting internal nodes and sampling ancestral reconstructions
	for(i = 0; i < nNodes; i++)
		{
		p = downPass[i];
		if(p->left != NULL && p->right != NULL && p->anc != NULL)
			{
			for(m = 0; m < nsamples; m++)
				{
				prob = 0.0;
				counter = 0;
				for(j = 0; j < nSites; j++)
					{
					str = SS + j;
					if(structure[j] == 0)
						{
						// Set fracs and posts to zeros
						for(n = 0; n < 20; n++)
							{
							fracs[n] = 0.0;
							posts[n] = 0.0;
							}
						// Get transition probabilities for hidden node
						SetHiddenSingeltonTransitionProbabilities(p, tprobs, 1.0);
						// Get fractional likelihoods
						for(n = 0; n < NSTATES; n++)
							{
							sum = 0.0;
							for(k = 0; k < NSTATES; k++)
								{
								sum += tprobs[k][n] * p->frac_like[(counter*20)+16+n];						
								}
							fracs[16+n] = sum;
							}
						// Get posterior probabilities
						sum = 0.0;
						for(k = 0; k < NSTATES; k++)
							{
							sumA = sumL = 0.0;					
							for(n = 0; n < NSTATES; n++)
								{
								sumA += tprobs[k][n] * p->frac_like[(counter*20)+16+n];
								sumL += root->t_probs[k][n] * fracs[16+n];
								}
							sum += posts[16+k] = (sumA * sumL) * parameters[28 + n];
							}
						// Normalize probabilities
						for(n = 0; n < NSTATES; n++)
							{
							posts[16+n] /= sum;
							//printf("%lf\n", posts[16+n]);
							}
						// Stochastically sample a state from the marginal posterior distribution
						sum = 0.0;
						ran = RandomNumber(seed);
						for(n = 0; n < NSTATES; n++)
							{
							sum += posts[16+n];
							if(ran < sum) // State has been selected
								{
								which = n;
								prob = log(posts[16+n]);
								n = NSTATES;
								}
							}
						// Assign the chosen state to the ancestral marix
						anc_seqs[(matrix_counter * nSites)+j] = which;
						counter++;
						}
					else if(structure[j] == 1)
						{
						// Set fracs and posts to zeros
						for(n = 0; n < 20; n++)
							{
							fracs[n] = 0.0;
							posts[n] = 0.0;
							}
						// Get transition probabilities for hidden node
						SetHiddenDoubletTransitionProbabilities(p, tprobs_rna, 1.0);
						// Get fractional likelihoods
						for(k = 0; k < NDBLSTATES; k++)
							{
							sum = 0.0;
							for(n = 0; n < NDBLSTATES; n++)
								{
								sum += tprobs_rna[k][n] * p->frac_like[(counter*20)+n];						
								}
							fracs[k] = sum;
							}
						// Get posterior probabilities
						sum = 0.0;
						for(k = 0; k < NDBLSTATES; k++)
							{
							sumA = sumL = 0.0;					
							for(n = 0; n < NDBLSTATES; n++)
								{
								sumA += tprobs_rna[k][n] * p->frac_like[(counter*20)+n];
								sumL += root->t_probs_rna[k][n] * fracs[n];
								}
							sum += posts[k] = (sumA * sumL) * parameters[12 + n];
							}
						// Normalize probabilities
						for(n = 0; n < NDBLSTATES; n++)
							{
							posts[n] /= sum;
							//printf("%lf\n", posts[n]);
							}
						// If user has requested that the reconstructions be conditional on the frequency of unpaired doublets in the real data at paired sites
						if(condition_pairs == YES && str->freq_canonical == 1.0)
							{
							// I think the easiest approach here is to get rid of non-pairing dinucs if the observed data has none for this column on the matrix
							sum = 0.0;
							for(n = 0; n < NSTATES; n++)
								{
								for(k = 0; k < NSTATES; k++)
									{
									if(IsCanonicalDoublet(n, k) == NO)
										posts[(n*NSTATES)+k] = 0.0;
									sum += posts[(n*NSTATES)+k];
									}
								}
							for(n = 0; n < NDBLSTATES; n++)
								posts[n] /= sum;
							}
						else if(condition_pairs == YES && str->freq_canonical < 1.0)
							{
							sum = 0.0;
							for(n = 0; n < NSTATES; n++)
								{
								for(k = 0; k < NSTATES; k++)
									{
									if(IsCanonicalDoublet(n, k) == YES)
										sum += posts[(n*NSTATES)+k] *= str->freq_canonical;
									else
										sum += posts[(n*NSTATES)+k] *= (1.0 - str->freq_canonical);
									}
								}
							for(n = 0; n < NDBLSTATES; n++)
								posts[n] /= sum;
							}
						// Stochastically sample a state from the marginal posterior distribution
						sum = 0.0;
						ran = RandomNumber(seed);
						for(n = 0; n < NDBLSTATES; n++)
							{
							sum += posts[n];
							if(ran < sum) // State has been selected
								{
								which = n;
								prob = log(posts[n]);
								n = NDBLSTATES;
								}
							}
						// Need to break doublet into singeltons and put them in the correct place in the matrix
						GetDoubletMembers(&pos1, &pos2, which);
						anc_seqs[(matrix_counter * nSites)+j] = pos1;
						anc_seqs[(matrix_counter * nSites)+(str->partner->site)] = pos2;
						counter++;
						}
					}
				log_P[matrix_counter] = prob;
				matrix_counter++;
				}
			}
		}
		
	free(fracs);
	free(posts);
}

void SetQMatrix(double *parameters)
{
	int	  i, j, k;
	double scaler, pi[NSTATES], rates[GTRPARAMS];
	
	for (i = 0; i < NSTATES; i++)
		{
		pi[i] = parameters[28 + i];
		//printf("%lf\t", pi[i]);
		}
	//printf("\n");

	for (i = 0; i < GTRPARAMS; i++)
		{
		rates[i] = parameters[6+i];
		//printf("%lf\t", rates[i]);
		}
	//printf("\n");

	/* set Q matrix to 0 */
	for (i = 0; i < NSTATES; i++)
		for (j = 0; j < NSTATES; j++)
			singlet_q_matrix[i][j] = 0.0;

	k = 0;
	scaler = 0.0;
	for (i = 0; i < NSTATES; i++)
		{
		for (j = i+1; j < NSTATES; j++)
			{
			singlet_q_matrix[i][i] -= (singlet_q_matrix[i][j] = pi[j] * rates[k]);
			singlet_q_matrix[j][j] -= (singlet_q_matrix[j][i] = pi[i] * rates[k]);
			scaler  += (singlet_q_matrix[i][j] * pi[i]);
			scaler  += (singlet_q_matrix[j][i] * pi[j]);
			k++;
			}
		}

	/* Rescale Q matrix */
	if(show_debug_info)
		printf("\nLoop Q-matrix:\n");
	scaler = 1.0 / scaler;
	for (i = 0; i < NSTATES; i++)
		{
		for (j = 0; j < NSTATES; j++)
			{
			singlet_q_matrix[i][j] *= scaler;
			if(show_debug_info)
				printf("%lf ", singlet_q_matrix[i][j]);
			}
		if(show_debug_info)
			printf("\n");
		}
	if(show_debug_info)
		printf("\n");
}

void SetDoubletQMatrix(double *parameters)
{
	int	  i, j, k;
	double scaler, pi[NDBLSTATES], rates[GTRPARAMS];
	
	for (i = 0; i < NDBLSTATES; i++)
		{
		pi[i] = parameters[12 + i];
		//printf("%lf\t", pi[i]);
		}
	//printf("\n");

	for (i = 0; i < GTRPARAMS; i++)
		{
		rates[i] = parameters[i];
		//printf("%lf\t", rates[i]);
		}
	//printf("\n");

	/* set Q matrix to 0 */
	for (i = 0; i < NDBLSTATES; i++)
		for (j = 0; j < NDBLSTATES; j++)
			doublet_q_matrix[i][j] = 0.0;

	scaler = 0.0;
	for (i = 0; i < NDBLSTATES; i++)
		{
		for (j = i+1; j < NDBLSTATES; j++)
			{
			k = r_ij[i][j];
			if(k != (-1))
				{
				doublet_q_matrix[i][i] -= (doublet_q_matrix[i][j] = pi[j] * rates[k]);
				doublet_q_matrix[j][j] -= (doublet_q_matrix[j][i] = pi[i] * rates[k]);
				scaler  += (doublet_q_matrix[i][j] * pi[i]);
				scaler  += (doublet_q_matrix[j][i] * pi[j]);
				}
			}
		}

	/* Rescale Q matrix */
	if(show_debug_info)
		printf("\nStem Q-matrix:\n");
	scaler = 1.0 / scaler;
	for (i = 0; i < NDBLSTATES; i++)
		{
		for (j = 0; j < NDBLSTATES; j++)
			{
			doublet_q_matrix[i][j] *= scaler;
			if(show_debug_info)
				printf("%lf ", doublet_q_matrix[i][j]);
			}
		if(show_debug_info)
			printf("\n");
		}
	if(show_debug_info)
		printf("\n");
}

void SetDoubletTransitionProbabilities(TreeNode **upPass, int nNodes)
{
	int      i, j, k, isComplex;
	double   rate, **q, *eigenValues, *eigvalsImag, **eigvecs, **inverseEigvecs, *c_ijk, **probs, length;
	complex  **Ceigvecs, **CinverseEigvecs;
	TreeNode *p;

	rate = 1.0;

	/* Allocate some memory */
	probs           = psdmatrix (NDBLSTATES);
	q               = psdmatrix (NDBLSTATES);
	eigenValues     = (double *)malloc((size_t) (NDBLSTATES * sizeof(double)));
	eigvalsImag     = (double *)malloc((size_t) (NDBLSTATES * sizeof(double)));
	eigvecs         = psdmatrix (NDBLSTATES);
	inverseEigvecs  = psdmatrix (NDBLSTATES);
	Ceigvecs        = pscmatrix (NDBLSTATES);
	CinverseEigvecs = pscmatrix (NDBLSTATES);

	for(i = 0; i < NDBLSTATES; i++)
		for(j = 0; j < NDBLSTATES; j++)
			q[i][j] = doublet_q_matrix[i][j];

	isComplex = GetEigens(NDBLSTATES, q, eigenValues, eigvalsImag, eigvecs, inverseEigvecs, Ceigvecs, CinverseEigvecs);
	if (isComplex == YES)
		{
		printf("Eigens are complex....\n");
		}
	else
		{
		c_ijk = (double *)malloc((size_t) (NDBLSTATES * NDBLSTATES * NDBLSTATES * sizeof(double)));
		CalcCijk (c_ijk, NDBLSTATES, eigvecs, inverseEigvecs);
		}

	for (k = 0; k < nNodes; k++)
		{
		p = upPass[k];
		if(p->anc != NULL)
			{
			if(p->anc->anc == NULL)
				{
				length = p->anc->left->current_length + p->anc->right->current_length;
				}
			else
				{
				length = p->current_length;				
				}
			if (isComplex == NO)
				{
				CalcPij(c_ijk, NDBLSTATES, eigenValues, p->current_length, rate, probs);
				}
			else
				{
				if(ComplexChangeMatrix (eigenValues, eigvalsImag, Ceigvecs, CinverseEigvecs, NDBLSTATES, probs, length, rate) == ERROR)
					printf ("Problem calculating transition probabilities for complex eigen values.\n");
				}
			for (i = 0; i < NDBLSTATES; i++)
				{
				for (j = 0; j < NDBLSTATES; j++)
					{
					p->t_probs_rna[i][j] = probs[i][j];
					}
				}
			}
		}

	if(show_debug_info)
		{
		printf("\nDoublet transition probabilities:\n");
		for (k = 0; k < nNodes; k++)
			{
			p = upPass[k];
			if(p->anc != NULL)
				{
				printf("Node %d:\n", p->index);
				for (i = 0; i < NDBLSTATES; i++)
					{
					printf("\t");
					for (j = 0; j < NDBLSTATES; j++)
						{
						printf("%lf  ", p->t_probs_rna[i][ j]);
						}
					printf("\n");
					}
				}
			}
		}

	free_psdmatrix (probs);
	free_psdmatrix (q);
	free (eigenValues);
	free (eigvalsImag);
	free_psdmatrix (eigvecs);
	free_psdmatrix (inverseEigvecs);
	free_pscmatrix (Ceigvecs);
	free_pscmatrix (CinverseEigvecs);
	if (isComplex == NO)
		free (c_ijk);
}

void SetHiddenDoubletTransitionProbabilities(TreeNode *p, double tprobs_rna[16][16], double lambda)
{
	int      i, j, isComplex;
	double   rate, **q, *eigenValues, *eigvalsImag, **eigvecs, **inverseEigvecs, *c_ijk, **probs;
	complex  **Ceigvecs, **CinverseEigvecs;

	rate = lambda;

	/* Allocate some memory */
	probs           = psdmatrix (NDBLSTATES);
	q               = psdmatrix (NDBLSTATES);
	eigenValues     = (double *)malloc((size_t) (NDBLSTATES * sizeof(double)));
	eigvalsImag     = (double *)malloc((size_t) (NDBLSTATES * sizeof(double)));
	eigvecs         = psdmatrix (NDBLSTATES);
	inverseEigvecs  = psdmatrix (NDBLSTATES);
	Ceigvecs        = pscmatrix (NDBLSTATES);
	CinverseEigvecs = pscmatrix (NDBLSTATES);

	for(i = 0; i < NDBLSTATES; i++)
		for(j = 0; j < NDBLSTATES; j++)
			q[i][j] = doublet_q_matrix[i][j];

	isComplex = GetEigens(NDBLSTATES, q, eigenValues, eigvalsImag, eigvecs, inverseEigvecs, Ceigvecs, CinverseEigvecs);
	if (isComplex == YES)
		{
		printf("Eigens are complex....\n");
		}
	else
		{
		c_ijk = (double *)malloc((size_t) (NDBLSTATES * NDBLSTATES * NDBLSTATES * sizeof(double)));
		CalcCijk (c_ijk, NDBLSTATES, eigvecs, inverseEigvecs);
		}

	if (isComplex == NO)
		{
		CalcPij(c_ijk, NDBLSTATES, eigenValues, p->avg_length_to_tips, rate, probs);
		}
	else
		{
		if(ComplexChangeMatrix (eigenValues, eigvalsImag, Ceigvecs, CinverseEigvecs, NDBLSTATES, probs, p->avg_length_to_tips, rate) == ERROR)
			printf ("Problem calculating transition probabilities for complex eigen values.\n");
		}
	for (i = 0; i < NDBLSTATES; i++)
		{
		for (j = 0; j < NDBLSTATES; j++)
			{
			tprobs_rna[i][j] = probs[i][j];
			}
		}

	if(show_debug_info)
		{
		printf("\nHidden doublet transition probabilities:\n");
		for (i = 0; i < NDBLSTATES; i++)
			{
			printf("\t");
			for (j = 0; j < NDBLSTATES; j++)
				{
				printf("%lf  ", tprobs_rna[i][ j]);
				}
			printf("\n");
			}
		}

	free_psdmatrix (probs);
	free_psdmatrix (q);
	free (eigenValues);
	free (eigvalsImag);
	free_psdmatrix (eigvecs);
	free_psdmatrix (inverseEigvecs);
	free_pscmatrix (Ceigvecs);
	free_pscmatrix (CinverseEigvecs);
	if (isComplex == NO)
		free (c_ijk);
}

void SetHiddenSingeltonTransitionProbabilities(TreeNode *p, double tprobs[4][4], double lambda)
{
	int      i, j, isComplex;
	double   rate, **q, *eigenValues, *eigvalsImag, **eigvecs, **inverseEigvecs, *c_ijk, **probs;
	complex  **Ceigvecs, **CinverseEigvecs;

	rate = lambda;

	/* Allocate some memory */
	probs           = psdmatrix (NSTATES);
	q               = psdmatrix (NSTATES);
	eigenValues     = (double *)malloc((size_t) (4 * sizeof(double)));
	eigvalsImag     = (double *)malloc((size_t) (4 * sizeof(double)));
	eigvecs         = psdmatrix (NSTATES);
	inverseEigvecs  = psdmatrix (NSTATES);
	Ceigvecs        = pscmatrix (NSTATES);
	CinverseEigvecs = pscmatrix (NSTATES);

	for(i = 0; i < NSTATES; i++)
		for(j = 0; j < NSTATES; j++)
			q[i][j] = singlet_q_matrix[i][j];

	isComplex = GetEigens(4, q, eigenValues, eigvalsImag, eigvecs, inverseEigvecs, Ceigvecs, CinverseEigvecs);
	if (isComplex == YES)
		{
		printf("Eigens are complex....\n");
		}
	else
		{
		c_ijk = (double *)malloc((size_t) (NSTATES * NSTATES * NSTATES * sizeof(double)));
		CalcCijk (c_ijk, 4, eigvecs, inverseEigvecs);
		}

	if (isComplex == NO)
		{
		CalcPij(c_ijk, 4, eigenValues, p->avg_length_to_tips, rate, probs);
		}
	else
		{
		if(ComplexChangeMatrix (eigenValues, eigvalsImag, Ceigvecs, CinverseEigvecs, 4, probs, p->avg_length_to_tips, rate) == ERROR)
			printf ("Problem calculating transition probabilities for complex eigen values.\n");
		}
	for (i = 0; i < NSTATES; i++)
		{
		for (j = 0; j < NSTATES; j++)
			{
			tprobs[i][j] = probs[i][j];
			}
		}

	if(show_debug_info)
		{
		printf("\nHidden singleton transition probabilities:\n");
		for (i = 0; i < NSTATES; i++)
			{
			for (j = 0; j < NSTATES; j++)
				{
				printf("\t\t%lf", tprobs[i][j]);
				}
			printf("\n");
			}
		}

	free_psdmatrix (probs);
	free_psdmatrix (q);
	free (eigenValues);
	free (eigvalsImag);
	free_psdmatrix (eigvecs);
	free_psdmatrix (inverseEigvecs);
	free_pscmatrix (Ceigvecs);
	free_pscmatrix (CinverseEigvecs);
	if (isComplex == NO)
		free (c_ijk);
}

void SetSingeltonTransitionProbabilities(TreeNode **upPass, int nNodes)
{
	int      i, j, k, isComplex;
	double   rate, **q, *eigenValues, *eigvalsImag, **eigvecs, **inverseEigvecs, *c_ijk, **probs, length;
	complex  **Ceigvecs, **CinverseEigvecs;
	TreeNode *p;

	rate = 1.0;

	/* Allocate some memory */
	probs           = psdmatrix (NSTATES);
	q               = psdmatrix (NSTATES);
	eigenValues     = (double *)malloc((size_t) (4 * sizeof(double)));
	eigvalsImag     = (double *)malloc((size_t) (4 * sizeof(double)));
	eigvecs         = psdmatrix (NSTATES);
	inverseEigvecs  = psdmatrix (NSTATES);
	Ceigvecs        = pscmatrix (NSTATES);
	CinverseEigvecs = pscmatrix (NSTATES);

	for(i = 0; i < NSTATES; i++)
		for(j = 0; j < NSTATES; j++)
			q[i][j] = singlet_q_matrix[i][j];

	isComplex = GetEigens(4, q, eigenValues, eigvalsImag, eigvecs, inverseEigvecs, Ceigvecs, CinverseEigvecs);
	if (isComplex == YES)
		{
		printf("Eigens are complex....\n");
		}
	else
		{
		c_ijk = (double *)malloc((size_t) (NSTATES * NSTATES * NSTATES * sizeof(double)));
		CalcCijk (c_ijk, 4, eigvecs, inverseEigvecs);
		}

	for (k = 0; k < nNodes; k++)
		{
		p = upPass[k];
		if(p->anc != NULL)
			{
			if(p->anc->anc == NULL)
				{
				length = p->anc->left->current_length + p->anc->right->current_length;
				}
			else
				{
				length = p->current_length;				
				}
			if (isComplex == NO)
				{
				CalcPij(c_ijk, 4, eigenValues, p->current_length, rate, probs);
				}
			else
				{
				if(ComplexChangeMatrix (eigenValues, eigvalsImag, Ceigvecs, CinverseEigvecs, 4, probs, length, rate) == ERROR)
					printf ("Problem calculating transition probabilities for complex eigen values.\n");
				}
			for (i = 0; i < NSTATES; i++)
				{
				for (j = 0; j < NSTATES; j++)
					{
					p->t_probs[i][j] = probs[i][j];
					}
				}
			}
		}

	if(show_debug_info)
		{
		printf("\nSingleton transition probabilities:\n");
		for (k = 0; k < nNodes; k++)
			{
			p = upPass[k];
			if(p->anc != NULL)
				{
				printf("\tNode index %d (time = %lf):\n", p->index, p->current_length);
				for (i = 0; i < NSTATES; i++)
					{
					for (j = 0; j < NSTATES; j++)
						{
						printf("\t\t%lf", p->t_probs[i][j]);
						}
					printf("\n");
					}
				}
			}
		}

	free_psdmatrix (probs);
	free_psdmatrix (q);
	free (eigenValues);
	free (eigvalsImag);
	free_psdmatrix (eigvecs);
	free_psdmatrix (inverseEigvecs);
	free_pscmatrix (Ceigvecs);
	free_pscmatrix (CinverseEigvecs);
	if (isComplex == NO)
		free (c_ijk);
}

void SetTreeStructureToDefaults(TreeNode *Nodes, TreeNode **upPass, TreeNode **downPass, TreeNode *root, int nNodes, int nSites)
{
	int i, j, k;
	TreeNode *p;

	root = NULL;
	for(i = 0; i < nNodes; i++)
		{
		upPass[i] = NULL;
		downPass[i] = NULL;
		p = Nodes + i;
		p->left = NULL;
		p->right = NULL;
		p->anc = NULL;
		p->current_state = NULLSTATE;
		p->index = NULLSTATE;
		p->d_matrix_pos = NULLSTATE;
		p->uncertain = NO;
		p->current_length = 0.0;
		p->avg_length_to_tips = 0.0;
		for(j = 0; j < 4; j++)
			for(k = 0; k < 4; k++)
				p->t_probs[j][k] = 0.0;
		for(j = 0; j < 16; j++)
			for(k = 0; k < 16; k++)
				p->t_probs_rna[j][k] = 0.0;
		for(j = 0; j < 20 * nSites; j++)
			{
			p->posterior[j] = 0.0;
			p->frac_like[j] = 0.0;
			}
		}
}

void SetUpTree(TreeNode *Node, TreeNode **upPass, TreeNode **downPass, TreeNode *root, TreeLabels *Names, char *s, int nTaxa)
{
	int    i = 0, k, j, counter = 0, convert = NO, numElementsToSkip = 0;
	char   temp_label[32], temp_len[32];
	double len;
	TreeNode *l, *r, *a, *p;
	TreeLabels *nm;

	int withBranches = 1;



	while (s[i] != NULLTERM)								/* BUILD TREE FROM FILE */
		{
		if (s[i] == LEFTPAR)									/* CREATE NODE - LEFT IF NULL ELSE RIGHT OR ROOT IF FIRST */
			{
			if (counter == 0)												/* root node */
				{
				p        = Node + numElementsToSkip +(counter);
				if(p==NULL)
					{
					printf("setUpStructureRepresentationOfTreeString err: problem when counter = 0 with ptr p");
					exit(1);
					}
				p->index = counter;
				root = p;
				counter++;
				}
			else if (p->left == NULL && p->right == NULL)					/* left node */
				{
				l        = Node + numElementsToSkip + (counter);
				l->index = counter;
				p->left  = 				l;
				l->anc   =				p;
				p        =				l;
				counter++;
				}
			else if (p->left != NULL && p->right == NULL)					/* right node */
				{
				r        = Node + numElementsToSkip + (counter);
				r->index = counter;
				p->right = 				r;
				r->anc   =				p;
				p        =				r;
				counter++;
				}
			else if (p->left != NULL && p->right != NULL)
				{
				a        = Node + numElementsToSkip + (counter);
				a->index = counter;
				a->left = p;
				p->anc = a;
				p = a;
				counter++;
				}
			}
		else if (s[i] == RIGHTPAR || s[i] == COMMA)								/* MOVE DOWN TO THE ANC NODE */
			{
			if(s[i] == RIGHTPAR && convert == YES)
				{
				if(withBranches)
					{
					if(p->left->current_length == 0.0)
						p = p->left;
					else if(p->right->current_length == 0.0)
						p = p->right;
					}
				}
			else
				{
				p = p->anc;
				}
			}
		else if (s[i] == COLON)									/* READ IN BRANCHLENGTH AND ASSIGN TO CURRENT NODE */
			{
			for (k = (i+1); k < ((i+1)+32); k++)
				{
				if (WhatChar(s[k]) == 0)
					temp_len[(k-(i+1))] = s[k];
				else if (WhatChar(s[k]) == 1 || WhatChar(s[k]) == 2)
					break;
				}
			temp_len[(k-(i+1))] = NULLTERM;
			i = (k-1);
			len = 0.0;
			if (p->anc == root && convert == YES)
				{
				sscanf(temp_len, "%lf", &len);
				p->current_length = len*(0.50);
				root->left->current_length = len*(0.50);
				p = root->left;
				}
			else
				{
				sscanf(temp_len, "%lf", &len);
				p->current_length = len;
				}
			}
		else if ( CompareChar(s[i]) == NO )				/* READ TAXON LABEL - IN THIS CASE WILL BE INTEGER AND ADD NEW NODE */
			{
			for (k = i; k < (k+32); k++)
				{
				if (WhatChar(s[k]) == 0)
					temp_label[k-i] = s[k];
				else if (WhatChar(s[k]) == 1)
					break;
				}
			temp_label[k-i] = NULLTERM;
			i = (k-1);
			if (p->left == NULL && p->right == NULL)						/* left node */
				{
				l        = Node +  numElementsToSkip + (counter);
				l->index = counter;
				p->left  = 				l;
				l->anc   =				p;
				p        =				l;
				counter++;
				}
			else if (p->left != NULL && p->right == NULL)					/* right node */
				{
				r        = Node +  numElementsToSkip + (counter);
				r->index = counter;
				p->right = 				r;
				r->anc   =				p;
				p        =				r;
				counter++;
				}
			else if (p->left != NULL && p->right != NULL && p->anc == NULL)	/* at root: tree being input is unrooted -> action = convert to rooted */
				{
				r        = Node +  numElementsToSkip + (counter);									/* New node for tip taxa */
				r->index = counter;
				counter++;
				a        = Node +  numElementsToSkip + (counter);
				a->index = counter;									/* New root node for tree, conversion from unrooted to rooted */
				root =      a;
				a->left  =              p;
				a->right =              r;
				p->anc   =              a;
				r->anc   =              a;
				p        =              r;
				convert  =            YES;
				}
			strcpy(p->label, temp_label);
			j = 0;
			while(j < nTaxa)
				{
				nm = Names + j;
				if(!strcmp(nm->nm_label, p->label))
					break;
				else
					j++;
				}
			p->d_matrix_pos = j;
			}
		else
			{
			/* Probably a white space - action = skip */
			printf("Shouldn't be here....(%c)\n", s[i]);
			exit(1);
			}
		i++;
		}

	// Re-index nodes
	i = 0;
	ReIndexNodesForTree(root, &i);

	// Get up and down pass sequences for tree
	i = 0;
	GetDownPassSequence(root, downPass, &i);
	i = 0;
	GetUpPassSequence(root, upPass, &i);
}

void ShowDataMatrix(TreeLabels *Names, int *matrix, int row, int col)
{
	int i, j;
	TreeLabels *nm;

	printf("Data matrix:\n");

	for(i = 0; i < row; i++)
		{
		nm = Names + i;
		printf("%s\t", nm->nm_label);
		for(j = 0; j < col; j++)
			{
			printf("%c", NucIntToChar(matrix[(i * col) + j]));
			}
		printf("\n");
		}
}

void ShowSecondaryStructure(int nSites, SecondStruct *SS)
{
	int i;
	SecondStruct *sst;
	
	printf("\n\nC-structure representation of secondary structure:\n");
	printf("\tsite\tpaired\tpartner\tfreq-canonical\tfreq-gaps\n");
	for(i = 0; i < nSites; i++)
		{
		sst = SS + i;
		if(sst->paired == YES)
			printf("\t%d\t%d\t%d\t%1.3lf\t%1.3lf\n", sst->site, sst->paired, sst->partner->site, sst->freq_canonical, sst->freq_gaps);
		else
			printf("\t%d\t%d\t*\t*\t%1.3lf\n", sst->site, sst->paired, sst->freq_gaps);
		}
}

void ShowAncestralSequences(double *log_P, int *structure, int *anc_seqs, int num_rows, int num_cols)
{
	int i, j;
	
	printf("\nAncestral sequences:\n");

	for(i = 0; i < num_cols; i++)
		{
		if(structure[i] == 0)
			printf(".");
		else if(structure[i] == 1)
			printf("(");
		else if(structure[i] == -1)
			printf(")");
		}
	printf("\n");
	for(i = 0; i < num_rows; i++)
		{
		for(j = 0; j < num_cols; j++)
			{
			WhichNuc(anc_seqs[(i*num_cols)+j]);
			}		
		printf("\t[P = %1.12lf]\n", exp(log_P[i]));
		}
	printf("\n");
}

void SimulateHiddenAncestralSequences(TreeNode **downPass, TreeNode *root, SecondStruct *SS, double *log_P,  double *parameters, int *structure, long int *seed, int *anc_seqs, int nsamples, int nSites, int nTaxa, int nNodes, double lambda, int condition_pairs)
{
	int i, j, n, m, k, counter=0, matrix_counter=0, anc_state, son_state, pos1, pos2;
	double ran, sum;
	double  tprobs[4][4];
	double tprobs_rna[16][16];
	TreeNode *p;
	SecondStruct *str;
	
	matrix_counter = nsamples*(nTaxa-2);
	
	// Pass down tree visiting internal nodes and sampling ancestral reconstructions
	for(i = 0; i < nNodes; i++)
		{
		p = downPass[i];
		if(p->left != NULL && p->right != NULL && p->anc != NULL)
			{
			for(m = 0; m < nsamples; m++)
				{
				counter = 0;
				for(j = 0; j < nSites; j++)
					{
					str = SS + j;
					if(structure[j] == 0)
						{
						// Get transition probabilities for hidden node, the branch length used here is the harmonic mean of lengths
						// across the tree multiplied by a user defined tuning parameter, lambda
						SetHiddenSingeltonTransitionProbabilities(p, tprobs, lambda);
						// Sample an ancestral node from the posterior
						sum = 0.0;
						ran = RandomNumber(seed);
						for(n = 0; n < NSTATES; n++)
							{
							sum += p->posterior[(counter*20)+16+n];
							if(ran < sum) // State has been selected
								{
								anc_state = n;
								n = NSTATES;
								}
							}
						// Determine simulated son state
						sum = 0.0;
						ran = RandomNumber(seed);
						for(n = 0; n < NSTATES; n++)
							{
							sum += tprobs[anc_state][n];
							if(ran < sum) // State has been selected
								{
								son_state = n;
								n = NSTATES;
								}
							}
						// Assign the chosen state to the ancestral marix
						anc_seqs[(matrix_counter * nSites)+j] = son_state;
						counter++;
						}
					else if(structure[j] == 1)
						{
						// Get transition probabilities for hidden node, the branch length used here is the harmonic mean of lengths
						// across the tree multiplied by a user defined tuning parameter, lambda
						SetHiddenDoubletTransitionProbabilities(p, tprobs_rna, lambda);
						sum = 0.0;
						ran = RandomNumber(seed);
						for(n = 0; n < NDBLSTATES; n++)
							{
							sum += p->posterior[(counter*20)+n];
							if(ran < sum) // State has been selected
								{
								anc_state = n;
								n = NDBLSTATES;
								}
							}
						// Make the transition probs conditional on the frequency of canonical pairs
						if(condition_pairs == YES && str->freq_canonical == 1.0)
							{
							for(n = 0; n < NDBLSTATES; n++)
								{
								sum = 0.0;
								for(k = 0; k < NDBLSTATES; k++)
									{
									if(IsDoubletChangeCanonicalChange(k) == NO)
										tprobs_rna[n][k] = 0.0;
									sum += tprobs_rna[n][k];
									}
								for(k = 0; k < NDBLSTATES; k++)
									tprobs_rna[n][k] /= sum;
								}
							}
						else if(condition_pairs == YES && str->freq_canonical < 0.0)
							{
							for(n = 0; n < NDBLSTATES; n++)
								{
								sum = 0.0;
								for(k = 0; k < NDBLSTATES; k++)
									{
									if(IsDoubletChangeCanonicalChange(k) == YES)
										sum += tprobs_rna[n][k] *= str->freq_canonical;
									else
										sum += tprobs_rna[n][k] *= (1.0 - str->freq_canonical);
									}
								for(k = 0; k < NDBLSTATES; k++)
									tprobs_rna[n][k] /= sum;
								}
							}
						// Determine simulated son state
						sum = 0.0;
						ran = RandomNumber(seed);
						for(n = 0; n < NDBLSTATES; n++)
							{
							sum += tprobs_rna[anc_state][n];
							if(ran < sum) // State has been selected
								{
								son_state = n;
								n = NDBLSTATES;
								}
							}
						// Need to break doublet into singeltons and put them in the correct place in the matrix
						GetDoubletMembers(&pos1, &pos2, son_state);
						anc_seqs[(matrix_counter * nSites)+j] = pos1;
						anc_seqs[(matrix_counter * nSites)+(str->partner->site)] = pos2;
						counter++;
						}
					}
				matrix_counter++;
				}
			}
		}
}

int WhatChar(char c)
{
	if (c == (char)COMMA)
		return (1);
	else if (c == (char)RIGHTPAR)
		return (1);
	else if (c == (char)LEFTPAR)
		return (1);
	else if (c == (char)COLON)
		return (1);
	else if (c == (char)SPACE)
		return (2);
	else
		return (0);
}

int WhichDoublet(int one, int two)
{
	if(one == 0 && two == 0)
		return 0;
	else if(one == 0 && two == 1)
		return 1;
	else if(one == 0 && two == 2)
		return 2;
	else if(one == 0 && two == 3)
		return 3;
	else if(one == 1 && two == 0)
		return 4;
	else if(one == 1 && two == 1)
		return 5;
	else if(one == 1 && two == 2)
		return 6;
	else if(one == 1 && two == 3)
		return 7;
	else if(one == 2 && two == 0)
		return 8;
	else if(one == 2 && two == 1)
		return 9;
	else if(one == 2 && two == 2)
		return 10;
	else if(one == 2 && two == 3)
		return 11;
	else if(one == 3 && two == 0)
		return 12;
	else if(one == 3 && two == 1)
		return 13;
	else if(one == 3 && two == 2)
		return 14;

	return 15;
}

int WhichDoubleIsSmallest(double s1, double s2)
{
	int which=NULLSTATE;
	double t=10000000.0;

	if(t > s1)
		{
		t = s1;
		which = 0;
		}
	if(t > s2)
		{
		t = s2;
		which = 1;
		}
	
	return which;
}

void WhichNuc(int x)
{		
	if (x == 0)
		printf("A");
	else if (x == 1)
		printf("C");
	else if (x == 2)
		printf("G");
	else if (x == 3)
		printf("U");
	else
		{
		printf("\n\nERROR: Problem printing character state.\n");		
		exit(1);
		}
}

void WriteResultsToFile_New(TreeNode *Nodes, TreeNode **downPass, char *outName, double *log_P, int *structure, int *anc_seqs, int write_original_plus, TreeLabels *Names, int *dmatrix, int nSites, int nTaxa, int nsamples, double lambda, int is_harmonic)
{
	int i, j, n, counter, matrix_counter=0, offset;
	char p_out[110];
	FILE *out;
	TreeLabels *nm;
	TreeNode *p;

	// Open the file in which data and structure is to be read from
	out = fopen(outName,"w+");
	if(!out)
		PrintErrorAndExitSafely("WriteResultsToFile Error: Unable to open output file ptr.");

	// Write the structure information
	fprintf(out, ">SS_cons\n");
	for(i = 0; i < nSites; i++)
		{
		if(structure[i] == 0)
			fprintf(out, "-");
		else if(structure[i] == 1)
			fprintf(out, "(");
		else if(structure[i] == -1)
			fprintf(out, ")");
		}
	fprintf(out, "\n");

	// If user has requested that the original data be written along with the ancestral sequences
	if(write_original_plus)
		{
		for(i = 0; i < nTaxa; i++)
			{
			nm = Names + i;
			fprintf(out, ">%s\n", nm->nm_label);
			for(j = 0; j < nSites; j++)
				{
				fprintf(out, "%c", NucIntToChar(dmatrix[(i*nSites)+j]));
				}
			fprintf(out, "\n");
			}
		}

	matrix_counter = 0;
	// Write anc sequences to file: naming convention is anc_node_rep
	for(i = 0; i < (2*(nTaxa-1)); i++)
		{
		p = downPass[i];
		if(p->left != NULL && p->right != NULL && p->anc != NULL)
			{
			counter = 0;
			for(n = 0; n < nsamples; n++)
				{
				fprintf(out, ">anc_%d_%d\n", p->index, counter+1);
				for(j = 0; j < nSites; j++)
					{
					fprintf(out, "%c", NucIntToChar(anc_seqs[matrix_counter+(n*nSites)+j]));
					}
				fprintf(out, "\n");
				counter++;
				}
			matrix_counter += (nSites * nsamples);
			}
		}
	
	// Write imaginary (hairy node) sequences to file: naming convention is anc_node_rep
	offset = (2*(nTaxa-1));
	for(i = 0; i < (2*(nTaxa-1))-1; i++)
		{
		p = downPass[i];
		if(p->left != NULL && p->right != NULL && p->anc != NULL)
			{
			counter = 0;
			for(n = 0; n < nsamples; n++)
				{
				fprintf(out, ">psr_%d_%d\n", p->index, counter+1);
				for(j = 0; j < nSites; j++)
					{
					fprintf(out, "%c", NucIntToChar(anc_seqs[matrix_counter+(n*nSites)+j]));
					}
				fprintf(out, "\n");
				counter++;
				}
			matrix_counter += (nSites * nsamples);
			}
		}

	fclose(out);

	sprintf(p_out, "%s.par", outName);
	
	// Open the file in which data and structure is to be read from
	out = fopen(p_out,"w+");
	if(!out)
		PrintErrorAndExitSafely("WriteResultsToFile Error: Unable to open p_out file ptr.");
	
	p = downPass[0];
	if(is_harmonic == 1)
		{
		fprintf(out, "h_mean\t%lf\n", p->avg_length_to_tips);
		fprintf(out, "lambda\t%lf\n", lambda);
		fprintf(out, "size\t%d\n", nsamples);
		}
	else
		{
		fprintf(out, "prd_brlen\t%lf\n", p->avg_length_to_tips);
		fprintf(out, "size\t%d\n", nsamples);
		}

	fclose(out);
}

void WriteResultsToFile_Old(char *outName, double *log_P, int *structure, int *anc_seqs, int write_original_plus, TreeLabels *Names, int *dmatrix, int nSites, int nTaxa, int nsamples)
{
	int i, j, counter;
	FILE *out;
	TreeLabels *nm;

	// Open the file in which data and structure is to be read from
	out = fopen(outName,"w+");
	if(!out)
		PrintErrorAndExitSafely("WriteResultsToFile Error: Unable to open output file ptr.");

	// Write the structure information
	fprintf(out, ">SS_cons\n");
	for(i = 0; i < nSites; i++)
		{
		if(structure[i] == 0)
			fprintf(out, "-");
		else if(structure[i] == 1)
			fprintf(out, "(");
		else if(structure[i] == -1)
			fprintf(out, ")");
		}
	fprintf(out, "\n");

	// If user has requested that the original dat be written along with the ancestral sequences
	if(write_original_plus)
		{
		for(i = 0; i < nTaxa; i++)
			{
			nm = Names + i;
			fprintf(out, ">%s\n", nm->nm_label);
			for(j = 0; j < nSites; j++)
				{
				fprintf(out, "%c", NucIntToChar(dmatrix[(i*nSites)+j]));
				}
			fprintf(out, "\n");
			}
		}

	counter = 0;
	for(i = 0; i < ((nTaxa-1)+(nTaxa-2))*nsamples; i++)
		{
		if(i < (nTaxa-1)*nsamples)
			fprintf(out, ">anc%d\n", i+1);
		else
			{
			fprintf(out, ">img%d\n", counter+1);
			counter++;
			}
		for(j = 0; j < nSites; j++)
			{
			fprintf(out, "%c", NucIntToChar(anc_seqs[(i*nSites)+j]));
			}
		fprintf(out, "\n");
		}

	fclose(out);
}








