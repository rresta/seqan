#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "errorexit.h"

#define YES		1
#define NO			0

typedef struct names
	{
	struct names       *next, *prev;
	char nm_label[150];
	int index;
	} TreeLabels;

typedef struct sec_struct
	{
	struct sec_struct       *next, *prev, *partner;
	int site, paired;
	} SecondStruct;
	
SecondStruct *SS;
TreeLabels *Names;

int isStructure;

void ConvertFile(char *dataName, char *outName);
int GetInputValues(const char *argv[], int commandline, char *dataName, char *outName);
int GetNumberOfTaxaSitesAndAllocate(char *fileName, int *nTaxa, int *nSites);
void InitializeStemsAndLoops(int nTaxa, int nSites, int *structure);
void ReadDataFile(char *fileName, int nTaxa, int nSites, int *nNodes, char *dmatrix, int *structure, int n);
void ShowSecondaryStructure(int nSites);
void WriteMrBayesDoubletFile(char *outName, int nTaxa, int nSites, char *dmatrix, int *structure);

int main (int argc, const char * argv[])
{
	int			commandline=NO, action;
	char			dataName[150], outName[150];
	
	isStructure = NO;
	
	if(argc > 1)
		commandline = YES;

	// Get input values - commandline or user input
	action = GetInputValues(argv, commandline, dataName, outName);
	if(action == 1)
		PrintErrorAndExitSafely("A problem with the input or commandline values was encountered. Check and retry.");

	ConvertFile(dataName, outName);

	return 0;
}

void ConvertFile(char *dataName, char *outName)
{
	int			 nTaxa, nSites, length, *structure;
	char			*dmatrix;
	
	// Get number of taxa and sites
	length = GetNumberOfTaxaSitesAndAllocate(dataName, &nTaxa, &nSites);
	// Allocate some things based on nTaxa and nSites
	structure = (int *)malloc(sizeof(int) * (nSites));
	if(!structure) { printf("ConvertFile Error: Problem allocating memory for structure.\n"); exit(1); }
	Names = (TreeLabels *)malloc(sizeof(TreeLabels) * (nTaxa));
	if(!Names) { printf("ConvertFile Error: Problem allocating memory for Names.\n"); exit(1); }
	SS = (SecondStruct *)malloc(sizeof(SecondStruct) * (nSites));
	if(!SS) { printf("ConvertFile Error: Problem allocating memory for Data.\n"); exit(1); }
	dmatrix = (char *)malloc(sizeof(char) * (nSites) * (nTaxa));
	if(!dmatrix) { printf("ConvertFile Error: Problem allocating memory for dmatrix.\n"); exit(1); }

	// Get data file contents -- data and structure
	ReadDataFile(dataName, nTaxa, nSites, 0, dmatrix, structure, length);
	
	// Intialize stems and loops in Data
	InitializeStemsAndLoops(nTaxa, nSites, structure);
		
	// Write MrBayes file for doublet analysis
	WriteMrBayesDoubletFile(outName, nTaxa, nSites, dmatrix, structure);

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
		//printf("%s\n", string);
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

int GetInputValues(const char *argv[], int commandline, char *dataName, char *outName)
{
	// Get user values by commandline or stdin
	if(commandline == 1)
		{
		strcpy (dataName,		argv[ 1]);
		strcpy (outName,		argv[ 2]);
		}
	else
		{
		printf("######################## INPUT OPTIONS ###########################\n");
		printf("\tEnter data/structure filename > ");
		scanf ("%s", dataName);
		printf("\tEnter output filename > ");
		scanf ("%s", outName);
		printf("##################################################################\n\n");
		}
		
	return 0;
}

void InitializeStemsAndLoops(int nTaxa, int nSites, int *structure)
{
	int i, j, tracker;
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
		s1->partner = NULL;
		}

	for(i = 0; i < nSites; i++)
		{
		s1 = SS + i;
		if(structure[i] == 1)
			{
			isStructure = YES;
			s1->paired = YES;
			tracker = 0;
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
			}
		}
		
	// DEBUG STATEMENT	
	//ShowSecondaryStructure(nSites);
}

void ReadDataFile(char *fileName, int nTaxa, int nSites, int *nNodes, char *dmatrix, int *structure, int n)
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
	
	while(fgets(string, n, in))
		{
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
				dmatrix[((nt-1) * ns) + (j++)] = x;
				}
			getSites=NO;
			getTaxa=YES;
			}
		isFirstSeq++;
		}
		
	free(string);
	fclose(in);
}

void ShowSecondaryStructure(int nSites)
{
	int i;
	SecondStruct *sst;
	
	printf("\nC-structure representation of secondary structure:\n");
	printf("\tsite\tpaired\tpartner\n");
	for(i = 0; i < nSites; i++)
		{
		sst = SS + i;
		if(sst->paired == YES)
			printf("\t%d\t%d\t%d\n", sst->site, sst->paired, sst->partner->site);
		else
			printf("\t%d\t%d\t-\n", sst->site, sst->paired);
		}
}

void WriteMrBayesDoubletFile(char *outName, int nTaxa, int nSites, char *dmatrix, int *structure)
{
	int i, j;
	FILE *out;
	SecondStruct *str;
	TreeLabels *nm;

	out = fopen(outName,"w+");
	if(!out)
		PrintErrorAndExitSafely("WriteMrBayesDoubletFile Error: Unable to open file.");

	// Write some begining info to file
	fprintf(out, "#NEXUS\n\nbegin data;\n\tdimensions ntax=%d nchar=%d;\n\tformat datatype=rna gap=- missing=? interleave=yes;\n\tmatrix\n", nTaxa, nSites);
	for(i = 0; i < nTaxa; i++)
		{
		nm = Names + i;
		fprintf(out, "\t%s\t\t", nm->nm_label);
		for(j = 0; j < nSites; j++)
			{
			fprintf(out, "%c", dmatrix[(i * nSites) + j]);
			}
		fprintf(out, "\n");
		}
	fprintf(out, ";\nend;\n\n");

	fprintf(out, "begin mrbayes;\n\tset autoclose=yes nowarn=yes;\n");
	
	if(isStructure)
		{
		fprintf(out, "\tpairs ");
		for(i = 0; i < nSites; i++)
			{
			str = SS + i;
			if(structure[i] == 1)
				fprintf(out, "%d:%d,", str->site+1, str->partner->site+1);
			}
		fseek(out, -1, SEEK_CUR);
		fprintf(out, ";\n\tcharset stems =");
		for(i = 0; i < nSites; i++)
			{
			str = SS + i;
			if(structure[i] == 1 || structure[i] == -1)
				fprintf(out, " %d", str->site+1);
			}
		fprintf(out, ";\n\tcharset loops =");
		for(i = 0; i < nSites; i++)
			{
			str = SS + i;
			if(structure[i] == 0)
				fprintf(out, " %d", str->site+1);
			}
		fprintf(out, ";\n");
		fprintf(out, "\tpartition by_stem_loop = 2:stems,loops;\n\tset partition = by_stem_loop;\n");
		fprintf(out, "\tlset applyto=(1)   nucmodel=doublet;\n\tlset applyto=(2)   nucmodel=4by4;\n\tlset applyto=(1) nst=6;\n\tlset applyto=(2) nst=6;\n\tunlink revmat=(all);\n\tmcmc ngen=1000000 savebrlens=yes nchains=1 nruns=1 printfreq=1000 samplefreq=1000;\n");
		}
	else
		{
		fprintf(out, "\tlset nucmodel=4by4 nst=6;\n\tmcmc ngen=1000000 savebrlens=yes nchains=1 nruns=1 printfreq=1000 samplefreq=1000;\n");
		}

	fprintf(out, "\tsump burnin=500 printtofile=yes;\n\tsumt burnin=500 contype=allcompat;\n");
	fprintf(out, "\tquit;\n");
	fprintf(out, "end;");

	fclose(out);
}
