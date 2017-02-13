#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "precomp_definitions.h"
#include "errorexit.h"
#include "main.h"

int main(int argc, char *argv[])
{
	int			commandline=NO, action, debug_status, nsamples, output_original, condition_pairs, is_harmonic;
	double		lambda;
	char			dataName[FILENAMESIZE], treeName[FILENAMESIZE], outName[FILENAMESIZE], paramName[FILENAMESIZE];

	if(argc > 1)
		commandline = YES;

	// Get input values - commandline or user input
	action = GetInputValues(argc, argv, commandline, &is_harmonic, &lambda, &condition_pairs, &debug_status, &nsamples, &output_original, dataName, treeName, outName, paramName);
	if(action == 1)
		PrintErrorAndExitSafely("A problem with the input or commandline values was encountered. Check and retry.");

	CalcAncestralStates(is_harmonic, lambda, debug_status, condition_pairs, output_original, nsamples, dataName, treeName, outName, paramName);

	return 0;
}

int GetInputValues(int argc, char *argv[], int commandline, int *is_harmonic, double *lambda, int *condition_pairs, int *debug_status, int *nsamples, int *output_original, char *dataName, char *treeName, char *outName, char *paramName)
{
	// Get user values by commandline or stdin
	if(commandline == 1)
		{
		strcpy (dataName,		argv[ 1]);
		strcpy (treeName,		argv[ 2]);
		strcpy(paramName,		argv[ 3]);
		strcpy (outName,		argv[ 4]);
		*(is_harmonic) = atoi(argv[5]);
		*(lambda) = atof(argv[6]);
		*(nsamples) = atoi(argv[7]);
		*(output_original) =atoi(argv[8]);
		*(condition_pairs) = atoi(argv[9]);
		if(argc == 11)
			*(debug_status) = atoi(argv[10]);
		else
			*(debug_status) = 0;
		}
	else
		{
		printf("############################### INPUT OPTIONS #############################\n");
		printf("\tEnter data/structure filename > ");
		scanf ("%s", dataName);
		printf("\tEnter tree filename > ");
		scanf ("%s", treeName);
		printf("\tEnter parameter filename > ");
		scanf ("%s", paramName);
		printf("\tEnter output filename > ");
		scanf ("%s", outName);
		double tharm;
		printf("\tEnter 1 to use harmonic mean, 0 otherwise > ");
		scanf ("%lf", &tharm);
		*(is_harmonic) = tharm;
		double tlambda;
		printf("\tEnter tuning parameter lambda > ");
		scanf ("%lf", &tlambda);
		*(lambda) = tlambda;
		int nsamps;
		printf("\tEnter number of ancestral samples > ");
		scanf ("%d", &nsamps);
		*(nsamples) = nsamps;
		int original;
		printf("\tTo include original data in output enter 1, otherwise enter 0 > ");
		scanf ("%d", &original);
		*(output_original) = original;
		int cnd_pr;
		printf("\tCondition on freq of pairs 1, otherwise unconditional 0 > ");
		scanf ("%d", &cnd_pr);
		*(condition_pairs) = cnd_pr;
		int db_stat;
		printf("\tTo show debug info enter 1, otherwise enter 0 > ");
		scanf ("%d", &db_stat);
		*(debug_status) = db_stat;
		printf("###########################################################################\n\n");
		}
		
	return 0;
}