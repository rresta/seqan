#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "precomp_definitions.h"
#include "errorexit.h"
#include "tokenize.h"

void StringToken(int *start, int *finish, int s_length, char *s, char *delimit, char *token)
{	
	int i;
		
	if((*start) < s_length)
		{
		// Move the start to first non-delimiter character
		while(strchr_b(delimit, (char)s[(*start)]) == 1 && (*start) < s_length)
			(*start)++;

		(*finish) = (*start);
		// Set finish to start and move it to next delimiter character
		while(strchr_b(delimit, (char)s[(*finish)]) == 0 && (*finish) < s_length)
			(*finish)++;
		
		for(i = 0; i < (*finish)-(*start); i++)
			token[i] = s[(*start)+i];
		token[i] = '\0';
		
		//printf("token (%d, %d): (%s)\n", (*start), (*finish), token);
		(*start) = (*finish);
		}
	else
		{
		token[0] = '\0';
		}
}

int strchr_b(char *delimit, char x)
{
	int n, i;
	
	n = strlen(delimit);

	for(i = 0; i < n; i++)
		{
		if(delimit[i] == x)
			return 1;
		}
	
	return 0;
}
