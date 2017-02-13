#include <stdlib.h>
#include <stdio.h>

#include "errorexit.h"

void PrintErrorAndExitSafely(char *string)
{
	printf("\n\nErr: %s\n", string);
	exit(1);
}
