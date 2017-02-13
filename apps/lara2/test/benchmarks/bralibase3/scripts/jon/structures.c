#include <stdlib.h>
#include <stdio.h>

#include "precomp_definitions.h"
#include "structures.h"

void AddAncNode(int X, TreeNode *nodes)
{
	MutationList L;
	Position P;
	TreeNode *p;
	
	p = nodes + X;
	if(p == NULL)
		{
		printf("AddAncNode Error: Problem trying to assign anc node.\n");
		exit(1);
		}

	p->anc = NULL;
	p->left = NULL;
	p->right = NULL;
	p->current_length = 0.0;
	p->proposed_length = 0.0;
	p->index = X;
	p->d_matrix_pos = NULLSTATE;
	p->uncertain = NO;
	p->current_state = NULLSTATE;
	p->proposed_state = NULLSTATE;

	L = MakeEmpty(NULL);
	P = Header(L);
	p->ptrToList = L;
	if(p->ptrToList == NULL)
		{
		printf("Problem in attaching 'ptrToList'.\n");
		exit(1);
		}

	//printf("Root node assigned (%d).\n", p->index);
}

void AddLeftNode(int X, double len, TreeNode *p, TreeNode *nodes)
{
	MutationList L;
	Position P;
	TreeNode *q;

	q = nodes + X;
	if(q == NULL)
		{
		printf("AddLeftNode Error: Problem trying to assign node.\n");
		exit(1);
		}

	q->anc = p;
	q->left = NULL;
	q->right = NULL;
	q->current_length = len;
	q->proposed_length = len;
	q->index = X;
	q->d_matrix_pos = NULLSTATE;
	q->uncertain = NO;
	q->current_state = NULLSTATE;
	q->proposed_state = NULLSTATE;

	L = MakeEmpty(NULL);
	P = Header(L);
	q->ptrToList = L;
	if(q->ptrToList == NULL)
		{
		printf("Problem in attaching 'ptrToList'.\n");
		exit(1);
		}
	
	p->left = q;

	//printf("Node (%d) assigned to the left (%d).\n", q->index, p->index);
}

void AddRightNode(int X, double len, TreeNode *p, TreeNode *nodes)
{
	MutationList L;
	Position P;
	TreeNode *q;

	q = nodes + X;
	if(q == NULL)
		{
		printf("AddRightNode Error: Problem trying to assign node.\n");
		exit(1);
		}

	q->anc = p;
	q->left = NULL;
	q->right = NULL;
	q->current_length = len;
	q->proposed_length = len;
	q->index = X;
	q->d_matrix_pos = NULLSTATE;
	q->uncertain = NO;
	q->current_state = NULLSTATE;
	q->proposed_state = NULLSTATE;

	L = MakeEmpty(NULL);
	P = Header(L);
	q->ptrToList = L;
	if(q->ptrToList == NULL)
		{
		printf("Problem in attaching 'ptrToList'.\n");
		exit(1);
		}
	
	p->right = q;

	//printf("Node (%d) assigned to the right (%d).\n", q->index, p->index);
}

Position Advance(Position P)
{
	return (P->next);
}

Position AdvanceNumberOfElements(int num_elements, Position P)
{
	int i=0;
	
	while(P != NULL && (i++) < num_elements)
		P = P->next;
		
	return P;
}

double* AllocateMemoryForFractionalLikelihoods(int X)
{
	double *a;
	
	a = malloc(sizeof(double) * X);
	if(!a)
		{
		printf("AllocateMemoryForFractionalLikelihoods Error(1): Problem allocating memory.\n");
		exit(1);
		}

	return a;
}

TreeNode* AllocateMemoryForNodeStruct(int num_nodes)
{
	TreeNode *nodes;

	nodes = (TreeNode *)malloc(sizeof(TreeNode)*num_nodes);
	if (nodes == NULL)
		{
		printf("AllocateMemoryForTree Error(1): Unable to allocate memory.\n");
		exit(1);
		}
		
	return nodes;
}

TreeNode** AllocateMemoryForPointerArray(int num_nodes)
{
    TreeNode **d;
	
	d = (TreeNode **)malloc(sizeof(TreeNode *) * num_nodes);
    if(!d)
        {
		printf("AllocateMemoryForPointerArray Error(1): Unable to allocate memory.\n");
		exit(1);
        }
		
	return d;
}

double* AllocateMemoryForPosteriorProbs(int X)
{
	double *a;
	
	a = malloc(sizeof(double) * X);
	if(!a)
		{
		printf("AllocateMemoryForPosteriorProbs Error(1): Problem allocating memory.\n");
		exit(1);
		}

	return a;
}

TreeStrings* AllocateMemoryForTreeStringStructure(int num_nodes)
{
	TreeStrings *d;

	d = (TreeStrings *)malloc(sizeof(TreeStrings) * num_nodes);
	if(!d)
		{
		printf("AllocateMemoryForTreeStringStructure Error(1): Unable to allocate memory.\n");
		exit(1);
		}
		
	return d;
}

TreeLabels* AllocateMemoryForTaxonNameStructure(int num_nodes)
{
    TreeLabels *d;
	
	d = (TreeLabels *)malloc(sizeof(TreeLabels) * num_nodes);
    if(!d)
        {
		printf("AllocateMemoryForTaxonNameStructure Error(1): Unable to allocate memory.\n");
		exit(1);
        }
		
	return d;
}

SecondStruct* AllocateMemoryForDataStructure(int num_sites)
{
	SecondStruct *d;

	d = (SecondStruct *)malloc(sizeof(SecondStruct) * num_sites);
	if(!d)
		{
		printf("AllocateMemoryForDataStructure Error(1): Unable to allocate memory.\n");
		exit(1);
		}
		
	return d;
}

void Delete(int X, MutationList L)
{
	Position P, TmpCell;

	P = FindPrevious(X, L);

	/* Assumption use of a header node */
	if (!IsLast(P, L))
		{
		/* Position (X) is found; delete it */
		TmpCell = P->next;

		P->next = TmpCell->next;
		free(TmpCell);
		}
}

void DeleteList(MutationList L)
{
  Position P, Tmp;

  P = L->next;
  L->next = NULL;
  
  while(P != NULL)
	{
	Tmp = P->next;
	free(P);
	P = Tmp;
	}
}

Position First(MutationList L)
{
	return (L->next);
}

Position Find(int X, MutationList L)
{
	Position P;

	P = L->next;
	while (P->next != NULL && P->next->index != X)
		{
		P = P->next;
		}

	return P;
}

Position FindPrevious(int X, MutationList L)
{
	Position P;

	P = L;
	while (P->next != NULL && P->next->index != X)
		{
		P = P->next;
		}

	return P;
}

Position FindSite(int X, Position C, TreeNode *p, int which)
{
	Position P;
	P = C;

	if(X > C->site)
		{
		while(P->site != X && P->next_site != NULL)
			P = P->next_site;
		//printf("A: Looking for site %d and found site %d\n", X, P->site);
		return P;
		}
	else if(X < C->site)
		{
		while(P->site != X && P->prev_site != NULL)
			P = P->prev_site;
		//printf("B: Looking for site %d and found site %d\n", X, P->site);
		return P;
		}
	else if(X == C->site)
		{
		if(C->prev_site != NULL)
			{
			if(C->prev_site->site == X)
				{
				while(P->site != X && P->prev_site != NULL)
					P = P->prev_site;
				//printf("C: Looking for site %d and found site %d\n", X, P->site);
				return P;
				}
			else
				{
				//printf("D: Looking for site %d and found site %d\n", X, P->site);
				return P;
				}
			}
		else
			{
			//printf("E: Looking for site %d and found site %d\n", X, P->site);
			return P;
			}
		}
	/*
		THIS IS FROM THE CPG PROJECT CODE. IT MAY NEED TO NE CHANGED. AS OF APRIL 18 2005.
	
	  if(X == C->site)
		{
		if(P->prev_site == NULL)
			{
			if(which == 0)
				P = p->start_current;
			else
				P = p->start_proposal;
			}
		else if(P->prev_site->next_site != P)
			{
			P = P->prev_site;
			P = P->next_site;
			}
		return P;
		}*/

	return NULL;
}

void GetDownPassSequence(TreeNode *p, TreeNode **s, int *i)
{
    if (p != NULL)
        {
        GetDownPassSequence(p->left,  s, i);
        GetDownPassSequence(p->right, s, i);
        s[(*i)++] = p;
        }
}

void GetUpPassSequence(TreeNode *p, TreeNode **s, int *i)
{
    if (p != NULL)
        {
        s[(*i)++] = p;
        GetUpPassSequence(p->left, s, i);
        GetUpPassSequence(p->right, s, i);
        }	
}

Position Header(MutationList L)
{
	return L;
}

void Insert(int X, int site_num, MutationList L, Position P)
{

	Position Temp;

	Temp = malloc( sizeof( struct List ) );
	if(Temp == NULL)
		{
		printf("ERROR: Can not allocate additional list element.\n");
		exit(1);
		}

	Temp->index = X;
	Temp->site = site_num;
	Temp->next = P->next;
	P->next = Temp;
}

int IsEmpty(MutationList L)
{
	return (L->next == NULL);
}

int IsLast(Position P, MutationList L)
{
	return (P->next == NULL);
}

MutationList MakeEmpty(MutationList L)
{
	if (L != NULL)
		{
		DeleteList(L);
		}

	L = malloc(sizeof(struct List));
	if (L == NULL)
		{
		printf("MakeEmpty Error: Unable to allocate more memory.\n");
		exit(1);
		}

	L->next     = NULL;
	L->next_mut = NULL;
	L->next_site = NULL;
	L->prev_site = NULL;
	L->isLast = YES;
	L->state    = 0;
	L->index	= 0;
	L->site     = 0;
	L->duration = 0.0;
	
	return L;
}

void PrintList(MutationList L)
{
	Position P = Header(L);

	if ( IsEmpty(L) )
		{
		printf("ERROR: List is empty.\n");
		}
	else
		{
		printf("Mutational list:\n");
		printf("\td = element index\n");
		printf("\ts = state\n");
		printf("\tt = time\n");
		printf("\tne = index of next element\n");
		printf("\tns = index of first element of next site\n");
		printf("\tnm = index of next mutation in history\n");
		printf("\t[0,1] = is last interval in current history\n");
		do
			{
			if(P->next_mut == NULL && P->next_site != NULL)
				printf("\t%d: s=%d t=%lf (ne=%d, ns=%d, nm=NULL) (%d) [%d]\n", RetrieveIndex(P), RetrieveState(P), RetrieveDuration(P), P->next->index, P->next_site->index, P->site, P->isLast);
			else if(P->next_mut != NULL && P->next_site == NULL)
				printf("\t%d: s=%d t=%lf (ne=%d, ns=NULL, nm=%d) (%d) [%d]\n", RetrieveIndex(P), RetrieveState(P), RetrieveDuration(P), P->next->index, P->next_mut->index, P->site, P->isLast);
			else if(P->next_mut == NULL && P->next_site == NULL)
				printf("\t%d: s=%d t=%lf (ne=end, ns=NULL, nm=NULL) (%d) [%d]\n", RetrieveIndex(P), RetrieveState(P), RetrieveDuration(P), P->site, P->isLast);
			else
				printf("\t%d: s=%d t=%lf (ne=%d, ns=%d, nm=%d) (%d) [%d]\n", RetrieveIndex(P), RetrieveState(P), RetrieveDuration(P), P->next->index, P->next_site->index, P->next_mut->index, P->site, P->isLast);
			P = Advance(P);
			} while (P != NULL);
		printf( "\n" );
		}
}

double RetrieveDuration(Position P)
{
	return (P->duration);
}

int RetrieveIndex(Position P)
{
	return (P->index);
}

int RetrieveState(Position P)
{
	return (P->state);
}

void ShowTreeStructure(int X, TreeNode *nodes, TreeNode **d, TreeNode **u)
{
	int i;
	TreeNode *p;

	printf("\nTree structure and ptr traversal sequences.\n\tTable representation of tree (index: left index,right index, anc index\n\t\t(length:avg) [matrix position, label]):\n");

	for(i = 0; i < X; i++)
		{
		p = nodes + i;
		if(p->left == NULL && p->right == NULL && p->anc != NULL)
			printf("\t%3d\t:nil,nil,%3d\t(%1.5lf:%1.5lf) [%d, %s]\n", p->index, p->anc->index, p->current_length, p->avg_length_to_tips, p->d_matrix_pos, p->label);
		else if(p->left != NULL && p->right != NULL && p->anc != NULL)
			printf("\t%3d\t:%3d,%3d,%3d\t(%1.5lf:%1.5lf)\n", p->index, p->left->index, p->right->index, p->anc->index, p->current_length, p->avg_length_to_tips);
		else if(p->left != NULL && p->right != NULL && p->anc == NULL)
			printf("\t%3d\t:%3d,%3d,nil\t(%1.5lf:%1.5lf)\n", p->index, p->left->index, p->right->index, p->current_length, p->avg_length_to_tips);
		else
			printf("ShowTreeStructure Error: Problem with pointers in tree structure.\n");
		}
	
	printf("\tPost-order and Pre-order traversal sequences by node index:\n");
	printf("\t\tPost-order:\t");
	for(i = 0; i < X; i++)
		{
		p = d[i];
		printf("%d\t", p->index);
		}
	printf("\n\t\tPre-order :\t");
	for(i = 0; i < X; i++)
		{
		p = u[i];
		printf("%d\t", p->index);
		}
	printf("\n\n");
}

