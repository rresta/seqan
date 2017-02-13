// Define structure
typedef struct Node
	{
	struct Node					*left, *right, *anc;
	double			current_length, proposed_length;
	double				avg_length_to_tips, cum_sum, harmonic;
	double					 *frac_like, *posterior;
	char					label[MAX_LABEL_LENGTH];
	int	             index, d_matrix_pos, uncertain;
	int				  current_state, proposed_state;
	int								   size_of_list;
	void				   *ptrToList, *current_pos;
	void			*start_current, *start_proposal;
	double  t_probs[4][4];
	double t_probs_rna[16][16];
	
	struct List
		{
		struct List *next, *next_site, *prev_site, *next_mut;
		int						  state, site, index, isLast;
		double										duration;
		} MutationList;
	} TreeNode;

typedef struct List *MutationList;
typedef MutationList Position;

typedef struct names
	{
	struct names       *next, *prev;
	char nm_label[MAX_LABEL_LENGTH];
	int index;
	} TreeLabels;

typedef struct sec_struct
	{
	struct sec_struct       *next, *prev, *partner;
	int site, paired;
	double freq_canonical, freq_gaps;
	} SecondStruct;

typedef struct process
	{
	int     n_ij[NUM_NUC_STATES][NUM_NUC_STATES];
	double			 dwell_times[NUM_NUC_STATES];
	} SubProcess;

typedef struct tree_strings
    {
    int       index;
    char  *t_string;
    } TreeStrings;

// Functions to manipulate structures
Position		Advance(Position P);
Position		AdvanceNumberOfElements(int num_elements, Position P);
void			AddAncNode(int X, TreeNode *nodes);
void			AddLeftNode(int X, double len, TreeNode *p, TreeNode *nodes);
void			AddRightNode(int X, double len, TreeNode *p, TreeNode *nodes);
double*		AllocateMemoryForFractionalLikelihoods(int X);
SecondStruct*	AllocateMemoryForDataStructure(int num_sites);
TreeNode*		AllocateMemoryForNodeStruct(int num_nodes);
double*		AllocateMemoryForPosteriorProbs(int X);
TreeNode**	AllocateMemoryForPointerArray(int num_nodes);
TreeStrings*	AllocateMemoryForTreeStringStructure(int num_nodes);
TreeLabels*	AllocateMemoryForTaxonNameStructure(int num_nodes);
void			Delete(int X, MutationList L);
void			DeleteList(MutationList L);
Position		Find(int X, MutationList L);
Position		FindPrevious(int X, MutationList L);
Position		FindSite(int X, Position C, TreeNode *p, int which);
Position		First(MutationList L);
void			GetDownPassSequence(TreeNode *p, TreeNode **s, int *i);
void			GetUpPassSequence(TreeNode *p, TreeNode **s, int *i);
Position		Header(MutationList L);
void			Insert(int X, int site_num, MutationList L, Position P);
int			IsEmpty(MutationList L);
int			IsLast(Position P, MutationList L);
MutationList	MakeEmpty(MutationList L);
void			PrintList(MutationList L);
double		RetrieveDuration(Position P);
int			RetrieveIndex(Position P);
int			RetrieveState(Position P);
void			ShowTreeStructure(int X, TreeNode *nodes, TreeNode **d, TreeNode **u);


