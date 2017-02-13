// The following are debugging switches. They turn on a variety of verbose information.
#undef		D_GETFILECONTENTS
#undef		SHOW_TI_PROBS_IN_FUNC
#undef		SHOW_TI_PROBS
#undef		SHOW_FRAC_LIKE
#undef		SHOW_ANC_STATE
#undef		SHOW_MAP_INFO
#undef		SHOW_TABULATION_OF_NIJ
#undef		SHOW_NULL_MATRIX

// These are defined here and are sometimes defined by the compiler on particular machines (in these case a compiler warning may be issued).
#undef		NO_ERROR
#undef		ERROR
#define		NO_ERROR	0
#define		ERROR		1

#undef		FALSE
#undef		TRUE
#define		FALSE		0
#define		TRUE		1

#define		NO			 0
#define		YES			 1

// This is used as a default value for character states, which are all positive integers.
#define		NULLSTATE   -9

// The following are some integer/double definitions. Caution should be made in changing these because some loops move over arrays of fixed size.
#define		NUM_CODONS				64
#define		NUM_NUC_STATES			4
#define		NUM_PROCESSES			144
#define		MAX_LABEL_LENGTH		100
#define		MAXELEMS				10000000
#define		SMALL_STRING_SIZE 		500
#define		FILENAMESIZE 			100
#define		NR_END					1
#define		FREE_ARG				char*
#define		NRATECATS				4
#define		NSTATES					4
#define		NDBLSTATES				16
#define		MODELPARAMS				14
#define		GTRPARAMS				6
#define		MODELVALUES				(GTRPARAMS + NSTATES + 1)
#define		MAXEVENTS				100
#define		MAXTREES				5000
#define		MAXSPECIES				1000
#define		BURNIN					1
#define		ORIGINAL				1
#define		CONTINUITY_CORRECTION 	1.0
#define		CONSTANT				1.0
#define		CRITICAL_CUTOFF			0.95
#define		COMMA					','
#define		SEMICOLON				';'
#define		RIGHTPAR    			')'
#define		LEFTPAR					'('
#define		COLON					':'
#define		SPACE	    			' '
#define		PERIOD					'.'
#define		EQUALSIGN				'='
#define		HTAB					'\t'
#define		VTAB					'\v'
#define		BKSPACE					'\b'
#define		FORMFEED				'\f' /* 12 */
#define		CRETURN					'\r' /* 13 Used by UNIX and DOS (2nd EOL sequence) */
#define		NEWLINE					'\n' /* 10 Used by Mac and DOS (1st EOL sequence) */
#define		NULLTERM				'\0'
#define		PIPE					'|'
#define		LSQRBRACKET				'['
#define		RSQRBRACKET 			']'
#define		NULLCHAR    			''
