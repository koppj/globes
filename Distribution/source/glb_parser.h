typedef union {
  double  val;  /* For returning numbers.                   */
  double *dpt;  /* for rules */
  glb_List *ptr; 
  glb_List **ptrq;
  glb_symrec  *tptr;  /* For returning symbol-table pointers      */
  char *name;
  char *iname;
  int in;
  glb_namerec *nameptr;
  
} YYSTYPE;
#define	NUM	257
#define	SFNCT	258
#define	VAR	259
#define	FNCT	260
#define	IDN	261
#define	CROSS	262
#define	FLUXP	263
#define	FLUXM	264
#define	GRP	265
#define	GID	266
#define	FNAME	267
#define	VERS	268
#define	SIGNAL	269
#define	BG	270
#define	GRPOPEN	271
#define	GRPCLOSE	272
#define	PM	273
#define	FLAVOR	274
#define	NOGLOBES	275
#define	CHANNEL	276
#define	RULESEP	277
#define	RULEMULT	278
#define	ENERGY	279
#define	NAME	280
#define	RDF	281
#define	NDEF	282
#define	NEG	283


extern YYSTYPE yylval;
