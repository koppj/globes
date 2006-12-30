/* A Bison parser, made by GNU Bison 1.875.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NUM = 258,
     SFNCT = 259,
     LVAR = 260,
     VAR = 261,
     FNCT = 262,
     IDN = 263,
     CROSS = 264,
     FLUXP = 265,
     FLUXM = 266,
     SYS_ON_FUNCTION = 267,
     SYS_OFF_FUNCTION = 268,
     GRP = 269,
     GID = 270,
     FNAME = 271,
     VERS = 272,
     SIGNAL = 273,
     BG = 274,
     GRPOPEN = 275,
     GRPCLOSE = 276,
     PM = 277,
     FLAVOR = 278,
     NOGLOBES = 279,
     CHANNEL = 280,
     RULESEP = 281,
     RULEMULT = 282,
     ENERGY = 283,
     NAME = 284,
     RDF = 285,
     NDEF = 286,
     NEG = 287
   };
#endif
#define NUM 258
#define SFNCT 259
#define LVAR 260
#define VAR 261
#define FNCT 262
#define IDN 263
#define CROSS 264
#define FLUXP 265
#define FLUXM 266
#define SYS_ON_FUNCTION 267
#define SYS_OFF_FUNCTION 268
#define GRP 269
#define GID 270
#define FNAME 271
#define VERS 272
#define SIGNAL 273
#define BG 274
#define GRPOPEN 275
#define GRPCLOSE 276
#define PM 277
#define FLAVOR 278
#define NOGLOBES 279
#define CHANNEL 280
#define RULESEP 281
#define RULEMULT 282
#define ENERGY 283
#define NAME 284
#define RDF 285
#define NDEF 286
#define NEG 287




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 986 "glb_parser.y"
typedef union YYSTYPE {
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
/* Line 1249 of yacc.c.  */
#line 112 "glb_parser.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



