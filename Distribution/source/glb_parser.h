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
     GRP = 267,
     GID = 268,
     FNAME = 269,
     VERS = 270,
     SIGNAL = 271,
     BG = 272,
     GRPOPEN = 273,
     GRPCLOSE = 274,
     PM = 275,
     FLAVOR = 276,
     NOGLOBES = 277,
     CHANNEL = 278,
     RULESEP = 279,
     RULEMULT = 280,
     ENERGY = 281,
     NAME = 282,
     RDF = 283,
     NDEF = 284,
     NEG = 285
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
#define GRP 267
#define GID 268
#define FNAME 269
#define VERS 270
#define SIGNAL 271
#define BG 272
#define GRPOPEN 273
#define GRPCLOSE 274
#define PM 275
#define FLAVOR 276
#define NOGLOBES 277
#define CHANNEL 278
#define RULESEP 279
#define RULEMULT 280
#define ENERGY 281
#define NAME 282
#define RDF 283
#define NDEF 284
#define NEG 285




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 942 "glb_parser.y"
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
#line 108 "glb_parser.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



