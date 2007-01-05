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
     NUFLUX = 267,
     SYS_ON_FUNCTION = 268,
     SYS_OFF_FUNCTION = 269,
     GRP = 270,
     GID = 271,
     FNAME = 272,
     VERS = 273,
     SIGNAL = 274,
     BG = 275,
     GRPOPEN = 276,
     GRPCLOSE = 277,
     PM = 278,
     FLAVOR = 279,
     NOGLOBES = 280,
     CHANNEL = 281,
     RULESEP = 282,
     RULEMULT = 283,
     ENERGY = 284,
     NAME = 285,
     RDF = 286,
     NDEF = 287,
     NEG = 288
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
#define NUFLUX 267
#define SYS_ON_FUNCTION 268
#define SYS_OFF_FUNCTION 269
#define GRP 270
#define GID 271
#define FNAME 272
#define VERS 273
#define SIGNAL 274
#define BG 275
#define GRPOPEN 276
#define GRPCLOSE 277
#define PM 278
#define FLAVOR 279
#define NOGLOBES 280
#define CHANNEL 281
#define RULESEP 282
#define RULEMULT 283
#define ENERGY 284
#define NAME 285
#define RDF 286
#define NDEF 287
#define NEG 288




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 1019 "glb_parser.y"
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
#line 114 "glb_parser.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



