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
     VAR = 260,
     FNCT = 261,
     IDN = 262,
     CROSS = 263,
     FLUXP = 264,
     FLUXM = 265,
     GRP = 266,
     GID = 267,
     FNAME = 268,
     VERS = 269,
     SIGNAL = 270,
     BG = 271,
     GRPOPEN = 272,
     GRPCLOSE = 273,
     PM = 274,
     FLAVOR = 275,
     NOGLOBES = 276,
     CHANNEL = 277,
     RULESEP = 278,
     RULEMULT = 279,
     ENERGY = 280,
     NAME = 281,
     RDF = 282,
     NDEF = 283,
     NEG = 284
   };
#endif
#define NUM 258
#define SFNCT 259
#define VAR 260
#define FNCT 261
#define IDN 262
#define CROSS 263
#define FLUXP 264
#define FLUXM 265
#define GRP 266
#define GID 267
#define FNAME 268
#define VERS 269
#define SIGNAL 270
#define BG 271
#define GRPOPEN 272
#define GRPCLOSE 273
#define PM 274
#define FLAVOR 275
#define NOGLOBES 276
#define CHANNEL 277
#define RULESEP 278
#define RULEMULT 279
#define ENERGY 280
#define NAME 281
#define RDF 282
#define NDEF 283
#define NEG 284




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 828 "glb_parser.y"
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
#line 107 "glb_parser.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



