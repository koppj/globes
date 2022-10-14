/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NUM = 258,
     SFNCT = 259,
     BOGUS = 260,
     LVAR = 261,
     VAR = 262,
     FNCT = 263,
     IDN = 264,
     CROSS = 265,
     FLUXP = 266,
     FLUXM = 267,
     NUFLUX = 268,
     EFT_FLUX_COEFF_FILE = 269,
     EFT_XSEC_COEFF_FILE = 270,
     SYS_ON_FUNCTION = 271,
     SYS_OFF_FUNCTION = 272,
     SYS_MULTIEX_ERRORS = 273,
     GRP = 274,
     GID = 275,
     FNAME = 276,
     VERS = 277,
     SIGNAL = 278,
     BG = 279,
     ENERGY = 280,
     CHANNEL = 281,
     NDEF = 282,
     GRPOPEN = 283,
     GRPCLOSE = 284,
     PM = 285,
     FLAVOR = 286,
     NOGLOBES = 287,
     RULESEP = 288,
     RULEMULT = 289,
     NAME = 290,
     RDF = 291,
     ENDEXP = 292,
     ENDDET = 293,
     NEG = 294
   };
#endif
/* Tokens.  */
#define NUM 258
#define SFNCT 259
#define BOGUS 260
#define LVAR 261
#define VAR 262
#define FNCT 263
#define IDN 264
#define CROSS 265
#define FLUXP 266
#define FLUXM 267
#define NUFLUX 268
#define EFT_FLUX_COEFF_FILE 269
#define EFT_XSEC_COEFF_FILE 270
#define SYS_ON_FUNCTION 271
#define SYS_OFF_FUNCTION 272
#define SYS_MULTIEX_ERRORS 273
#define GRP 274
#define GID 275
#define FNAME 276
#define VERS 277
#define SIGNAL 278
#define BG 279
#define ENERGY 280
#define CHANNEL 281
#define NDEF 282
#define GRPOPEN 283
#define GRPCLOSE 284
#define PM 285
#define FLAVOR 286
#define NOGLOBES 287
#define RULESEP 288
#define RULEMULT 289
#define NAME 290
#define RDF 291
#define ENDEXP 292
#define ENDDET 293
#define NEG 294




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 1195 "glb_parser.y"
{
  double  val;  /* For returning numbers.                   */
  double *dpt;  /* for rules */
  glb_List *ptr;
  glb_List **ptrq;
  glb_symrec  *tptr;  /* For returning symbol-table pointers      */
  char *name;
  char *iname;
  int in;
  glb_namerec *nameptr;
}
/* Line 1529 of yacc.c.  */
#line 139 "glb_parser.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

