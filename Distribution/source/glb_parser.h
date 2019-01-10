/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

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

#ifndef YY_YY_GLB_PARSER_H_INCLUDED
# define YY_YY_GLB_PARSER_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
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
    SYS_ON_FUNCTION = 269,
    SYS_OFF_FUNCTION = 270,
    SYS_MULTIEX_ERRORS = 271,
    GRP = 272,
    GID = 273,
    FNAME = 274,
    VERS = 275,
    SIGNAL = 276,
    BG = 277,
    ENERGY = 278,
    CHANNEL = 279,
    NDEF = 280,
    GRPOPEN = 281,
    GRPCLOSE = 282,
    PM = 283,
    FLAVOR = 284,
    NOGLOBES = 285,
    RULESEP = 286,
    RULEMULT = 287,
    NAME = 288,
    RDF = 289,
    ENDEXP = 290,
    ENDDET = 291,
    NEG = 292
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
#define SYS_ON_FUNCTION 269
#define SYS_OFF_FUNCTION 270
#define SYS_MULTIEX_ERRORS 271
#define GRP 272
#define GID 273
#define FNAME 274
#define VERS 275
#define SIGNAL 276
#define BG 277
#define ENERGY 278
#define CHANNEL 279
#define NDEF 280
#define GRPOPEN 281
#define GRPCLOSE 282
#define PM 283
#define FLAVOR 284
#define NOGLOBES 285
#define RULESEP 286
#define RULEMULT 287
#define NAME 288
#define RDF 289
#define ENDEXP 290
#define ENDDET 291
#define NEG 292

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 1191 "glb_parser.y" /* yacc.c:1909  */

  double  val;  /* For returning numbers.                   */
  double *dpt;  /* for rules */
  glb_List *ptr;
  glb_List **ptrq;
  glb_symrec  *tptr;  /* For returning symbol-table pointers      */
  char *name;
  char *iname;
  int in;
  glb_namerec *nameptr;

#line 140 "glb_parser.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_GLB_PARSER_H_INCLUDED  */
