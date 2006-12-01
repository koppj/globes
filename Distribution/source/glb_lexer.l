/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2004,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */




%pointer
%option noyywrap


%x incl
%s grp
%s rule
%s chl
%s energy
%s cross
%s flux
%s defr
%s defc
%s defe
%s defx
%s deff
%x comment
%x str


%{
/* need this for the call to atof() below */
#include <math.h>
#include <string.h>

#include "glb_parser_type.h"
#include "glb_parser.h"
#include "glb_error.h"


 #define MAX_INCLUDE_DEPTH 10
static YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
static int include_stack_ptr = 0;
int glb_line_num=0;
char *glb_file_id;
static int line_buf[MAX_INCLUDE_DEPTH];
static char *file_name_buf[MAX_INCLUDE_DEPTH];
static int cond_buf; 
/* I am not sure wether this works with nested contexts */
static int cond_buf_i; 
static int cond_buf_str; 
static char *fnbuf,*fnbuf2;
static size_t len,i;

void glb_lexer_cleanup();


%}

EOE       #NEXT#
FN       [a-zA-Z0-9./][a-zA-Z0-9_./-]*

GLOB     GLoBES
RL       @signal
DIGIT    [0-9]
IVAL     $[a-zA-Z][a-zA-Z0-9_]*
ID       [a-zA-Z][a-zA-Z0-9_]*
LID      ![a-zA-Z][a-zA-Z0-9_]*

%%

  

"/*"         {cond_buf=YY_START;BEGIN(comment);}

<comment>[^*\n]*        /* eat anything that's not a '*' */
<comment>"*"+[^*/\n]*   /* eat up '*'s not followed by '/'s */
<comment>\n             ++glb_line_num;
<comment>"*"+"/"        BEGIN(cond_buf);

{EOE}   {
glb_copy_buff();
glbReset();
}


<INITIAL>{
rule    {BEGIN(defr);yylval.name=yytext;return GID;}
channel {BEGIN(defc);yylval.name=yytext;return GID;}
energy  {BEGIN(defe);yylval.name=yytext;return GID;}
cross   {BEGIN(defx);yylval.name=yytext;return GID;}
flux    {BEGIN(deff);yylval.name=yytext;return GID;}
}

include             {cond_buf_i=YY_START;BEGIN(incl);}

<incl>[ \t]*      /* eat the whitespace */
<incl>[^ \t\n]+   { /* got the include file name */
if ( include_stack_ptr >= MAX_INCLUDE_DEPTH )
{
fprintf( stderr, "FATAL: Includes nested too deeply\n" );
exit( 1 );
}

if(include_stack_ptr==0) 
{
file_name_buf[include_stack_ptr]=glb_file_id;
line_buf[include_stack_ptr]=glb_line_num;
}



include_stack[include_stack_ptr++] =
YY_CURRENT_BUFFER;

fnbuf=strdup(yytext);
len=strlen(fnbuf);
if(fnbuf[len-1]!='\"') {fprintf(stderr,"FATAL:" 
" Unterminated filename\n");exit(1) ;}
if(fnbuf[0]!='\"') {fprintf(stderr,"FATAL:" 
" Unquoted filename\n");exit(1) ;}
fnbuf2=(char *) glb_malloc((len-1)*sizeof(char));
for(i=0;i<len-1;i++) fnbuf2[i]=fnbuf[i+1];
fnbuf2[len-2]='\0';
yyin = glb_fopen( fnbuf2, "r" );

if ( ! yyin ){
fprintf(stderr,"FATAL: Could not open file %s\n",fnbuf2);exit(1);}

file_name_buf[include_stack_ptr]=(char*) strdup(fnbuf2);
glb_file_id=file_name_buf[include_stack_ptr];
line_buf[include_stack_ptr]=glb_line_num;
glb_line_num=0;
glb_free(fnbuf);
glb_free(fnbuf2);

yy_switch_to_buffer(
		    yy_create_buffer( yyin, YY_BUF_SIZE ) );


BEGIN(cond_buf_i);
}



<<EOF>> { if ( --include_stack_ptr < 0 )
{
glb_copy_buff();
glbReset();
/* BUGFIX #9 -- This lines fixes a tricky bug in the parser */
include_stack_ptr=0;
yyterminate();
}
else
{
yy_delete_buffer( YY_CURRENT_BUFFER );
yy_switch_to_buffer(
		    include_stack[include_stack_ptr] );
glb_line_num=line_buf[include_stack_ptr];
glb_file_id=file_name_buf[include_stack_ptr];
}
}


<rule,chl,energy,cross,flux>[<]  {
yylval.name=yytext;

return GRPOPEN;}

<chl>@channel          {
yylval.name=yytext;

return CHANNEL;}

<chl>"+" {yylval.in=+1;return PM;}
<chl>"-" {yylval.in=-1;return PM;}

[eE][lectronLECTRON]* {yylval.in=1;return FLAVOR;}
[mM][uonUON]*  {yylval.in=2;return FLAVOR;}
[tT][auAU]* {yylval.in=3;return FLAVOR;}

NOSC_[eE][lectronLECTRON]* {yylval.in=11;return FLAVOR;}
NOSC_[mM][uonUON]*  {yylval.in=12;return FLAVOR;}
NOSC_[tT][auAU]* {yylval.in=13;return FLAVOR;}


<energy>@energy         {
yylval.name=yytext;

return ENERGY;}

<rule>@signal           {
yylval.name=yytext;

return SIGNAL;}

<rule>@background          {
yylval.name=yytext;

return BG;}

<rule>@sys_on_function   {
yylval.name=yytext;
return SYS_ON_FUNCTION;
}

<rule>@sys_off_function   {
yylval.name=yytext;
return SYS_OFF_FUNCTION;
}



\" {cond_buf_str=YY_START;BEGIN(str);}

<str>\" {BEGIN(cond_buf_str);}

<str>\n {glb_fatal("Lexer error: Unterminated string");exit(1);}

<str>{FN} {
yylval.name=yytext;
return FNAME;}

$version {
yylval.name=yytext;
return VERS;
}

<cross>@cross_file           {
yylval.name=yytext;
return CROSS;}

<flux>@flux_file          {
yylval.name=yytext;
return FLUXP;}


<flux>@flux_file_m           {
yylval.name=yytext;
return FLUXM;}


:        {return RULESEP;}
@        {return RULEMULT;} 



<defr>#{ID}        {glb_namerec *s;
s = glb_getname (yytext,"rule");
if (s == 0) 
s = glb_putname (yytext,"rule",NAME); 
BEGIN(INITIAL);
BEGIN(rule); 
if(strcmp(s->context,"rule")!=0) 
{yylval.nameptr=s;return RDF;}

yylval.nameptr = s; 

return s->type;
}

<defc>#{ID}        {glb_namerec *s;
s = glb_getname (yytext,"channel");
if (s == 0) 
s = glb_putname (yytext,"channel",NAME);
BEGIN(INITIAL);
BEGIN(chl); 
if(strcmp(s->context,"channel")!=0) 
{yylval.nameptr=s;return RDF;}
 
yylval.nameptr = s;  

return s->type;
}

<defe>#{ID}        {glb_namerec *s;
s = glb_getname (yytext,"energy");
if (s == 0) 
s = glb_putname (yytext,"energy",NAME);
BEGIN(INITIAL);
BEGIN(energy);
if(strcmp(s->context,"energy")!=0) 
{yylval.nameptr=s;return RDF;}
 
yylval.nameptr = s; 
 
return s->type;
}

<defx>#{ID}        {glb_namerec *s;
s = glb_getname (yytext,"cross");
if (s == 0) 
s = glb_putname (yytext,"cross",NAME);
BEGIN(INITIAL);
BEGIN(cross);
if(strcmp(s->context,"cross")!=0) 
{yylval.nameptr=s;return RDF;}
 
yylval.nameptr = s; 
 
return s->type;
}


<deff>#{ID}        {glb_namerec *s;
s = glb_getname (yytext,"flux");
if (s == 0) 
s = glb_putname (yytext,"flux",NAME);
BEGIN(INITIAL);
BEGIN(flux);
if(strcmp(s->context,"flux")!=0) 
{yylval.nameptr=s;return RDF;}
 
yylval.nameptr = s; 
 
return s->type;
}

#{ID}        {glb_namerec *s;
s = glb_getname (yytext,"rule");
if (s == 0) return NDEF;
yylval.nameptr = s; 
return s->type; 
}

<rule,chl,energy,cross,flux>[>]        {
BEGIN(INITIAL);
return GRPCLOSE;}


<chl,rule,energy,cross,flux>@{ID} {
yylval.name=yytext;
return IDN;
}




{DIGIT}+    {yylval.val=atof(yytext); 
return NUM;
}

{IVAL}  {yylval.name=yytext;return IDN;}

{DIGIT}+"."{DIGIT}*        {
yylval.val=atof(yytext); 
return NUM;
}

 "."{DIGIT}+       {
 yylval.val=atof(yytext); 
 return NUM;
}

{DIGIT}+"."{DIGIT}*[eE][-+]?{DIGIT}+     {
yylval.val=atof(yytext); 
return NUM;
}

"."{DIGIT}+[eE][-+]?{DIGIT}+     {
yylval.val=atof(yytext); 
return NUM;
}

{DIGIT}+[eE][-+]?{DIGIT}+     {
yylval.val=atof(yytext); 
return NUM;
}




{ID}        {glb_symrec *s;
s = glb_getsym (yytext);
if (s == 0) 
s = glb_putsym (yytext, VAR); 
yylval.tptr = s; 
return s->type; 
}


{LID}        {glb_symrec *s;
s = glb_getsym (yytext);
if (s == 0) 
s = glb_putsym (yytext, LVAR); 
yylval.tptr = s; 
return s->type; 
}





[ \t]*          /* eat up whitespace */

"%!GLoBES" /* empty */

.           return (int) yytext[0];


\n         /* yylloc.first_line++;*/ glb_line_num++; return '\n';
%%


/* OLD CODE, DOES NOT RUN WITH NEWER FLEX VERSIONS ANYMORE [WW, JK]:

void glb_lexer_cleanup() {
     if( yy_current_buffer != NULL )
     {
         yy_delete_buffer( yy_current_buffer );
         yy_current_buffer = 0;
     }
}

*/


void glb_lexer_cleanup() {
      if( YY_CURRENT_BUFFER != NULL )
          yy_delete_buffer( YY_CURRENT_BUFFER );
}


