/*
 * Simple Script
 * @author tweber
 */


%{
	#define YYDEBUG 1
	#define YYERROR_VERBOSE 1
	#define YYTOKEN_TABLE 1

	#include "parseobj.h"
	
	#define YYLEX_PARAM pParseObj
	#define YYPARSE_PARAM pParseObj

	#include <iostream>
	#include <math.h>
	#include "yylexer.h"

	//#define YYSTYPE Node*
%}

//%glr-parser

%union
{
	void *pNode;
}


%token<pNode> TOK_DOUBLE
%token<pNode> TOK_INT
%token<pNode> TOK_STRING
%token<pNode> TOK_IDENT
%type<pNode> prog
%type<pNode> funcs
%type<pNode> func
%type<pNode> stmts
%type<pNode> stmt
%type<pNode> expr
%type<pNode> args
%type<pNode> idents
%type<pNode> ident


%right '='
%left '+' '-'
%left '*' '/'
%left UMINUS UPLUS
%right '^'


%%

prog:	funcs			{ $$ = $1; ((ParseObj*)pParseObj)->pRoot = (Node*)$1; }
	;
	


funcs: func funcs		{ $$ = new NodeBinaryOp($1, $2, NODE_FUNCS); }
	| /* eps */			{ $$ = 0; }
	;

func: ident '(' idents ')' '{' stmts '}'		{ $$ = new NodeFunction($1, $3, $6); }
	; 




stmts:	stmt stmts		{ $$ = new NodeBinaryOp($1, $2, NODE_STMTS); }
	| '{' stmts '}'		{ $$ = new NodeUnaryOp($2, NODE_STMTS); }
	| /*eps*/			{ $$ = 0; } 
	;
	
stmt:	expr ';'		{ $$ = $1; }
	;
	



args:	expr ',' args		{ $$ = new NodeBinaryOp($1, $3, NODE_ARGS); }
	|   expr				{ $$ = $1; }
	| 	/*eps*/				{ $$ = 0; }
	;

idents:	ident ',' idents	{ $$ = new NodeBinaryOp($1, $3, NODE_IDENTS); }
	|   ident				{ $$ = $1; }
	| 	/*eps*/				{ $$ = 0; }
	;



expr:	'(' expr ')'		{ $$ = $2; }
	| ident '=' expr		{ $$ = new NodeBinaryOp($1, $3, NODE_ASSIGN); }
	| expr '+' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_PLUS); }
	| expr '-' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_MINUS); }
	| expr '*' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_MULT); }
	| expr '/' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_DIV); }
	| expr '^' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_POW); }
	| '-' expr %prec UMINUS	{ $$ = new NodeUnaryOp($2, NODE_UMINUS); }
	| '+' expr %prec UPLUS	{ $$ = $2; }
	| TOK_DOUBLE			{ $$ = $1; }
	| TOK_INT				{ $$ = $1; }
	| TOK_STRING			{ $$ = $1; }
	| ident					{ $$ = $1; }
	| ident '(' args ')'	{ $$ = new NodeCall($1, $3);  }
	;

	
ident: TOK_IDENT			{ $$ = $1; }

%%
