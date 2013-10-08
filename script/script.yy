/*
 * Simple Script grammar
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

%token<pNode> TOK_IF
%token<pNode> TOK_ELSE
%token<pNode> TOK_FOR
%token<pNode> TOK_WHILE
%token<pNode> TOK_RETURN
%token<pNode> TOK_BREAK
%token<pNode> TOK_CONTINUE

%token<pNode> TOK_LOG_AND
%token<pNode> TOK_LOG_OR
%token<pNode> TOK_LOG_NOT
%token<pNode> TOK_LOG_EQ
%token<pNode> TOK_LOG_NEQ
%token<pNode> TOK_LOG_LESS
%token<pNode> TOK_LOG_GREATER
%token<pNode> TOK_LOG_LEQ
%token<pNode> TOK_LOG_GEQ


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
%left TOK_LOG_OR
%left TOK_LOG_AND
%left TOK_LOG_EQ
%left TOK_LOG_NEQ
%left TOK_LOG_LESS
%left TOK_LOG_LEQ
%left TOK_LOG_GREATER
%left TOK_LOG_GEQ
%left '+' '-'
%left '*' '/'
%right TOK_LOG_NOT
%right UMINUS UPLUS
%right '^'


%%

prog:	funcs			{ $$ = $1; ((ParseObj*)pParseObj)->pRoot = (Node*)$1; }
	;
	


funcs: func funcs		{ $$ = new NodeBinaryOp($1, $2, NODE_FUNCS); }
	| /* eps */			{ $$ = 0; }
	;

func: ident '(' idents ')' stmt 		{ $$ = new NodeFunction($1, $3, $5); }
	; 




stmts:	stmt stmts		{ $$ = new NodeBinaryOp($1, $2, NODE_STMTS); }
	| /*eps*/			{ $$ = 0; } 
	;
	
stmt:	expr ';'		{ $$ = $1; }
	| 	'{' stmts '}'	{ $$ = $2; }
	
	|	TOK_IF '(' expr ')' stmt 					{ $$ = new NodeIf($3, $5); }
	|	TOK_IF '(' expr ')' stmt TOK_ELSE stmt 		{ $$ = new NodeIf($3, $5, $7); }
	|	TOK_WHILE '(' expr ')' stmt					{ $$ = new NodeWhile($3, $5); }
	|	TOK_FOR '(' ident ':' expr ')' stmt			{ $$ = new NodeRangedFor($3, $5, $7); }
	|	TOK_RETURN expr ';'							{ $$ = new NodeReturn($2); }
	|	TOK_RETURN ';'								{ $$ = new NodeReturn(); }
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
	| expr '=' expr		{ $$ = new NodeBinaryOp($1, $3, NODE_ASSIGN); }

	| expr '+' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_PLUS); }
	| expr '-' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_MINUS); }
	| expr '*' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_MULT); }
	| expr '/' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_DIV); }
	| expr '^' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_POW); }
	| '-' expr %prec UMINUS	{ $$ = new NodeUnaryOp($2, NODE_UMINUS); }
	| '+' expr %prec UPLUS	{ $$ = $2; }
	
	| expr TOK_LOG_AND expr		{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_AND); }
	| expr TOK_LOG_OR expr 		{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_OR); }
	| expr TOK_LOG_EQ expr 		{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_EQ); }
	| expr TOK_LOG_NEQ expr		{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_NEQ); }
	| expr TOK_LOG_LESS expr 	{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_LESS); }
	| expr TOK_LOG_GREATER expr	{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_GREATER); }
	| expr TOK_LOG_LEQ expr 	{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_LEQ); }
	| expr TOK_LOG_GEQ expr		{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_GEQ); }
	| TOK_LOG_NOT expr			{ $$ = new NodeUnaryOp($2, NODE_LOG_NOT); }		
		
	| TOK_DOUBLE			{ $$ = $1; }
	| TOK_INT				{ $$ = $1; }
	| TOK_STRING			{ $$ = $1; }
	| ident					{ $$ = $1; }
	| ident '(' args ')'	{ $$ = new NodeCall($1, $3);  }
	
	| '[' args ']'				{ $$ = new NodeArray($2); }
	| ident '[' expr ']'		{ $$ = new NodeArrayAccess($1, $3); }
	;

	
ident: TOK_IDENT			{ $$ = $1; }

%%
