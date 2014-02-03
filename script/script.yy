/*
 * Simple Script grammar
 * @author tweber
 * @date 2013
 */


%{
	#include "parseobj.h"
	
	#define YYLEX_PARAM pParseObj
	#define YYPARSE_PARAM pParseObj
	//#define YYLTYPE int

	#include <iostream>
	#include <math.h>
	#include "yylexer.h"

	//#define YYSTYPE Node*
	
	static inline void set_linenr(void* dollardollar, void* pParseObj)
	{
		if(dollardollar && pParseObj)
			((Node*)dollardollar)->m_iLine = ((ParseObj*)pParseObj)->iCurLine;
	}
%}

%language "C"
%pure-parser
//%define api.pure full
//%glr-parser
%token-table

//%verbose
%error-verbose
//%define parse.error verbose
//%debug

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
%token<pNode> TOK_GLOBAL

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
%type<pNode> pair
%type<pNode> ranged_expr
%type<pNode> arr_args
%type<pNode> map_args
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
%left '*' '/' '%'
%right TOK_LOG_NOT
%right UMINUS UPLUS
%right '^'
%right '['

%right TOK_IF TOK_ELSE

%%

prog:	funcs			{ $$ = $1; ((ParseObj*)pParseObj)->pRoot = (Node*)$1; set_linenr($$, pParseObj); }
	;
	


funcs: func funcs		{ $$ = new NodeBinaryOp($1, $2, NODE_FUNCS); set_linenr($$, pParseObj); }
	| /* eps */			{ $$ = 0; }
	;

func: ident '(' idents ')' stmt 		{ $$ = new NodeFunction($1, $3, $5); set_linenr($$, pParseObj); }
	; 




stmts:	stmt stmts		{ $$ = new NodeBinaryOp($1, $2, NODE_STMTS); set_linenr($$, pParseObj); }
	| /*eps*/			{ $$ = 0; } 
	;
	
stmt:	expr ';'		{ $$ = $1; }
	| 	'{' stmts '}'	{ $$ = $2; }
	
	|	TOK_IF '(' expr ')' stmt 								%prec TOK_IF 			{ $$ = new NodeIf($3, $5); set_linenr($$, pParseObj); }
	|	TOK_IF '(' expr ')' stmt TOK_ELSE stmt 	%prec TOK_ELSE 	{ $$ = new NodeIf($3, $5, $7); set_linenr($$, pParseObj); }
	|	TOK_WHILE '(' expr ')' stmt							{ $$ = new NodeWhile($3, $5); set_linenr($$, pParseObj); }
	|	TOK_FOR '(' ident ':' expr ')' stmt					{ $$ = new NodeRangedFor($3, $5, $7); set_linenr($$, pParseObj); }
	|	TOK_RETURN expr ';'										{ $$ = new NodeReturn($2); set_linenr($$, pParseObj); }
	|	TOK_RETURN ';'												{ $$ = new NodeReturn(); set_linenr($$, pParseObj); }
	|	TOK_BREAK ';'												{ $$ = new NodeBreak(); set_linenr($$, pParseObj); }
	|	TOK_CONTINUE ';'											{ $$ = new NodeContinue(); set_linenr($$, pParseObj); }
	;


args:	expr ',' args		{ $$ = new NodeBinaryOp($1, $3, NODE_ARGS); set_linenr($$, pParseObj); }
	|	expr				{ $$ = $1; }
	|	/*eps*/				{ $$ = 0; }
	;

ranged_expr: expr				{ $$ = $1; }
	|	':'						{ $$ = new NodeRange(RANGE_FULL); set_linenr($$, pParseObj);}
	|	expr ':' expr			{ $$ = new NodeRange($1, $3); set_linenr($$, pParseObj);}
	;

arr_args: ranged_expr ',' arr_args	{ $$ = new NodeBinaryOp($1, $3, NODE_ARGS); set_linenr($$, pParseObj); }
	|	ranged_expr					{ $$ = $1; }
	|	/*eps*/						{ $$ = 0; }
	;

pair: expr ':' expr		{ $$ = new NodePair($1, $3); set_linenr($$, pParseObj); }
	;

map_args: pair ',' map_args	{ $$ = new NodeBinaryOp($1, $3, NODE_ARGS); set_linenr($$, pParseObj); }
	|	pair			{ $$ = $1; }
//	|	/*eps*/		{ $$ = 0; }
	;

idents:	ident ',' idents	{ $$ = new NodeBinaryOp($1, $3, NODE_IDENTS); set_linenr($$, pParseObj); }
	|   ident				{ $$ = $1; } 
	| 	/*eps*/				{ $$ = 0; }
	;



expr:	'(' expr ')'		{ $$ = $2; }
	| TOK_GLOBAL expr '=' expr	{ $$ = new NodeBinaryOp($2, $4, NODE_ASSIGN); 
												((NodeBinaryOp*)$$)->m_bGlobal = 1; 
												set_linenr($$, pParseObj); }
	| expr '=' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_ASSIGN); set_linenr($$, pParseObj); }

	| expr '+' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_PLUS); set_linenr($$, pParseObj); }
	| expr '-' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_MINUS); set_linenr($$, pParseObj); }
	| expr '*' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_MULT); set_linenr($$, pParseObj); }
	| expr '/' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_DIV); set_linenr($$, pParseObj); }
	| expr '%' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_MOD); set_linenr($$, pParseObj); }
	| expr '^' expr			{ $$ = new NodeBinaryOp($1, $3, NODE_POW); set_linenr($$, pParseObj); }
	| '-' expr %prec UMINUS		{ $$ = new NodeUnaryOp($2, NODE_UMINUS); set_linenr($$, pParseObj); }
	| '+' expr %prec UPLUS		{ $$ = $2; }
	
	| expr TOK_LOG_AND expr		{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_AND); set_linenr($$, pParseObj); }
	| expr TOK_LOG_OR expr 		{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_OR); set_linenr($$, pParseObj); }
	| expr TOK_LOG_EQ expr 		{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_EQ); set_linenr($$, pParseObj); }
	| expr TOK_LOG_NEQ expr		{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_NEQ); set_linenr($$, pParseObj); }
	| expr TOK_LOG_LESS expr 	{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_LESS); set_linenr($$, pParseObj); }
	| expr TOK_LOG_GREATER expr	{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_GREATER); set_linenr($$, pParseObj); }
	| expr TOK_LOG_LEQ expr 	{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_LEQ); set_linenr($$, pParseObj); }
	| expr TOK_LOG_GEQ expr		{ $$ = new NodeBinaryOp($1, $3, NODE_LOG_GEQ); set_linenr($$, pParseObj); }
	| TOK_LOG_NOT expr			{ $$ = new NodeUnaryOp($2, NODE_LOG_NOT); set_linenr($$, pParseObj); }
	
	| '*' expr					{ $$ = new NodeUnaryOp($2, NODE_UNPACK); set_linenr($$, pParseObj); }
		
	| TOK_DOUBLE			{ $$ = $1; }
	| TOK_INT				{ $$ = $1; }
	| TOK_STRING			{ $$ = $1; }
	| ident					{ $$ = $1; }
	| ident '(' args ')'	{ $$ = new NodeCall($1, $3);  set_linenr($$, pParseObj); }
	
	| '[' args ']'			{ $$ = new NodeArray($2); set_linenr($$, pParseObj); }
	| '[' map_args ']'		{ $$ = new NodeMap($2); set_linenr($$, pParseObj); }
	| expr '[' arr_args ']'	{ $$ = new NodeArrayAccess($1, $3); set_linenr($$, pParseObj); }
	;

	
ident: TOK_IDENT			{ $$ = $1; }

%%
