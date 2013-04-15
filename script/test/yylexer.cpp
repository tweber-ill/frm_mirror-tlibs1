/*
 * Simple Script
 * @author tweber
 */

#include <iostream>
#include <string.h>

#include "tokens.h"
#include "yylexer.h"
#include "parseobj.h"

extern "C" int yylex(void* _pObj)
{
	Lexer* pLexer = ((ParseObj*)_pObj)->pLexer;
	//pLexer->print();
	
	const Token& tok = pLexer->lex();
	if(tok.type == LEX_TOKEN_END)
		return 0;
	else if(tok.type == LEX_TOKEN_INVALID)
		return 1;
	else if(tok.type == LEX_TOKEN_DOUBLE)
	{
		Node *pNode = new NodeDouble(tok.dVal);
		yylval.pNode = pNode;
		return TOK_DOUBLE;
	}
	else if(tok.type == LEX_TOKEN_CHAROP)
		return tok.cOp;
	else if(tok.type == LEX_TOKEN_IDENT)
	{
		Node *pNode = new NodeIdent(tok.strVal);
		yylval.pNode = pNode;
		return TOK_IDENT;
	}
	else if(tok.type == LEX_TOKEN_STRING)
	{
		// TODO
	}
}


extern "C" void yyerror(const char* pc)
{
	std::cerr << "Error: " << pc << std::endl;
}
