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
	else if(tok.type == LEX_TOKEN_INT)
	{
		Node *pNode = new NodeInt(tok.iVal);
		yylval.pNode = pNode;
		return TOK_INT;
	}
	else if(tok.type == LEX_TOKEN_STRING)
	{
		Node *pNode = new NodeString(tok.strVal);
		yylval.pNode = pNode;
		return TOK_STRING;
	}
	else if(tok.type == LEX_TOKEN_CHAROP)
		return tok.cOp;
	else if(tok.type == LEX_TOKEN_IDENT)
	{
		Node *pNode = new NodeIdent(tok.strVal);
		yylval.pNode = pNode;
		return TOK_IDENT;
	}

	else if(tok.type == LEX_TOKEN_IF)
		return TOK_IF;
	else if(tok.type == LEX_TOKEN_ELSE)
		return TOK_ELSE;
	else if(tok.type == LEX_TOKEN_FOR)
		return TOK_FOR;
	else if(tok.type == LEX_TOKEN_WHILE)
		return TOK_WHILE;

	else if(tok.type == LEX_TOKEN_RETURN)
		return TOK_RETURN;
	else if(tok.type == LEX_TOKEN_BREAK)
		return TOK_BREAK;
	else if(tok.type == LEX_TOKEN_CONTINUE)
		return TOK_CONTINUE;

	else if(tok.type == LEX_TOKEN_LOG_AND)
		return TOK_LOG_AND;
	else if(tok.type == LEX_TOKEN_LOG_OR)
		return TOK_LOG_OR;
	else if(tok.type == LEX_TOKEN_LOG_NOT)
		return TOK_LOG_NOT;
	else if(tok.type == LEX_TOKEN_LOG_EQ)
		return TOK_LOG_EQ;
	else if(tok.type == LEX_TOKEN_LOG_NEQ)
		return TOK_LOG_NEQ;
	else if(tok.type == LEX_TOKEN_LOG_LESS)
		return TOK_LOG_LESS;
	else if(tok.type == LEX_TOKEN_LOG_GREATER)
		return TOK_LOG_GREATER;
	else if(tok.type == LEX_TOKEN_LOG_LEQ)
		return TOK_LOG_LEQ;
	else if(tok.type == LEX_TOKEN_LOG_GEQ)
		return TOK_LOG_GEQ;

	return 0;
}


extern "C" void yyerror(const char* pc)
{
	std::cerr << "Error: " << pc << std::endl;
}
