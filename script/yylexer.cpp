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
	else
	{
		switch(tok.type)
		{
			case LEX_TOKEN_IF: 	return TOK_IF;
			case LEX_TOKEN_ELSE:	return TOK_ELSE;
			case LEX_TOKEN_FOR:	return TOK_FOR;
			case LEX_TOKEN_WHILE:	return TOK_WHILE;
			case LEX_TOKEN_RETURN:	return TOK_RETURN;
			case LEX_TOKEN_BREAK:	return TOK_BREAK;
			case LEX_TOKEN_CONTINUE:	return TOK_CONTINUE;
			case LEX_TOKEN_LOG_AND:	return TOK_LOG_AND;
			case LEX_TOKEN_LOG_OR:	return TOK_LOG_OR;
			case LEX_TOKEN_LOG_NOT:	return TOK_LOG_NOT;
			case LEX_TOKEN_LOG_EQ:	return TOK_LOG_EQ;
			case LEX_TOKEN_LOG_NEQ:	return TOK_LOG_NEQ;
			case LEX_TOKEN_LOG_LESS:	return TOK_LOG_LESS;
			case LEX_TOKEN_LOG_GREATER:	return TOK_LOG_GREATER;
			case LEX_TOKEN_LOG_LEQ:	return TOK_LOG_LEQ;
			case LEX_TOKEN_LOG_GEQ:	return TOK_LOG_GEQ;
			case LEX_TOKEN_GLOBAL:	return TOK_GLOBAL;
		}
	}

	std::cerr << "Error: Invalid token: " << tok.type << std::endl;
	return 0;
}


extern "C" void yyerror(const char* pc)
{
	std::cerr << "Error: " << pc << std::endl;
}
