/*
 * Simple Script
 * @author tweber
 * @date 2013
 */

#ifndef __MIEZE_SCRIPT_GLOBAL__
#define __MIEZE_SCRIPT_GLOBAL__

#include "types.h"

#include "lexer.h"
#include "node.h"
#include "symbol.h"
#include <string>

struct ParseObj
{
	Lexer* pLexer;
	Node* pRoot;

	// only used during parsing/lexing for yyerror(), NOT during exec
	unsigned int iCurLine;
	t_string strCurFile;
};

#endif
