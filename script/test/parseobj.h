/*
 * Simple Script
 * @author tweber
 */

#ifndef __MIEZE_SCRIPT_GLOBAL__
#define __MIEZE_SCRIPT_GLOBAL__

#include "lexer.h"
#include "node.h"
#include "symbol.h"

struct ParseObj
{
	Lexer* pLexer;
	SymbolTable* pSym;
	Node* pRoot;
};

#endif
