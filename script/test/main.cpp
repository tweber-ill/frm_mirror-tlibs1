/*
 * Simple Script
 * @author tweber
 */
 
#include <iostream>
#include "parseobj.h"


extern int yydebug;
int yyparse(void* pLexer);


int main(int argc, char** argv)
{
	ParseObj par;
	par.pLexer = new Lexer("a = -(2+4)*2;\nb = (75+25)^2;\nc=a+b;");
	par.pSym = new SymbolTable();

	yydebug = 0;
	yyparse(&par);

	par.pRoot->eval(par.pSym);

	par.pSym->print();

	delete par.pLexer;
	delete par.pSym;
	delete par.pRoot;
	return 0;
}
