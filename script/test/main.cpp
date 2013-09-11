/*
 * Simple Script
 * @author tweber
 */
 
#include <iostream>
#include <fstream>
#include "parseobj.h"


extern int yydebug;
int yyparse(void* pLexer);


int main(int argc, char** argv)
{
	if(argc<=1)
	{
		std::cout << "Usage: " << argv[0] << " <file>" << std::endl;
		return -1;
	}

	std::ifstream ifstr(argv[1]);
	if(!ifstr.is_open())
	{
		std::cerr << "Error: Cannot open \"" << argv[1] << "\"." << std::endl;
		return -2;
	}

	ifstr.seekg(0, std::ios::end);
	std::size_t iFileLen = ifstr.tellg();
	ifstr.seekg(0, std::ios::beg);

	char* pcInput = new char[iFileLen+1];
	ifstr.read(pcInput, iFileLen);
	pcInput[iFileLen] = 0;


	ParseObj par;
	par.pLexer = new Lexer(pcInput);
	par.pSym = new SymbolTable();

	yydebug = 0;
	yyparse(&par);

	par.pRoot->eval(par.pSym);

	par.pSym->print();


	delete[] pcInput;

	delete par.pLexer;
	delete par.pSym;
	delete par.pRoot;

	return 0;
}
