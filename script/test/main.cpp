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
		std::cout << "This is the Hermelin script interpreter." << std::endl;
		std::cout << "\tUsage: " << argv[0] << " <script file>" << std::endl;
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
	ParseInfo info;
	par.pLexer = new Lexer(pcInput);
	info.pGlobalSyms = new SymbolTable();

	yydebug = 0;
	yyparse(&par);

	par.pRoot->eval(info);
	//info.pGlobalSyms->print();


	delete[] pcInput;

	delete par.pLexer;
	delete info.pGlobalSyms;
	delete par.pRoot;

	return 0;
}
