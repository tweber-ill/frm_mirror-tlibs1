/*
 * Simple Script
 * @author tweber
 */
 
#include <iostream>
#include <fstream>
#include "parseobj.h"


extern int yydebug;
int yyparse(void*);


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
	delete[] pcInput;
	pcInput = 0;

	if(!par.pLexer->IsOk())
	{
		std::cerr << "Error: Lexer returned with errors." << std::endl;
		return -3;
	}


	info.pGlobalSyms = new SymbolTable();

	yydebug = 0;
	int iParseRet = yyparse(&par);

	if(iParseRet != 0)
	{
		std::cerr << "Error: Parser returned with error code " << iParseRet << "." << std::endl;
		return -4;
	}

	par.pRoot->eval(info);
	//info.pGlobalSyms->print();


	delete par.pLexer;
	delete info.pGlobalSyms;
	delete par.pRoot;

	return 0;
}
