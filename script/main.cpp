/*
 * Simple Script
 * @author tweber
 * @date 2013
 */

#include "helper/flags.h"
#include <iostream>
#include <fstream>
#include <cmath>

#include "parseobj.h"
#include "script_helper.h"

#include "globals.h"


//extern int yydebug;
extern int yyparse(void*);

static inline int script_main(int argc, char** argv)
{
	if(argc<=1)
	{
		std::cout << "This is the " << g_pcVersion << "." << std::endl;
		std::cout << "\tUsage: " << argv[0] << " <script file> [arguments]" << std::endl;
		return -1;
	}

	const char* pcFile = argv[1];
	char* pcInput = load_file(pcFile);
	if(!pcInput)
		return -2;

	ParseObj par;
	ParseInfo info;

	par.strCurFile = pcFile;
	par.pLexer = new Lexer(pcInput, pcFile);

	delete[] pcInput;
	pcInput = 0;

	if(!par.pLexer->IsOk())
	{
		std::cerr << "Error: Lexer returned with errors." << std::endl;
		return -3;
	}

	init_global_syms(info.pGlobalSyms);

	//yydebug = 0;
	int iParseRet = yyparse(&par);

	delete par.pLexer;
	par.pLexer = 0;

	if(iParseRet != 0)
	{
		std::cerr << "Error: Parser returned with error code " << iParseRet << "." << std::endl;
		return -4;
	}

	SymbolArray *parrMainArgs = new SymbolArray();
	for(int iArg=1; iArg<argc; ++iArg)
	{
		SymbolString *pSymArg = new SymbolString();
		pSymArg->m_strVal = argv[iArg];
		parrMainArgs->m_arr.push_back(pSymArg);
	}
	//std::vector<Symbol*> vecMainArgs = { &arrMainArgs };

	SymbolTable *pTableSup = new SymbolTable();

	info.pmapModules->insert(ParseInfo::t_mods::value_type(pcFile, par.pRoot));
	info.strExecFkt = "main";
	//info.pvecExecArg = &vecMainArgs;
	info.strInitScrFile = pcFile;

	SymbolArray arrMainArgs;
	arrMainArgs.m_arr.push_back(parrMainArgs);
	pTableSup->InsertSymbol("<args>", &arrMainArgs);
	par.pRoot->eval(info, pTableSup);
	pTableSup->RemoveSymbolNoDelete("<args>");

	//info.pGlobalSyms->print();

	delete pTableSup;

	return 0;
}


int main(int argc, char** argv)
{
	int iRet = -99;

	try
	{
		iRet = script_main(argc, argv);
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Critical error in script interpreter: " << ex.what() << std::endl;
	}

	return iRet;
}
