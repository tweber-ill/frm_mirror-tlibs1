/*
 * Simple Script
 * @author tweber
 * @date 2013
 */

#include "flags.h" 
#include <iostream>
#include <fstream>
#include <cmath>

#include "parseobj.h"
#include "script_helper.h"

#include "globals.h"
#include "calls_plot.h"
#include "calls_math.h"
#include "calls_file.h"
#include "calls_thread.h"

const char* g_pcVersion = "Hermelin script interpreter, version 0.4";

//extern int yydebug;
extern int yyparse(void*);

void init_all_externals(SymbolTable* pGlobals)
{
	init_global_syms(pGlobals);

	init_ext_thread_calls();
	init_ext_file_calls();
	init_ext_math_calls();
	init_ext_plot_calls();
}

int main(int argc, char** argv)
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

	init_all_externals(info.pGlobalSyms);

	//yydebug = 0;
	int iParseRet = yyparse(&par);

	delete par.pLexer;
	par.pLexer = 0;

	if(iParseRet != 0)
	{
		std::cerr << "Error: Parser returned with error code " << iParseRet << "." << std::endl;
		return -4;
	}

	SymbolArray arrMainArgs;
	for(int iArg=1; iArg<argc; ++iArg)
	{
		SymbolString *pSymArg = new SymbolString();
		pSymArg->m_strVal = argv[iArg];
		arrMainArgs.m_arr.push_back(pSymArg);
	}
	std::vector<Symbol*> vecMainArgs = { &arrMainArgs };

	info.pmapModules->insert(ParseInfo::t_mods::value_type(pcFile, par.pRoot));
	info.strExecFkt = "main";
	info.pvecExecArg = &vecMainArgs;
	info.strInitScrFile = pcFile;
	par.pRoot->eval(info);
	//info.pGlobalSyms->print();

	return 0;
}
