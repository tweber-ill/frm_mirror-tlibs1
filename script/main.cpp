/*
 * Simple Script
 * @author tweber
 * @date 2013
 */

#include "types.h"
#include "helper/flags.h"
#include "helper/string.h"
#include "helper/spec_char.h"
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
		G_COUT << "This is the " << g_pcVersion << "." << std::endl;
		G_COUT << "\tUsage: " << argv[0] << " <script file> [arguments]" << std::endl;
		return -1;
	}

	const char* pcFile = argv[1];
	t_string strFile = STR_TO_WSTR(pcFile);

	t_char* pcInput = load_file(pcFile);
	if(!pcInput)
		return -2;

	ParseObj par;
	ParseInfo info;

	par.strCurFile = strFile;
	par.pLexer = new Lexer(pcInput, strFile.c_str());

	delete[] pcInput;
	pcInput = 0;

	if(!par.pLexer->IsOk())
	{
		G_CERR << "Error: Lexer returned with errors." << std::endl;
		return -3;
	}

	init_global_syms(info.pGlobalSyms);

	//yydebug = 0;
	int iParseRet = yyparse(&par);

	delete par.pLexer;
	par.pLexer = 0;

	if(iParseRet != 0)
	{
		G_CERR << "Error: Parser returned with error code " << iParseRet << "." << std::endl;
		return -4;
	}

	SymbolArray *parrMainArgs = new SymbolArray();
	for(int iArg=1; iArg<argc; ++iArg)
	{
		SymbolString *pSymArg = new SymbolString();
		pSymArg->SetVal(STR_TO_WSTR(argv[iArg]));
		parrMainArgs->GetArr().push_back(pSymArg);
	}
	//std::vector<Symbol*> vecMainArgs = { &arrMainArgs };

	SymbolTable *pTableSup = new SymbolTable();

	info.pmapModules->insert(ParseInfo::t_mods::value_type(strFile, par.pRoot));
	info.strExecFkt = T_STR"main";
	//info.pvecExecArg = &vecMainArgs;
	info.strInitScrFile = strFile;

	SymbolArray arrMainArgs;
	arrMainArgs.GetArr().push_back(parrMainArgs);
	pTableSup->InsertSymbol(T_STR"<args>", &arrMainArgs);
	par.pRoot->eval(info, pTableSup);
	pTableSup->RemoveSymbolNoDelete(T_STR"<args>");

	//info.pGlobalSyms->print();

	delete pTableSup;

	return 0;
}


int main(int argc, char** argv)
{
	int iRet = -99;

	try
	{
		init_spec_chars();
		iRet = script_main(argc, argv);
	}
	catch(const std::exception& ex)
	{
		G_CERR << "Critical error in script interpreter: " << ex.what() << std::endl;
	}

	return iRet;
}
