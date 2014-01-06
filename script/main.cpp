/*
 * Simple Script
 * @author tweber
 */

#include "flags.h" 
#include <iostream>
#include <fstream>
#include <cmath>
#include "parseobj.h"
#include "script_helper.h"
#include "helper/neutrons.hpp"

const char* g_pcVersion = "Hermelin script interpreter, version 0.4";

//extern int yydebug;
extern int yyparse(void*);

static void init_global_syms(SymbolTable *pSymTab)
{
	pSymTab->InsertSymbol("pi", new SymbolDouble(M_PI));

	// hbar in eVs
	pSymTab->InsertSymbol("hbar_eVs", new SymbolDouble(co::hbar / one_eV / units::si::second));

	// hbar in Js
	pSymTab->InsertSymbol("hbar", new SymbolDouble(co::hbar / units::si::joule / units::si::second));

	// neutron mass
	pSymTab->InsertSymbol("m_n", new SymbolDouble(co::m_n / units::si::kilogram));
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


	info.pGlobalSyms = new SymbolTable();
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
