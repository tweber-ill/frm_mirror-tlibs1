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
#include "calls.h"

extern int yyparse(void*);
static bool g_bShowTiming = 0;

static inline int script_main(int argc, char** argv)
{
	if(argc<=1)
	{
		G_COUT << "This is the " << g_pcVersion << "." << "\n";
		G_COUT << "Built on " << __DATE__ << ", " << __TIME__;
		G_COUT << " with CC version " << __VERSION__ << ".\n\n";
		G_COUT << "Usage: " << argv[0] << " [arguments to hermelin]" 
					<< " <script file> [arguments to script]" 
					<< "\n";
		G_COUT << "\nArguments to hermelin:" << "\n";
		G_COUT << "\t-t, --timing\t\tShow timing information.\n";
		G_COUT << "\t-s, --symbols\t\tShow symbol tables.\n";
		//G_COUT << "\t-d, --debug\t\tEnable debug output.\n";
		G_COUT << std::endl;
		return -1;
	}

	bool bShowSymbols = 0;
	unsigned int iStartArg = 1;
	for(iStartArg=1; iStartArg<argc; ++iStartArg)
	{
		t_string strArg = STR_TO_WSTR(argv[iStartArg]);
		trim(strArg);

		// end of arguments to hermelin
		if(strArg[0] != T_STR'-')
			break;

		if(strArg=="-t" || strArg == "--timing")
			g_bShowTiming = 1;
		else if(strArg=="-s" || strArg == "--symbols")
			bShowSymbols = 1;
		else if(strArg=="-d" || strArg == "--debug")
			g_bDebug = 1;
	}

	// debug in script.yy needs to be set
	yydebug = g_bDebug;

	if(iStartArg >= argc)
	{
		G_CERR << "Error: No input file given." << std::endl;
		return -1;
	}

	const char* pcFile = argv[iStartArg];
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
	for(int iArg=iStartArg; iArg<argc; ++iArg)
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

	if(bShowSymbols)
	{
		G_COUT << "\n";
		G_COUT << "================================================================================\n";	
		G_COUT << "Global symbols:\n";
		info.pGlobalSyms->print();


		G_COUT << "\n";
		G_COUT << "Script functions:\n";
		for(const NodeFunction* pFunc : info.vecFuncs)
			G_COUT << pFunc->GetName() << ", ";
		G_COUT << "\n";


		const t_mapFkts* pExtFkts = get_ext_calls();

		G_COUT << "\n";
		G_COUT << "System functions:\n";
		for(const auto& fktpair : *pExtFkts)
			G_COUT << fktpair.first << ", ";
		G_COUT << "\n";

		G_COUT << "================================================================================";
		G_COUT << std::endl;
	}

	delete pTableSup;

	return 0;
}


#include <chrono>
#include <ctime>
typedef std::chrono::system_clock::time_point t_tp;
typedef std::chrono::system_clock::duration t_dur;

int main(int argc, char** argv)
{
	int iRet = -99;
	t_tp timeStart = std::chrono::system_clock::now();

	try
	{
		init_spec_chars();
		iRet = script_main(argc, argv);
	}
	catch(const std::exception& ex)
	{
		G_CERR << "Critical error in script interpreter: " << ex.what() << std::endl;
	}

	if(g_bShowTiming)
	{
		t_tp timeStop = std::chrono::system_clock::now();
		t_dur dur = timeStop-timeStart;
		double dDur = double(t_dur::period::num)/double(t_dur::period::den) * double(dur.count());

		std::time_t tStart = std::chrono::system_clock::to_time_t(timeStart);
		std::time_t tStop = std::chrono::system_clock::to_time_t(timeStop);

		std::tm tmStart = *std::localtime(&tStart);
		std::tm tmStop = *std::localtime(&tStop);

		char cStart[128], cStop[128];
		std::strftime(cStart, sizeof cStart, "%Y-%b-%d %H:%M:%S", &tmStart);
		std::strftime(cStop, sizeof cStop, "%Y-%b-%d %H:%M:%S", &tmStop);

		G_COUT << "\n";
		G_COUT << "================================================================================\n";	
		G_COUT << "Script start time:     " << cStart << "\n";
		G_COUT << "Script stop time:      " << cStop << "\n";
		G_COUT << "Script execution time: " << dDur << " s.\n";
		G_COUT << "================================================================================";
		G_COUT << std::endl;
	}

	return iRet;
}
