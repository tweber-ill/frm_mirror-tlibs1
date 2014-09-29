/*
 * Simple Script
 * @author tweber
 * @date 2013
 */

#include "types.h"
#include "helper/flags.h"
#include "helper/string.h"
#include "helper/spec_char.h"
#include "helper/log.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <array>

#include "parseobj.h"
#include "script_helper.h"
#include "globals.h"
#include "calls.h"
#include "node_opt.h"

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
		G_COUT << "\t-d[0-4] \t\tVerbosity (0=none, 1=errors, 2=warnings, 3=infos, 4=debug).\n";
		G_COUT << std::endl;
		return -1;
	}

	bool bShowSymbols = 0;
	unsigned int uiDebugLevel = 3;
#ifndef NDEBUG
	uiDebugLevel = 4;
#endif
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
		else if(strArg=="-d0") uiDebugLevel = 0;
		else if(strArg=="-d1") uiDebugLevel = 1;
		else if(strArg=="-d2") uiDebugLevel = 2;
		else if(strArg=="-d3") uiDebugLevel = 3;
		else if(strArg=="-d4") uiDebugLevel = 4;
	}

	const std::array<Log*, 5> arrLogs{&log_crit, &log_err, &log_warn, &log_info, &log_debug};
	for(unsigned int iLog=0; iLog<arrLogs.size(); ++iLog)
		arrLogs[iLog]->SetEnabled(uiDebugLevel>=iLog);

	// debug in script.yy needs to be set
	yydebug = (uiDebugLevel>=4);

	if(iStartArg >= argc)
	{
		log_err("No input file given.");
		return -1;
	}



	// loading of input file
	const char* pcFile = argv[iStartArg];
	t_string strFile = STR_TO_WSTR(pcFile);

	t_char* pcInput = load_file(pcFile);
	if(!pcInput)
		return -2;


	ParseObj par;
	ParseInfo info;
	info.bEnableDebug = (uiDebugLevel>=4);


	// lexing
	par.strCurFile = strFile;
	par.pLexer = new Lexer(pcInput, strFile.c_str());

	delete[] pcInput;
	pcInput = 0;

	if(!par.pLexer->IsOk())
	{
		log_err("Lexer returned with errors.");
		return -3;
	}

	init_global_syms(info.pGlobalSyms);


	// parsing
	int iParseRet = yyparse(&par);

	delete par.pLexer;
	par.pLexer = 0;

	if(iParseRet != 0)
	{
		log_err("Parser returned with error code ", iParseRet, ".");
		return -4;
	}


	// optimizing
	par.pRoot = par.pRoot->optimize();



	// executing
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
	delete pTableSup;


	if(bShowSymbols)
	{
		log_info("================================================================================");
		log_info("Global symbols:");
		info.pGlobalSyms->print();

		std::ostringstream ostrFkts;
		for(const NodeFunction* pFunc : info.vecFuncs)
			ostrFkts << pFunc->GetName() << ", ";
		log_info("Script functions: ", ostrFkts.str());


		const t_mapFkts* pExtFkts = get_ext_calls();

		std::ostringstream ostrSysFkts;
		for(const auto& fktpair : *pExtFkts)
			ostrSysFkts << fktpair.first << ", ";
		log_info("System functions: ", ostrSysFkts.str());
		log_info("================================================================================");
	}

	return 0;
}


#include <chrono>
#include <ctime>
typedef std::chrono::system_clock::time_point t_tp;
typedef std::chrono::system_clock::duration t_dur;

int main(int argc, char** argv)
{
	const std::array<Log*, 5> arrLogs{&log_crit, &log_err, &log_warn, &log_info, &log_debug};
	for(Log* pLog : arrLogs)
	{
		pLog->SetShowDate(0);
		pLog->SetShowThread(0);
	}

	int iRet = -99;
	t_tp timeStart = std::chrono::system_clock::now();

	try
	{
		init_spec_chars();
		iRet = script_main(argc, argv);
	}
	catch(const std::exception& ex)
	{
		log_crit(ex.what());
	}

	if(g_bShowTiming && iRet==0)
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

		log_info("================================================================================");
		log_info("Script start time:     ", cStart);
		log_info("Script stop time:      ", cStop);
		log_info("Script execution time: ", dDur, " s");
		log_info("================================================================================");
	}

	return iRet;
}
