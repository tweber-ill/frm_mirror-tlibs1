/*
 * External Functions
 * @author tweber
 */

#include "flags.h"
#include "calls.h"
#include "parseobj.h"
#include "script_helper.h"
#include "helper/gnuplot.h"
#include "helper/fourier.h"
#include "helper/string.h"
#include "loader/loadtxt.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <cstdio>
#include <cmath>
#include <wait.h>

static std::string linenr(const std::string& strErr, const ParseInfo &info)
{
	if(info.pCurCaller)
		return info.pCurCaller->linenr(strErr, info);
	return strErr + ": ";
}

static Symbol* fkt_version(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	extern const char* g_pcVersion;
	return new SymbolString(g_pcVersion);
}

static Symbol* fkt_int(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() == 0)
		return new SymbolInt(0);

	return vecSyms[0]->ToType(SYMBOL_INT);
}

static Symbol* fkt_double(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() == 0)
		return new SymbolDouble(0.);

	return vecSyms[0]->ToType(SYMBOL_DOUBLE);
}

static Symbol* fkt_str(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	SymbolString *pSymRet = new SymbolString;
	for(const Symbol *pSym : vecSyms)
		pSymRet->m_strVal += pSym->print();
	return pSymRet;
}

static Symbol* fkt_output(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	std::ostream& ostr = std::cout;
	
	for(Symbol *pSym : vecSyms)
		if(pSym)
			ostr << pSym->print();

	ostr.flush();
	return 0;
}

static Symbol* fkt_print(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	fkt_output(vecSyms, info, pSymTab);
	std::cout << std::endl;
	return 0;
}

static Symbol* fkt_input(const std::vector<Symbol*>& vecSyms,
		ParseInfo& info, SymbolTable* pSymTab)
{
	fkt_output(vecSyms, info, pSymTab);

	SymbolString* pSymStr = new SymbolString();
	std::getline(std::cin, pSymStr->m_strVal);
	return pSymStr;
}

static Symbol* fkt_exec(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	std::string strExec;
	
	for(Symbol *pSym : vecSyms)
		if(pSym)
		{
			strExec += pSym->print();
			strExec += " ";
		}
	
	bool bOk = 0;
	FILE *pPipe = ::popen(strExec.c_str(), "w");
	if(pPipe)
	{
		bOk = 1;
		int iRet = ::pclose(pPipe);
		if(iRet == -1)
		{
			bOk = 0;
		}
		else
		{
			//int bHasExited = WIFEXITED(iRet);
			int iExitCode = int(char(WEXITSTATUS(iRet)));
			//std::cout << "Exit code: " << iExitCode << std::endl;
			bOk = (iExitCode==0);
		}
	}

	return new SymbolInt(bOk);
}

extern int yyparse(void*);
static bool _import_file(const std::string& strFile, ParseInfo& info, SymbolTable* pSymTab)
{
	ParseInfo::t_mods::iterator iterMod = info.pmapModules->find(strFile);
	if(iterMod != info.pmapModules->end())
	{
		//std::cerr << "Warning: Module \"" << strFile << "\" already loaded." << std::endl;
		return 0;
	}

	char* pcInput = load_file(strFile.c_str());
	if(!pcInput)
		return 0;

	ParseObj par;
	par.strCurFile = strFile;
	par.pLexer = new Lexer(pcInput, strFile.c_str());

	delete[] pcInput;
	pcInput = 0;

	if(!par.pLexer->IsOk())
	{
		std::cerr << linenr("Error", info) << "Lexer returned with errors." << std::endl;
		return 0;
	}

	int iParseRet = yyparse(&par);

	delete par.pLexer;
	par.pLexer = 0;

	if(iParseRet != 0)
	{
		std::cerr << linenr("Error", info) << "Parser returned with error code " << iParseRet << "." << std::endl;
		return 0;
	}

	Node *pRoot = par.pRoot;
	info.pmapModules->insert(ParseInfo::t_mods::value_type(strFile, pRoot));
	info.strExecFkt = "";
	info.strInitScrFile = strFile;
	pRoot->eval(info);

	return 1;
}

static Symbol* fkt_import(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	for(Symbol *pSym : vecSyms)
		if(pSym && pSym->GetType()==SYMBOL_STRING)
		{
			const std::string& strFile = ((SymbolString*)pSym)->m_strVal;
			bool bOk = _import_file(strFile, info, pSymTab);
		}

	return 0;
}


static Symbol* fkt_typeof(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=1)
	{
		std::cerr << linenr("Error", info) << "typeof takes exactly one argument." << std::endl;
		return 0;
	}

	Symbol *pSymbol = vecSyms[0];
	if(!pSymbol)
	{
		std::cerr << linenr("Error", info) << "Invalid argument for typename." << std::endl;
		return 0;
	}

	SymbolString *pType = new SymbolString(pSymbol->GetTypeName().c_str());
	return pType;
}


// --------------------------------------------------------------------------------
// map

static Symbol* fkt_map(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	return new SymbolMap();
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// array

static Symbol* fkt_array(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
		return new SymbolArray();
	
	Symbol *pSymSize = vecSyms[0];
	if(pSymSize->GetType() != SYMBOL_INT)
	{
		std::cerr << linenr("Error", info) << "\"num\" in vec(num, val=0) has to be integer." << std::endl;
		return 0;
	}
	
	int iVal = ((SymbolInt*)pSymSize)->m_iVal;
	if(iVal < 0) iVal = 0;
	


	bool bOwnVal = 0;
	Symbol *pSymVal = 0;
	if(vecSyms.size()>1)
	{
		pSymVal = vecSyms[1];
	}
	else
	{
		pSymVal = new SymbolDouble(0.);
		bOwnVal = 1;
	}


	SymbolArray* pSymRet = new SymbolArray;
	pSymRet->m_arr.reserve(iVal);
	for(int i=0; i<iVal; ++i)
		pSymRet->m_arr.push_back(pSymVal->clone());

	pSymRet->UpdateIndices();

	if(bOwnVal)
		delete pSymVal;

	return pSymRet;
}

static Symbol* fkt_array_size(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		std::cerr << linenr("Error", info) << "vec_size(vec) needs one argument." << std::endl;
		return 0;
	}
	
	Symbol *pSymArr = vecSyms[0];
	SymbolInt *pSymRet = new SymbolInt(0);
	
	if(pSymArr->GetType() != SYMBOL_ARRAY)
	{
		std::cerr << linenr("Warning", info) << "vec_size needs an array type argument." << std::endl;
		return pSymRet;
	}

	
	pSymRet->m_iVal = ((SymbolArray*)pSymArr)->m_arr.size();
	return pSymRet;
}

static Symbol* fkt_cur_iter(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		std::cerr << linenr("Error", info) << "cur_iter takes exactly one argument." << std::endl;
		return 0;
	}

	const std::string& strIdent = vecSyms[0]->m_strIdent;
	if(strIdent == "")
	{
		std::cerr << linenr("Error", info) << "No identifier given for cur_iter." << std::endl;
		return 0;
	}


	std::string strIter = "<cur_iter_" + strIdent + ">";

	Symbol* pSymIter = pSymTab->GetSymbol(strIter);
	//pSymTab->print();
	if(!pSymIter || pSymIter->GetType()!=SYMBOL_INT)
	{
		std::cerr << linenr("Error", info) << "cur_iter could not determine iteration index \""
					<< strIter << "\"." << std::endl;
		return 0;
	}

	return pSymIter;
}

static inline Symbol* _fkt_linlogspace(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab, bool bLog)
{
	if(vecSyms.size()<3)
	{
		std::cerr << linenr("Error", info) << "Invalid call to linspace(start, end, count)." << std::endl;
		return 0;
	}

	SymbolDouble *pSymStart = (SymbolDouble*)vecSyms[0]->ToType(SYMBOL_DOUBLE);
	SymbolDouble *pSymEnd = (SymbolDouble*)vecSyms[1]->ToType(SYMBOL_DOUBLE);
	SymbolInt *pSymCount = (SymbolInt*)vecSyms[2]->ToType(SYMBOL_INT);

	int iNum = pSymCount->m_iVal;
	double dmin = pSymStart->m_dVal;
	double dmax = pSymEnd->m_dVal;

	delete pSymStart;
	delete pSymEnd;
	delete pSymCount;


	SymbolArray* pSymRet = new SymbolArray;
	pSymRet->m_arr.reserve(iNum);
	for(int i=0; i<iNum; ++i)
	{
		SymbolDouble *pSymD = new SymbolDouble();
		double dDiv = (iNum!=1 ? double(iNum-1) : 1);
		pSymD->m_dVal = double(i)*(dmax-dmin)/dDiv + dmin;
		if(bLog)
		{
			const double dBase = 10.;
			pSymD->m_dVal = std::pow(dBase, pSymD->m_dVal);
		}
		pSymRet->m_arr.push_back(pSymD);
	}

	pSymRet->UpdateIndices();

	return pSymRet;
}

static Symbol* fkt_linspace(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	return _fkt_linlogspace(vecSyms, info, pSymTab, false);
}

static Symbol* fkt_logspace(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	return _fkt_linlogspace(vecSyms, info, pSymTab, true);
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// string operations
static Symbol* fkt_trim(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=1)
	{
		std::cerr << linenr("Error", info)
				<< "Need one argument for trim." << std::endl;
		return 0;
	}

	if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray *pArr = new SymbolArray();
		pArr->m_arr.reserve(((SymbolArray*)vecSyms[0])->m_arr.size());

		for(Symbol* pSymArr : ((SymbolArray*)vecSyms[0])->m_arr)
		{
			std::vector<Symbol*> vecDummy = { pSymArr };
			pArr->m_arr.push_back(fkt_trim(vecDummy, info, pSymTab));
		}

		return pArr;
	}
	else if(vecSyms[0]->GetType() == SYMBOL_STRING)
	{
		std::string str = ((SymbolString*)vecSyms[0])->m_strVal;
		::trim(str);
		return new SymbolString(str);
	}

	// simply copy non-string arguments
	//std::cerr << linenr("Warning", info)
	//		<< "Called trim with invalid argument." << std::endl;
	return vecSyms[0]->clone();
}

static Symbol* fkt_split(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	std::string *pstrInput = 0;
	std::string strDelim = " \t\n";

	// split("Test 123")
	if(vecSyms.size() >= 1 && vecSyms[0]->GetType()==SYMBOL_STRING)
		pstrInput = &((SymbolString*)vecSyms[0])->m_strVal;

	// split("Test 123", " \t\n")
	if(vecSyms.size() >= 2 && vecSyms[1]->GetType()==SYMBOL_STRING)
		strDelim = ((SymbolString*)vecSyms[1])->m_strVal;

	if(!pstrInput)
	{
		std::cerr << linenr("Error", info)
				<< "Called split with invalid arguments." << std::endl;
		return 0;
	}

	std::vector<std::string> vecTokens;
	::get_tokens<std::string>(*pstrInput, strDelim, vecTokens);

	SymbolArray* pArr = new SymbolArray;
	for(const std::string& strTok : vecTokens)
		pArr->m_arr.push_back(new SymbolString(strTok));

	return pArr;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// thread

std::vector<Symbol*>* clone_symbols(const std::vector<Symbol*>* pvecSyms,
								unsigned int iBegin=0)
{
	if(!pvecSyms)
		return 0;

	std::vector<Symbol*> *pvec = new std::vector<Symbol*>;
	pvec->reserve(pvecSyms->size());

	for(unsigned int i=iBegin; i<pvecSyms->size(); ++i)
	{
		Symbol *pSym = (*pvecSyms)[i];
		if(pSym) pSym = pSym->clone();
		pvec->push_back(pSym);
	}

	return pvec;
}

void delete_symbols(std::vector<Symbol*>* pvecSyms)
{
	for(Symbol *pSym : *pvecSyms)
		delete pSym;
	delete pvecSyms;
}

static void thread_proc(NodeFunction* pFunc, ParseInfo* pinfo, std::vector<Symbol*>* pvecSyms)
{
	if(!pFunc || !pinfo) return;

	pinfo->pmutexInterpreter->lock();
		NodeFunction* pThreadFunc = (NodeFunction*)pFunc->clone();
		pThreadFunc->SetArgSyms(pvecSyms);
		ParseInfo *pinfo2 = new ParseInfo(*pinfo);	// threads cannot share the same bWantReturn etc.
		pinfo2->bDestroyParseInfo = 0;
	pinfo->pmutexInterpreter->unlock();

	Symbol* pRet = pThreadFunc->eval(*pinfo2, 0);

	if(pRet) delete pRet;
	if(pvecSyms) delete_symbols(pvecSyms);
	if(pThreadFunc) delete pThreadFunc;
	if(pinfo2) delete pinfo2;
}

static Symbol* fkt_thread(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		std::cerr << linenr("Error", info) << "Need thread proc identifier." << std::endl;
		return 0;
	}

	Symbol* _pSymIdent = vecSyms[0];
	if(_pSymIdent->GetType() != SYMBOL_STRING)
	{
		std::cerr << linenr("Error", info) << "Thread proc identifier needs to be a string." << std::endl;
		return 0;
	}

	SymbolString *pSymIdent = (SymbolString*)_pSymIdent;
	const std::string& strIdent = pSymIdent->m_strVal;


	NodeFunction* pFunc = info.GetFunction(strIdent);
	if(pFunc == 0)
	{
		std::cerr << linenr("Error", info) << "Thread proc \"" << strIdent << "\" not defined." << std::endl;
		return 0;
	}

	std::vector<Symbol*>* vecThreadSymsClone = clone_symbols(&vecSyms, 1);
	std::thread* pThread = new std::thread(::thread_proc, pFunc, &info, vecThreadSymsClone);
	unsigned int iHandle = info.phandles->AddHandle(new HandleThread(pThread));

	return new SymbolInt(iHandle);
}

static Symbol* fkt_thread_hwcount(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	unsigned int iNumThreads = std::thread::hardware_concurrency();
	if(iNumThreads == 0)
		iNumThreads = 1;
	
	return new SymbolInt(iNumThreads);
}	

static Symbol* fkt_begin_critical(const std::vector<Symbol*>& vecSyms,
								ParseInfo& info, SymbolTable* pSymTab)
{
	info.pmutexGlobal->lock();
	return 0;
}

static Symbol* fkt_end_critical(const std::vector<Symbol*>& vecSyms,
								ParseInfo& info, SymbolTable* pSymTab)
{
	info.pmutexGlobal->unlock();
	return 0;
}


static Symbol* fkt_thread_join(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		std::cerr << linenr("Error", info) << "join needs at least one argument." << std::endl;
		return 0;
	}

	for(Symbol* pSym : vecSyms)
	{
		if(pSym == 0) continue;

		if(pSym->GetType() == SYMBOL_ARRAY)
		{
			return fkt_thread_join(((SymbolArray*)pSym)->m_arr, info, pSymTab);
		}

		if(pSym->GetType() != SYMBOL_INT)
		{
			std::cerr << linenr("Error", info) << "join needs thread handles." << std::endl;
			continue;
		}

		int iHandle = ((SymbolInt*)pSym)->m_iVal;
		Handle *pHandle = info.phandles->GetHandle(iHandle);

		if(pHandle==0 || pHandle->GetType()!=HANDLE_THREAD)
		{
			std::cerr << linenr("Error", info) << "Handle (" << iHandle << ") does not exist"
					 << " or is not a thread handle." << std::endl;
			continue;
		}

		HandleThread *pThreadHandle = (HandleThread*)pHandle;
		std::thread *pThread = pThreadHandle->GetInternalHandle();

		pThread->join();
	}
	return 0;
}


// --------------------------------------------------------------------------------
// nthread

// nthread(iNumThreads, strFunc, vecArgs, ...)
static Symbol* fkt_nthread(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()<3)
	{
		std::cerr << linenr("Error", info) << "nthread needs at least 3 arguments: N, func, arg." << std::endl;
		return 0;
	}

	Symbol* _pSymN = vecSyms[0];
	if(_pSymN->GetType() != SYMBOL_INT)
	{
		std::cerr << linenr("Error", info) << "Number of threads has to be integer." << std::endl;
		return 0;
	}

	SymbolInt *pSymN = (SymbolInt*)_pSymN;
	int iNumThreads = pSymN->m_iVal;



	Symbol* _pSymIdent = vecSyms[1];
	if(_pSymIdent->GetType() != SYMBOL_STRING)
	{
		std::cerr << linenr("Error", info) << "Thread proc identifier needs to be a string." << std::endl;
		return 0;
	}

	SymbolString *pSymIdent = (SymbolString*)_pSymIdent;
	const std::string& strIdent = pSymIdent->m_strVal;



	Symbol* _pSymArr = vecSyms[2];
	if(_pSymArr->GetType() != SYMBOL_ARRAY)
	{
		std::cerr << linenr("Error", info) << "Thread arg has to be an array." << std::endl;
		return 0;
	}

	SymbolArray *pSymArr = (SymbolArray*)_pSymArr;
	const std::vector<Symbol*>& vecArr = pSymArr->m_arr;



	NodeFunction* pFunc = info.GetFunction(strIdent);
	if(pFunc == 0)
	{
		std::cerr << linenr("Error", info) << "Thread proc \"" << strIdent << "\" not defined." << std::endl;
		return 0;
	}





	if(iNumThreads > vecArr.size())
	{
		iNumThreads = vecArr.size();
		std::cerr << linenr("Warning", info) << "More threads requested in nthread than necessary, "
						  << "reducing to array size (" << iNumThreads << ")."
						  << std::endl;
	}


	std::vector<SymbolArray*> vecSymArrays;
	vecSymArrays.resize(iNumThreads);

	int iCurTh = 0;
	for(Symbol* pThisSym : vecArr)
	{
		if(!vecSymArrays[iCurTh])
			vecSymArrays[iCurTh] = new SymbolArray();

		vecSymArrays[iCurTh]->m_arr.push_back(pThisSym->clone());

		++iCurTh;
		if(iCurTh == iNumThreads)
			iCurTh = 0;
	}



	std::vector<std::thread*> vecThreads;
	vecThreads.reserve(iNumThreads);

	for(iCurTh=0; iCurTh<iNumThreads; ++iCurTh)
	{
		std::vector<Symbol*>* vecThreadSyms = new std::vector<Symbol*>;
		vecThreadSyms->reserve(vecSyms.size()-3+1);

		vecThreadSyms->push_back(vecSymArrays[iCurTh]);

		for(unsigned int iSym=3; iSym<vecSyms.size(); ++iSym)
			vecThreadSyms->push_back(vecSyms[iSym]->clone());

		std::thread *pth = new std::thread(::thread_proc, pFunc, &info, vecThreadSyms);
		vecThreads.push_back(pth);
	}

	/*
	// automatically join
	for(iCurTh=0; iCurTh<iNumThreads; ++iCurTh)
	{
		vecThreads[iCurTh]->join();
		delete vecThreads[iCurTh];
		vecThreads[iCurTh] = 0;
	}*/


	SymbolArray* pArrThreads = new SymbolArray();

	for(iCurTh=0; iCurTh<iNumThreads; ++iCurTh)
	{
		std::thread* pCurThread = vecThreads[iCurTh];
		unsigned int iHandle = info.phandles->AddHandle(new HandleThread(pCurThread));
		SymbolInt *pSymThreadHandle = new SymbolInt(iHandle);

		pArrThreads->m_arr.push_back(pSymThreadHandle);
	}

	return pArrThreads;
}

// --------------------------------------------------------------------------------





// --------------------------------------------------------------------------------
// math
enum MathFkts
{
	MATH_MAX,
	MATH_MIN
};

template<MathFkts fkt, typename T>
const T& math_fkt(const T& t1, const T& t2)
{
	if(fkt == MATH_MAX)
		return std::max<T>(t1, t2);
	else if(fkt == MATH_MIN)
		return std::min<T>(t1, t2);

	static const T tErr = T(0);
	std::cerr << "Error: Invalid function selected in math_fkt." << std::endl;
	return tErr;
}

template<MathFkts fkt>
static Symbol* fkt_math_for_every(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() < 1)
	{
		std::cerr << linenr("Error", info) << "fkt_math_for_every needs at least one argument" << std::endl;
		return 0;
	}

	double dRes;
	int iRes;

	bool bHadInt = 0, 
		bHadDouble = 0;

	for(Symbol* pSym : vecSyms)
	{
		Symbol *pThisSym = pSym;
		bool bCleanSym = 0;
		if(pSym->GetType() == SYMBOL_ARRAY)
		{
			pThisSym = fkt_math_for_every<fkt>(
					((SymbolArray*)pSym)->m_arr, 
					info, pSymTab);

			bCleanSym = 1;
		}

		if(pThisSym->GetType() == SYMBOL_INT)
		{
			if(!bHadInt)
				iRes = ((SymbolInt*)pThisSym)->m_iVal;
			else
				iRes = math_fkt<fkt, int>(iRes, ((SymbolInt*)pThisSym)->m_iVal);

			bHadInt = 1;
		}
		else if(pThisSym->GetType() == SYMBOL_DOUBLE)
		{
			if(!bHadDouble)
				dRes = ((SymbolDouble*)pThisSym)->m_dVal;
			else
				dRes = math_fkt<fkt, double>(dRes, ((SymbolDouble*)pThisSym)->m_dVal);

			bHadDouble = 1;
		}

		if(bCleanSym)
			delete pThisSym;
	}

	if(bHadInt && !bHadDouble)
		return new SymbolInt(iRes);
	else if(bHadInt && bHadDouble)
	{
		dRes = math_fkt<fkt, double>(dRes, double(iRes));
		return new SymbolDouble(dRes);
	}
	else if(!bHadInt && bHadDouble)
		return new SymbolDouble(dRes);

	std::cerr << linenr("Error", info) << "No valid arguments given for fkt_math_for_every." << std::endl;
	return 0;
}


template<double (*FKT)(double)>
static Symbol* fkt_math_1arg(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		std::cerr << linenr("Error", info) << "fkt_math_1arg takes exactly one argument." << std::endl;
		return 0;
	}
	
	if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray* pArrRet = new SymbolArray();
		
		SymbolArray* pSymArr = (SymbolArray*)vecSyms[0];
		for(Symbol* pArrElem : pSymArr->m_arr)
		{
			std::vector<Symbol*> vecDummy;
			vecDummy.push_back(pArrElem);
			
			pArrRet->m_arr.push_back(fkt_math_1arg<FKT>(vecDummy, info, pSymTab));
		}
		
		pArrRet->UpdateIndices();
		return pArrRet;
	}
	else
	{
		SymbolDouble* pSym = (SymbolDouble*)vecSyms[0]->ToType(SYMBOL_DOUBLE);
		double dResult = FKT(pSym->m_dVal);
		
		// recycle this SymbolDouble
		pSym->m_dVal = dResult;
		return pSym;
	}
	
	return 0;
}

// TODO: for arrays
template<double (*FKT)(double, double)>
static Symbol* fkt_math_2args(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2)
	{
		std::cerr << linenr("Error", info) << "fkt_math_2args takes exactly two arguments." << std::endl;
		return 0;
	}
	
	SymbolDouble* pSym1 = (SymbolDouble*)vecSyms[0]->ToType(SYMBOL_DOUBLE);
	SymbolDouble* pSym2 = (SymbolDouble*)vecSyms[1]->ToType(SYMBOL_DOUBLE);
	double dResult = FKT(pSym1->m_dVal, pSym2->m_dVal);
	delete pSym2;
	
	// recycle this SymbolDouble
	pSym1->m_dVal = dResult;
	return pSym1;
}


template<typename T> T myabs(T t)
{
	if(t < T(0))
		return -t;
	return t;
}

template<> double myabs(double t) { return ::fabs(t); }


// TODO: integrate with fkt_math_1arg
static Symbol* fkt_math_abs(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		std::cerr << linenr("Error", info) << "abs takes exactly one argument." << std::endl;
		return 0;
	}
	
	if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray* pArrRet = new SymbolArray();
		
		SymbolArray* pSymArr = (SymbolArray*)vecSyms[0];
		for(Symbol* pArrElem : pSymArr->m_arr)
		{
			std::vector<Symbol*> vecDummy;
			vecDummy.push_back(pArrElem);
			
			pArrRet->m_arr.push_back(fkt_math_abs(vecDummy, info, pSymTab));
		}
		
		pArrRet->UpdateIndices();
		return pArrRet;
	}
	else if(vecSyms[0]->GetType() == SYMBOL_INT)
	{
		SymbolInt* pSymInt = (SymbolInt*)vecSyms[0];
		return new SymbolInt(myabs(pSymInt->m_iVal));
	}
	else if(vecSyms[0]->GetType() == SYMBOL_DOUBLE)
	{
		SymbolDouble* pSymD = (SymbolDouble*)vecSyms[0];
		return new SymbolDouble(myabs(pSymD->m_dVal));
	}
	

	std::cerr << linenr("Error", info) << "abs received unsupported symbol type." << std::endl;
	return 0;
}


// --------------------------------------------------------------------------------
// FFT

static Symbol* _fkt_fft(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab,
						bool bInv)
{
	bool (Fourier::*pFkt)(const double*, const double*, double*, double*)
								= (bInv ? &Fourier::ifft : &Fourier::fft);

	bool bArgsOk=1;
	std::vector<double> vecRealIn, vecImagIn;

	// real and imag part as separate arguments
	if(vecSyms.size()==2 && vecSyms[0]->GetType()==SYMBOL_ARRAY 
			&& vecSyms[1]->GetType()==SYMBOL_ARRAY)
	{
		vecRealIn = ((SymbolArray*)vecSyms[0])->ToDoubleArray();
		vecImagIn = ((SymbolArray*)vecSyms[1])->ToDoubleArray();
	}
	// arrays in one array
	else if(vecSyms.size()==1 && vecSyms[0]->GetType()==SYMBOL_ARRAY)
	{
		SymbolArray* pSymArr = (SymbolArray*)vecSyms[0];
		unsigned int iSymArrSize = pSymArr->m_arr.size();

		if(iSymArrSize==0)
			bArgsOk = 0;
		if(iSymArrSize>=1 && pSymArr->m_arr[0]->GetType()==SYMBOL_ARRAY)
			vecRealIn = ((SymbolArray*)pSymArr->m_arr[0])->ToDoubleArray();
		if(iSymArrSize>=2 && pSymArr->m_arr[1]->GetType()==SYMBOL_ARRAY)
			vecImagIn = ((SymbolArray*)pSymArr->m_arr[1])->ToDoubleArray();

		// simple array containing real data
		if(pSymArr->m_arr[0]->GetType()!=SYMBOL_ARRAY)
			vecRealIn = pSymArr->ToDoubleArray();
	}
	
	if(!bArgsOk)
	{
		std::cerr << linenr("Error", info)
				<< "fft received invalid arguments."
				<< std::endl;
		return 0;
	}

	if(vecRealIn.size() != vecImagIn.size())
	{
		unsigned int iSize = std::max(vecRealIn.size(), vecImagIn.size());
		vecRealIn.resize(iSize);
		vecImagIn.resize(iSize);
	}

	std::vector<double> vecRealOut, vecImagOut;
	vecRealOut.resize(vecRealIn.size());
	vecImagOut.resize(vecImagIn.size());

	Fourier fourier(vecRealIn.size());
	(fourier.*pFkt)(vecRealIn.data(), vecImagIn.data(), 
			vecRealOut.data(), vecImagOut.data());

	SymbolArray* pArrReal = new SymbolArray();
	SymbolArray* pArrImag = new SymbolArray();
	
	pArrReal->FromDoubleArray(vecRealOut);
	pArrImag->FromDoubleArray(vecImagOut);

	SymbolArray* pRet = new SymbolArray();
	pRet->m_arr.push_back(pArrReal);
	pRet->m_arr.push_back(pArrImag);

	return pRet;
}

static Symbol* fkt_fft(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{ return _fkt_fft(vecSyms, info, pSymTab, false); }
static Symbol* fkt_ifft(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{ return _fkt_fft(vecSyms, info, pSymTab, true); }

// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// linalg stuff

static Symbol* fkt_dot(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2)
	{
		std::cerr << linenr("Error", info) << "dot needs two arguments." << std::endl;
		return 0;
	}

	if(vecSyms[0]->GetType()!=SYMBOL_ARRAY || vecSyms[1]->GetType()!=SYMBOL_ARRAY)
	{
		std::cerr << linenr("Error", info) << "dot needs vector arguments." << std::endl;
		return 0;
	}

	SymbolArray* pLeft = (SymbolArray*)vecSyms[0];
	SymbolArray* pRight = (SymbolArray*)vecSyms[1];

	// TODO: avoid new/deletes with NODE_PLUSEQUAL
	Symbol *pSum = 0;
	unsigned int iSize = std::min(pLeft->m_arr.size(), pRight->m_arr.size());
	for(unsigned int i=0; i<iSize; ++i)
	{
		Symbol *pProd = Node::Op(pLeft->m_arr[i], pRight->m_arr[i], NODE_MULT);
		if(i==0)
			pSum = pProd->clone();
		else
		{
			Symbol *pSumNew = Node::Op(pSum, pProd, NODE_PLUS);
			delete pSum;
			pSum = pSumNew;
		}
		delete pProd;
	}

	return pSum;
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// plotting

GnuPlot g_plot;

static inline bool is_array_of_arrays(const Symbol* pSym)
{
	if(!pSym) return 0;
	if(pSym->GetType()!=SYMBOL_ARRAY) return 0;

	SymbolArray* pSymArr = (SymbolArray*)pSym;
	if(pSymArr->m_arr.size()==0) return 0;

	Symbol* pSymInArr = pSymArr->m_arr[0];
	return (pSymInArr->GetType()==SYMBOL_ARRAY);
}

static inline bool is_array_of_array_of_arrays(const Symbol* pSym)
{
	if(!pSym) return 0;
	if(pSym->GetType() != SYMBOL_ARRAY)
		return 0;

	if(((SymbolArray*)pSym)->m_arr.size() == 0)
		return 0;

	return is_array_of_arrays(((SymbolArray*)pSym)->m_arr[0]);
}

struct XYLimits
{
	bool bHasX, bHasY, bHasCB, bCBCyclic;
	double dMinX, dMaxX;
	double dMinY, dMaxY;
	double dMinCB, dMaxCB;

	XYLimits() : bHasX(0), bHasY(0), bHasCB(0), bCBCyclic(0)
	{}
};

static XYLimits get_plot_limits(SymbolMap* pParamMap)
{
	XYLimits lim;

	bool bHasVal=0;
	std::string strVal = pParamMap->GetStringVal("xylimits", &bHasVal);
	if(bHasVal)
	{
		std::istringstream istr(strVal);
		istr >> lim.dMinX >> lim.dMaxX >> lim.dMinY >> lim.dMaxY;

		lim.bHasX = 1;
		lim.bHasY = 1;
	}
	else
	{
		std::string strX = pParamMap->GetStringVal("xlimits", &lim.bHasX);
		std::string strY = pParamMap->GetStringVal("ylimits", &lim.bHasY);

		if(lim.bHasX)
		{
			std::istringstream istrX(strX);
			istrX >> lim.dMinX >> lim.dMaxX;
		}

		if(lim.bHasY)
		{
			std::istringstream istrY(strY);
			istrY >> lim.dMinY >> lim.dMaxY;
		}
	}
	//std::cout << "xlimits: " << lim.dMinX << ", " << lim.dMaxX << std::endl;
	//std::cout << "ylimits: " << lim.dMinY << ", " << lim.dMaxY << std::endl;


	std::string strValCB = pParamMap->GetStringVal("cblimits", &lim.bHasCB);
	std::istringstream istrCB(strValCB);
	istrCB >> lim.dMinCB >> lim.dMaxCB;
	//std::cout << "colorbar: " << lim.dMinCB << ", " << lim.dMaxCB << std::endl;

	bool bHasCyc = 0;
	std::string strValCBCyc = pParamMap->GetStringVal("cbcyclic", &bHasCyc);
	if(bHasCyc)
	{
		std::istringstream istrCyc(strValCBCyc);
		istrCyc >> lim.bCBCyclic;
	}

	return lim;
}


static Symbol* fkt_plot(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	g_plot.Init();

	unsigned int iNumSyms = vecSyms.size();

	// plot([[x, y, yerr, xerr, mapParams], ...]);
	if(iNumSyms==1 && is_array_of_array_of_arrays(vecSyms[0]))
		return fkt_plot(((SymbolArray*)vecSyms[0])->m_arr, info, pSymTab);
	// plot([x, y, yerr, xerr, mapParams], [x2, y2, yerr2, xerr2, mapParams2], ...)
	else if(iNumSyms>=1 && is_array_of_arrays(vecSyms[0]))
	{
		g_plot.StartPlot();
		for(Symbol *pArr : vecSyms)
		{
			// ignore non-array arguments
			if(pArr->GetType() != SYMBOL_ARRAY)
				continue;

			fkt_plot(((SymbolArray*)pArr)->m_arr, info, pSymTab);
		}
		g_plot.FinishPlot();
	}
	// plot(x, y, yerr, xerr, mapParams)
	else if(iNumSyms >= 2 && vecSyms[0]->GetType()==SYMBOL_ARRAY && vecSyms[1]->GetType()==SYMBOL_ARRAY)
	{
		std::vector<double> vecX = ((SymbolArray*)vecSyms[0])->ToDoubleArray();
		std::vector<double> vecY = ((SymbolArray*)vecSyms[1])->ToDoubleArray();
		std::vector<double> vecYErr, vecXErr;

		if(iNumSyms >= 3 && vecSyms[2]->GetType()==SYMBOL_ARRAY)
			vecYErr = ((SymbolArray*)vecSyms[2])->ToDoubleArray();
		if(iNumSyms >= 4 && vecSyms[3]->GetType()==SYMBOL_ARRAY)
			vecXErr = ((SymbolArray*)vecSyms[3])->ToDoubleArray();


		PlotObj obj;
		obj.vecX = vecX;
		obj.vecY = vecY;
		obj.vecErrX = vecXErr;
		obj.vecErrY = vecYErr;

		// parameter map given as last argument
		if(vecSyms[iNumSyms-1]->GetType()==SYMBOL_MAP)
		{
			SymbolMap *pParamMap = (SymbolMap*)vecSyms[iNumSyms-1];

			bool bHasVal = 0;
			std::string strTitle = pParamMap->GetStringVal("title", &bHasVal);
			if(bHasVal) g_plot.SetTitle(strTitle.c_str());

			std::string strXLab = pParamMap->GetStringVal("xlabel", &bHasVal);
			if(bHasVal) g_plot.SetXLabel(strXLab.c_str());

			std::string strYLab = pParamMap->GetStringVal("ylabel", &bHasVal);
			if(bHasVal) g_plot.SetYLabel(strYLab.c_str());

			std::string strStyle = pParamMap->GetStringVal("style", &bHasVal);
			if(bHasVal) obj.bConnectLines = (strStyle=="line");

			std::string strLegend = pParamMap->GetStringVal("legend", &bHasVal);
			if(bHasVal) obj.strLegend = strLegend;

			XYLimits lim = get_plot_limits(pParamMap);
			if(lim.bHasX) g_plot.SetXRange(lim.dMinX, lim.dMaxX);
			if(lim.bHasY) g_plot.SetYRange(lim.dMinY, lim.dMaxY);
		}

		g_plot.StartPlot();
		g_plot.AddLine(obj);
		g_plot.FinishPlot();
	}
	else
	{
		std::cerr << linenr("Error", info) << "Invalid call to plot." << std::endl;
		return 0;
	}

	return 0;
}

static Symbol* fkt_plot2d(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	g_plot.Init();

	unsigned int iNumSyms = vecSyms.size();

	// e.g. plot2d([[1,2],[3,4]], params)
	if(iNumSyms >= 1 && is_array_of_arrays(vecSyms[0]))
	{
		SymbolArray* pArr = (SymbolArray*)vecSyms[0];
		SymbolMap* pMapParam = 0;

		// has parameter map
		if(iNumSyms == 2 && vecSyms[1]->GetType()==SYMBOL_MAP)
			pMapParam = (SymbolMap*)vecSyms[1];

		std::vector<std::vector<double> > vecXY;
		vecXY.reserve(pArr->m_arr.size());

		for(const Symbol* _pX : pArr->m_arr)
		{
			SymbolArray* pX = (SymbolArray*)_pX;
			std::vector<double> vecX = pX->ToDoubleArray();

			vecXY.push_back(vecX);
		}

		double dRMinX=1., dRMaxX=-1., dRMinY=1., dRMaxY=-1.;
		if(pMapParam)
		{
			bool bHasVal = 0;
			std::string strTitle = pMapParam->GetStringVal("title", &bHasVal);
			if(bHasVal) g_plot.SetTitle(strTitle.c_str());

			std::string strXLab = pMapParam->GetStringVal("xlabel", &bHasVal);
			if(bHasVal) g_plot.SetXLabel(strXLab.c_str());

			std::string strYLab = pMapParam->GetStringVal("ylabel", &bHasVal);
			if(bHasVal) g_plot.SetYLabel(strYLab.c_str());

			XYLimits lim = get_plot_limits(pMapParam);
			if(lim.bHasX) { dRMinX = lim.dMinX; dRMaxX = lim.dMaxX; }
			if(lim.bHasY) { dRMinY = lim.dMinY; dRMaxY = lim.dMaxY; }
			if(lim.bHasCB) g_plot.SetColorBarRange(lim.dMinCB, lim.dMaxCB, lim.bCBCyclic);
		}

		g_plot.SimplePlot2d(vecXY, dRMinX, dRMaxX, dRMinY, dRMaxY);
	}
	else
	{
		std::cerr << linenr("Error", info) << "Invalid call to plot2d." << std::endl;
		return 0;
	}

	return 0;
}

// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// loading and saving .dat files

static Symbol* fkt_loadtxt(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		std::cerr << linenr("Error", info) << "loadtxt takes exactly one argument." << std::endl;
		return 0;
	}

	if(vecSyms[0]->GetType() != SYMBOL_STRING)
	{
		std::cerr << linenr("Error", info) << "loadtxt needs a string argument." << std::endl;
		return 0;
	}

	const std::string& strFile = ((SymbolString*)vecSyms[0])->m_strVal;


	SymbolArray *pArr = new SymbolArray();

	LoadTxt dat;

	bool bLoaded = dat.Load(strFile.c_str());
	if(!bLoaded)
	{
		std::cerr << linenr("Error", info) << "loadtxt could not open \"" << strFile << "\"." << std::endl;
		return pArr;
	}

	pArr->m_arr.reserve(dat.GetColCnt());
	for(unsigned int iCol=0; iCol<dat.GetColCnt(); ++iCol)
	{
		const unsigned int iColLen = dat.GetColLen();
		const double *pCol = dat.GetColumn(iCol);

		SymbolArray *pArrCol = new SymbolArray;
		pArrCol->m_arr.reserve(iColLen);

		for(unsigned int iRow=0; iRow<iColLen; ++iRow)
		{
			SymbolDouble* pSymD = new SymbolDouble();
			pSymD->m_dVal = pCol[iRow];

			pArrCol->m_arr.push_back(pSymD);
		}

		pArr->m_arr.push_back(pArrCol);
	}

	// load the parameter map
	typedef std::map<std::string, std::string> tmapcomm;
	tmapcomm mapComm = dat.GetCommMapSingle();

	SymbolMap *pSymMap = new SymbolMap();
	for(const tmapcomm::value_type &val : mapComm)
	{
		SymbolString *pSymStrVal = new SymbolString;
		pSymStrVal->m_strVal = val.second;

		pSymMap->m_map.insert(SymbolMap::t_map::value_type(val.first, pSymStrVal));
	}
	pArr->m_arr.push_back(pSymMap);

	return pArr;
}

static void get_2darr_size(const SymbolArray* pArr,
				unsigned int& iColLen, unsigned int& iRowLen)
{
	iColLen = pArr->m_arr.size();
	iRowLen = 0;

	if(iColLen)
	{
		// look for first real array (not the parameter map)
		for(unsigned int iCol=0; iCol<iColLen; ++iCol)
		{
			Symbol* pSym = pArr->m_arr[iCol];
			if(pSym->GetType() == SYMBOL_ARRAY)
			{
				iRowLen = ((SymbolArray*)pSym)->m_arr.size();
				break;
			}
		}
	}

	unsigned int iNonArray = 0;
	// don't count the parameter map
	for(unsigned int iCol=0; iCol<iColLen; ++iCol)
	{
		Symbol* pSym = pArr->m_arr[iCol];
		if(pSym->GetType() != SYMBOL_ARRAY)
			++iNonArray;
	}

	iColLen -= iNonArray;
}

static double get_2darr_val(const SymbolArray* pArr,
				unsigned int iCol, unsigned int iRow)
{
	unsigned int iColLen = pArr->m_arr.size();
	if(iCol >= iColLen)
		return 0.;

	bool bFoundCol = 0;
	unsigned int iColRealArray = 0;
	for(unsigned int iCurCol=0; iCurCol<iColLen; ++iCurCol)
	{
		Symbol *pSym = pArr->m_arr[iCurCol];
		if(pSym->GetType() == SYMBOL_ARRAY)
		{
			if(iColRealArray == iCol)
			{
				bFoundCol = 1;
				iCol = iCurCol;
				break;
			}

			++iColRealArray;
		}
	}

	if(!bFoundCol)
	{
		std::cerr << "Error: Invalid column index: " << iCol << "." << std::endl;
		return 0.;
	}


	Symbol *pSym = pArr->m_arr[iCol];

	const std::vector<Symbol*>& veccol = ((SymbolArray*)pSym)->m_arr;
	if(iRow >= veccol.size())
		return 0.;

	const Symbol *pVal = veccol[iRow];
	if(pVal->GetType() == SYMBOL_DOUBLE)
	{
		return ((SymbolDouble*)pVal)->m_dVal;
	}
	else
	{
		SymbolDouble* pVal2 = (SymbolDouble*)pVal->ToType(SYMBOL_DOUBLE);
		double dVal = pVal2->m_dVal;
		delete pVal2;

		return dVal;
	}

	return 0.;
}

static Symbol* fkt_savetxt(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=2 || 
		(vecSyms[0]->GetType()!=SYMBOL_STRING || vecSyms[1]->GetType()!=SYMBOL_ARRAY))
	{
		std::cerr << linenr("Error", info) << "savetxt takes two arguments (file name, 2d array)." << std::endl;
		return 0;
	}

	const std::string& strFile = ((SymbolString*)vecSyms[0])->m_strVal;
	SymbolArray* pArr = (SymbolArray*)vecSyms[1];

	
	std::ofstream ofstr(strFile.c_str());
	if(!ofstr.is_open())
	{
		std::cerr << linenr("Error", info) << "Cannot open \"" << strFile << "\"." << std::endl;
		return 0;
	}

	// save parameter map
	for(unsigned int iCol=0; iCol<pArr->m_arr.size(); ++iCol)
	{
		Symbol *_pSymMap = pArr->m_arr[iCol];
		if(_pSymMap && _pSymMap->GetType() == SYMBOL_MAP)
		{
			SymbolMap *pSymMap = (SymbolMap*)_pSymMap;

			for(const SymbolMap::t_map::value_type& val : pSymMap->m_map)
			{
				ofstr << "# " << val.first << " : "
							<< (val.second?val.second->print():"") << "\n";
			}
		}
	}

	unsigned int iColLen=0, iRowLen=0;
	get_2darr_size(pArr, iColLen, iRowLen);
	//std::cout << "col len: " << iColLen << ", row len: " << iRowLen << std::endl;

	for(unsigned int iRow=0; iRow<iRowLen; ++iRow)
	{
		for(unsigned int iCol=0; iCol<iColLen; ++iCol)
			ofstr << get_2darr_val(pArr, iCol, iRow) << " ";
		ofstr << "\n";
	}

	ofstr.flush();
	ofstr.close();
	return 0;
}

// --------------------------------------------------------------------------------



typedef std::map<std::string, Symbol*(*)(const std::vector<Symbol*>&, ParseInfo&, SymbolTable*)> t_mapFkts;
static t_mapFkts g_mapFkts =
{
	// basic stuff
	t_mapFkts::value_type("ver", fkt_version),

	// input/output
	t_mapFkts::value_type("output", fkt_output),
	t_mapFkts::value_type("input", fkt_input),
	t_mapFkts::value_type("print", fkt_print),	// output with "\n" at the end

	// modules
	t_mapFkts::value_type("exec", fkt_exec),
	t_mapFkts::value_type("import", fkt_import),
	
	// casts
	t_mapFkts::value_type("int", fkt_int),
	t_mapFkts::value_type("real", fkt_double),
	t_mapFkts::value_type("str", fkt_str),
	t_mapFkts::value_type("map", fkt_map),
	t_mapFkts::value_type("vec", fkt_array),
	t_mapFkts::value_type("typeof", fkt_typeof),

	// string operations
	t_mapFkts::value_type("trim", fkt_trim),
	t_mapFkts::value_type("split", fkt_split),

	// array operations
	t_mapFkts::value_type("vec_size", fkt_array_size),
	t_mapFkts::value_type("cur_iter", fkt_cur_iter),

	// vector / matrix operations
	t_mapFkts::value_type("dot", fkt_dot),

	// plotting
	t_mapFkts::value_type("plot", fkt_plot),
	t_mapFkts::value_type("plot2d", fkt_plot2d),

	// threads
	t_mapFkts::value_type("thread", fkt_thread),
	t_mapFkts::value_type("nthread", fkt_nthread),
	t_mapFkts::value_type("thread_hwcount", fkt_thread_hwcount),
	t_mapFkts::value_type("join", fkt_thread_join),
	t_mapFkts::value_type("begin_critical", fkt_begin_critical),
	t_mapFkts::value_type("end_critical", fkt_end_critical),
	
	// math stuff
	t_mapFkts::value_type("sqrt", fkt_math_1arg< ::sqrt >),
	t_mapFkts::value_type("cbrt", fkt_math_1arg< ::cbrt >),
	t_mapFkts::value_type("exp", fkt_math_1arg< ::exp >),
	t_mapFkts::value_type("exp2", fkt_math_1arg< ::exp2 >),
	t_mapFkts::value_type("expm1", fkt_math_1arg< ::expm1 >),
	t_mapFkts::value_type("log", fkt_math_1arg< ::log >),
	t_mapFkts::value_type("log1p", fkt_math_1arg< ::log1p >),
	t_mapFkts::value_type("log10", fkt_math_1arg< ::log10 >),
	t_mapFkts::value_type("log2", fkt_math_1arg< ::log2 >),
	t_mapFkts::value_type("logb", fkt_math_1arg< ::logb >),
	t_mapFkts::value_type("pow", fkt_math_2args< ::pow >),
		
	t_mapFkts::value_type("sin", fkt_math_1arg< ::sin >),
	t_mapFkts::value_type("cos", fkt_math_1arg< ::cos >),
	t_mapFkts::value_type("tan", fkt_math_1arg< ::tan >),
	t_mapFkts::value_type("asin", fkt_math_1arg< ::asin >),
	t_mapFkts::value_type("acos", fkt_math_1arg< ::acos >),
	t_mapFkts::value_type("atan", fkt_math_1arg< ::atan >),
	t_mapFkts::value_type("atan2", fkt_math_2args< ::atan2 >),
	t_mapFkts::value_type("hypot", fkt_math_2args< ::hypot >),
	
	t_mapFkts::value_type("sinh", fkt_math_1arg< ::sinh >),
	t_mapFkts::value_type("cosh", fkt_math_1arg< ::cosh >),
	t_mapFkts::value_type("tanh", fkt_math_1arg< ::tanh >),	
	t_mapFkts::value_type("asinh", fkt_math_1arg< ::asinh >),
	t_mapFkts::value_type("acosh", fkt_math_1arg< ::acosh >),
	t_mapFkts::value_type("atanh", fkt_math_1arg< ::atanh >),		
	
	t_mapFkts::value_type("erf", fkt_math_1arg< ::erf >),
	t_mapFkts::value_type("erfc", fkt_math_1arg< ::erfc >),
	t_mapFkts::value_type("tgamma", fkt_math_1arg< ::tgamma >),
	t_mapFkts::value_type("lgamma", fkt_math_1arg< ::lgamma >),
	
	t_mapFkts::value_type("round", fkt_math_1arg< ::round >),
	t_mapFkts::value_type("trunc", fkt_math_1arg< ::trunc >),
	t_mapFkts::value_type("rint", fkt_math_1arg< ::rint >),
	t_mapFkts::value_type("nearbyint", fkt_math_1arg< ::nearbyint >),
	t_mapFkts::value_type("fmod", fkt_math_2args< ::fmod >),
	t_mapFkts::value_type("nextafter", fkt_math_2args< ::nextafter >),
	//t_mapFkts::value_type("nexttoward", fkt_math_2args< ::nexttoward >),
	t_mapFkts::value_type("ceil", fkt_math_1arg< ::ceil >),
	t_mapFkts::value_type("floor", fkt_math_1arg< ::floor >),
	t_mapFkts::value_type("abs", fkt_math_abs),
	t_mapFkts::value_type("max", fkt_math_for_every<MATH_MAX>),
	t_mapFkts::value_type("min", fkt_math_for_every<MATH_MIN>),
	t_mapFkts::value_type("fdim", fkt_math_2args< ::fdim >),
	t_mapFkts::value_type("remainder", fkt_math_2args< ::remainder >),

	// fft
	t_mapFkts::value_type("fft", fkt_fft),
	t_mapFkts::value_type("ifft", fkt_ifft),

	// arrays
	t_mapFkts::value_type("linspace", fkt_linspace),
	t_mapFkts::value_type("logspace", fkt_logspace),

	// dat files
	t_mapFkts::value_type("loadtxt", fkt_loadtxt),
	t_mapFkts::value_type("savetxt", fkt_savetxt),
};

extern Symbol* ext_call(const std::string& strFkt,
						const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	t_mapFkts::iterator iter = g_mapFkts.find(strFkt);
	if(iter == g_mapFkts.end())
	{
		std::cerr << linenr("Error", info)
					<< "Tried to call unknown function \""
					<< strFkt << "\"."
					<< std::endl;
		return 0;
	}
	
	return (*iter).second(vecSyms, info, pSymTab);
}
