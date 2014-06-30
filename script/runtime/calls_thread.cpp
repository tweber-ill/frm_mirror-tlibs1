/*
 * external thread functions
 * @author tweber
 * @date dec 2013
 */

#include "../types.h"
#include "../helper/flags.h"
#include "../helper/string.h"
#include "calls_thread.h"
#include "../calls.h"
#include <thread>
#include <wait.h>
#include <cstdlib>

static inline Symbol* fkt_exec(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	t_string strExec;

	for(Symbol *pSym : vecSyms)
		if(pSym)
		{
			strExec += pSym->print();
			strExec += T_STR" ";
		}

	bool bOk = 0;
	//G_COUT << "Executing " << strExec << std::endl;
	//int iRet = system(strExec.c_str());

	std::string _strExec = WSTR_TO_STR(strExec);
	FILE *pPipe = (FILE*)::my_popen(_strExec.c_str(), "w");

//	fflush(pPipe);
	if(pPipe)
	{
		bOk = 1;
		int iRet = ::my_pclose(pPipe);
		if(iRet == -1)
		{
			bOk = 0;
		}
		else
		{
			int iExitCode = int(char(WEXITSTATUS(iRet)));
			//G_COUT << "Exit code: " << iExitCode << std::endl;
			bOk = (iExitCode==0);
		}
	}

	return new SymbolInt(bOk);
}

static inline Symbol* fkt_exit(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	int iStatus = 0;
	if(vecSyms.size()>=1 && vecSyms[0])
		iStatus = vecSyms[0]->GetValInt();

	std::exit(iStatus);		// TODO: change to less brutal way to exit
	return 0;
}

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

	SymbolTable *pTable = new SymbolTable();
	SymbolArray arrArgs;
	arrArgs.SetDontDel(1);

	pinfo->pmutexInterpreter->lock();
		const NodeFunction* pThreadFunc = (NodeFunction*)pFunc/*->clone()*/;
		ParseInfo *pinfo2 = new ParseInfo(*pinfo);	// threads cannot share the same bWantReturn etc.
		pinfo2->bDestroyParseInfo = 0;
		if(pvecSyms) arrArgs.GetArr() = *pvecSyms;
	pinfo->pmutexInterpreter->unlock();

	pTable->InsertSymbol(T_STR"<args>", &arrArgs);
	Symbol* pRet = pThreadFunc->eval(*pinfo2, pTable);
	pTable->RemoveSymbolNoDelete(T_STR"<args>");

	if(pTable) delete pTable;
	if(pRet) delete pRet;
	if(pvecSyms) delete_symbols(pvecSyms);
	//if(pThreadFunc) delete pThreadFunc;
	if(pinfo2) delete pinfo2;
}

static Symbol* fkt_thread(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		G_CERR << linenr(T_STR"Error", info) << "Need thread proc identifier." << std::endl;
		return 0;
	}

	Symbol* _pSymIdent = vecSyms[0];
	if(_pSymIdent->GetType() != SYMBOL_STRING)
	{
		G_CERR << linenr(T_STR"Error", info) << "Thread proc identifier needs to be a string." << std::endl;
		return 0;
	}

	SymbolString *pSymIdent = (SymbolString*)_pSymIdent;
	const t_string& strIdent = pSymIdent->GetVal();


	NodeFunction* pFunc = info.GetFunction(strIdent);
	if(pFunc == 0)
	{
		G_CERR << linenr(T_STR"Error", info) << "Thread proc \"" << strIdent << "\" not defined." << std::endl;
		return 0;
	}

	std::vector<Symbol*>* vecThreadSymsClone = clone_symbols(&vecSyms, 1);
	std::thread* pThread = new std::thread(::thread_proc, pFunc, &info, vecThreadSymsClone);
	t_int iHandle = info.phandles->AddHandle(new HandleThread(pThread));

	return new SymbolInt(iHandle);
}

static Symbol* fkt_thread_hwcount(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	unsigned int iNumThreads = std::thread::hardware_concurrency();
	if(iNumThreads == 0)
		iNumThreads = 1;

	return new SymbolInt(t_int(iNumThreads));
}

static Symbol* fkt_mutex(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	std::mutex* pMutex = new std::mutex;
	t_int iHandle = info.phandles->AddHandle(new HandleMutex(pMutex));

	return new SymbolInt(iHandle);
}

static Symbol* fkt_begin_critical(const std::vector<Symbol*>& vecSyms,
								ParseInfo& info, SymbolTable* pSymTab)
{
	// no argument given: lock global mutex
	if(vecSyms.size() == 0)
	{
		info.pmutexGlobal->lock();
	}
	else
	{
		for(Symbol* pSym : vecSyms)
		{
			if(pSym->GetType() == SYMBOL_ARRAY)
				fkt_begin_critical(((SymbolArray*)pSym)->GetArr(), info, pSymTab);
			else if(pSym->GetType() == SYMBOL_INT)
			{
				int iHandle = pSym->GetValInt();
				Handle *pHandle = info.phandles->GetHandle(iHandle);

				if(pHandle==0 || pHandle->GetType()!=HANDLE_MUTEX)
				{
					G_CERR << linenr(T_STR"Error", info) << "Handle (" << iHandle << ") does not exist"
							 << " or is not a mutex handle." << std::endl;
					continue;
				}

				std::mutex *pMutex = ((HandleMutex*)pHandle)->GetInternalHandle();
				pMutex->lock();
			}
			else
			{
				G_CERR << linenr(T_STR"Error", info) << "Invalid mutex handle: "
						<< pSym->print() << std::endl;
			}
		}
	}

	return 0;
}

static Symbol* fkt_end_critical(const std::vector<Symbol*>& vecSyms,
								ParseInfo& info, SymbolTable* pSymTab)
{
	// no argument given: unlock global mutex
	if(vecSyms.size() == 0)
	{
		info.pmutexGlobal->unlock();
	}
	else
	{
		for(Symbol* pSym : vecSyms)
		{
			if(pSym->GetType() == SYMBOL_ARRAY)
				fkt_begin_critical(((SymbolArray*)pSym)->GetArr(), info, pSymTab);
			else if(pSym->GetType() == SYMBOL_INT)
			{
				int iHandle = pSym->GetValInt();
				Handle *pHandle = info.phandles->GetHandle(iHandle);

				if(pHandle==0 || pHandle->GetType()!=HANDLE_MUTEX)
				{
					G_CERR << linenr(T_STR"Error", info) << "Handle (" << iHandle << ") does not exist"
							 << " or is not a mutex handle." << std::endl;
					continue;
				}

				std::mutex *pMutex = ((HandleMutex*)pHandle)->GetInternalHandle();
				pMutex->unlock();
			}
			else
			{
				G_CERR << linenr(T_STR"Error", info) << "Invalid mutex handle: "
						<< pSym->print() << std::endl;
			}
		}
	}

	return 0;
}


static Symbol* fkt_thread_join(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		G_CERR << linenr(T_STR"Error", info) << "join needs at least one argument." << std::endl;
		return 0;
	}

	for(Symbol* pSym : vecSyms)
	{
		if(pSym == 0) continue;

		if(pSym->GetType() == SYMBOL_ARRAY)
		{
			return fkt_thread_join(((SymbolArray*)pSym)->GetArr(), info, pSymTab);
		}

		if(pSym->GetType() != SYMBOL_INT)
		{
			G_CERR << linenr(T_STR"Error", info) << "join needs thread handles." << std::endl;
			continue;
		}

		t_int iHandle = ((SymbolInt*)pSym)->GetVal();
		Handle *pHandle = info.phandles->GetHandle(iHandle);

		if(pHandle==0 || pHandle->GetType()!=HANDLE_THREAD)
		{
			G_CERR << linenr(T_STR"Error", info) << "Handle (" << iHandle << ") does not exist"
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
		G_CERR << linenr(T_STR"Error", info) << "nthread needs at least 3 arguments: N, func, arg." << std::endl;
		return 0;
	}

	Symbol* _pSymN = vecSyms[0];
	if(_pSymN->GetType() != SYMBOL_INT)
	{
		G_CERR << linenr(T_STR"Error", info) << "Number of threads has to be integer." << std::endl;
		return 0;
	}

	SymbolInt *pSymN = (SymbolInt*)_pSymN;
	t_int iNumThreads = pSymN->GetVal();



	Symbol* _pSymIdent = vecSyms[1];
	if(_pSymIdent->GetType() != SYMBOL_STRING)
	{
		G_CERR << linenr(T_STR"Error", info) << "Thread proc identifier needs to be a string." << std::endl;
		return 0;
	}

	SymbolString *pSymIdent = (SymbolString*)_pSymIdent;
	const t_string& strIdent = pSymIdent->GetVal();



	Symbol* _pSymArr = vecSyms[2];
	if(_pSymArr->GetType() != SYMBOL_ARRAY)
	{
		G_CERR << linenr(T_STR"Error", info) << "Thread arg has to be an array." << std::endl;
		return 0;
	}

	SymbolArray *pSymArr = (SymbolArray*)_pSymArr;
	const std::vector<Symbol*>& vecArr = pSymArr->GetArr();



	NodeFunction* pFunc = info.GetFunction(strIdent);
	if(pFunc == 0)
	{
		G_CERR << linenr(T_STR"Error", info) << "Thread proc \"" << strIdent << "\" not defined." << std::endl;
		return 0;
	}





	if(iNumThreads > vecArr.size())
	{
		iNumThreads = vecArr.size();
		G_CERR << linenr(T_STR"Warning", info) << "More threads requested in nthread than necessary, "
						  << "reducing to array size (" << iNumThreads << ")."
						  << std::endl;
	}


	std::vector<SymbolArray*> vecSymArrays;
	vecSymArrays.resize(iNumThreads);

	t_int iCurTh = 0;
	for(Symbol* pThisSym : vecArr)
	{
		if(!vecSymArrays[iCurTh])
			vecSymArrays[iCurTh] = new SymbolArray();

		vecSymArrays[iCurTh]->GetArr().push_back(pThisSym->clone());

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
		t_int iHandle = info.phandles->AddHandle(new HandleThread(pCurThread));
		SymbolInt *pSymThreadHandle = new SymbolInt(iHandle);

		pArrThreads->GetArr().push_back(pSymThreadHandle);
	}

	return pArrThreads;
}

// --------------------------------------------------------------------------------


extern void init_ext_thread_calls()
{
	t_mapFkts mapFkts =
	{
		// threads
		t_mapFkts::value_type(T_STR"thread", fkt_thread),
		t_mapFkts::value_type(T_STR"nthread", fkt_nthread),
		t_mapFkts::value_type(T_STR"thread_hwcount", fkt_thread_hwcount),
		t_mapFkts::value_type(T_STR"join", fkt_thread_join),
		t_mapFkts::value_type(T_STR"mutex", fkt_mutex),
		t_mapFkts::value_type(T_STR"begin_critical", fkt_begin_critical),
		t_mapFkts::value_type(T_STR"end_critical", fkt_end_critical),

		// processes
		t_mapFkts::value_type(T_STR"exec", fkt_exec),
		t_mapFkts::value_type(T_STR"exit", fkt_exit),
	};

	add_ext_calls(mapFkts);
}
