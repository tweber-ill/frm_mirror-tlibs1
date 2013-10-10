/*
 * External Functions
 * @author tweber
 */

#include "calls.h"
#include <sstream>
#include <iostream>
#include <map>


static Symbol* fkt_version(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info,
							SymbolTable* pSymTab)
{
	return new SymbolString("Hermelin Interpreter, Version 0.1");
}


static Symbol* fkt_print(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	std::ostream& ostr = std::cout;
	
	for(Symbol *pSym : vecSyms)
		if(pSym)
			ostr << pSym->print();

	ostr << std::endl;
	return 0;
}

static Symbol* fkt_array(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		std::cerr << "Error: vec(num, val=0) needs at least one argument." << std::endl;
		return 0;
	}
	
	Symbol *pSymSize = vecSyms[0];
	if(pSymSize->GetType() != SYMBOL_INT)
	{
		std::cerr << "Error: \"num\" in vec(num, val=0) has to be integer." << std::endl;
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
		pSymVal->m_strName = "<const>";
		bOwnVal = 1;
	}
	
	
	SymbolArray* pSymRet = new SymbolArray;
	pSymRet->m_arr.reserve(iVal);
	for(int i=0; i<iVal; ++i)
		pSymRet->m_arr.push_back(pSymVal->clone());
	
	if(bOwnVal)
		delete pSymVal;
	
	return pSymRet;
}

static Symbol* fkt_array_size(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info,
							SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		std::cerr << "Error: vec_size(vec) needs one argument." << std::endl;
		return 0;
	}
	
	Symbol *pSymArr = vecSyms[0];
	SymbolInt *pSymRet = new SymbolInt(0);
	
	if(pSymArr->GetType() != SYMBOL_ARRAY)
	{
		std::cerr << "Warning: vec_size needs an array type argument." << std::endl;
		return pSymRet;
	}

	
	pSymRet->m_iVal = ((SymbolArray*)pSymArr)->m_arr.size();
	return pSymRet;
}

static Symbol* fkt_cur_iter(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info,
							SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		std::cerr << "Error: cur_iter takes exactly one argument." << std::endl;
		return 0;
	}

	const std::string& strIdent = vecSyms[0]->m_strIdent;
	std::string strIter = "<cur_iter_" + strIdent + ">";

	Symbol* pSymIter = pSymTab->GetSymbol(strIter);
	if(!pSymIter || pSymIter->GetType()!=SYMBOL_INT)
	{
		std::cerr << "Error: cur_iter could not determine iteration index."
					<< std::endl;
		return 0;
	}

	return pSymIter;
}

static Symbol* fkt_typeof(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info,
							SymbolTable* pSymTab)
{
	if(vecSyms.size()!=1)
	{
		std::cerr << "Error: typeof takes exactly one argument." << std::endl;
		return 0;
	}

	Symbol *pSymbol = vecSyms[0];
	if(!pSymbol)
	{
		std::cerr << "Error: Invalid argument for typename." << std::endl;
		return 0;
	}

	SymbolString *pType = new SymbolString(pSymbol->GetTypeName().c_str());
	return pType;
}


static void thread_proc(NodeFunction* pThreadFunc, ParseInfo* pinfo, std::vector<Symbol*>* pvecSyms)
{
	pThreadFunc->SetArgSyms(pvecSyms);
	Symbol* pRet = pThreadFunc->eval(*pinfo, 0);

	if(pRet) delete pRet;
	if(pvecSyms) delete pvecSyms;
	if(pThreadFunc) delete pThreadFunc;
}

static Symbol* fkt_thread(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		std::cerr << "Error: Need thread proc identifier." << std::endl;
		return 0;
	}

	Symbol* _pSymIdent = vecSyms[0];
	if(_pSymIdent->GetType() != SYMBOL_STRING)
	{
		std::cerr << "Error: Thread proc identifier needs to be a string." << std::endl;
		return 0;
	}

	SymbolString *pSymIdent = (SymbolString*)_pSymIdent;
	const std::string& strIdent = pSymIdent->m_strVal;


	NodeFunction* pFunc = info.GetFunction(strIdent);
	if(pFunc == 0)
	{
		std::cerr << "Error: Thread proc \"" << strIdent << "\" not defined." << std::endl;
		return 0;
	}

	NodeFunction* pFunc_clone = (NodeFunction*)pFunc->clone();



	std::vector<Symbol*>* vecThreadSyms = new std::vector<Symbol*>(vecSyms.size()-1);
	std::copy(vecSyms.begin()+1, vecSyms.end(), vecThreadSyms->begin());

	std::thread* pThread = new std::thread(::thread_proc, pFunc_clone, &info, vecThreadSyms);
	unsigned int iHandle = info.handles.AddHandle(new HandleThread(pThread));

	return new SymbolInt(iHandle);
}

static Symbol* fkt_thread_join(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		std::cerr << "Error: join needs at least one argument." << std::endl;
		return 0;
	}

	for(Symbol* pSym : vecSyms)
	{
		if(pSym == 0) continue;

		if(pSym->GetType() != SYMBOL_INT)
		{
			std::cerr << "Error: join needs thread handles." << std::endl;
			continue;
		}

		int iHandle = ((SymbolInt*)pSym)->m_iVal;
		Handle *pHandle = info.handles.GetHandle(iHandle);

		if(pHandle==0 || pHandle->GetType()!=HANDLE_THREAD)
		{
			std::cerr << "Error: Handle (" << iHandle << ") does not exist"
					 << " or is not a thread handle." << std::endl;
			continue;
		}

		HandleThread *pThreadHandle = (HandleThread*)pHandle;
		std::thread *pThread = pThreadHandle->GetInternalHandle();

		pThread->join();
	}
	return 0;
}




typedef std::map<std::string, Symbol*(*)(const std::vector<Symbol*>&, ParseInfo&, SymbolTable*)> t_mapFkts;
static t_mapFkts g_mapFkts =
{
	t_mapFkts::value_type("typeof", fkt_typeof),

	t_mapFkts::value_type("ver", fkt_version),
	t_mapFkts::value_type("print", fkt_print),
	
	t_mapFkts::value_type("vec", fkt_array),
	t_mapFkts::value_type("vec_size", fkt_array_size),

	t_mapFkts::value_type("cur_iter", fkt_cur_iter),

	t_mapFkts::value_type("thread", fkt_thread),
	t_mapFkts::value_type("join", fkt_thread_join)
};

extern Symbol* ext_call(const std::string& strFkt,
						const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	t_mapFkts::iterator iter = g_mapFkts.find(strFkt);
	if(iter == g_mapFkts.end())
	{
		std::cerr << "Error: Tried to call unknown function \"" 
					<< strFkt << "\"."
					<< std::endl;
		return 0;
	}
	
	return (*iter).second(vecSyms, info, pSymTab);
}
