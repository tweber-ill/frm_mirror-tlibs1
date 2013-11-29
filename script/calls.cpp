/*
 * External Functions
 * @author tweber
 */

#ifdef __CYGWIN__
	#undef __STRICT_ANSI__
#endif

#include "calls.h"
#include <sstream>
#include <iostream>
#include <map>
#include <cstdio>
#include <cmath>

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

static Symbol* fkt_exec(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
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
		if(::pclose(pPipe) == -1)
			bOk = 0;
	}
		
	return new SymbolInt(bOk);
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


// --------------------------------------------------------------------------------
// array

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


// --------------------------------------------------------------------------------




// --------------------------------------------------------------------------------
// thread

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

static Symbol* fkt_thread_hwcount(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	unsigned int iNumThreads = std::thread::hardware_concurrency();
	if(iNumThreads == 0)
		iNumThreads = 1;
	
	return new SymbolInt(iNumThreads);
}	

static Symbol* fkt_begin_critical(const std::vector<Symbol*>& vecSyms,
								ParseInfo& info,
								SymbolTable* pSymTab)
{
	info.mutexGlobal.lock();
	return 0;
}

static Symbol* fkt_end_critical(const std::vector<Symbol*>& vecSyms,
								ParseInfo& info,
								SymbolTable* pSymTab)
{
	info.mutexGlobal.unlock();
	return 0;
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
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	if(vecSyms.size() < 1)
	{
		std::cerr << "Error: fkt_math_for_every needs at least one argument" << std::endl;
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

	std::cerr << "Error: No valid arguments given for fkt_math_for_every." << std::endl;
	return 0;
}


template<double (*FKT)(double)>
static Symbol* fkt_math_1arg(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		std::cerr << "Error: fkt_math_1arg takes exactly one argument." << std::endl;
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
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2)
	{
		std::cerr << "Error: fkt_math_2args takes exactly two arguments." << std::endl;
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
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		std::cerr << "Error: abs takes exactly one argument." << std::endl;
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
	

	std::cerr << "Error: abs received unsupported symbol type." << std::endl;
	return 0;
}

// --------------------------------------------------------------------------------



typedef std::map<std::string, Symbol*(*)(const std::vector<Symbol*>&, ParseInfo&, SymbolTable*)> t_mapFkts;
static t_mapFkts g_mapFkts =
{
	t_mapFkts::value_type("typeof", fkt_typeof),
	t_mapFkts::value_type("ver", fkt_version),
	t_mapFkts::value_type("print", fkt_print),
	t_mapFkts::value_type("exec", fkt_exec),
	
	t_mapFkts::value_type("vec", fkt_array),
	t_mapFkts::value_type("vec_size", fkt_array_size),
	t_mapFkts::value_type("cur_iter", fkt_cur_iter),

	t_mapFkts::value_type("thread", fkt_thread),
	t_mapFkts::value_type("thread_hwcount", fkt_thread_hwcount),
	t_mapFkts::value_type("join", fkt_thread_join),
	t_mapFkts::value_type("begin_critical", fkt_begin_critical),
	t_mapFkts::value_type("end_critical", fkt_end_critical),
	
	
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
