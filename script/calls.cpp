/*
 * external functions
 * @author tweber
 * @date 2013-2014
 */

#include "helper/flags.h"
#include "calls.h"
#include "parseobj.h"
#include "script_helper.h"
#include "helper/string.h"

#include <sstream>
#include <iostream>
#include <map>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <utility>
#include <initializer_list>

#include <ratio>
#include <chrono>
#include <thread>

//#include <regex>

extern t_string linenr(const t_string& strErr, const ParseInfo &info)
{
	if(info.pCurCaller)
		return info.pCurCaller->linenr(strErr, info);
	return strErr + T_STR": ";
}

extern bool check_args(ParseInfo& info, 
			const std::vector<Symbol*>& vecSyms, 
			const std::initializer_list<unsigned int>& lstTypes, 
			const std::initializer_list<bool> &lstOptional, 
			const char* pcFkt, const char* pcErr=0)
{
	const std::size_t iSyms = vecSyms.size();
	const std::size_t iTypes = lstTypes.size();
	const std::size_t iTotalOpt = lstOptional.size();
	std::size_t iCompulsory = 0;
	for(bool bOpt : lstOptional)
		if(!bOpt) ++iCompulsory;

	if(iSyms < iCompulsory)
	{
		G_CERR << linenr(T_STR"Error", info) 
				<< "Function \"" << pcFkt << "\""
				<< " requires " << iCompulsory
				<< " arguments, but only "
				<< iSyms << " were given.";
		if(pcErr) G_CERR << " " << pcErr;
		G_CERR << std::endl;
		return 0;
	}


	std::vector<Symbol*>::const_iterator iterSym = vecSyms.begin();
	std::initializer_list<unsigned int>::iterator iterTypes = lstTypes.begin();
	std::initializer_list<bool>::iterator iterOptional = lstOptional.begin();

	std::size_t iCurSym = 0;
	for(; iterSym!=vecSyms.end(); ++iterSym, ++iterTypes, ++iterOptional, ++iCurSym)
	{
		// ignore remaining symbols
		if(iCurSym >= iTypes || iCurSym >= iTotalOpt)
			break;

		if(!*iterSym)
		{
			G_CERR << linenr(T_STR"Error", info)
					<< "Function \"" << pcFkt << "\""
					<< " argument " << (iCurSym+1) 
					<< " is invalid.";
			if(pcErr) G_CERR << " " << pcErr;
			G_CERR << std::endl;

			return 0;
		}

		if(!(*iterTypes & (*iterSym)->GetType()))
		{
			G_CERR << linenr(T_STR"Error", info)
					<< "Function \"" << pcFkt << "\""
					<< " argument " << (iCurSym+1)
					<< " has wrong type.";
			if(pcErr) G_CERR << " " << pcErr;
			G_CERR << std::endl;

			return 0;
		}
	}

	return 1;
}

static Symbol* fkt_version(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	extern const t_char* g_pcVersion;
	return new SymbolString(g_pcVersion);
}

static Symbol* fkt_is_valid(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	for(const Symbol* pSym : vecSyms)
	{
		if(!pSym)
			return new SymbolInt(0);
	}
	return new SymbolInt(1);
}

static Symbol* fkt_null(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	return 0;
}

// call("fkt", [arg1, arg2, ...]);
static Symbol* fkt_call(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(!check_args(info, vecSyms, {SYMBOL_STRING, SYMBOL_ARRAY}, {0,1}, "call", "Function name needed."))
		return 0;

	const t_string& strFkt = ((SymbolString*)vecSyms[0])->GetVal();
	NodeFunction *pFkt = info.GetFunction(strFkt);

	SymbolArray *pArr = 0;
	if(vecSyms.size()>=2 && vecSyms[1]->GetType()==SYMBOL_ARRAY)
		pArr = (SymbolArray*)vecSyms[1];

	Symbol *pSymRet = 0;
	if(pFkt)	// user function
	{
		pSymTab->InsertSymbol(T_STR"<args>", pArr);
		pSymRet = pFkt->eval(info, pSymTab);
		info.bWantReturn = 0;
		pSymTab->RemoveSymbolNoDelete(T_STR"<args>");
	}
	else		// system function
	{
		pSymRet = ext_call(strFkt, pArr ? pArr->GetArr() : std::vector<Symbol*>(), info, pSymTab);
	}

	return pSymRet;
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

static Symbol* fkt_double_vec(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(!check_args(info, vecSyms, {SYMBOL_ARRAY}, {0}, "real_vec"))
		return 0;

	SymbolArray* pSymRet = new SymbolArray();
	pSymRet->GetArr().reserve(vecSyms.size());

	for(Symbol *pSym : ((SymbolArray*)vecSyms[0])->GetArr())
	{
		Symbol *pSymCast = 0;

		if(pSym->GetType() == SYMBOL_ARRAY)
			pSymCast = fkt_double_vec(((SymbolArray*)pSym)->GetArr(), info, pSymTab);
		else
		{
			std::vector<Symbol*> vecSymTmp{pSym};
			pSymCast = fkt_double(vecSymTmp, info, pSymTab);
		}

		if(pSymCast)
			pSymRet->GetArr().push_back(pSymCast);
	}

	return pSymRet;
}

static Symbol* fkt_str(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	SymbolString *pSymRet = new SymbolString;
	for(const Symbol *pSym : vecSyms)
	{
		if(!pSym)
			continue;
		pSymRet->GetVal() += pSym->print();
	}
	return pSymRet;
}

static Symbol* fkt_output(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	t_ostream& ostr = G_COUT;

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
	G_COUT << std::endl;
	return 0;
}

static Symbol* fkt_input(const std::vector<Symbol*>& vecSyms,
		ParseInfo& info, SymbolTable* pSymTab)
{
	fkt_output(vecSyms, info, pSymTab);

	SymbolString* pSymStr = new SymbolString();
	std::getline(G_CIN, pSymStr->GetVal());
	return pSymStr;
}

static Symbol* fkt_sleep(const std::vector<Symbol*>& vecSyms,
		ParseInfo& info, SymbolTable* pSymTab)
{
	if(!check_args(info, vecSyms, {SYMBOL_DOUBLE|SYMBOL_INT}, {0,1}, "sleep", "Number of ms needed."))
		return 0;

	t_real dMillis = vecSyms[0]->GetValDouble();

	typedef std::chrono::duration<t_real, std::milli> t_millis;
	t_millis delay(dMillis);
	std::this_thread::sleep_for(delay);
	return 0;
}


extern int yyparse(void*);
static bool _import_file(const t_string& strFile, ParseInfo& info, SymbolTable* pSymTab)
{
	ParseInfo::t_mods::iterator iterMod = info.pmapModules->find(strFile);
	if(iterMod != info.pmapModules->end())
	{
		//G_CERR << "Warning: Module \"" << strFile << "\" already loaded." << std::endl;
		return 0;
	}

	t_string _strFile = WSTR_TO_STR(strFile);
	t_char* pcInput = load_file(_strFile.c_str());
	if(!pcInput)
		return 0;

	ParseObj par;
	par.strCurFile = strFile;
	par.pLexer = new Lexer(pcInput, strFile.c_str());

	delete[] pcInput;
	pcInput = 0;

	if(!par.pLexer->IsOk())
	{
		G_CERR << linenr(T_STR"Error", info) << "Lexer returned with errors." << std::endl;
		return 0;
	}

	int iParseRet = yyparse(&par);

	delete par.pLexer;
	par.pLexer = 0;

	if(iParseRet != 0)
	{
		G_CERR << linenr(T_STR"Error", info) << "Parser returned with error code " << iParseRet << "." << std::endl;
		return 0;
	}

	Node *pRoot = par.pRoot;
	info.pmapModules->insert(ParseInfo::t_mods::value_type(strFile, pRoot));
	info.strExecFkt = T_STR"";
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
			const t_string& strFile = ((SymbolString*)pSym)->GetVal();
			bool bOk = _import_file(strFile, info, pSymTab);
		}

	return 0;
}

static Symbol* fkt_has_var(const std::vector<Symbol*>& vecSyms,
			ParseInfo& info, SymbolTable* pSymTab)
{
	if(!check_args(info, vecSyms, {SYMBOL_STRING}, {0}, "has_var", "Symbol name needed."))
		return 0;

	const t_string& strVar = ((SymbolString*)vecSyms[0])->GetVal();
	bool bHasVar = 0;

	// check local variables
	if(pSymTab->GetSymbol(strVar))
		bHasVar = 1;

	// check global variables
	if(info.pGlobalSyms->GetSymbol(strVar))
		bHasVar = 1;

	return new SymbolInt(bHasVar);
}

// register a variable in the symbol table
static Symbol* fkt_register_var(const std::vector<Symbol*>& vecSyms,
			ParseInfo& info, SymbolTable* pSymTab)
{
	bool bUseGlobal = 0;

	if(vecSyms.size()<1 || vecSyms.size()>2)
	{
		G_CERR << linenr(T_STR"Error", info)
			<< "Need either a symbol name and a symbol or a map for register_var."
			<< std::endl;
		return 0;
	}

	if(vecSyms.size() == 2)
	{
		const t_string strVar = vecSyms[0]->print();
		Symbol* pVar = vecSyms[1]->clone();

		SymbolTable *pThisSymTab = pSymTab;
		if(bUseGlobal)
			pThisSymTab = info.pGlobalSyms;
		pThisSymTab->InsertSymbol(strVar, pVar);
	}
	else if(vecSyms.size() == 1)
	{
		if(vecSyms[0]==0 || vecSyms[0]->GetType() != SYMBOL_MAP)
		{
			G_CERR << linenr(T_STR"Error", info)
				<< "register_var needs a map."
				<< std::endl;
			return 0;
		}


		SymbolMap::t_map& varmap = ((SymbolMap*)vecSyms[0])->GetMap();
		for(SymbolMap::t_map::value_type& val : varmap)
		{
			SymbolString symKey(val.first);
			Symbol *pSymVal = val.second;

			std::vector<Symbol*> vecDummy;
			vecDummy.resize(2);

			vecDummy[0] = &symKey;
			vecDummy[1] = pSymVal;

			fkt_register_var(vecDummy, info, pSymTab);
		}
	}

	return 0;
}

static Symbol* fkt_typeof(const std::vector<Symbol*>& vecSyms,
			ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=1)
	{
		G_CERR << linenr(T_STR"Error", info) << "typeof takes exactly one argument." << std::endl;
		return 0;
	}

	Symbol *pSymbol = vecSyms[0];
	if(!pSymbol)
	{
		G_CERR << linenr(T_STR"Error", info) << "Invalid argument for typeof." << std::endl;
		return 0;
	}

	SymbolString *pType = new SymbolString(pSymbol->GetTypeName().c_str());
	return pType;
}

static Symbol* fkt_setprec(const std::vector<Symbol*>& vecSyms,
			ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=1)
	{
		G_CERR << linenr(T_STR"Error", info) << "set_prec takes exactly one argument." << std::endl;
		return 0;
	}

	Symbol *pSymbol = vecSyms[0];
	if(!pSymbol)
	{
		G_CERR << linenr(T_STR"Error", info) << "Invalid argument for set_prec." << std::endl;
		return 0;
	}

	int iPrec = pSymbol->GetValInt();

	if(iPrec >= 0)
		SymbolDouble::SetPrec(iPrec);
	else
		SymbolDouble::SetPrec(SymbolDouble::GetDefPrec());

	return 0;
}


// --------------------------------------------------------------------------------
// map

static Symbol* fkt_map(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()>=1 && vecSyms[0]->GetType()==SYMBOL_MAP)
		return vecSyms[0]->clone();

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

	if(vecSyms.size()>=1 && vecSyms[0]->GetType()==SYMBOL_ARRAY)
		return vecSyms[0]->clone();


	Symbol *pSymSize = vecSyms[0];
	if(pSymSize->GetType() != SYMBOL_INT)
	{
		G_CERR << linenr(T_STR"Error", info) << "\"num\" in vec(num, val=0) has to be integer." << std::endl;
		return 0;
	}

	t_int iVal = ((SymbolInt*)pSymSize)->GetVal();
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
	pSymRet->GetArr().reserve(iVal);
	for(t_int i=0; i<iVal; ++i)
		pSymRet->GetArr().push_back(pSymVal->clone());

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
		G_CERR << linenr(T_STR"Error", info) << "vec_size needs one argument." << std::endl;
		return 0;
	}

	Symbol *pSymArr = vecSyms[0];
	SymbolInt *pSymRet = new SymbolInt(0);

	if(pSymArr->GetType() == SYMBOL_ARRAY)
		pSymRet->SetVal(((SymbolArray*)pSymArr)->GetArr().size());
	else if(pSymArr->GetType() == SYMBOL_STRING)
		pSymRet->SetVal(((SymbolString*)pSymArr)->GetVal().length());
	else if(pSymArr->GetType() == SYMBOL_MAP)
		pSymRet->SetVal(((SymbolMap*)pSymArr)->GetMap().size());
	else
		G_CERR << linenr(T_STR"Error", info) << "vec_size needs a vector type argument." << std::endl;

	return pSymRet;
}

static int pos_in_string(const SymbolString* pSymStr, const SymbolString *pSym)
{
	const t_string& str = pSymStr->GetVal();
	const t_string& strToFind = pSym->GetVal();

	std::size_t pos = str.find(strToFind);
	if(pos == t_string::npos)
		return -1;

	return int(pos);
}

static int pos_in_array(const SymbolArray* pSymArr, const Symbol *pSym)
{
	unsigned int iPos;
	bool bFound = 0;

	for(iPos=0; iPos<pSymArr->GetArr().size(); ++iPos)
	{
		const Symbol& sym = *pSymArr->GetArr()[iPos];
		bool bSym0Scalar = sym.IsScalar() || sym.GetType()==SYMBOL_STRING;
		bool bSym1Scalar = pSym->IsScalar() || pSym->GetType()==SYMBOL_STRING;

		if(bSym0Scalar != bSym1Scalar)
			continue;
		else if(bSym0Scalar && bSym1Scalar)
		{
			Symbol *pEqu = Node::Op(&sym, pSym, NODE_LOG_EQ);
			if(!pEqu) continue;
			int iEqu = pEqu->GetValInt();
			delete pEqu;

			if(iEqu)
			{
				bFound = 1;
				break;
			}
		}
		else
		{
			G_CERR << "Error: Array compare not yet implemented." << std::endl;
			continue;

			// TODO: array/map compare
		}
	}

	if(!bFound) return -1;
	return int(iPos);
}

static Symbol* fkt_find(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()<2)
	{
		G_CERR << linenr(T_STR"Error", info) 
			<< "Find needs two arguments." << std::endl;
		return 0;
	}

	const Symbol *pContainer = vecSyms[0];
	const Symbol *pContainee = vecSyms[1];

	if(pContainer->GetType() == SYMBOL_ARRAY)
	{
		int iPos = pos_in_array((SymbolArray*)pContainer, pContainee);
		return new SymbolInt(iPos);
	}
	else if(pContainer->GetType() == SYMBOL_MAP)
	{
		G_CERR << linenr(T_STR"Error", info) 
			<< "Find not yet implemented for map." << std::endl;

		// TODO: Implement and also adapt fkt_contains
	}
	else if(pContainer->GetType() == SYMBOL_STRING)
	{
		if(pContainee->GetType() != SYMBOL_STRING)
		{
			G_CERR << linenr(T_STR"Error", info)
				<< "Second argument to find has to be of string type."
				<< std::endl;
		}

		int iPos = pos_in_string((SymbolString*)pContainer, (SymbolString*)pContainee);
		return new SymbolInt(iPos);
	}

	return 0;
}

static Symbol* fkt_contains(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	int iIdx = -1;
	Symbol *pFind = fkt_find(vecSyms, info, pSymTab);
	if(pFind)
	{
		iIdx = pFind->GetValInt();
		delete pFind;
	}

	return new SymbolInt(iIdx>=0);
}

static Symbol* fkt_cur_iter(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		G_CERR << linenr(T_STR"Error", info) << "cur_iter takes exactly one argument." << std::endl;
		return 0;
	}

	const t_string& strIdent = vecSyms[0]->GetIdent();
	//std::cout << "Ident: " << strIdent << std::endl;
	if(strIdent == T_STR"")
	{
		G_CERR << linenr(T_STR"Error", info) << "No identifier given for cur_iter." << std::endl;
		return 0;
	}


	t_string strIter = T_STR"<cur_iter_" + strIdent + T_STR">";

	Symbol* pSymIter = pSymTab->GetSymbol(strIter);
	//pSymTab->print();
	if(!pSymIter || pSymIter->GetType()!=SYMBOL_INT)
	{
		G_CERR << linenr(T_STR"Error", info) << "cur_iter could not determine iteration index \""
					<< strIter << "\"." << std::endl;
		return 0;
	}

	return pSymIter;
}


typedef std::tuple<const Symbol*, unsigned int> t_symtup;

static void _sortarr(std::vector<t_symtup>& vec)
{
	auto comp = [](const t_symtup& tup1, const t_symtup& tup2) -> bool
	{
		const Symbol* pSym0 = std::get<0>(tup1);
		const Symbol* pSym1 = std::get<0>(tup2);

		return pSym0->IsLessThan(*pSym1);
	};

	std::sort(vec.begin(), vec.end(), comp);
}

static void _rearrangearr(std::vector<Symbol*>& vecSyms, const std::vector<unsigned int>& vecIdx)
{
	std::vector<Symbol*> vecTmp = vecSyms;
	for(unsigned int i=0; i<vecSyms.size(); ++i)
		vecSyms[i] = vecTmp[vecIdx[i]];
}

static Symbol* fkt_sort(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1 || vecSyms[0]->GetType()!=SYMBOL_ARRAY)
	{
		G_CERR << linenr(T_STR"Error", info) 
			<< "Arguments to sort have to be arrays." 
			<< std::endl;
		return 0;
	}

	const SymbolArray* pArr = (SymbolArray*)vecSyms[0];
	const unsigned int iArrSize = pArr->GetArr().size();

	std::vector<t_symtup> vecTups;
	vecTups.reserve(iArrSize);
	for(unsigned int iElem=0; iElem<iArrSize; ++iElem)
		vecTups.push_back(t_symtup(pArr->GetArr()[iElem], iElem));

	_sortarr(vecTups);

	SymbolArray* pArrRet = new SymbolArray();
	pArrRet->GetArr().reserve(iArrSize);

	std::vector<unsigned int> vecSortedIndices;
	vecSortedIndices.reserve(iArrSize);

	for(unsigned int iElem=0; iElem<iArrSize; ++iElem)
	{
		pArrRet->GetArr().push_back(std::get<0>(vecTups[iElem])->clone());
		vecSortedIndices.push_back(std::get<1>(vecTups[iElem]));
	}

	// no other arguments to sort
	if(vecSyms.size() == 1)
		return pArrRet;



	// sort other arrays in the same way
	SymbolArray *pArrArr = new SymbolArray();
	pArrArr->GetArr().reserve(vecSyms.size());
	pArrArr->GetArr().push_back(pArrRet);

	for(unsigned int iElem=1; iElem<vecSyms.size(); ++iElem)
	{
		const Symbol* pSym = vecSyms[iElem];
		if(pSym->GetType() != SYMBOL_ARRAY)
		{
			G_CERR << linenr(T_STR"Error", info) 
				<< "Arguments to sort have to be arrays." 
				<< std::endl;
			continue;
		}

		if(((SymbolArray*)pSym)->GetArr().size() != vecSortedIndices.size())
		{
			G_CERR << linenr(T_STR"Error", info)
				<< "Array size mismatch in sort."
				<< std::endl;
			continue;
		}

		SymbolArray *pNextArr = (SymbolArray*)pSym->clone();
		_rearrangearr(pNextArr->GetArr(), vecSortedIndices);

		pArrArr->GetArr().push_back(pNextArr);
	}

	return pArrArr;
}


static Symbol* fkt_zip(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	SymbolArray* pArrRet = new SymbolArray();
	std::vector<Symbol*>* pVec = &pArrRet->GetArr();

	bool bFirst = 1;
	unsigned int iSize = 0;
	for(Symbol* pSym : vecSyms)
	{
		if(pSym->GetType() != SYMBOL_ARRAY)
		{
			G_CERR << linenr(T_STR"Error", info) << "Symbol \"" << pSym->GetName()
					<< "\" is of type " << pSym->GetTypeName()
					<< ", but zip needs vectors." << std::endl;
			continue;
		}

		std::vector<Symbol*>& curSym = ((SymbolArray*)pSym)->GetArr();
		if(bFirst)
		{
			iSize = curSym.size();
			pVec->reserve(iSize);
			for(unsigned int iCurSym=0; iCurSym<iSize; ++iCurSym)
				pVec->push_back(new SymbolArray());

			bFirst = 0;
		}

		if(curSym.size() < iSize)
		{
			for(unsigned int i=0; i<iSize-curSym.size(); ++i)
			{
				delete *pVec->rbegin();
				pVec->pop_back();
			}

			iSize = curSym.size();
		}

		for(unsigned int i=0; i<iSize; ++i)
			((SymbolArray*)(*pVec)[i])->GetArr().push_back(curSym[i]->clone());
	}

	return pArrRet;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// string operations
static Symbol* fkt_trim(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=1)
	{
		G_CERR << linenr(T_STR"Error", info)
				<< "Need one argument for trim." << std::endl;
		return 0;
	}

	if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray *pArr = new SymbolArray();
		pArr->GetArr().reserve(((SymbolArray*)vecSyms[0])->GetArr().size());

		for(Symbol* pSymArr : ((SymbolArray*)vecSyms[0])->GetArr())
		{
			std::vector<Symbol*> vecDummy = { pSymArr };
			pArr->GetArr().push_back(fkt_trim(vecDummy, info, pSymTab));
		}

		return pArr;
	}
	else if(vecSyms[0]->GetType() == SYMBOL_STRING)
	{
		t_string str = ((SymbolString*)vecSyms[0])->GetVal();
		::trim(str);
		return new SymbolString(str);
	}

	// simply copy non-string arguments
	//G_CERR << linenr(T_STR"Warning", info)
	//		<< "Called trim with invalid argument." << std::endl;
	return vecSyms[0]->clone();
}

// split(string, delim)
static Symbol* fkt_split(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2 || vecSyms[0]->GetType()!=SYMBOL_STRING || vecSyms[1]->GetType()!=SYMBOL_STRING)
	{
		G_CERR << linenr(T_STR"Error", info)
			<< "Split needs two string arguments." << std::endl;
		return 0;
	}

	const t_string& str = ((SymbolString*)vecSyms[0])->GetVal();
	const t_string& strDelim = ((SymbolString*)vecSyms[1])->GetVal();
	std::size_t iPos = str.find(strDelim);

	t_string str0 = str.substr(0, iPos);
	t_string str1;

	if(iPos != t_string::npos)
		str1 = str.substr(iPos + strDelim.length(), t_string::npos);

	SymbolArray* pSymRet = new SymbolArray();
	pSymRet->GetArr().push_back(new SymbolString(str0));
	pSymRet->GetArr().push_back(new SymbolString(str1));
	return pSymRet;
}

// tokens(string, delim)
static Symbol* fkt_tokens(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	t_string *pstrInput = 0;
	t_string strDelim = T_STR" \t\n";

	// split("Test 123")
	if(vecSyms.size() >= 1 && vecSyms[0]->GetType()==SYMBOL_STRING)
		pstrInput = &((SymbolString*)vecSyms[0])->GetVal();

	// split("Test 123", " \t\n")
	if(vecSyms.size() >= 2 && vecSyms[1]->GetType()==SYMBOL_STRING)
		strDelim = ((SymbolString*)vecSyms[1])->GetVal();

	if(!pstrInput)
	{
		G_CERR << linenr(T_STR"Error", info)
				<< "Called tokens with invalid arguments." << std::endl;
		return 0;
	}

	std::vector<t_string> vecTokens;
	::get_tokens<t_string>(*pstrInput, strDelim, vecTokens);

	SymbolArray* pArr = new SymbolArray;
	for(const t_string& strTok : vecTokens)
		pArr->GetArr().push_back(new SymbolString(strTok));

	return pArr;
}


static Symbol* fkt_replace(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=3 || vecSyms[0]->GetType()!=SYMBOL_STRING
						|| vecSyms[1]->GetType()!=SYMBOL_STRING
						|| vecSyms[2]->GetType()!=SYMBOL_STRING)
	{
		G_CERR << linenr(T_STR"Error", info)
				<< "Replace needs three string arguments." << std::endl;
		return 0;
	}

	t_string str = ((SymbolString*)vecSyms[0])->GetVal();
	const t_string& strOld = ((SymbolString*)vecSyms[1])->GetVal();
	const t_string& strNew = ((SymbolString*)vecSyms[2])->GetVal();

	find_all_and_replace(str, strOld, strNew);

	return new SymbolString(str);
}

/*
// replace_regex(str, regex, repl);
static Symbol* fkt_replace_regex(const std::vector<Symbol*>& vecSyms,
				ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=3 || vecSyms[0]->GetType()!=SYMBOL_STRING
				|| vecSyms[1]->GetType()!=SYMBOL_STRING
				|| vecSyms[2]->GetType()!=SYMBOL_STRING)
	{
		G_CERR << linenr(T_STR"Error", info)
			<< "Replace needs three string arguments." << std::endl;
		return 0;
	}

	t_string str = ((SymbolString*)vecSyms[0])->GetVal();
	const t_string& strRegex = ((SymbolString*)vecSyms[1])->GetVal();
	const t_string& strRepl = ((SymbolString*)vecSyms[2])->GetVal();

	std::string strRet;
	try
	{
		//std::cout << "Regex: " << strRegex << std::endl;
		std::basic_regex<t_char> rex(strRegex, std::regex::ECMAScript);
		strRet = std::regex_replace(str, rex, strRepl);
	}
	catch(const std::exception& ex)
	{
		G_CERR << linenr(T_STR"Error", info)
			<< "Regex evaluation failed with error: " << ex.what() 
			<< std::endl;
		return 0;
	}

	return new SymbolString(std::move(strRet));
}
*/
// --------------------------------------------------------------------------------




// --------------------------------------------------------------------------------
// map operations

static Symbol* fkt_has_key(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=2 || vecSyms[0]->GetType()!=SYMBOL_MAP
				|| vecSyms[1]->GetType()!=SYMBOL_STRING)
	{
		G_CERR << linenr(T_STR"Error", info)
			<< "Has_key needs a map and a string argument."
			<< std::endl;
		return 0;
	}

	const SymbolMap* pMap = (SymbolMap*)vecSyms[0];
	const std::string& strKey = ((SymbolString*)vecSyms[1])->GetVal();

	int bHasKey = (pMap->GetMap().find(strKey) != pMap->GetMap().end());
	return new SymbolInt(bHasKey);
}

// --------------------------------------------------------------------------------





// --------------------------------------------------------------------------------
// basic external functions
static t_mapFkts g_mapFkts =
{
	// basic stuff
	t_mapFkts::value_type(T_STR"interp_ver", fkt_version),
	t_mapFkts::value_type(T_STR"register_var", fkt_register_var),
	t_mapFkts::value_type(T_STR"is_valid", fkt_is_valid),
	t_mapFkts::value_type(T_STR"null", fkt_null),
	t_mapFkts::value_type(T_STR"call", fkt_call),

	// input/output
	t_mapFkts::value_type(T_STR"input", fkt_input),
	t_mapFkts::value_type(T_STR"output", fkt_output),
	t_mapFkts::value_type(T_STR"print", fkt_print),	// output with "\n" at the end
	t_mapFkts::value_type(T_STR"sleep", fkt_sleep),

	// modules
	t_mapFkts::value_type(T_STR"import", fkt_import),

	// symbols & casts
	t_mapFkts::value_type(T_STR"int", fkt_int),
	t_mapFkts::value_type(T_STR"real", fkt_double),
	t_mapFkts::value_type(T_STR"real_vec", fkt_double_vec),
	t_mapFkts::value_type(T_STR"str", fkt_str),
	t_mapFkts::value_type(T_STR"map", fkt_map),
	t_mapFkts::value_type(T_STR"vec", fkt_array),
	t_mapFkts::value_type(T_STR"has_var", fkt_has_var),
	t_mapFkts::value_type(T_STR"typeof", fkt_typeof),
	t_mapFkts::value_type(T_STR"set_prec", fkt_setprec),

	// string operations
	t_mapFkts::value_type(T_STR"trim", fkt_trim),
	t_mapFkts::value_type(T_STR"split", fkt_split),
	t_mapFkts::value_type(T_STR"tokens", fkt_tokens),
	t_mapFkts::value_type(T_STR"replace", fkt_replace),
	t_mapFkts::value_type(T_STR"length", fkt_array_size),
	//t_mapFkts::value_type(T_STR"replace_regex", fkt_replace_regex),


	// array operations
	t_mapFkts::value_type(T_STR"vec_size", fkt_array_size),	// deprecated, use "length" instead
	t_mapFkts::value_type(T_STR"cur_iter", fkt_cur_iter),
	t_mapFkts::value_type(T_STR"zip", fkt_zip),
	t_mapFkts::value_type(T_STR"sort", fkt_sort),

	// map/array operations
	t_mapFkts::value_type(T_STR"contains", fkt_contains),
	t_mapFkts::value_type(T_STR"has_key", fkt_has_key),
	t_mapFkts::value_type(T_STR"find", fkt_find),
};

// --------------------------------------------------------------------------------

extern Symbol* ext_call(const t_string& strFkt,
						const std::vector<Symbol*>& vecSyms,
						ParseInfo& info,
						SymbolTable* pSymTab)
{
	t_mapFkts::iterator iter = g_mapFkts.find(strFkt);
	if(iter == g_mapFkts.end())
	{
		G_CERR << linenr(T_STR"Error", info)
					<< "Tried to call unknown function \""
					<< strFkt << "\"."
					<< std::endl;
		return 0;
	}

	return (*iter).second(vecSyms, info, pSymTab);
}

// --------------------------------------------------------------------------------

extern bool add_ext_call(const t_string& strFkt, t_extcall pExtCall)
{
	bool bInserted = g_mapFkts.insert(t_mapFkts::value_type(strFkt, pExtCall)).second;
	return bInserted;
}

extern void add_ext_calls(t_mapFkts& mapFkt)
{
	g_mapFkts.insert(mapFkt.begin(), mapFkt.end());
}

extern bool has_ext_call(const t_string& strFkt)
{
	return g_mapFkts.find(strFkt) != g_mapFkts.end();
}

extern const t_mapFkts* get_ext_calls()
{
	return &g_mapFkts;
}
