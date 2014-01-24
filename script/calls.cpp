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

extern std::string linenr(const std::string& strErr, const ParseInfo &info)
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
	{
		if(!pSym)
			continue;
		pSymRet->m_strVal += pSym->print();
	}
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

static Symbol* fkt_has_var(const std::vector<Symbol*>& vecSyms,
			ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=1)
	{
		std::cerr << linenr("Error", info)
			<< "Need a symbol name for has_var."
			<< std::endl;
		return 0;
	}

	const std::string strVar = vecSyms[0]->print();
	bool bHasVar = 0;

	// check local variables
	if(pSymTab->GetSymbol(strVar))
		bHasVar = 1;

	// check global variables
	if(info.pGlobalSyms->GetSymbol(strVar))
		bHasVar = 1;

	return new SymbolInt(bHasVar);
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
		std::cerr << linenr("Error", info) << "vec_size needs one argument." << std::endl;
		return 0;
	}
	
	Symbol *pSymArr = vecSyms[0];
	SymbolInt *pSymRet = new SymbolInt(0);
	
	if(pSymArr->GetType() == SYMBOL_ARRAY)
		pSymRet->m_iVal = ((SymbolArray*)pSymArr)->m_arr.size();
	else if(pSymArr->GetType() == SYMBOL_STRING)
		pSymRet->m_iVal = ((SymbolString*)pSymArr)->m_strVal.length();
	else
		std::cerr << linenr("Error", info) << "vec_size needs a vector type argument." << std::endl;

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


static Symbol* fkt_zip(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	SymbolArray* pArrRet = new SymbolArray();
	std::vector<Symbol*>* pVec = &pArrRet->m_arr;

	bool bFirst = 1;
	unsigned int iSize = 0;
	for(Symbol* pSym : vecSyms)
	{
		if(pSym->GetType() != SYMBOL_ARRAY)
		{
			std::cerr << linenr("Error", info) << "Symbol \"" << pSym->m_strName
					<< "\" is of type " << pSym->GetTypeName()
					<< ", but zip needs vectors." << std::endl;
			continue;
		}

		std::vector<Symbol*>& curSym = ((SymbolArray*)pSym)->m_arr;
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
			((SymbolArray*)(*pVec)[i])->m_arr.push_back(curSym[i]->clone());
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

// split(string, delim)
static Symbol* fkt_split(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2 || vecSyms[0]->GetType()!=SYMBOL_STRING || vecSyms[1]->GetType()!=SYMBOL_STRING)
	{
		std::cerr << linenr("Error", info)
			<< "Split needs two string arguments." << std::endl;
		return 0;
	}

	const std::string& str = ((SymbolString*)vecSyms[0])->m_strVal;
	const std::string& strDelim = ((SymbolString*)vecSyms[1])->m_strVal;
	std::size_t iPos = str.find(strDelim);

	std::string str0 = str.substr(0, iPos);
	std::string str1;

	if(iPos != std::string::npos)
		str1 = str.substr(iPos + strDelim.length(), std::string::npos);

	SymbolArray* pSymRet = new SymbolArray();
	pSymRet->m_arr.push_back(new SymbolString(str0));
	pSymRet->m_arr.push_back(new SymbolString(str1));
	return pSymRet;
}

// tokens(string, delim)
static Symbol* fkt_tokens(const std::vector<Symbol*>& vecSyms,
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
				<< "Called tokens with invalid arguments." << std::endl;
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
// basic external functions
static t_mapFkts g_mapFkts =
{
	// basic stuff
	t_mapFkts::value_type("ver", fkt_version),

	// input/output
	t_mapFkts::value_type("output", fkt_output),
	t_mapFkts::value_type("input", fkt_input),
	t_mapFkts::value_type("print", fkt_print),	// output with "\n" at the end

	// modules
	t_mapFkts::value_type("import", fkt_import),

	// symbols & casts
	t_mapFkts::value_type("int", fkt_int),
	t_mapFkts::value_type("real", fkt_double),
	t_mapFkts::value_type("str", fkt_str),
	t_mapFkts::value_type("map", fkt_map),
	t_mapFkts::value_type("vec", fkt_array),
	t_mapFkts::value_type("has_var", fkt_has_var),
	t_mapFkts::value_type("typeof", fkt_typeof),

	// string operations
	t_mapFkts::value_type("trim", fkt_trim),
	t_mapFkts::value_type("split", fkt_split),
	t_mapFkts::value_type("tokens", fkt_tokens),
	t_mapFkts::value_type("length", fkt_array_size),

	// array operations
	t_mapFkts::value_type("vec_size", fkt_array_size),	// deprecated, use "length" instead
	t_mapFkts::value_type("cur_iter", fkt_cur_iter),
	t_mapFkts::value_type("zip", fkt_zip),
};

// --------------------------------------------------------------------------------

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

// --------------------------------------------------------------------------------

extern bool add_ext_call(const std::string& strFkt, t_extcall pExtCall)
{
	bool bInserted = g_mapFkts.insert(t_mapFkts::value_type(strFkt, pExtCall)).second;
	return bInserted;
}

extern void add_ext_calls(t_mapFkts& mapFkt)
{
	g_mapFkts.insert(mapFkt.begin(), mapFkt.end());
}
