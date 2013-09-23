/*
 * External Functions
 * @author tweber
 */

#include "calls.h"
#include <sstream>
#include <iostream>
#include <map>


static Symbol* fkt_version(const std::vector<Symbol*>& vecSyms)
{
	return new SymbolString("Hermelin Interpreter, Version 0.1");
}


static Symbol* fkt_print(const std::vector<Symbol*>& vecSyms)
{
	std::ostream& ostr = std::cout;
	
	for(Symbol *pSym : vecSyms)
		if(pSym)
			ostr << pSym->print();

	ostr << std::endl;
	return 0;
}

static Symbol* fkt_array(const std::vector<Symbol*>& vecSyms)
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

static Symbol* fkt_array_size(const std::vector<Symbol*>& vecSyms)
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


typedef std::map<std::string, Symbol*(*)(const std::vector<Symbol*>&)> t_mapFkts;
static t_mapFkts g_mapFkts =
{
	t_mapFkts::value_type("ver", fkt_version),
	t_mapFkts::value_type("print", fkt_print),
	
	t_mapFkts::value_type("vec", fkt_array),
	t_mapFkts::value_type("vec_size", fkt_array_size),
};

extern Symbol* ext_call(const std::string& strFkt, const std::vector<Symbol*>& vecSyms)
{
	t_mapFkts::iterator iter = g_mapFkts.find(strFkt);
	if(iter == g_mapFkts.end())
	{
		std::cerr << "Error: Tried to call unknown function \"" 
					<< strFkt << "\"."
					<< std::endl;
		return 0;
	}
	
	return (*iter).second(vecSyms);
}
