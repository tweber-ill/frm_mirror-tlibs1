/*
 * Simple Script
 * @author tweber
 */

#include "symbol.h"

SymbolTable::SymbolTable()
{
}

SymbolTable::~SymbolTable()
{
}


Symbol* SymbolDouble::ToType(SymbolType stype) const
{
	Symbol *pNewSym = 0;

	if(stype == SYMBOL_DOUBLE)
	{
		pNewSym = new SymbolDouble();
		*pNewSym = *this;
	}

	return pNewSym;
}

std::string SymbolDouble::print() const
{
	std::ostringstream ostr;
	ostr << dVal;
	return ostr.str();
}

Symbol* SymbolDouble::clone() const
{
	SymbolDouble *pSym = new SymbolDouble;
	*pSym = *this;
	return pSym;
}





void SymbolTable::print() const
{
	for(const auto& val : m_syms)
	{
		Symbol *pSym = val.second;
		std::string strSym;
		if(pSym)
			strSym = pSym->print();
		else
			strSym = "<null>";

		std::cout << val.first << " = " << strSym << std::endl;
	}
}

Symbol* SymbolTable::GetSymbol(const std::string& strKey)
{
	auto sym = m_syms.find(strKey);
	if(sym == m_syms.end())
		return 0;
	return sym->second;
}

void SymbolTable::InsertSymbol(const std::string& strKey, Symbol *pSym)
{
	m_syms[strKey] = pSym;
}

bool SymbolTable::IsPtrInMap(const Symbol* pSym) const
{
	for(auto entry : m_syms)
	{
		const Symbol *pSymInMap = entry.second;
		if(pSymInMap == pSym)
			return 1;
	}

	return 0;
}
