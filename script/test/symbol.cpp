/*
 * Simple Script
 * @author tweber
 */

#include "symbol.h"
#include <set>


Symbol* SymbolDouble::ToType(SymbolType stype) const
{
	Symbol *pNewSym = 0;

	if(stype == SYMBOL_DOUBLE)
	{
		pNewSym = this->clone();
	}
	else if(stype == SYMBOL_INT)
	{
		SymbolInt *pNewSymI = new SymbolInt();
		pNewSymI->strName = this->strName;
		pNewSymI->iVal = int(this->dVal);

		pNewSym = pNewSymI;
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



Symbol* SymbolInt::ToType(SymbolType stype) const
{
	Symbol *pNewSym = 0;

	if(stype == SYMBOL_INT)
	{
		pNewSym = this->clone();
	}
	else if(stype == SYMBOL_DOUBLE)
	{
		SymbolDouble *pNewSymD = new SymbolDouble();
		pNewSymD->strName = this->strName;
		pNewSymD->dVal = double(this->iVal);

		pNewSym = pNewSymD;
	}

	return pNewSym;
}

std::string SymbolInt::print() const
{
	std::ostringstream ostr;
	ostr << iVal;
	return ostr.str();
}

Symbol* SymbolInt::clone() const
{
	SymbolInt *pSym = new SymbolInt;
	*pSym = *this;
	return pSym;
}


//--------------------------------------------------------------------------------

SymbolTable::SymbolTable()
{
}

SymbolTable::~SymbolTable()
{
	std::set<Symbol*> m_setDeleted;

	for(t_syms::iterator iter=m_syms.begin(); iter!=m_syms.end(); ++iter)
	{
		Symbol *pSym = iter->second;
		if(pSym && m_setDeleted.find(pSym)!=m_setDeleted.end())
		{
			m_setDeleted.insert(pSym);

			delete pSym;
			iter->second = 0;
		}
	}
}


void SymbolTable::print() const
{
	for(const auto& val : m_syms)
	{
		Symbol *pSym = val.second;
		std::string strSym;
		if(pSym)
		{
			std::string strType;
			if(pSym->GetType() == SYMBOL_DOUBLE)
				strType = "double";
			else if(pSym->GetType() == SYMBOL_INT)
				strType = "int";

			strSym = pSym->print() + " (" + strType + ")";
		}
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
