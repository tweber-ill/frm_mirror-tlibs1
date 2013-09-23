/*
 * Symbol Table
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
		pNewSymI->m_strName = this->m_strName;
		pNewSymI->m_iVal = int(this->m_dVal);

		pNewSym = pNewSymI;
	}
	else if(stype == SYMBOL_STRING)
	{
		SymbolString *pNewSymS = new SymbolString();
		pNewSymS->m_strName = this->m_strName;
		pNewSymS->m_strVal = print();

		pNewSym = pNewSymS;
	}

	return pNewSym;
}

std::string SymbolDouble::print() const
{
	std::ostringstream ostr;
	ostr << m_dVal;
	return ostr.str();
}

Symbol* SymbolDouble::clone() const
{
	SymbolDouble *pSym = new SymbolDouble;
	*pSym = *this;
	return pSym;
}

void SymbolDouble::assign(Symbol *pSym)
{
	SymbolDouble *pOther = (SymbolDouble*)pSym->ToType(GetType());
	this->m_dVal = pOther->m_dVal;
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
		pNewSymD->m_strName = this->m_strName;
		pNewSymD->m_dVal = double(this->m_iVal);

		pNewSym = pNewSymD;
	}
	else if(stype == SYMBOL_STRING)
	{
		SymbolString *pNewSymS = new SymbolString();
		pNewSymS->m_strName = this->m_strName;
		pNewSymS->m_strVal = print();

		pNewSym = pNewSymS;
	}

	return pNewSym;
}

std::string SymbolInt::print() const
{
	std::ostringstream ostr;
	ostr << m_iVal;
	return ostr.str();
}

Symbol* SymbolInt::clone() const
{
	SymbolInt *pSym = new SymbolInt;
	*pSym = *this;
	return pSym;
}

void SymbolInt::assign(Symbol *pSym)
{
	SymbolInt *pOther = (SymbolInt*)pSym->ToType(GetType());
	this->m_iVal = pOther->m_iVal;
}





Symbol* SymbolString::ToType(SymbolType stype) const
{
	Symbol *pNewSym = 0;

	if(stype == SYMBOL_STRING)
	{
		pNewSym = this->clone();
	}
	else if(stype == SYMBOL_INT)
	{
		std::istringstream istr(m_strVal);

		SymbolInt *pNewSymI = new SymbolInt();
		pNewSymI->m_strName = this->m_strName;
		istr >> pNewSymI->m_iVal;

		pNewSym = pNewSymI;
	}
	else if(stype == SYMBOL_DOUBLE)
	{
		std::istringstream istr(m_strVal);

		SymbolDouble *pNewSymD = new SymbolDouble();
		pNewSymD->m_strName = this->m_strName;
		istr >> pNewSymD->m_dVal;

		pNewSym = pNewSymD;
	}

	return pNewSym;
}

std::string SymbolString::print() const
{
	std::ostringstream ostr;
	ostr << m_strVal;
	return ostr.str();
}

Symbol* SymbolString::clone() const
{
	SymbolString *pSym = new SymbolString;
	*pSym = *this;
	return pSym;
}

void SymbolString::assign(Symbol *pSym)
{
	SymbolString *pOther = (SymbolString*)pSym->ToType(GetType());
	this->m_strVal = pOther->m_strVal;
}




SymbolArray::~SymbolArray()
{
	for(Symbol *pSym : m_arr)
		delete pSym;
	
	m_arr.clear();
}

Symbol* SymbolArray::ToType(SymbolType stype) const
{
	Symbol *pNewSym = 0;
	
	if(stype == SYMBOL_ARRAY)
	{
		pNewSym = this->clone();
	}
	if(stype == SYMBOL_STRING)
	{
		SymbolString *pNewSymS = new SymbolString();
		pNewSymS->m_strName = this->m_strName;
		pNewSymS->m_strVal = print();

		pNewSym = pNewSymS;
	}
	else
		std::cerr << "Error: Cannot convert array to other type."
					<< std::endl;

	return pNewSym;
}

std::string SymbolArray::print() const
{
	std::ostringstream ostr;

	ostr << "[";
	for(unsigned int i=0; i<m_arr.size(); ++i)
	{
		const Symbol* pSym = m_arr[i];
		ostr << pSym->print();
		
		if(i<m_arr.size()-1)
			ostr << ", ";
	}
	ostr << "]";
	
	return ostr.str();
}

Symbol* SymbolArray::clone() const
{
	SymbolArray *pSym = new SymbolArray;
	//*pSym = *this;
	
	pSym->m_arr.reserve(m_arr.size());
	for(Symbol *pArrSym : m_arr)
		pSym->m_arr.push_back(pArrSym->clone());
	
	return pSym;
}

void SymbolArray::assign(Symbol *pSym)
{
	SymbolArray *pOther = (SymbolArray*)pSym->ToType(GetType());
	this->m_arr = pOther->m_arr;
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
			else if(pSym->GetType() == SYMBOL_STRING)
				strType = "string";
			else
				strType = "<unknown>";

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
