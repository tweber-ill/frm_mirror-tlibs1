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
		if(pSym) delete pSym;
	
	m_arr.clear();

	//std::cout << "symarr -> del" << std::endl;
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
	
	pSym->m_arr.reserve(m_arr.size());

	for(Symbol *pArrSym : m_arr)
		pSym->m_arr.push_back(pArrSym->clone());

	pSym->UpdateIndices();
	return pSym;
}

void SymbolArray::assign(Symbol *pSym)
{
	SymbolArray *pOther = (SymbolArray*)pSym->ToType(GetType());
	this->m_arr = pOther->m_arr;
}

void SymbolArray::UpdateIndex(unsigned int iIdx)
{
	m_arr[iIdx]->m_pArr = this;
	m_arr[iIdx]->m_iArrIdx = iIdx;
}

void SymbolArray::UpdateIndices()
{
	for(unsigned int iIdx=0; iIdx<m_arr.size(); ++iIdx)
		UpdateIndex(iIdx);
}

std::vector<double> SymbolArray::ToDoubleArray() const
{
	std::vector<double> vec;
	vec.reserve(m_arr.size());

	for(const Symbol *pSym : m_arr)
	{
		double dVal = ((SymbolDouble*)pSym->ToType(SYMBOL_DOUBLE))->m_dVal;
		vec.push_back(dVal);
	}

	return vec;
}

void SymbolArray::FromDoubleArray(const std::vector<double>& vec)
{
	m_arr.reserve(m_arr.size() + vec.size());

	for(double d : vec)
	{
		SymbolDouble *pSym = new SymbolDouble;
		pSym->m_dVal = d;
		m_arr.push_back(pSym);
	}
}



SymbolMap::~SymbolMap()
{
	for(t_map::value_type& val : m_map)
		if(val.second) delete val.second;

	m_map.clear();
}

Symbol* SymbolMap::ToType(SymbolType stype) const
{
	Symbol *pNewSym = 0;

	if(stype == SYMBOL_MAP)
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
		std::cerr << "Error: Cannot convert map to other type."
					<< std::endl;

	return pNewSym;
}

std::string SymbolMap::print() const
{
	std::ostringstream ostr;

	ostr << "[";
	unsigned int iIter = 0;
	for(const t_map::value_type& val : m_map)
	{
		const Symbol* pSym = val.second;
		ostr << val.first << " : " << (val.second ? val.second->print() : "");

		if(iIter < m_map.size()-1)
			ostr << ", ";

		++iIter;
	}
	ostr << "]";

	return ostr.str();
}

Symbol* SymbolMap::clone() const
{
	//std::cout << "SymbolMap::clone" << std::endl;
	SymbolMap *pSym = new SymbolMap;

	for(const t_map::value_type& val : m_map)
		pSym->m_map.insert(t_map::value_type(val.first, val.second->clone()));

	pSym->UpdateIndices();
	return pSym;
}

void SymbolMap::assign(Symbol *pSym)
{
	//std::cout << "SymbolMap::assign" << std::endl;

	SymbolMap *pOther = (SymbolMap*)pSym->ToType(GetType());
	this->m_map = pOther->m_map;
}

void SymbolMap::UpdateIndex(const t_map::key_type& strKey)
{
	t_map::iterator iter = m_map.find(strKey);
	if(iter != m_map.end() && iter->second)
	{
		//std::cout << "updating index for " << strKey << std::endl;

		iter->second->m_pMap = this;
		iter->second->m_strMapKey = iter->first;
	}
}

void SymbolMap::UpdateIndices()
{
	for(const t_map::value_type& val : m_map)
	{
		//std::cout << "updating index for " << val.first << std::endl;

		val.second->m_pMap = this;
		val.second->m_strMapKey = val.first;
	}
}

std::string SymbolMap::GetStringVal(const std::string& strKey, bool *pbHasVal) const
{
	if(pbHasVal) *pbHasVal = 0;

	t_map::const_iterator iter = m_map.find(strKey);
	if(iter == m_map.end())
		return "";

	if(pbHasVal) *pbHasVal = 1;

	Symbol *pSym = iter->second;
	if(!pSym)
		return "";

	if(pbHasVal) *pbHasVal = 1;
	return pSym->print();
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
	RemoveSymbol(strKey);
	m_syms[strKey] = pSym;
}

void SymbolTable::RemoveSymbol(const std::string& strKey)
{
	t_syms::iterator iter = m_syms.find(strKey);
	if(iter != m_syms.end())
	{
		if(iter->second) delete iter->second;
		m_syms.erase(iter);
	}
}

void SymbolTable::RemoveSymbolNoDelete(const std::string& strKey)
{
	m_syms.erase(strKey);
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
