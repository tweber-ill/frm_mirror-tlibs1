/*
 * Symbol Table
 * @author tweber
 * @date 2013
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

t_string SymbolDouble::print() const
{
	t_ostringstream ostr;
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

t_string SymbolInt::print() const
{
	t_ostringstream ostr;
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
		t_istringstream istr(m_strVal);

		SymbolInt *pNewSymI = new SymbolInt();
		pNewSymI->m_strName = this->m_strName;
		istr >> pNewSymI->m_iVal;

		pNewSym = pNewSymI;
	}
	else if(stype == SYMBOL_DOUBLE)
	{
		t_istringstream istr(m_strVal);

		SymbolDouble *pNewSymD = new SymbolDouble();
		pNewSymD->m_strName = this->m_strName;
		istr >> pNewSymD->m_dVal;

		pNewSym = pNewSymD;
	}

	return pNewSym;
}

t_string SymbolString::print() const
{
	t_ostringstream ostr;
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
	if(!m_bDontDel)
		for(Symbol *pSym : m_arr)
			if(pSym) delete pSym;

	m_arr.clear();
}

Symbol* SymbolArray::ToType(SymbolType stype) const
{
	Symbol *pNewSym = 0;

	if(stype == SYMBOL_ARRAY)
	{
		pNewSym = this->clone();
	}
	else if(stype == SYMBOL_STRING)
	{
		SymbolString *pNewSymS = new SymbolString();
		pNewSymS->m_strName = this->m_strName;
		pNewSymS->m_strVal = print();

		pNewSym = pNewSymS;
	}
	else
		G_CERR << "Error: Cannot convert array to type "
			<< stype << "." << std::endl;

	return pNewSym;
}

t_string SymbolArray::print() const
{
	t_ostringstream ostr;

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
		G_CERR << "Error: Cannot convert map to other type."
					<< std::endl;

	return pNewSym;
}

t_string SymbolMap::print() const
{
	t_ostringstream ostr;

	ostr << "[";
	unsigned int iIter = 0;
	for(const t_map::value_type& val : m_map)
	{
		const Symbol* pSym = val.second;
		ostr << val.first << " : " << (val.second ? val.second->print() : T_STR"");

		if(iIter < m_map.size()-1)
			ostr << ", ";

		++iIter;
	}
	ostr << "]";

	return ostr.str();
}

Symbol* SymbolMap::clone() const
{
	//G_COUT << "SymbolMap::clone" << std::endl;
	SymbolMap *pSym = new SymbolMap;

	for(const t_map::value_type& val : m_map)
		pSym->m_map.insert(t_map::value_type(val.first, val.second->clone()));

	pSym->UpdateIndices();
	return pSym;
}

void SymbolMap::assign(Symbol *pSym)
{
	//G_COUT << "SymbolMap::assign" << std::endl;

	SymbolMap *pOther = (SymbolMap*)pSym->ToType(GetType());
	this->m_map = pOther->m_map;
}

void SymbolMap::UpdateIndex(const t_map::key_type& strKey)
{
	t_map::iterator iter = m_map.find(strKey);
	if(iter != m_map.end() && iter->second)
	{
		//G_COUT << "updating index for " << strKey << std::endl;

		iter->second->m_pMap = this;
		iter->second->m_strMapKey = iter->first;
	}
}

void SymbolMap::UpdateIndices()
{
	for(const t_map::value_type& val : m_map)
	{
		//G_COUT << "updating index for " << val.first << std::endl;

		val.second->m_pMap = this;
		val.second->m_strMapKey = val.first;
	}
}

t_string SymbolMap::GetStringVal(const t_string& strKey, bool *pbHasVal) const
{
	if(pbHasVal) *pbHasVal = 0;

	t_map::const_iterator iter = m_map.find(strKey);
	if(iter == m_map.end())
		return T_STR"";

	if(pbHasVal) *pbHasVal = 1;

	Symbol *pSym = iter->second;
	if(!pSym)
		return T_STR"";

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
		t_string strSym;
		if(pSym)
		{
			t_string strType;
			if(pSym->GetType() == SYMBOL_DOUBLE)
				strType = T_STR"double";
			else if(pSym->GetType() == SYMBOL_INT)
				strType = T_STR"int";
			else if(pSym->GetType() == SYMBOL_STRING)
				strType = T_STR"string";
			else
				strType = T_STR"<unknown>";

			strSym = pSym->print() + T_STR" (" + strType + T_STR")";
		}
		else
			strSym = T_STR"<null>";

		G_COUT << val.first << " = " << strSym << std::endl;
	}
}

Symbol* SymbolTable::GetSymbol(const t_string& strKey)
{
	auto sym = m_syms.find(strKey);
	if(sym == m_syms.end())
		return 0;
	return sym->second;
}

void SymbolTable::InsertSymbol(const t_string& strKey, Symbol *pSym)
{
	RemoveSymbol(strKey);
	m_syms[strKey] = pSym;
}

void SymbolTable::RemoveSymbol(const t_string& strKey)
{
	t_syms::iterator iter = m_syms.find(strKey);
	if(iter != m_syms.end())
	{
		if(iter->second) delete iter->second;
		m_syms.erase(iter);
	}
}

void SymbolTable::RemoveSymbolNoDelete(const t_string& strKey)
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


// --------------------------------------------------------------------------------


bool is_vec(const Symbol* pSym)
{
	if(!pSym)
		return false;
	if(pSym->GetType() != SYMBOL_ARRAY)
		return false;

	const SymbolArray* pSymArr = (SymbolArray*)pSym;
	for(const Symbol* pSymInArr : pSymArr->m_arr)
	{
		if(!pSymInArr->IsScalar())
			return false;
	}
	return true;
}

bool is_mat(const Symbol* pSym, unsigned int *piNumCols, unsigned int *piNumRows)
{
	if(!pSym)
		return false;
	if(pSym->GetType() != SYMBOL_ARRAY)
		return false;

	const SymbolArray* pSymArr = (SymbolArray*)pSym;
	if(piNumRows) *piNumRows = pSymArr->m_arr.size();

	unsigned int iVecSize = 0;
	bool bHasSize = 0;
	for(const Symbol* pSymInArr : pSymArr->m_arr)
	{
		if(!is_vec(pSymInArr))
			return false;

		unsigned int iSize = ((SymbolArray*)pSymInArr)->m_arr.size();
		if(!bHasSize)
		{
			iVecSize = iSize;
			bHasSize = 1;

			if(piNumCols) *piNumCols = iVecSize;
		}

		// element vectors have to be of the same size
		if(iSize != iVecSize)
			return false;
	}

	return true;
}
