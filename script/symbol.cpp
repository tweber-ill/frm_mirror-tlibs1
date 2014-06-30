/*
 * Symbol Table
 * @author tweber
 * @date 2013
 */

#include "symbol.h"
#include <set>
#include <limits>

const int SymbolDouble::m_defprec = std::numeric_limits<t_real>::digits10;
int SymbolDouble::m_prec = m_defprec;

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
		pNewSymI->SetName(this->GetName());
		pNewSymI->SetVal(t_int(this->GetVal()));

		pNewSym = pNewSymI;
	}
	else if(stype == SYMBOL_STRING)
	{
		SymbolString *pNewSymS = new SymbolString();
		pNewSymS->SetName(this->GetName());
		pNewSymS->SetVal(print());

		pNewSym = pNewSymS;
	}

	return pNewSym;
}

t_string SymbolDouble::print() const
{
	t_ostringstream ostr;
	ostr.precision(m_prec);

	ostr << GetVal();
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
	this->SetVal(pOther->GetVal());
}

bool SymbolDouble::IsLessThan(const Symbol& sym) const
{
	return GetValDouble() < sym.GetValDouble();
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
		pNewSymD->SetName(this->GetName());
		pNewSymD->SetVal(t_real(this->GetVal()));

		pNewSym = pNewSymD;
	}
	else if(stype == SYMBOL_STRING)
	{
		SymbolString *pNewSymS = new SymbolString();
		pNewSymS->SetName(this->GetName());
		pNewSymS->SetVal(print());

		pNewSym = pNewSymS;
	}

	return pNewSym;
}

t_string SymbolInt::print() const
{
	t_ostringstream ostr;
	ostr << GetVal();
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
	this->SetVal(pOther->GetVal());
}

bool SymbolInt::IsLessThan(const Symbol& sym) const
{
       	return GetValInt() < sym.GetValInt();
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
		t_istringstream istr(GetVal());

		SymbolInt *pNewSymI = new SymbolInt();
		pNewSymI->SetName(this->GetName());

		t_int iVal;
		istr >> iVal;
		pNewSymI->SetVal(iVal);

		pNewSym = pNewSymI;
	}
	else if(stype == SYMBOL_DOUBLE)
	{
		t_istringstream istr(GetVal());

		SymbolDouble *pNewSymD = new SymbolDouble();
		pNewSymD->SetName(this->GetName());
		t_real dVal;
		istr >> dVal;
		pNewSymD->SetVal(dVal);

		pNewSym = pNewSymD;
	}

	return pNewSym;
}

t_string SymbolString::print() const
{
	t_ostringstream ostr;
	ostr << GetVal();
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
	this->SetVal(pOther->GetVal());
}

bool SymbolString::IsLessThan(const Symbol& sym) const
{
       	return GetVal() < sym.print();
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
		pNewSymS->SetName(this->GetName());
		pNewSymS->SetVal(print());

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

	pSym->GetArr().reserve(m_arr.size());

	for(Symbol *pArrSym : m_arr)
		pSym->GetArr().push_back(pArrSym->clone());

	pSym->UpdateIndices();
	return pSym;
}

void SymbolArray::assign(Symbol *pSym)
{
	SymbolArray *pOther = (SymbolArray*)pSym->ToType(GetType());
	this->m_arr = pOther->GetArr();
}

void SymbolArray::UpdateIndex(unsigned int iIdx)
{
	m_arr[iIdx]->SetArrPtr(this);
	m_arr[iIdx]->SetArrIdx(iIdx);
}

void SymbolArray::UpdateIndices()
{
	for(unsigned int iIdx=0; iIdx<m_arr.size(); ++iIdx)
		UpdateIndex(iIdx);
}

std::vector<t_real> SymbolArray::ToDoubleArray() const
{
	std::vector<t_real> vec;
	vec.reserve(m_arr.size());

	for(const Symbol *pSym : m_arr)
	{
		t_real dVal = ((SymbolDouble*)pSym->ToType(SYMBOL_DOUBLE))->GetVal();
		vec.push_back(dVal);
	}

	return vec;
}

void SymbolArray::FromDoubleArray(const std::vector<t_real>& vec)
{
	m_arr.reserve(m_arr.size() + vec.size());

	for(t_real d : vec)
	{
		SymbolDouble *pSym = new SymbolDouble;
		pSym->SetVal(d);
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
		pNewSymS->SetName(this->GetName());
		pNewSymS->SetVal(print());

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

		iter->second->SetMapPtr(this);
		iter->second->SetMapKey(iter->first);
	}
}

void SymbolMap::UpdateIndices()
{
	for(const t_map::value_type& val : m_map)
	{
		//G_COUT << "updating index for " << val.first << std::endl;

		val.second->SetMapPtr(this);
		val.second->SetMapKey(val.first);
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

t_int SymbolMap::GetIntVal(const t_string& strKey, bool *pbHasVal) const
{
	if(pbHasVal) *pbHasVal = 0;

	t_map::const_iterator iter = m_map.find(strKey);
	if(iter == m_map.end())
		return 0;

	if(pbHasVal) *pbHasVal = 1;

	const Symbol *pSym = iter->second;
	if(!pSym)
		return 0;

	if(pbHasVal) *pbHasVal = 1;
		return pSym->GetValInt();
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
			t_string strType = pSym->GetTypeName();
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
	for(const Symbol* pSymInArr : pSymArr->GetArr())
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
	if(piNumRows) *piNumRows = pSymArr->GetArr().size();

	unsigned int iVecSize = 0;
	bool bHasSize = 0;
	for(const Symbol* pSymInArr : pSymArr->GetArr())
	{
		if(!is_vec(pSymInArr))
			return false;

		unsigned int iSize = ((SymbolArray*)pSymInArr)->GetArr().size();
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
