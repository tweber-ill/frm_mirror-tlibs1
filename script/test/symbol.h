/*
 * Simple Script
 * @author tweber
 */

#ifndef __MIEZE_SYM__
#define __MIEZE_SYM__

#include <map>
#include <string>
#include <iostream>
#include <sstream>

enum SymbolType
{
	SYMBOL_DOUBLE
};

struct Symbol
{
	std::string strName;

	virtual SymbolType GetType() const = 0;
	virtual Symbol* ToType(SymbolType stype) const = 0;
	virtual std::string print() const = 0;
	virtual Symbol* clone() const = 0;
};

struct SymbolDouble : public Symbol
{
	double dVal;

	virtual SymbolType GetType() const { return SYMBOL_DOUBLE; }
	virtual Symbol* ToType(SymbolType stype) const
	{
		Symbol *pNewSym = 0;

		if(stype == SYMBOL_DOUBLE)
		{
			pNewSym = new SymbolDouble();
			*pNewSym = *this;
		}

		return pNewSym;
	}

	virtual std::string print() const
	{
		std::ostringstream ostr;
		ostr << dVal;
		return ostr.str();
	}

	virtual Symbol* clone() const
	{
		SymbolDouble *pSym = new SymbolDouble;
		*pSym = *this;
		return pSym;
	}
};

static std::ostream& operator<<(std::ostream& ostr, const Symbol& sym)
{
	std::string str = sym.print();
	ostr << str;
	return ostr;
}


class SymbolTable
{
protected:
	typedef std::map<std::string, Symbol*> t_syms;
	t_syms m_syms;

public:
	SymbolTable();
	virtual ~SymbolTable();
	
	void print() const
	{
		for(auto val : m_syms)
			std::cout << val.first << " = " << (*val.second) << std::endl;
	}
	
	Symbol* GetSymbol(const std::string& strKey)
	{
		auto sym = m_syms.find(strKey);
		if(sym == m_syms.end())
			return 0;
		return sym->second;
	}
	
	void InsertSymbol(const std::string& strKey, Symbol *pSym)
	{
		m_syms[strKey] = pSym;
	}
	
	bool IsPtrInMap(const Symbol* pSym) const
	{
		for(auto entry : m_syms)
		{
			const Symbol *pSymInMap = entry.second;
			if(pSymInMap == pSym)
				return 1;
		}
		
		return 0;
	}
};

#endif
