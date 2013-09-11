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
	virtual Symbol* ToType(SymbolType stype) const;
	virtual std::string print() const;
	virtual Symbol* clone() const;
};

/*
static std::ostream& operator<<(std::ostream& ostr, const Symbol& sym)
{
	std::string str = sym.print();
	ostr << str;
	return ostr;
}
*/

class SymbolTable
{
protected:
	typedef std::map<std::string, Symbol*> t_syms;
	t_syms m_syms;

public:
	SymbolTable();
	virtual ~SymbolTable();
	
	void print() const;
	
	Symbol* GetSymbol(const std::string& strKey);
	void InsertSymbol(const std::string& strKey, Symbol *pSym);
	bool IsPtrInMap(const Symbol* pSym) const;
};

#endif
