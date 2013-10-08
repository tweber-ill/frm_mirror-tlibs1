/*
 * Symbol Table
 * @author tweber
 */

#ifndef __MIEZE_SYM__
#define __MIEZE_SYM__

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

enum SymbolType
{
	SYMBOL_DOUBLE,
	SYMBOL_INT,
	SYMBOL_STRING,
	
	SYMBOL_ARRAY
};

struct Symbol
{
	std::string m_strName;
	std::string m_strIdent;		// last seen identifier
	
	virtual SymbolType GetType() const = 0;
	virtual std::string GetTypeName() const = 0;
	virtual Symbol* ToType(SymbolType stype) const = 0;

	virtual std::string print() const = 0;
	virtual Symbol* clone() const = 0;
	virtual void assign(Symbol *pSym) = 0;

	virtual bool IsNotZero() const = 0;
};

struct SymbolDouble : public Symbol
{
	double m_dVal;
	
	SymbolDouble() {}
	SymbolDouble(double dVal) : m_dVal(dVal) {}

	virtual SymbolType GetType() const { return SYMBOL_DOUBLE; }
	virtual std::string GetTypeName() const { return "double"; }
	virtual Symbol* ToType(SymbolType stype) const;

	virtual std::string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);

	virtual bool IsNotZero() const { return m_dVal != 0.; }
};

struct SymbolInt : public Symbol
{
	int m_iVal;
	
	SymbolInt() {}
	SymbolInt(int iVal) : m_iVal(iVal) {}

	virtual SymbolType GetType() const { return SYMBOL_INT; }
	virtual std::string GetTypeName() const { return "int"; }
	virtual Symbol* ToType(SymbolType stype) const;

	virtual std::string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);

	virtual bool IsNotZero() const { return m_iVal != 0.; }
};

struct SymbolString : public Symbol
{
	std::string m_strVal;
	
	SymbolString() {}
	SymbolString(const char* pcStr) : m_strVal(pcStr) {}

	virtual SymbolType GetType() const { return SYMBOL_STRING; }
	virtual std::string GetTypeName() const { return "string"; }
	virtual Symbol* ToType(SymbolType stype) const;

	virtual std::string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);

	virtual bool IsNotZero() const { return 0; }
};




struct SymbolArray : public Symbol
{
	std::vector<Symbol*> m_arr;
	
	virtual ~SymbolArray();
	
	virtual SymbolType GetType() const { return SYMBOL_ARRAY; }
	virtual std::string GetTypeName() const { return "vector"; }
	virtual Symbol* ToType(SymbolType stype) const;

	virtual std::string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);
	
	virtual bool IsNotZero() const { return 0; }
};




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
	void RemoveSymbol(const std::string& strKey);
	bool IsPtrInMap(const Symbol* pSym) const;
};

#endif
