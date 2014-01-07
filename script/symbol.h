/*
 * Symbol Table
 * @author tweber
 * @date 2013
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
	
	SYMBOL_ARRAY,
	SYMBOL_MAP
};

struct SymbolArray;
struct SymbolMap;
struct Symbol
{
	std::string m_strName;
	std::string m_strIdent;			// last seen identifier

	unsigned int m_iArrIdx;			// if symbol is contained in an array
	SymbolArray *m_pArr;

	std::string m_strMapKey;		// if symbol is contained in a map
	SymbolMap *m_pMap;

	Symbol() : m_iArrIdx(0), m_pArr(0), m_pMap(0) {}
	virtual ~Symbol() {}
	
	virtual SymbolType GetType() const = 0;
	virtual std::string GetTypeName() const = 0;

	// cast and clone symbol
	virtual Symbol* ToType(SymbolType stype) const = 0;

	virtual std::string print() const = 0;
	virtual Symbol* clone() const = 0;
	virtual void assign(Symbol *pSym) = 0;

	virtual bool IsNotZero() const = 0;

	// cast value without converting or cloning symbol
	virtual int GetValInt() const { return 0; }
	virtual double GetValDouble() const { return 0.; }
};

struct SymbolDouble : public Symbol
{
	double m_dVal;
	
	SymbolDouble() : Symbol() {}
	SymbolDouble(double dVal) : m_dVal(dVal) {}

	virtual SymbolType GetType() const { return SYMBOL_DOUBLE; }
	virtual std::string GetTypeName() const { return "real"; }
	virtual Symbol* ToType(SymbolType stype) const;

	virtual std::string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);

	virtual bool IsNotZero() const { return m_dVal != 0.; }

	virtual int GetValInt() const { return int(m_dVal); }
	virtual double GetValDouble() const { return m_dVal; }
};

struct SymbolInt : public Symbol
{
	int m_iVal;
	
	SymbolInt() : Symbol() {}
	SymbolInt(int iVal) : m_iVal(iVal) {}

	virtual SymbolType GetType() const { return SYMBOL_INT; }
	virtual std::string GetTypeName() const { return "int"; }
	virtual Symbol* ToType(SymbolType stype) const;

	virtual std::string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);

	virtual bool IsNotZero() const { return m_iVal != 0; }

	virtual int GetValInt() const { return m_iVal; }
	virtual double GetValDouble() const { return double(m_iVal); }
};

struct SymbolString : public Symbol
{
	std::string m_strVal;
	
	SymbolString() : Symbol() {}
	SymbolString(const char* pcStr) : m_strVal(pcStr) {}
	SymbolString(const std::string& str) : m_strVal(str) {}

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
	
	SymbolArray() : Symbol() { /*std::cout << "symarr -> new" << std::endl;*/ }
	virtual ~SymbolArray();
	
	virtual SymbolType GetType() const { return SYMBOL_ARRAY; }
	virtual std::string GetTypeName() const { return "vector"; }
	virtual Symbol* ToType(SymbolType stype) const;

	std::vector<double> ToDoubleArray() const;
	void FromDoubleArray(const std::vector<double>& vec);

	virtual std::string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);
	
	virtual bool IsNotZero() const { return 0; }

	void UpdateIndex(unsigned int);
	void UpdateIndices();
};

struct SymbolMap : public Symbol
{
	typedef std::map<std::string, Symbol*> t_map;
	t_map m_map;

	SymbolMap() : Symbol() {}
	virtual ~SymbolMap();

	virtual SymbolType GetType() const { return SYMBOL_MAP; }	
	virtual std::string GetTypeName() const { return "map"; }
	virtual Symbol* ToType(SymbolType stype) const;

	virtual std::string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);

	virtual bool IsNotZero() const { return 0; }

	void UpdateIndex(const t_map::key_type& strKey);
	void UpdateIndices();

	std::string GetStringVal(const std::string& strKey, bool *pbHasVal=0) const;
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
	void RemoveSymbolNoDelete(const std::string& strKey);
	bool IsPtrInMap(const Symbol* pSym) const;
};

#endif
