/*
 * Symbol Table
 * @author tweber
 * @date 2013
 */

#ifndef __HERMELIN_SYM__
#define __HERMELIN_SYM__

#include "types.h"
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "helper/exception.h"

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
protected:
	t_string m_strName;
	t_string m_strIdent;			// last seen identifier

	unsigned int m_iArrIdx;			// if symbol is contained in an array
	SymbolArray *m_pArr;

	t_string m_strMapKey;			// if symbol is contained in a map
	SymbolMap *m_pMap;

public:
	Symbol() : m_iArrIdx(0), m_pArr(0), m_pMap(0) {}
	virtual ~Symbol() {}

	virtual SymbolType GetType() const = 0;
	virtual t_string GetTypeName() const = 0;

	virtual const t_string& GetName() const { return m_strName; }
	virtual const t_string& GetIdent() const { return m_strIdent; }
	virtual void SetName(const t_string& strName) { m_strName = strName; }
	virtual void SetIdent(const t_string& strIdent) { m_strIdent = strIdent; }

	// cast and clone symbol
	virtual Symbol* ToType(SymbolType stype) const = 0;

	virtual t_string print() const = 0;
	virtual Symbol* clone() const = 0;
	virtual void assign(Symbol *pSym) = 0;
	//virtual bool equals(Symbol *pSym) const = 0;
	virtual bool IsNotZero() const = 0;
	virtual bool IsLessThan(const Symbol&) const { return 0; }

	// cast value without converting or cloning symbol
	virtual t_int GetValInt() const { return 0; }
	virtual t_real GetValDouble() const { return 0.; }

	virtual bool IsScalar() const { return 0; }

public:
	SymbolMap* GetMapPtr() { return m_pMap; }
	void SetMapPtr(SymbolMap* pMap) { m_pMap = pMap; }

	const t_string& GetMapKey() const { return m_strMapKey; }
	void SetMapKey(const t_string& str) { m_strMapKey = str; }


	SymbolArray* GetArrPtr() const { return m_pArr; }
	void SetArrPtr(SymbolArray* pArr) { m_pArr = pArr; }

	unsigned int GetArrIdx() const { return m_iArrIdx; }
	void SetArrIdx(unsigned int iIdx) { m_iArrIdx = iIdx; }
};

struct SymbolDouble : public Symbol
{
protected:
	t_real m_dVal;
	static const int m_defprec;
	static int m_prec;

public:
	SymbolDouble() : Symbol(), m_dVal(0.) {}
	SymbolDouble(t_real dVal) : m_dVal(dVal) {}
	SymbolDouble(const t_string&) { throw Err("Invalid SymbolDouble constructor."); }

	virtual SymbolType GetType() const { return SYMBOL_DOUBLE; }
	virtual t_string GetTypeName() const { return T_STR"real"; }
	virtual Symbol* ToType(SymbolType stype) const;

	virtual t_string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);
	//virtual bool equals(Symbol *pSym) const;
	virtual bool IsLessThan(const Symbol&) const;
	virtual bool IsNotZero() const { return m_dVal != 0.; }

	virtual t_int GetValInt() const { return t_int(m_dVal); }
	virtual t_real GetValDouble() const { return m_dVal; }

	static const int GetDefPrec() { return m_defprec; }
	static const int GetPrec() { return m_prec; }
	static void SetPrec(int iPrec) { m_prec = iPrec; }

	void SetVal(double dVal) { m_dVal = dVal; }
	const t_real& GetVal() const { return m_dVal; }
	t_real& GetVal() { return m_dVal; }

	virtual bool IsScalar() const { return 1; }
};

struct SymbolInt : public Symbol
{
protected:
	t_int m_iVal;

public:
	SymbolInt() : Symbol(), m_iVal(0) {}
	SymbolInt(t_int iVal) : m_iVal(iVal) {}

	virtual SymbolType GetType() const { return SYMBOL_INT; }
	virtual t_string GetTypeName() const { return T_STR"int"; }
	virtual Symbol* ToType(SymbolType stype) const;

	virtual t_string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);
	//virtual bool equals(Symbol *pSym) const;
	virtual bool IsLessThan(const Symbol&) const;
	virtual bool IsNotZero() const { return m_iVal != 0; }

	virtual t_int GetValInt() const { return m_iVal; }
	virtual t_real GetValDouble() const { return t_real(m_iVal); }

	void SetVal(int iVal) { m_iVal = iVal; }
	const t_int& GetVal() const { return m_iVal; }
	t_int& GetVal() { return m_iVal; }

	virtual bool IsScalar() const { return 1; }
};

struct SymbolString : public Symbol
{
protected:
	t_string m_strVal;

public:
	SymbolString() : Symbol() {}
	SymbolString(const t_char* pcStr) : m_strVal(pcStr) {}
	SymbolString(const t_string& str) : m_strVal(str) {}
	SymbolString(t_real dVal) { throw Err("Invalid SymbolString constructor."); }

	virtual SymbolType GetType() const { return SYMBOL_STRING; }
	virtual t_string GetTypeName() const { return T_STR"string"; }
	virtual Symbol* ToType(SymbolType stype) const;
	//virtual bool equals(Symbol *pSym) const;
	virtual bool IsLessThan(const Symbol&) const;

	virtual t_string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);

	void SetVal(const t_string& str) { m_strVal = str; }
	const t_string& GetVal() const { return m_strVal; }
	t_string& GetVal() { return m_strVal; }

	virtual bool IsNotZero() const { return 0; }
};


struct SymbolArray : public Symbol
{
protected:
	bool m_bDontDel;
	std::vector<Symbol*> m_arr;

public:
	SymbolArray() : Symbol(), m_bDontDel(0) { /*std::cout << "symarr -> new" << std::endl;*/ }
	virtual ~SymbolArray();

	virtual SymbolType GetType() const { return SYMBOL_ARRAY; }
	virtual t_string GetTypeName() const { return T_STR"vector"; }
	virtual Symbol* ToType(SymbolType stype) const;

	std::vector<t_real> ToDoubleArray() const;
	void FromDoubleArray(const std::vector<t_real>& vec);

	virtual t_string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);
	//virtual bool equals(Symbol *pSym) const;

	virtual bool IsNotZero() const { return 0; }

	void UpdateIndex(unsigned int);
	void UpdateIndices();

	const std::vector<Symbol*>& GetArr() const { return m_arr; }
	std::vector<Symbol*>& GetArr() { return m_arr; }

	virtual bool IsScalar() const { return 0; }


	bool GetDontDel() const { return m_bDontDel; }
	void SetDontDel(bool b) { m_bDontDel = b; }
};

struct SymbolMap : public Symbol
{
public:
	typedef std::map<t_string, Symbol*> t_map;
protected:
	t_map m_map;

public:
	SymbolMap() : Symbol() {}
	virtual ~SymbolMap();

	// no copy/move constructors -> use clone
	//SymbolMap(const SymbolMap& map);
	//SymbolMap(SymbolMap&& map);

	virtual SymbolType GetType() const { return SYMBOL_MAP; }
	virtual t_string GetTypeName() const { return T_STR"map"; }
	virtual Symbol* ToType(SymbolType stype) const;

	virtual t_string print() const;
	virtual Symbol* clone() const;
	virtual void assign(Symbol *pSym);
	//virtual bool equals(Symbol *pSym) const;

	virtual bool IsNotZero() const { return 0; }

	void UpdateIndex(const t_map::key_type& strKey);
	void UpdateIndices();

	t_string GetStringVal(const t_string& strKey, bool *pbHasVal=0) const;
	t_int GetIntVal(const t_string& strKey, bool *pbHasVal=0) const;

	virtual bool IsScalar() const { return 0; }

	const t_map& GetMap() const { return m_map; }
	t_map& GetMap() { return m_map; }
};


// --------------------------------------------------------------------------------


class SymbolTable
{
protected:
	typedef std::map<t_string, Symbol*> t_syms;
	t_syms m_syms;

public:
	SymbolTable();
	virtual ~SymbolTable();

	void print() const;

	Symbol* GetSymbol(const t_string& strKey);
	void InsertSymbol(const t_string& strKey, Symbol *pSym);
	void RemoveSymbol(const t_string& strKey);
	void RemoveSymbolNoDelete(const t_string& strKey);
	bool IsPtrInMap(const Symbol* pSym) const;
};



// --------------------------------------------------------------------------------
// conversions

#include <type_traits>

template<typename T> static T convert_symbol(const Symbol* pSym)
{ std::cerr << "Error: Invalid symbol conversion." << std::endl; return T(); }

template<> t_string convert_symbol<t_string>(const Symbol* pSym)
{ return pSym->print(); }
template<>  t_real convert_symbol<t_real>(const Symbol* pSym)
{ return pSym->GetValDouble(); }
template<> t_int convert_symbol<t_int>(const Symbol* pSym)
{ return pSym->GetValInt(); }


template<typename T> static Symbol* create_symbol(const T& t)
{ G_CERR << "Error: Invalid symbol creation." << std::endl; return 0; }

template<> Symbol* create_symbol<t_real>(const t_real& t)
{ return new SymbolDouble(t); }
template<> Symbol* create_symbol<t_int>(const t_int& t)
{ return new SymbolInt(t); }
template<> Symbol* create_symbol<t_string>(const t_string& t)
{ return new SymbolString(t); }


template<typename T1=t_string, typename T2=t_real>
static std::map<T1, T2> sym_to_map(const Symbol* pSym)
{
	if(pSym->GetType() != SYMBOL_MAP)
		return std::map<T1, T2>();

	SymbolMap* pSymMap = (SymbolMap*)pSym;
	std::map<T1, T2> _map;

	unsigned int iIdx = 0;
	for(const typename SymbolMap::t_map::value_type& pair : pSymMap->GetMap())
		_map[pair.first] = convert_symbol<T2>(pair.second);

	return _map;
}

extern bool is_vec(const Symbol* pSym);
extern bool is_mat(const Symbol* pSym, unsigned int *piNumCols=0, unsigned int *piNumRows=0);

template<template<class> class t_vec, typename T=t_real>
static t_vec<T> sym_to_vec(const Symbol* pSym)
{
	if(!pSym || pSym->GetType() != SYMBOL_ARRAY)
		return t_vec<T>();

	SymbolArray* pSymArr = (SymbolArray*)pSym;
	t_vec<T> vec(pSymArr->GetArr().size());

	unsigned int iIdx = 0;
	for(const Symbol* pSymInArr : pSymArr->GetArr())
	{
		vec[iIdx] = convert_symbol<T>(pSymInArr);
		++iIdx;
	}

	return vec;
}

template<template<class> class t_vec, typename T=t_real>
static Symbol* vec_to_sym(const t_vec<T>& vec)
{
	SymbolArray* pSym = new SymbolArray();
	pSym->GetArr().reserve(vec.size());

	for(const T& t : vec)
		pSym->GetArr().push_back(create_symbol<T>(t));

	return pSym;
}

template<template<class> class t_mat, template<class> class t_vec, typename T=t_real>
static t_mat<T> sym_to_mat(const Symbol* pSym, bool* pbIsMat=0)
{
	unsigned int iNumCols=0, iNumRows=0;
	if(!is_mat(pSym, &iNumCols, &iNumRows))
	{
		if(pbIsMat) *pbIsMat = 0;
		return t_mat<T>();
	}
	if(pbIsMat) *pbIsMat = 1;

	t_mat<T> mat(iNumRows, iNumCols);
	const SymbolArray* pSymArr = (SymbolArray*)pSym;

	unsigned int iRow=0;
	for(const Symbol* pSymInArr : pSymArr->GetArr())
	{
		t_vec<T> vecRow = sym_to_vec<t_vec>(pSymInArr);
		unsigned int iNumActCols = std::min<unsigned int>(vecRow.size(), iNumCols);

		for(unsigned int iCol=0; iCol<iNumActCols; ++iCol)
			mat(iRow, iCol) = vecRow[iCol];

		// fill rest with 0
		for(unsigned int iCol=iNumActCols; iCol<iNumCols; ++iCol)
			mat(iRow, iCol) = 0.;

		++iRow;
	}

	return mat;
}

template<template<class> class t_mat, typename T=t_real>
static Symbol* mat_to_sym(const t_mat<T>& mat)
{
	unsigned int iNumRows = mat.size1();
	unsigned int iNumCols = mat.size2();

	SymbolArray* pSym = new SymbolArray();
	pSym->GetArr().reserve(iNumRows);

	for(unsigned int iRow=0; iRow<iNumRows; ++iRow)
	{
		SymbolArray* pRow = new SymbolArray();
		pRow->GetArr().reserve(iNumCols);

		for(unsigned int iCol=0; iCol<iNumCols; ++iCol)
		{
			Symbol *pSymVal = create_symbol<T>(mat(iRow, iCol));
			pRow->GetArr().push_back(pSymVal);
		}

		pSym->GetArr().push_back(pRow);
	}

	return pSym;
}

#endif
