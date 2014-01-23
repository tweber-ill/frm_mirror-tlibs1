/*
 * external math functions
 * @author tweber
 * @date 2013
 */

#ifndef __SCRIPT_CALLS_MATH_H__
#define __SCRIPT_CALLS_MATH_H__

#include "../symbol.h"

extern bool is_vec(const Symbol* pSym);
extern bool is_mat(const Symbol* pSym, unsigned int *piNumCols=0, unsigned int *piNumRows=0);

template<template<class> class t_vec, typename T=double>
static t_vec<T> sym_to_vec(const Symbol* pSym)
{
	if(pSym->GetType() != SYMBOL_ARRAY)
		return t_vec<T>();

	SymbolArray* pSymArr = (SymbolArray*)pSym;
	t_vec<T> vec(pSymArr->m_arr.size());

	unsigned int iIdx = 0;
	for(const Symbol* pSymInArr : pSymArr->m_arr)
	{
		vec[iIdx] = pSymInArr->GetValDouble();
		++iIdx;
	}

	return vec;
}

template<template<class> class t_vec, typename T=double>
static Symbol* vec_to_sym(const t_vec<T>& vec)
{
	SymbolArray* pSym = new SymbolArray();
	pSym->m_arr.reserve(vec.size());

	for(const T& t : vec)
		pSym->m_arr.push_back(new SymbolDouble(t));

	return pSym;
}

template<template<class> class t_mat, template<class> class t_vec, typename T=double>
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
	for(const Symbol* pSymInArr : pSymArr->m_arr)
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

template<template<class> class t_mat, typename T=double>
static Symbol* mat_to_sym(const t_mat<T>& mat)
{
	unsigned int iNumRows = mat.size1();
	unsigned int iNumCols = mat.size2();

	SymbolArray* pSym = new SymbolArray();
	pSym->m_arr.reserve(iNumRows);

	for(unsigned int iRow=0; iRow<iNumRows; ++iRow)
	{
		SymbolArray* pRow = new SymbolArray();
		pRow->m_arr.reserve(iNumCols);

		for(unsigned int iCol=0; iCol<iNumCols; ++iCol)
		{
			SymbolDouble *pSymVal = new SymbolDouble(mat(iRow, iCol));
			pRow->m_arr.push_back(pSymVal);
		}

		pSym->m_arr.push_back(pRow);
	}

	return pSym;
}

extern void init_ext_math_calls();

#endif
