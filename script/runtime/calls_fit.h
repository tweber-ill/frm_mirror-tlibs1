/*
 * external fitter functions
 * @author tweber
 * @date jan 2014
 */

#ifndef __SCRIPT_CALLS_FIT_H__
#define __SCRIPT_CALLS_FIT_H__

#include "../symbol.h"


extern void init_ext_fit_calls();


template<typename T1=std::string, typename T2=double>
static std::map<T1, T2> sym_to_map(const Symbol* pSym)
{
	if(pSym->GetType() != SYMBOL_MAP)
		return std::map<T1, T2>();

	SymbolMap* pSymMap = (SymbolMap*)pSym;
	std::map<T1, T2> _map;

	unsigned int iIdx = 0;
	for(const typename SymbolMap::t_map::value_type& pair : pSymMap->m_map)
		_map[pair.first] = pair.second->GetValDouble();

	return _map;
}

#endif
