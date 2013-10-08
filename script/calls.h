/*
 * External Functions
 * @author tweber
 */
#ifndef __EXT_CALLS__
#define __EXT_CALLS__

#include <vector>
#include "symbol.h"

extern Symbol* ext_call(const std::string& strFkt,
						const std::vector<Symbol*>& vecSyms,
						SymbolTable* pSymTab);

#endif
