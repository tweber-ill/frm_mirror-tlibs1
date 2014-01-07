/*
 * external functions
 * @author tweber
 * @date 2013-2014
 */
#ifndef __EXT_CALLS__
#define __EXT_CALLS__

#include <vector>
#include "symbol.h"
#include "node.h"

// helper functions
extern std::string linenr(const std::string& strErr, const ParseInfo &info);


// typedefs
typedef Symbol*(*t_extcall)(const std::vector<Symbol*>&, ParseInfo&, SymbolTable*);
typedef std::map<std::string, t_extcall> t_mapFkts;


// adding new external calls
extern bool add_ext_call(const std::string& strFkt, t_extcall pExtCall);
extern void add_ext_calls(t_mapFkts&);


// calling external functions from interpreter
extern Symbol* ext_call(const std::string& strFkt,
						const std::vector<Symbol*>& vecSyms,
						ParseInfo &info,
						SymbolTable* pSymTab);

#endif
