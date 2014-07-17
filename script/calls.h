/*
 * external functions
 * @author tweber
 * @date 2013-2014
 */
#ifndef __EXT_CALLS__
#define __EXT_CALLS__

#include "types.h"

#include <vector>
#include "symbol.h"
#include "node.h"

// helper functions
extern t_string linenr(const t_string& strErr, const ParseInfo &info);


// typedefs
typedef Symbol*(*t_extcall)(const std::vector<Symbol*>&, ParseInfo&, SymbolTable*);
typedef std::map<t_string, t_extcall> t_mapFkts;


// adding new external calls
extern bool add_ext_call(const t_string& strFkt, t_extcall pExtCall);
extern void add_ext_calls(t_mapFkts&);

extern bool has_ext_call(const t_string& strFkt);
extern const t_mapFkts* get_ext_calls();

// calling external functions from interpreter
extern Symbol* ext_call(const t_string& strFkt,
						const std::vector<Symbol*>& vecSyms,
						ParseInfo &info,
						SymbolTable* pSymTab);

#endif
