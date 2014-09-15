/*
 * external functions
 * @author tweber
 * @date 2013-2014
 */

#include "helper/flags.h"
#include "calls.h"
#include "parseobj.h"
#include "script_helper.h"



// external functions
static t_mapFkts g_mapFkts;


extern t_string linenr(const ParseInfo &info)
{
	if(info.pCurCaller)
		return info.pCurCaller->linenr(info);
	return T_STR"";
}


extern Symbol* ext_call(const t_string& strFkt,
			const std::vector<Symbol*>& vecSyms,
			ParseInfo& info,
			SymbolTable* pSymTab)
{
	t_mapFkts::iterator iter = g_mapFkts.find(strFkt);
	if(iter == g_mapFkts.end())
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(info)
			<< "Tried to call unknown function \""
			<< strFkt << "\"."
			<< std::endl;
		throw Err(ostrErr.str(),0);
	}

	return (*iter).second(vecSyms, info, pSymTab);
}

// --------------------------------------------------------------------------------

extern bool add_ext_call(const t_string& strFkt, t_extcall pExtCall)
{
	bool bInserted = g_mapFkts.insert(t_mapFkts::value_type(strFkt, pExtCall)).second;
	return bInserted;
}

extern void add_ext_calls(t_mapFkts& mapFkt)
{
	g_mapFkts.insert(mapFkt.begin(), mapFkt.end());
}

extern bool has_ext_call(const t_string& strFkt)
{
	return g_mapFkts.find(strFkt) != g_mapFkts.end();
}

extern const t_mapFkts* get_ext_calls()
{
	return &g_mapFkts;
}
