/*
 * Parse & Runtime Info
 * @author tweber
 * @date 12-oct-14
 */

#include "info.h"

ParseInfo::ParseInfo()
{
	pmapModules = new t_mods();
	pGlobalSyms = new SymbolTable();
	phandles = new HandleManager();
	pmutexGlobal = new std::mutex();
	pmutexInterpreter = new std::mutex();
	//pmutexCrit = new std::mutex();
}

ParseInfo::~ParseInfo()
{
	if(!bDestroyParseInfo)
		return;

	if(phandles) { delete phandles; phandles=0; }
	if(pGlobalSyms) { delete pGlobalSyms; pGlobalSyms=0; }
	if(pmutexGlobal) { delete pmutexGlobal; pmutexGlobal=0; }
	if(pmutexInterpreter) { delete pmutexInterpreter; pmutexInterpreter=0; }

	if(pmapModules)
	{
		for(ParseInfo::t_mods::value_type vals : *pmapModules)
		{
			if(vals.second)
				delete vals.second;
		}

		pmapModules->clear();
		delete pmapModules;
		pmapModules = 0;
	}
}

void ParseInfo::PushTraceback(std::string&& strTrace)
{
	if(!bEnableDebug) return;

	t_oneTraceback *pStck = 0;
	{
		std::lock_guard<std::mutex> lck(*pmutexInterpreter);
		pStck = &stckTraceback[std::this_thread::get_id()];
	}

	pStck->push_front(std::move(strTrace));
}

void ParseInfo::PopTraceback()
{
	if(!bEnableDebug) return;

	t_oneTraceback *pStck = 0;
	{
		std::lock_guard<std::mutex> lck(*pmutexInterpreter);
		pStck = &stckTraceback[std::this_thread::get_id()];
	}

	pStck->pop_front();
}
