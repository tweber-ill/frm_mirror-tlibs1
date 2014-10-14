/*
 * external file functions
 * @author tweber
 * @date dec-2013
 */

#include "../types.h"
#include "../helper/string.h"
#include "../helper/file.h"
#include "../helper/log.h"
#include "calls_file.h"
#include "../calls.h"
#include "../loader/loadtxt.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>

// --------------------------------------------------------------------------------
// file operations
static Symbol* fkt_file_exists(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING}, {0}, "file_exists"))
		return 0;

	const t_string& strFile = ((SymbolString*)vecSyms[0])->GetVal();
	t_ifstream ifstr(strFile);

	bool bFileExists = ifstr.is_open();
	return new SymbolInt(bFileExists);
}

static Symbol* fkt_read_file(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING}, {0}, "read_file"))
		return 0;

	const t_string& strFile = ((SymbolString*)vecSyms[0])->GetVal();

	t_ifstream ifstr(strFile);
	if(!ifstr.is_open())
		return 0;

	//std::streampos iFileSize = get_file_size<t_char>(ifstr);
	t_ostringstream ostr;

	std::copy(std::istreambuf_iterator<t_char>(ifstr),
				std::istreambuf_iterator<t_char>(),
				std::ostreambuf_iterator<t_char>(ostr));

	return new SymbolString(ostr.str());
}

static Symbol* fkt_write_file(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_ANY}, {0,0}, "write_file"))
		return 0;

	const t_string& strFile = ((SymbolString*)vecSyms[0])->GetVal();
	t_ofstream ofstr(strFile);
	if(!ofstr.is_open())
		return new SymbolInt(0);


	t_string* pStr = 0;
	bool bAllocatedStr = 0;
	if(vecSyms[1]->GetType() == SYMBOL_STRING)
		pStr = &((SymbolString*)vecSyms[1])->GetVal();
	else
	{
		pStr = new std::string;
		bAllocatedStr = 1;
		*pStr = vecSyms[1]->print();
	}

	std::copy(pStr->begin(), pStr->end(),
				std::ostreambuf_iterator<t_char>(ofstr));

	if(bAllocatedStr)
		delete pStr;
	return new SymbolInt(1);
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// loading and saving of .dat files

static Symbol* fkt_loadtxt(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING}, {0}, "loadtxt"))
		return 0;

	const t_string& strFile = ((SymbolString*)vecSyms[0])->GetVal();
	SymbolArray *pArr = new SymbolArray();
	LoadTxt dat;

	bool bLoaded = dat.Load(WSTR_TO_STR(strFile).c_str());
	if(!bLoaded)
	{
		log_err(linenr(runinfo), "loadtxt could not open \"", strFile, "\".");
		return pArr;
	}

	pArr->GetArr().reserve(dat.GetColCnt());
	for(unsigned int iCol=0; iCol<dat.GetColCnt(); ++iCol)
	{
		const unsigned int iColLen = dat.GetColLen();
		const t_real *pCol = dat.GetColumn(iCol);

		SymbolArray *pArrCol = new SymbolArray;
		pArrCol->GetArr().reserve(iColLen);

		for(unsigned int iRow=0; iRow<iColLen; ++iRow)
		{
			SymbolReal* pSymD = new SymbolReal();
			pSymD->SetVal(pCol[iRow]);

			pArrCol->GetArr().push_back(pSymD);
		}

		pArrCol->UpdateIndices();
		pArr->GetArr().push_back(pArrCol);
	}

	// load the parameter map
	typedef std::map<std::string, std::string> tmapcomm;
	tmapcomm mapComm = dat.GetCommMapSingle();

	SymbolMap *pSymMap = new SymbolMap();
	for(const tmapcomm::value_type &val : mapComm)
	{
		t_string strKey = STR_TO_WSTR(val.first);

		SymbolString *pSymStrVal = new SymbolString;
		pSymStrVal->SetVal(STR_TO_WSTR(val.second));

		pSymMap->GetMap().insert(SymbolMap::t_map::value_type(strKey, pSymStrVal));
	}
	pArr->GetArr().push_back(pSymMap);

	pArr->UpdateIndices();
	return pArr;
}

static void get_2darr_size(const SymbolArray* pArr,
				unsigned int& iColLen, unsigned int& iRowLen)
{
	iColLen = pArr->GetArr().size();
	iRowLen = 0;

	if(iColLen)
	{
		// look for first real array (not the parameter map)
		for(unsigned int iCol=0; iCol<iColLen; ++iCol)
		{
			Symbol* pSym = pArr->GetArr()[iCol];
			if(pSym->GetType() == SYMBOL_ARRAY)
			{
				iRowLen = ((SymbolArray*)pSym)->GetArr().size();
				break;
			}
		}
	}

	unsigned int iNonArray = 0;
	// don't count the parameter map
	for(unsigned int iCol=0; iCol<iColLen; ++iCol)
	{
		Symbol* pSym = pArr->GetArr()[iCol];
		if(pSym->GetType() != SYMBOL_ARRAY)
			++iNonArray;
	}

	iColLen -= iNonArray;
}

static std::string get_2darr_strval(const SymbolArray* pArr,
				unsigned int iCol, unsigned int iRow)
{
	unsigned int iColLen = pArr->GetArr().size();
	if(iCol >= iColLen)
		return "0";

	bool bFoundCol = 0;
	unsigned int iColRealArray = 0;
	for(unsigned int iCurCol=0; iCurCol<iColLen; ++iCurCol)
	{
		Symbol *pSym = pArr->GetArr()[iCurCol];
		if(pSym->GetType() == SYMBOL_ARRAY)
		{
			if(iColRealArray == iCol)
			{
				bFoundCol = 1;
				iCol = iCurCol;
				break;
			}

			++iColRealArray;
		}
	}

	if(!bFoundCol)
	{
		log_err("Invalid column index: ", iCol, ".");
		return "0";
	}


	Symbol *pSym = pArr->GetArr()[iCol];

	const std::vector<Symbol*>& veccol = ((SymbolArray*)pSym)->GetArr();
	if(iRow >= veccol.size())
		return "0";

	return veccol[iRow]->print();
}

static Symbol* fkt_savetxt(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, RuntimeInfo &runinfo, SymbolTable* pSymTab)
{
	if(!check_args(runinfo, vecSyms, {SYMBOL_STRING, SYMBOL_ARRAY}, {0,0}, "savetxt"))
		return 0;

	const t_string& strFile = ((SymbolString*)vecSyms[0])->GetVal();
	SymbolArray* pArr = (SymbolArray*)vecSyms[1];

	t_ofstream ofstr(WSTR_TO_STR(strFile).c_str());
	if(!ofstr.is_open())
	{
		log_err(linenr(runinfo), "Cannot open \"", strFile, "\".");
		return 0;
	}

	// save parameter map
	for(unsigned int iCol=0; iCol<pArr->GetArr().size(); ++iCol)
	{
		Symbol *_pSymMap = pArr->GetArr()[iCol];
		if(_pSymMap && _pSymMap->GetType() == SYMBOL_MAP)
		{
			SymbolMap *pSymMap = (SymbolMap*)_pSymMap;

			for(const SymbolMap::t_map::value_type& val : pSymMap->GetMap())
			{
				ofstr << "# " << val.first.strKey << " : "
						<< (val.second?val.second->print():T_STR"") << "\n";
			}
		}
	}

	unsigned int iColLen=0, iRowLen=0;
	get_2darr_size(pArr, iColLen, iRowLen);
	//G_COUT << "col len: " << iColLen << ", row len: " << iRowLen << std::endl;

	for(unsigned int iRow=0; iRow<iRowLen; ++iRow)
	{
		for(unsigned int iCol=0; iCol<iColLen; ++iCol)
			ofstr << std::setw(20) << std::left << get_2darr_strval(pArr, iCol, iRow) << " ";
		ofstr << "\n";
	}

	ofstr.flush();
	ofstr.close();
	return 0;
}

// --------------------------------------------------------------------------------



extern void init_ext_file_calls()
{
	t_mapFkts mapFkts =
	{
		// general files
		//t_mapFkts::value_type(T_STR"file_in", fkt_file_in),
		//t_mapFkts::value_type(T_STR"file_out", fkt_file_out),
		//t_mapFkts::value_type(T_STR"file_close", fkt_file_close),
		t_mapFkts::value_type(T_STR"read_file", fkt_read_file),
		t_mapFkts::value_type(T_STR"write_file", fkt_write_file),

		t_mapFkts::value_type(T_STR"file_exists", fkt_file_exists),

		// dat files
		t_mapFkts::value_type(T_STR"loadtxt", fkt_loadtxt),
		t_mapFkts::value_type(T_STR"savetxt", fkt_savetxt),
	};

	add_ext_calls(mapFkts);
}
