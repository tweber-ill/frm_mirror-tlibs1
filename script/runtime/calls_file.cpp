/*
 * external file functions
 * @author tweber
 * @date dec-2013
 */

#include "../types.h"
#include "../helper/string.h"
#include "calls_file.h"
#include "../calls.h"
#include "../loader/loadtxt.h"
#include <sstream>
#include <fstream>

// --------------------------------------------------------------------------------
// loading and saving of .dat files

static Symbol* fkt_loadtxt(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		G_CERR << linenr(T_STR"Error", info) << "loadtxt takes exactly one argument." << std::endl;
		return 0;
	}

	if(vecSyms[0]->GetType() != SYMBOL_STRING)
	{
		G_CERR << linenr(T_STR"Error", info) << "loadtxt needs a string argument." << std::endl;
		return 0;
	}

	const t_string& strFile = ((SymbolString*)vecSyms[0])->m_strVal;


	SymbolArray *pArr = new SymbolArray();

	LoadTxt dat;

	bool bLoaded = dat.Load(WSTR_TO_STR(strFile).c_str());
	if(!bLoaded)
	{
		G_CERR << linenr(T_STR"Error", info) << "loadtxt could not open \"" << strFile << "\"." << std::endl;
		return pArr;
	}

	pArr->m_arr.reserve(dat.GetColCnt());
	for(unsigned int iCol=0; iCol<dat.GetColCnt(); ++iCol)
	{
		const unsigned int iColLen = dat.GetColLen();
		const double *pCol = dat.GetColumn(iCol);

		SymbolArray *pArrCol = new SymbolArray;
		pArrCol->m_arr.reserve(iColLen);

		for(unsigned int iRow=0; iRow<iColLen; ++iRow)
		{
			SymbolDouble* pSymD = new SymbolDouble();
			pSymD->m_dVal = pCol[iRow];

			pArrCol->m_arr.push_back(pSymD);
		}

		pArr->m_arr.push_back(pArrCol);
	}

	// load the parameter map
	typedef std::map<std::string, std::string> tmapcomm;
	tmapcomm mapComm = dat.GetCommMapSingle();

	SymbolMap *pSymMap = new SymbolMap();
	for(const tmapcomm::value_type &val : mapComm)
	{
		t_string strKey = STR_TO_WSTR(val.first);

		SymbolString *pSymStrVal = new SymbolString;
		pSymStrVal->m_strVal = STR_TO_WSTR(val.second);

		pSymMap->m_map.insert(SymbolMap::t_map::value_type(strKey, pSymStrVal));
	}
	pArr->m_arr.push_back(pSymMap);

	return pArr;
}

static void get_2darr_size(const SymbolArray* pArr,
				unsigned int& iColLen, unsigned int& iRowLen)
{
	iColLen = pArr->m_arr.size();
	iRowLen = 0;

	if(iColLen)
	{
		// look for first real array (not the parameter map)
		for(unsigned int iCol=0; iCol<iColLen; ++iCol)
		{
			Symbol* pSym = pArr->m_arr[iCol];
			if(pSym->GetType() == SYMBOL_ARRAY)
			{
				iRowLen = ((SymbolArray*)pSym)->m_arr.size();
				break;
			}
		}
	}

	unsigned int iNonArray = 0;
	// don't count the parameter map
	for(unsigned int iCol=0; iCol<iColLen; ++iCol)
	{
		Symbol* pSym = pArr->m_arr[iCol];
		if(pSym->GetType() != SYMBOL_ARRAY)
			++iNonArray;
	}

	iColLen -= iNonArray;
}

static double get_2darr_val(const SymbolArray* pArr,
				unsigned int iCol, unsigned int iRow)
{
	unsigned int iColLen = pArr->m_arr.size();
	if(iCol >= iColLen)
		return 0.;

	bool bFoundCol = 0;
	unsigned int iColRealArray = 0;
	for(unsigned int iCurCol=0; iCurCol<iColLen; ++iCurCol)
	{
		Symbol *pSym = pArr->m_arr[iCurCol];
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
		G_CERR << "Error: Invalid column index: " << iCol << "." << std::endl;
		return 0.;
	}


	Symbol *pSym = pArr->m_arr[iCol];

	const std::vector<Symbol*>& veccol = ((SymbolArray*)pSym)->m_arr;
	if(iRow >= veccol.size())
		return 0.;

	return veccol[iRow]->GetValDouble();
}

static Symbol* fkt_savetxt(const std::vector<Symbol*>& vecSyms,
							ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=2 ||
		(vecSyms[0]->GetType()!=SYMBOL_STRING || vecSyms[1]->GetType()!=SYMBOL_ARRAY))
	{
		G_CERR << linenr(T_STR"Error", info) << "savetxt takes two arguments (file name, 2d array)." << std::endl;
		return 0;
	}

	const t_string& strFile = ((SymbolString*)vecSyms[0])->m_strVal;
	SymbolArray* pArr = (SymbolArray*)vecSyms[1];


	t_ofstream ofstr(WSTR_TO_STR(strFile).c_str());
	if(!ofstr.is_open())
	{
		G_CERR << linenr(T_STR"Error", info) << "Cannot open \"" << strFile << "\"." << std::endl;
		return 0;
	}

	// save parameter map
	for(unsigned int iCol=0; iCol<pArr->m_arr.size(); ++iCol)
	{
		Symbol *_pSymMap = pArr->m_arr[iCol];
		if(_pSymMap && _pSymMap->GetType() == SYMBOL_MAP)
		{
			SymbolMap *pSymMap = (SymbolMap*)_pSymMap;

			for(const SymbolMap::t_map::value_type& val : pSymMap->m_map)
			{
				ofstr << "# " << val.first << " : "
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
			ofstr << get_2darr_val(pArr, iCol, iRow) << " ";
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
			// dat files
			t_mapFkts::value_type(T_STR"loadtxt", fkt_loadtxt),
			t_mapFkts::value_type(T_STR"savetxt", fkt_savetxt),
	};

	add_ext_calls(mapFkts);
}
