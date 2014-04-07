/*
 * external plot functions
 * @author tweber
 * @date dec 2013
 */

#include "../types.h"
#include "calls_plot.h"
#include "../calls.h"
#include "../helper/gnuplot.h"
#include "../helper/string.h"


// --------------------------------------------------------------------------------
// plotting

#define DEFAULT_TERM T_STR"qt";
//#define DEFAULT_TERM T_STR"x11";
static GnuPlot g_plot;

static inline bool is_array_of_arrays(const Symbol* pSym)
{
	if(!pSym) return 0;
	if(pSym->GetType()!=SYMBOL_ARRAY) return 0;

	SymbolArray* pSymArr = (SymbolArray*)pSym;
	if(pSymArr->m_arr.size()==0) return 0;

	Symbol* pSymInArr = pSymArr->m_arr[0];
	return (pSymInArr->GetType()==SYMBOL_ARRAY);
}

static inline bool is_array_of_array_of_arrays(const Symbol* pSym)
{
	if(!pSym) return 0;
	if(pSym->GetType() != SYMBOL_ARRAY)
		return 0;

	if(((SymbolArray*)pSym)->m_arr.size() == 0)
		return 0;

	return is_array_of_arrays(((SymbolArray*)pSym)->m_arr[0]);
}

struct XYLimits
{
	bool bHasX, bHasY, bHasCB, bCBCyclic;
	double dMinX, dMaxX;
	double dMinY, dMaxY;
	double dMinCB, dMaxCB;

	XYLimits() : bHasX(0), bHasY(0), bHasCB(0), bCBCyclic(0)
	{}
};

static XYLimits get_plot_limits(SymbolMap* pParamMap)
{
	XYLimits lim;

	bool bHasVal=0;
	t_string strVal = pParamMap->GetStringVal(T_STR"xylimits", &bHasVal);
	if(bHasVal)
	{
		t_istringstream istr(strVal);
		istr >> lim.dMinX >> lim.dMaxX >> lim.dMinY >> lim.dMaxY;

		lim.bHasX = 1;
		lim.bHasY = 1;
	}
	else
	{
		t_string strX = pParamMap->GetStringVal(T_STR"xlimits", &lim.bHasX);
		t_string strY = pParamMap->GetStringVal(T_STR"ylimits", &lim.bHasY);

		if(lim.bHasX)
		{
			t_istringstream istrX(strX);
			istrX >> lim.dMinX >> lim.dMaxX;
		}

		if(lim.bHasY)
		{
			t_istringstream istrY(strY);
			istrY >> lim.dMinY >> lim.dMaxY;
		}
	}
	//G_COUT << "xlimits: " << lim.dMinX << ", " << lim.dMaxX << std::endl;
	//G_COUT << "ylimits: " << lim.dMinY << ", " << lim.dMaxY << std::endl;


	t_string strValCB = pParamMap->GetStringVal(T_STR"cblimits", &lim.bHasCB);
	t_istringstream istrCB(strValCB);
	istrCB >> lim.dMinCB >> lim.dMaxCB;
	//G_COUT << "colorbar: " << lim.dMinCB << ", " << lim.dMaxCB << std::endl;

	bool bHasCyc = 0;
	t_string strValCBCyc = pParamMap->GetStringVal(T_STR"cbcyclic", &bHasCyc);
	if(bHasCyc)
	{
		t_istringstream istrCyc(strValCBCyc);
		istrCyc >> lim.bCBCyclic;
	}

	return lim;
}


static void set_plot_params(GnuPlot& plot, SymbolMap* pParamMap, PlotObj* pCurPlotObj=0, XYLimits* pLimits=0)
{
	bool bHasVal = 0;
	t_string strTitle = pParamMap->GetStringVal(T_STR"title", &bHasVal);
	if(bHasVal) plot.SetTitle(WSTR_TO_STR(strTitle).c_str());

	t_string strXLab = pParamMap->GetStringVal(T_STR"xlabel", &bHasVal);
	if(bHasVal) plot.SetXLabel(WSTR_TO_STR(strXLab).c_str());

	t_string strYLab = pParamMap->GetStringVal(T_STR"ylabel", &bHasVal);
	if(bHasVal) plot.SetYLabel(WSTR_TO_STR(strYLab).c_str());

	if(pCurPlotObj)
	{
		t_string strStyle = pParamMap->GetStringVal(T_STR"style", &bHasVal);
		if(bHasVal) pCurPlotObj->bConnectLines = (strStyle==T_STR"line");

		t_string strLegend = pParamMap->GetStringVal(T_STR"legend", &bHasVal);
		if(bHasVal) pCurPlotObj->strLegend = WSTR_TO_STR(strLegend);

		t_string strSize = pParamMap->GetStringVal(T_STR"size", &pCurPlotObj->bHasSize);
		if(pCurPlotObj->bHasSize)
		{
			t_istringstream istrSize(strSize);
			istrSize >> pCurPlotObj->dSize;
		}
	}

	XYLimits lim = get_plot_limits(pParamMap);
	if(lim.bHasX) plot.SetXRange(lim.dMinX, lim.dMaxX);
	if(lim.bHasY) plot.SetYRange(lim.dMinY, lim.dMaxY);


	// for 2D plot
	if(pLimits)
	{
		*pLimits = lim;
		if(pLimits->bHasCB)
			plot.SetColorBarRange(pLimits->dMinCB, pLimits->dMaxCB, pLimits->bCBCyclic);
	}


	// terminal
	t_string strTerm = DEFAULT_TERM;
	t_string strUserTerm = pParamMap->GetStringVal(T_STR"term", &bHasVal);
	if(bHasVal) strTerm = strUserTerm;

	int iPlotWnd = 0;
	int iUserPlotWnd = atoi(WSTR_TO_STR(pParamMap->GetStringVal(T_STR"window", &bHasVal)).c_str());
	if(bHasVal) iPlotWnd = iUserPlotWnd;

	plot.SetTerminal(iPlotWnd, WSTR_TO_STR(strTerm).c_str());
}

static Symbol* fkt_plot(const std::vector<Symbol*>& vecSyms,
			ParseInfo& info, SymbolTable* pSymTab)
{
	g_plot.Init();

	unsigned int iNumSyms = vecSyms.size();

	// plot([[x, y, yerr, xerr, mapParams], ...]);
	if(iNumSyms==1 && is_array_of_array_of_arrays(vecSyms[0]))
		return fkt_plot(((SymbolArray*)vecSyms[0])->m_arr, info, pSymTab);
	// plot([x, y, yerr, xerr, mapParams], [x2, y2, yerr2, xerr2, mapParams2], ...)
	else if(iNumSyms>=1 && is_array_of_arrays(vecSyms[0]))
	{
		g_plot.StartPlot();
		for(Symbol *pArr : vecSyms)
		{
			SymbolType symType = pArr->GetType();

			// ignore non-array arguments
			if(symType == SYMBOL_ARRAY)
				fkt_plot(((SymbolArray*)pArr)->m_arr, info, pSymTab);
			else if(symType == SYMBOL_MAP)
				set_plot_params(g_plot, (SymbolMap*)pArr);
		}
		g_plot.FinishPlot();
	}
	// plot(x, y, yerr, xerr, mapParams)
	else if(iNumSyms >= 2 && vecSyms[0]->GetType()==SYMBOL_ARRAY && vecSyms[1]->GetType()==SYMBOL_ARRAY)
	{
		std::vector<double> vecX = ((SymbolArray*)vecSyms[0])->ToDoubleArray();
		std::vector<double> vecY = ((SymbolArray*)vecSyms[1])->ToDoubleArray();
		std::vector<double> vecYErr, vecXErr;

		if(iNumSyms >= 3 && vecSyms[2]->GetType()==SYMBOL_ARRAY)
			vecYErr = ((SymbolArray*)vecSyms[2])->ToDoubleArray();
		if(iNumSyms >= 4 && vecSyms[3]->GetType()==SYMBOL_ARRAY)
			vecXErr = ((SymbolArray*)vecSyms[3])->ToDoubleArray();


		PlotObj obj;
		obj.vecX = vecX;
		obj.vecY = vecY;
		obj.vecErrX = vecXErr;
		obj.vecErrY = vecYErr;

		// parameter map given as last argument
		if(vecSyms[iNumSyms-1]->GetType()==SYMBOL_MAP)
			set_plot_params(g_plot, (SymbolMap*)vecSyms[iNumSyms-1], &obj);

		g_plot.StartPlot();
		g_plot.AddLine(obj);
		g_plot.FinishPlot();
	}
	else
	{
		G_CERR << linenr(T_STR"Error", info) << "Invalid call to plot." << std::endl;
		return 0;
	}

	return 0;
}

static Symbol* fkt_plot2d(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	g_plot.Init();

	unsigned int iNumSyms = vecSyms.size();

	// e.g. plot2d([[1,2],[3,4]], params)
	if(iNumSyms >= 1 && is_array_of_arrays(vecSyms[0]))
	{
		SymbolArray* pArr = (SymbolArray*)vecSyms[0];
		SymbolMap* pMapParam = 0;

		// has parameter map
		if(iNumSyms == 2 && vecSyms[1]->GetType()==SYMBOL_MAP)
			pMapParam = (SymbolMap*)vecSyms[1];

		std::vector<std::vector<double> > vecXY;
		vecXY.reserve(pArr->m_arr.size());

		for(const Symbol* _pX : pArr->m_arr)
		{
			SymbolArray* pX = (SymbolArray*)_pX;
			std::vector<double> vecX = pX->ToDoubleArray();

			vecXY.push_back(vecX);
		}

		double dRMinX=1., dRMaxX=-1., dRMinY=1., dRMaxY=-1.;
		if(pMapParam)
		{
			XYLimits lim;
			set_plot_params(g_plot, (SymbolMap*)pMapParam, 0, &lim);
        	        if(lim.bHasX) { dRMinX = lim.dMinX; dRMaxX = lim.dMaxX; }
	                if(lim.bHasY) { dRMinY = lim.dMinY; dRMaxY = lim.dMaxY; }
		}

		g_plot.SimplePlot2d(vecXY, dRMinX, dRMaxX, dRMinY, dRMaxY);
	}
	else
	{
		G_CERR << linenr(T_STR"Error", info) << "Invalid call to plot2d." << std::endl;
		return 0;
	}

	return 0;
}


static Symbol* _fkt_fileplot(const std::vector<Symbol*>& vecSyms, 
				ParseInfo& info, SymbolTable* pSymTab,
				Symbol* (*pPltFkt)(const::std::vector<Symbol*>&, ParseInfo&, SymbolTable*))
{
	g_plot.Init();
	if(vecSyms.size() < 1 || vecSyms[0]->GetType()!=SYMBOL_STRING)
	{
		G_CERR << linenr(T_STR"Error", info)
			<< "First argument to fileplot has to be the file name." << std::endl;
		return 0;
	}

	const t_string& strFile = ((SymbolString*)vecSyms[0])->m_strVal;
	g_plot.SetFileTerminal(WSTR_TO_STR(strFile).c_str());
	g_plot.LockTerminal();

	std::vector<Symbol*> vecPlot;
	vecPlot.reserve(vecSyms.size()-1);
	for(unsigned int iSym=1; iSym<vecSyms.size(); ++iSym)
		vecPlot.push_back(vecSyms[iSym]);

	Symbol* pSymRet = pPltFkt(vecPlot, info, pSymTab);

	g_plot.UnlockTerminal();
	return pSymRet;
}

static Symbol* fkt_fileplot(const std::vector<Symbol*>& vecSyms,
				ParseInfo& info, SymbolTable* pSymTab)
{
	_fkt_fileplot(vecSyms, info, pSymTab, fkt_plot);
} 

static Symbol* fkt_fileplot2d(const std::vector<Symbol*>& vecSyms,
				ParseInfo& info, SymbolTable* pSymTab)
{
	_fkt_fileplot(vecSyms, info, pSymTab, fkt_plot2d);
} 


// --------------------------------------------------------------------------------



extern void init_ext_plot_calls()
{
	t_mapFkts mapFkts =
	{
		t_mapFkts::value_type(T_STR"plot", fkt_plot),
		t_mapFkts::value_type(T_STR"plot2d", fkt_plot2d),

		t_mapFkts::value_type(T_STR"fileplot", fkt_fileplot),
		t_mapFkts::value_type(T_STR"fileplot2d", fkt_fileplot2d),
	};

	add_ext_calls(mapFkts);
}
