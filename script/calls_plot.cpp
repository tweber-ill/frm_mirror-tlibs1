/*
 * external plot functions
 * @author tweber
 * @date dec 2013
 */

#include "calls_plot.h"
#include "calls.h"
#include "helper/gnuplot.h"


// --------------------------------------------------------------------------------
// plotting

#define DEFAULT_TERM "qt";
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
	std::string strVal = pParamMap->GetStringVal("xylimits", &bHasVal);
	if(bHasVal)
	{
		std::istringstream istr(strVal);
		istr >> lim.dMinX >> lim.dMaxX >> lim.dMinY >> lim.dMaxY;

		lim.bHasX = 1;
		lim.bHasY = 1;
	}
	else
	{
		std::string strX = pParamMap->GetStringVal("xlimits", &lim.bHasX);
		std::string strY = pParamMap->GetStringVal("ylimits", &lim.bHasY);

		if(lim.bHasX)
		{
			std::istringstream istrX(strX);
			istrX >> lim.dMinX >> lim.dMaxX;
		}

		if(lim.bHasY)
		{
			std::istringstream istrY(strY);
			istrY >> lim.dMinY >> lim.dMaxY;
		}
	}
	//std::cout << "xlimits: " << lim.dMinX << ", " << lim.dMaxX << std::endl;
	//std::cout << "ylimits: " << lim.dMinY << ", " << lim.dMaxY << std::endl;


	std::string strValCB = pParamMap->GetStringVal("cblimits", &lim.bHasCB);
	std::istringstream istrCB(strValCB);
	istrCB >> lim.dMinCB >> lim.dMaxCB;
	//std::cout << "colorbar: " << lim.dMinCB << ", " << lim.dMaxCB << std::endl;

	bool bHasCyc = 0;
	std::string strValCBCyc = pParamMap->GetStringVal("cbcyclic", &bHasCyc);
	if(bHasCyc)
	{
		std::istringstream istrCyc(strValCBCyc);
		istrCyc >> lim.bCBCyclic;
	}

	return lim;
}


static void set_plot_params(GnuPlot& plot, SymbolMap* pParamMap, PlotObj* pCurPlotObj=0, XYLimits* pLimits=0)
{
	bool bHasVal = 0;
	std::string strTitle = pParamMap->GetStringVal("title", &bHasVal);
	if(bHasVal) plot.SetTitle(strTitle.c_str());

	std::string strXLab = pParamMap->GetStringVal("xlabel", &bHasVal);
	if(bHasVal) plot.SetXLabel(strXLab.c_str());

	std::string strYLab = pParamMap->GetStringVal("ylabel", &bHasVal);
	if(bHasVal) plot.SetYLabel(strYLab.c_str());

	if(pCurPlotObj)
	{
		std::string strStyle = pParamMap->GetStringVal("style", &bHasVal);
		if(bHasVal) pCurPlotObj->bConnectLines = (strStyle=="line");

		std::string strLegend = pParamMap->GetStringVal("legend", &bHasVal);
		if(bHasVal) pCurPlotObj->strLegend = strLegend;
	}

	XYLimits lim = get_plot_limits(pParamMap);
	if(lim.bHasX) plot.SetXRange(lim.dMinX, lim.dMaxX);
	if(lim.bHasY) plot.SetYRange(lim.dMinY, lim.dMaxY);


	// for 2D plot
	if(pLimits)
	{
		*pLimits = lim;
		if(pLimits->bHasCB) plot.SetColorBarRange(pLimits->dMinCB, pLimits->dMaxCB, pLimits->bCBCyclic);
	}


	// terminal
	std::string strTerm = DEFAULT_TERM;
	std::string strUserTerm = pParamMap->GetStringVal("term", &bHasVal);
	if(bHasVal) strTerm = strUserTerm;

	int iPlotWnd = 0;
	int iUserPlotWnd = atoi(pParamMap->GetStringVal("window", &bHasVal).c_str());
	if(bHasVal) iPlotWnd = iUserPlotWnd;

	plot.SetTerminal(iPlotWnd, strTerm.c_str());
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
		std::cerr << linenr("Error", info) << "Invalid call to plot." << std::endl;
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
		std::cerr << linenr("Error", info) << "Invalid call to plot2d." << std::endl;
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
		std::cerr << linenr("Error", info) 
			<< "First argument to fileplot has to be the file name." << std::endl;
		return 0;
	}

	const std::string& strFile = ((SymbolString*)vecSyms[0])->m_strVal;
	g_plot.SetFileTerminal(strFile.c_str());
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
		t_mapFkts::value_type("plot", fkt_plot),
		t_mapFkts::value_type("plot2d", fkt_plot2d),

		t_mapFkts::value_type("fileplot", fkt_fileplot),
		t_mapFkts::value_type("fileplot2d", fkt_fileplot2d),
	};

	add_ext_calls(mapFkts);
}
