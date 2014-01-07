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
			// ignore non-array arguments
			if(pArr->GetType() != SYMBOL_ARRAY)
				continue;

			fkt_plot(((SymbolArray*)pArr)->m_arr, info, pSymTab);
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
		{
			SymbolMap *pParamMap = (SymbolMap*)vecSyms[iNumSyms-1];

			bool bHasVal = 0;
			std::string strTitle = pParamMap->GetStringVal("title", &bHasVal);
			if(bHasVal) g_plot.SetTitle(strTitle.c_str());

			std::string strXLab = pParamMap->GetStringVal("xlabel", &bHasVal);
			if(bHasVal) g_plot.SetXLabel(strXLab.c_str());

			std::string strYLab = pParamMap->GetStringVal("ylabel", &bHasVal);
			if(bHasVal) g_plot.SetYLabel(strYLab.c_str());

			std::string strStyle = pParamMap->GetStringVal("style", &bHasVal);
			if(bHasVal) obj.bConnectLines = (strStyle=="line");

			std::string strLegend = pParamMap->GetStringVal("legend", &bHasVal);
			if(bHasVal) obj.strLegend = strLegend;

			XYLimits lim = get_plot_limits(pParamMap);
			if(lim.bHasX) g_plot.SetXRange(lim.dMinX, lim.dMaxX);
			if(lim.bHasY) g_plot.SetYRange(lim.dMinY, lim.dMaxY);
		}

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
			bool bHasVal = 0;
			std::string strTitle = pMapParam->GetStringVal("title", &bHasVal);
			if(bHasVal) g_plot.SetTitle(strTitle.c_str());

			std::string strXLab = pMapParam->GetStringVal("xlabel", &bHasVal);
			if(bHasVal) g_plot.SetXLabel(strXLab.c_str());

			std::string strYLab = pMapParam->GetStringVal("ylabel", &bHasVal);
			if(bHasVal) g_plot.SetYLabel(strYLab.c_str());

			XYLimits lim = get_plot_limits(pMapParam);
			if(lim.bHasX) { dRMinX = lim.dMinX; dRMaxX = lim.dMaxX; }
			if(lim.bHasY) { dRMinY = lim.dMinY; dRMaxY = lim.dMaxY; }
			if(lim.bHasCB) g_plot.SetColorBarRange(lim.dMinCB, lim.dMaxCB, lim.bCBCyclic);
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

// --------------------------------------------------------------------------------



extern void init_ext_plot_calls()
{
	t_mapFkts mapFkts =
	{
		t_mapFkts::value_type("plot", fkt_plot),
		t_mapFkts::value_type("plot2d", fkt_plot2d),
	};

	add_ext_calls(mapFkts);
}
