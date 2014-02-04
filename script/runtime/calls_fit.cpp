/*
 * external fit functions
 * @author tweber
 * @date jan 2014
 */

#include "calls_fit.h"
#include "calls_math.h"
#include "../calls.h"
#include "../node.h"


template<typename T=double> using t_stdvec = std::vector<T>;

// --------------------------------------------------------------------------------
// fitting

#include "../fitter/fitter.h"
#include "../fitter/chi2.h"
#include <algorithm>
#include <exception>

#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>


class GenericModel : public FunctionModel
{
protected:
	NodeFunction *m_pFkt;
	ParseInfo *m_pinfo;
	SymbolTable *m_pCallerSymTab;

	std::string m_strFreeParam;
	std::vector<std::string> m_vecParamNames;
	std::vector<Symbol*> m_vecSyms;

public:
	GenericModel(const GenericModel& mod)
			: m_pinfo(new ParseInfo(*mod.m_pinfo)),
			  m_pCallerSymTab(mod.m_pCallerSymTab),
			  m_pFkt((NodeFunction*)mod.m_pFkt/*->clone()*/)
	{
		m_pinfo->bDestroyParseInfo = 0;
		m_strFreeParam = mod.m_strFreeParam;
		m_vecParamNames = mod.m_vecParamNames;

		m_vecSyms.reserve(m_vecParamNames.size()+1);
		for(unsigned int i=0; i<m_vecParamNames.size()+1; ++i)
			m_vecSyms.push_back(new SymbolDouble(mod.m_vecSyms[i]->GetValDouble()));
	}

	GenericModel(const NodeFunction *pFkt, ParseInfo& info, SymbolTable *pCallerSymTab)
				: m_pinfo(new ParseInfo(info)),
				  m_pCallerSymTab(pCallerSymTab),
				  m_pFkt((NodeFunction*)pFkt/*->clone()*/)
	{
		m_pinfo->bDestroyParseInfo = 0;
		std::vector<std::string> vecParams = m_pFkt->GetParamNames();
		m_strFreeParam = vecParams[0];

		m_vecParamNames.resize(vecParams.size()-1);
		std::copy(vecParams.begin()+1, vecParams.end(), m_vecParamNames.begin());

		m_vecSyms.reserve(m_vecParamNames.size()+1);
		for(unsigned int i=0; i<m_vecParamNames.size()+1; ++i)
			m_vecSyms.push_back(new SymbolDouble(0.));

		/*std::cout << "free param: " << m_strFreeParam << std::endl;
		std::cout << "args: ";
		for(const std::string& strName : m_vecParamNames)
			std::cout << strName << ", ";
		std::cout << std::endl;*/
	}

	virtual ~GenericModel()
	{
		//if(m_pFkt) { delete m_pFkt; m_pFkt=0; }
		if(m_pinfo) { delete m_pinfo; m_pinfo=0; }
	}

	virtual bool SetParams(const std::vector<double>& vecParams)
	{
		if(vecParams.size() != m_vecParamNames.size())
		{
			std::cerr << "Error: Parameter array length mismatch."
					<< std::endl;
			return 0;
		}

		for(unsigned int iParam=0; iParam<vecParams.size(); ++iParam)
		{
			//const std::string& strName = m_vecParamNames[iParam];
			double dVal = vecParams[iParam];

			((SymbolDouble*)m_vecSyms[iParam+1])->m_dVal = dVal;
		}
	}

	virtual double operator()(double x) const
	{
		((SymbolDouble*)m_vecSyms[0])->m_dVal = x;
		m_pFkt->SetArgSyms(&m_vecSyms);

		double dRetVal = 0.;
		Symbol *pSymRet = m_pFkt->eval(*m_pinfo, 0);
		m_pinfo->bWantReturn = 0;
		if(pSymRet)
			dRetVal = pSymRet->GetValDouble();
		safe_delete(pSymRet, m_pCallerSymTab, m_pinfo->pGlobalSyms);

		return dRetVal;
	}

	virtual GenericModel* copy() const
	{ return new GenericModel(*this); }
	virtual std::string print(bool bFillInSyms=true) const
	{ return "<not implemented>"; }
	virtual const char* GetModelName() const
	{ return "generic fitter model"; }
	virtual std::vector<std::string> GetParamNames() const
	{ return m_vecParamNames; }

	virtual std::vector<double> GetParamValues() const
	{ throw "Called invalid function in generic fitter model"; }
	virtual std::vector<double> GetParamErrors() const
	{ throw "Called invalid function in generic fitter model"; }
};


static void get_values(const std::vector<std::string>& vecParamNames,
						const Symbol* pSym,
						std::vector<double>& vec, std::vector<bool>& vecActive)
{
	vecActive.resize(vecParamNames.size());

	if(pSym->GetType() == SYMBOL_ARRAY)
	{
		vec = sym_to_vec<t_stdvec>(pSym);

		for(unsigned int i=0; i<vecActive.size(); ++i)
			vecActive[i] = 1;
	}
	else if(pSym->GetType() == SYMBOL_MAP)
	{
		vec.resize(vecParamNames.size());
		std::map<std::string, double> mymap = sym_to_map<std::string, double>(pSym);

		for(unsigned int iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			const std::string& strKey = vecParamNames[iParam];

			std::map<std::string, double>::iterator iter = mymap.find(strKey);
			if(iter == mymap.end())
			{
				vecActive[iParam] = 0;
			}
			else
			{
				vecActive[iParam] = 1;
				vec[iParam] = iter->second;
			}
		}
	}
}

// fit("function", x, y, yerr, params)
static Symbol* fkt_fit(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()<4 || !is_vec(vecSyms[1]) || !is_vec(vecSyms[2]) || !is_vec(vecSyms[3]))
	{
		std::cerr << linenr("Error", info)
			<< "Invalid arguments for fit."
			<< std::endl;
		return 0;
	}

	if(vecSyms[0]->GetType() != SYMBOL_STRING)
	{
		std::cerr << linenr("Error", info)
			<< "Need a fit function name."
			<< std::endl;
		return 0;
	}

	const std::string& strFkt = ((SymbolString*)vecSyms[0])->m_strVal;
	NodeFunction *pFkt = info.GetFunction(strFkt);
	if(!pFkt)
	{
		std::cerr << linenr("Error", info) << "Invalid function \""
			<< strFkt << "\"."
			<< std::endl;
		return 0;
	}


	GenericModel mod(pFkt, info, pSymTab);
	std::vector<std::string> vecParamNames = mod.GetParamNames();
	const unsigned int iParamSize = vecParamNames.size();


	std::vector<double> vecHints, vecHintsErr;
	std::vector<bool> vecHintsActive, vecHintsErrActive;

	std::vector<double> vecLimMin, vecLimMax;
	std::vector<bool> vecLimMinActive, vecLimMaxActive;

	vecLimMinActive.resize(iParamSize);
	vecLimMaxActive.resize(iParamSize);

	// parameter map
	if(vecSyms.size()==5 && vecSyms[4]->GetType()==SYMBOL_MAP)
	{
		SymbolMap::t_map& mapSym = ((SymbolMap*)vecSyms[4])->m_map;

		SymbolMap::t_map::iterator iterHints = mapSym.find("hints");
		SymbolMap::t_map::iterator iterHintsErr = mapSym.find("hints_errors");

		SymbolMap::t_map::iterator iterLimitsMin = mapSym.find("lower_limits");
		SymbolMap::t_map::iterator iterLimitsMax = mapSym.find("upper_limits");

		if(iterHints != mapSym.end())
			get_values(vecParamNames, iterHints->second, vecHints, vecHintsActive);
		if(iterHintsErr != mapSym.end())
			get_values(vecParamNames, iterHintsErr->second, vecHintsErr, vecHintsErrActive);

		if(iterLimitsMin != mapSym.end())
			get_values(vecParamNames, iterLimitsMin->second, vecLimMin, vecLimMinActive);
		if(iterLimitsMax != mapSym.end())
			get_values(vecParamNames, iterLimitsMax->second, vecLimMax, vecLimMaxActive);
	}


	std::vector<double> vecX = sym_to_vec<t_stdvec>(vecSyms[1]);
	std::vector<double> vecY = sym_to_vec<t_stdvec>(vecSyms[2]);
	std::vector<double> vecYErr = sym_to_vec<t_stdvec>(vecSyms[3]);

	unsigned int iSize = std::min<unsigned int>(vecX.size(), vecY.size());
	iSize = std::min<unsigned int>(iSize, vecYErr.size());

	Chi2Function chi2fkt(&mod, iSize, vecX.data(), vecY.data(), vecYErr.data());


	ROOT::Minuit2::MnUserParameters params;
	for(unsigned int iParam=0; iParam<iParamSize; ++iParam)
	{
		double dHint = 0.;
		double dErr = 0.;

		if(iParam < vecHints.size())
			dHint = vecHints[iParam];
		if(iParam < vecHintsErr.size())
			dErr = vecHintsErr[iParam];

		//std::cout << "hints for " << vecParamNames[iParam] << ": "
		//		<< dHint << " +- " << dErr << std::endl;
		params.Add(vecParamNames[iParam], dHint, dErr);



		double dLimMin = 0.;
		double dLimMax = 0.;

		if(iParam < vecLimMin.size())
			dLimMin = vecLimMin[iParam];
		if(iParam < vecLimMax.size())
			dLimMax = vecLimMax[iParam];

/*		if(vecLimMinActive[iParam])
			std::cout << "lower limit for " << vecParamNames[iParam] << ": "
						<< dLimMin << std::endl;
		if(vecLimMaxActive[iParam])
			std::cout << "upper limit for " << vecParamNames[iParam] << ": "
						<< dLimMax << std::endl;*/

		if(vecLimMinActive[iParam] && vecLimMaxActive[iParam])
			params.SetLimits(vecParamNames[iParam], dLimMin, dLimMax);
		else if(vecLimMinActive[iParam] && vecLimMaxActive[iParam]==0)
			params.SetLowerLimit(vecParamNames[iParam], dLimMin);
		else if(vecLimMinActive[iParam]==0 && vecLimMaxActive[iParam])
			params.SetUpperLimit(vecParamNames[iParam], dLimMax);

		//params.Fix(vecParamNames[iParam]);
	}


	bool bValidFit = 1;

	std::vector<ROOT::Minuit2::FunctionMinimum> minis;
	minis.reserve(2);

	{
		// step 1: free fit (limited)

		ROOT::Minuit2::MnMigrad migrad(chi2fkt, params, 1);
		ROOT::Minuit2::FunctionMinimum mini = migrad();
		bValidFit = mini.IsValid() && mini.HasValidParameters();

		for(const std::string& strSym : vecParamNames)
		{
			params.SetValue(strSym, mini.UserState().Value(strSym));
			params.SetError(strSym, mini.UserState().Error(strSym));
		}

		minis.push_back(mini);
	}


	{
		// step 2: free fit (unlimited)

		for(const std::string& strSym : vecParamNames)
			params.RemoveLimits(strSym);

		ROOT::Minuit2::MnMigrad migrad2(chi2fkt, params, 2);
		ROOT::Minuit2::FunctionMinimum mini = migrad2();
		bValidFit = mini.IsValid() && mini.HasValidParameters();

		minis.push_back(mini);
	}


	const ROOT::Minuit2::FunctionMinimum& lastmini = *minis.rbegin();


	SymbolMap *pSymMap = new SymbolMap();
	for(const std::string& strSym : vecParamNames)
	{
		double dVal = lastmini.UserState().Value(strSym);
		double dErr = lastmini.UserState().Error(strSym);
		dErr = fabs(dErr);

		SymbolArray* pArr = new SymbolArray();
		pArr->m_arr.push_back(new SymbolDouble(dVal));
		pArr->m_arr.push_back(new SymbolDouble(dErr));

		pSymMap->m_map.insert(SymbolMap::t_map::value_type(strSym, pArr));
	}

	bool bFitterDebug = 0;
	if(bFitterDebug)
	{
		std::cerr << "--------------------------------------------------------------------------------" << std::endl;
		unsigned int uiMini=0;
		for(const auto& mini : minis)
		{
			std::cerr << "result of user-defined fit step " << (++uiMini) << std::endl;
			std::cerr << "=================================" << std::endl;
			std::cerr << mini << std::endl;
		}
		std::cerr << "--------------------------------------------------------------------------------" << std::endl;
	}

	if(!bValidFit)
		std::cerr << "Error: Fit invalid!" << std::endl;
	return pSymMap;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// interpolation

#include "../fitter/models/interpolation.h"
namespace ublas = boost::numeric::ublas;


enum FktParam
{
	FKT_SPLINE,
	FKT_BEZIER
};


// bezier(x, y, 128)
// spline(x, y, 128, degree)
static Symbol* _fkt_param(FktParam whichfkt, const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() < 2 || !is_vec(vecSyms[0]) || !is_vec(vecSyms[1]))
	{
		std::cerr << linenr("Error", info) << "Function needs x and y vector arguments."
			<< std::endl;
		return 0;
	}

	unsigned int ideg = 3;
	unsigned int N = 128;

	if(vecSyms.size() > 2)
		N = vecSyms[2]->GetValInt();
	if(vecSyms.size() > 3)
		ideg = vecSyms[3]->GetValInt();

	std::vector<double> vecX = sym_to_vec<t_stdvec>(vecSyms[0]);
	std::vector<double> vecY = sym_to_vec<t_stdvec>(vecSyms[1]);

	unsigned int iSize = std::min(vecX.size(), vecY.size());

	FunctionModel_param* pfkt = 0;
	if(whichfkt == FKT_SPLINE)
		pfkt = new BSpline(iSize, vecX.data(), vecY.data(), ideg);
	else if(whichfkt == FKT_BEZIER)
		pfkt = new Bezier(iSize, vecX.data(), vecY.data());
	else
	{
		std::cerr << linenr("Error", info) << "Unknown parametric function selected."
					<< std::endl;
		return 0;
	}

	SymbolArray* pArrX = new SymbolArray();
	SymbolArray* pArrY = new SymbolArray();
	pArrX->m_arr.reserve(N);
	pArrY->m_arr.reserve(N);

	for(unsigned int i=0; i<N; ++i)
	{
		double t = double(i)/double(N-1);
		ublas::vector<double> vecSpl = (*pfkt)(t);

		if(vecSpl.size() < 2)
			continue;

		pArrX->m_arr.push_back(new SymbolDouble(vecSpl[0]));
		pArrY->m_arr.push_back(new SymbolDouble(vecSpl[1]));
	}

	delete pfkt;

	SymbolArray *pArr = new SymbolArray();
	pArr->m_arr.push_back(pArrX);
	pArr->m_arr.push_back(pArrY);
	return pArr;
}

static Symbol* fkt_bezier(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	return _fkt_param(FKT_BEZIER, vecSyms, info, pSymTab);
}

static Symbol* fkt_spline(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	return _fkt_param(FKT_SPLINE, vecSyms, info, pSymTab);
}
// --------------------------------------------------------------------------------


// ["a" : [val0, val1]]  =>  ["a" : val0]
static Symbol* fkt_map_vec_to_val(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1 || vecSyms[0]->GetType()!=SYMBOL_MAP)
	{
		std::cerr << linenr("Error", info)
					<< "Need a map of vectors."
					<< std::endl;

		return 0;
	}


	// index parameter
	int iIdx = 0;
	if(vecSyms.size()>1)
	{
		iIdx = vecSyms[1]->GetValInt();
		if(iIdx < 0)
		{
			std::cerr << linenr("Warning", info)
						<< "Ignoring negative index."
						<< std::endl;
			iIdx = 0;
		}
	}


	SymbolMap* pMapRet = new SymbolMap();

	for(const SymbolMap::t_map::value_type& pair : ((SymbolMap*)vecSyms[0])->m_map)
	{
		const std::string& strKey = pair.first;
		const Symbol* pSym = pair.second;

		if(pSym->GetType() != SYMBOL_ARRAY)
		{
			std::cerr << linenr("Warning", info)
						<< "Ignoring non-vector variable in map."
						<< std::endl;
			continue;
		}

		const std::vector<Symbol*>& arr = ((SymbolArray*)pSym)->m_arr;
		if(iIdx >= arr.size())
		{
			std::cerr << linenr("Warning", info)
						<< "Ignoring invalid index."
						<< std::endl;
			continue;
		}

		const Symbol* pElem = arr[iIdx];

		pMapRet->m_map.insert(SymbolMap::t_map::value_type(strKey, pElem->clone()));
	}

	return pMapRet;
}

extern void init_ext_fit_calls()
{
	t_mapFkts mapFkts =
	{
		t_mapFkts::value_type("fit", fkt_fit),
//		t_mapFkts::value_type("fit_sin", fkt_fit_sin),
//		t_mapFkts::value_type("fit_gauss", fkt_fit_sin),

		t_mapFkts::value_type("map_vec_to_val", fkt_map_vec_to_val),

		t_mapFkts::value_type("spline", fkt_spline),
		t_mapFkts::value_type("bezier", fkt_bezier),
	};

	add_ext_calls(mapFkts);
}
