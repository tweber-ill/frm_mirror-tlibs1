/*
 * external fit functions
 * @author tweber
 * @date jan 2014
 */

#include "../types.h"
#include "../helper/string.h"
#include "calls_fit.h"
#include "../calls.h"
#include "../node.h"


template<typename T> using t_stdvec = std::vector<T>;

// --------------------------------------------------------------------------------
// fitting

#include "../fitter/fitter.h"
#include "../fitter/models/interpolation.h"
#include "../fitter/chi2.h"
#include <algorithm>
#include <exception>

#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

/*
static inline std::vector<std::string>
convert_string_vector(const std::vector<t_string>& vecStr)
{
	std::vector<std::string> vecRet;
	vecRet.reserve(vecStr.size());

	for(const t_string& str : vecStr)
		vecRet.push_back(WSTR_TO_STR(str));

	return vecRet;
}

static inline std::vector<t_string>
convert_string_vector(const std::vector<std::string>& vecStr)
{
	std::vector<t_string> vecRet;
	vecRet.reserve(vecStr.size());

	for(const std::string& str : vecStr)
		vecRet.push_back(STR_TO_WSTR(str));

	return vecRet;
}*/

class GenericModel : public FunctionModel
{
protected:
	const NodeFunction *m_pFkt;
	ParseInfo *m_pinfo;
	const SymbolTable *m_pCallerSymTab;

	SymbolTable *m_pTable;

	std::string m_strFreeParam;
	std::vector<std::string> m_vecParamNames;
	std::vector<Symbol*> m_vecSyms;

	bool m_bUseVecParams = 0;

public:
	GenericModel(const GenericModel& mod)
			: m_pinfo(new ParseInfo(*mod.m_pinfo)),
			  m_pCallerSymTab(mod.m_pCallerSymTab),
			  m_pFkt((NodeFunction*)mod.m_pFkt/*->clone()*/),
			  m_pTable(new SymbolTable())
	{
		m_pinfo->bDestroyParseInfo = 0;
		m_strFreeParam = mod.m_strFreeParam;
		m_vecParamNames = mod.m_vecParamNames;
		m_bUseVecParams = mod.m_bUseVecParams;

		if(m_bUseVecParams)
		{
			m_vecSyms.resize(2);
			m_vecSyms[0] = mod.m_vecSyms[0]->clone();
			m_vecSyms[1] = mod.m_vecSyms[1]->clone();
		}
		else
		{
			m_vecSyms.reserve(m_vecParamNames.size()+1);
			for(unsigned int i=0; i<m_vecParamNames.size()+1; ++i)
				m_vecSyms.push_back(new SymbolDouble(mod.m_vecSyms[i]->GetValDouble()));
		}
	}

	GenericModel(const NodeFunction *pFkt, ParseInfo& info, SymbolTable *pCallerSymTab,
			const std::vector<t_string>* pvecParamNames=0)
				: m_pinfo(new ParseInfo(info)),
				  m_pCallerSymTab(pCallerSymTab),
				  m_pFkt((NodeFunction*)pFkt/*->clone()*/),
				  m_pTable(new SymbolTable())
	{
		m_pinfo->bDestroyParseInfo = 0;
		std::vector<std::string> vecParams = /*convert_string_vector*/(m_pFkt->GetParamNames());
		m_strFreeParam = vecParams[0];

		if(pvecParamNames)
		{
			m_bUseVecParams = 1;
			m_vecParamNames = *pvecParamNames;

			SymbolArray* pSymParams = new SymbolArray();
			for(unsigned int i=0; i<m_vecParamNames.size(); ++i)
				pSymParams->GetArr().push_back(new SymbolDouble(0.));

			m_vecSyms.resize(2);
			m_vecSyms[0] = new SymbolDouble(0.);	// free parameter
			m_vecSyms[1] = pSymParams;		// parameter vector
		}
		else
		{
			m_bUseVecParams = 0;
			m_vecParamNames.resize(vecParams.size()-1);
			std::copy(vecParams.begin()+1, vecParams.end(), m_vecParamNames.begin());

			m_vecSyms.reserve(m_vecParamNames.size()+1);
			for(unsigned int i=0; i<m_vecParamNames.size()+1; ++i)
				m_vecSyms.push_back(new SymbolDouble(0.));
		}

		/*G_COUT << "free param: " << m_strFreeParam << std::endl;
		G_COUT << "args: ";
		for(const std::string& strName : m_vecParamNames)
			G_COUT << strName << ", ";
		G_COUT << std::endl;*/
	}

	virtual ~GenericModel()
	{
		//if(m_pFkt) { delete m_pFkt; m_pFkt=0; }
		if(m_pinfo) { delete m_pinfo; m_pinfo=0; }
		if(m_pTable) { delete m_pTable; m_pTable=0; }
	}

	virtual bool SetParams(const std::vector<t_real>& vecParams)
	{
		if(vecParams.size() != m_vecParamNames.size())
		{
			G_CERR << "Error: Parameter array length mismatch."
					<< std::endl;
			return 0;
		}

		for(unsigned int iParam=0; iParam<vecParams.size(); ++iParam)
		{
			//const std::string& strName = m_vecParamNames[iParam];
			t_real dVal = vecParams[iParam];

			if(m_bUseVecParams)
				((SymbolDouble*)((SymbolArray*)m_vecSyms[1])->GetArr()[iParam])->SetVal(dVal);
			else
				((SymbolDouble*)m_vecSyms[iParam+1])->SetVal(dVal);
		}

		return 1;
	}

	virtual t_real operator()(t_real x) const
	{
		((SymbolDouble*)m_vecSyms[0])->SetVal(x);

		SymbolArray arrArgs;
		arrArgs.SetDontDel(1);
		arrArgs.GetArr() = m_vecSyms;
		//m_pFkt->SetArgSyms(&m_vecSyms);

		m_pTable->InsertSymbol(T_STR"<args>", &arrArgs);
		Symbol *pSymRet = m_pFkt->eval(*m_pinfo, m_pTable);
		m_pinfo->bWantReturn = 0;
		m_pTable->RemoveSymbolNoDelete(T_STR"<args>");

		t_real dRetVal = 0.;
		if(pSymRet)
		{
			dRetVal = pSymRet->GetValDouble();
			//if(std::isnan(dRetVal) || std::isinf(dRetVal))
			//	dRetVal = std::numeric_limits<double>::max();

			/*std::cout << "parameters: ";
			for(Symbol* pSym : m_vecSyms)
				std::cout << pSym->GetValDouble() << ", ";

			std::cout << "evaluated: f(" << x << ") = " << dRetVal << std::endl;*/
		}
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

	virtual std::vector<t_real> GetParamValues() const
	{ throw Err("Called invalid function in generic fitter model"); }
	virtual std::vector<t_real> GetParamErrors() const
	{ throw Err("Called invalid function in generic fitter model"); }
};


static void get_values(const std::vector<t_string>& vecParamNames,
						const Symbol* pSym,
						std::vector<t_real>& vec, std::vector<bool>& vecActive)
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
		std::map<t_string, t_real> mymap = sym_to_map<t_string, t_real>(pSym);

		for(unsigned int iParam=0; iParam<vecParamNames.size(); ++iParam)
		{
			const t_string& strKey = vecParamNames[iParam];

			std::map<t_string, t_real>::iterator iter = mymap.find(strKey);
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
	int iDebug = 0;

	if(vecSyms.size()<4 || !is_vec(vecSyms[1]) || !is_vec(vecSyms[2]) || !is_vec(vecSyms[3]))
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(T_STR"Error", info)
			<< "Invalid arguments for fit."
			<< std::endl;
		throw Err(ostrErr.str(),0);
	}

	if(vecSyms[0]->GetType() != SYMBOL_STRING)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(T_STR"Error", info)
			<< "Need a fit function name."
			<< std::endl;
		throw Err(ostrErr.str(),0);
	}

	const t_string& strFkt = ((SymbolString*)vecSyms[0])->GetVal();
	NodeFunction *pFkt = info.GetFunction(strFkt);
	if(!pFkt)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(T_STR"Error", info) << "Invalid function \""
			<< strFkt << "\"."
			<< std::endl;
		throw Err(ostrErr.str(),0);
	}


	std::vector<t_string> vecParamNames;

	bool bUseParamVec = 0;
	bool bHasParamMap = 0;
	if(vecSyms.size()==5 && vecSyms[4]->GetType()==SYMBOL_MAP)
		bHasParamMap = 1;

	if(bHasParamMap)
	{
		SymbolMap::t_map& mapSym = ((SymbolMap*)vecSyms[4])->GetMap();

		SymbolMap::t_map::iterator iterParamNames = mapSym.find(T_STR"use_param_vec");
		if(iterParamNames != mapSym.end())
		{
			vecParamNames = sym_to_vec<t_stdvec, t_string>(iterParamNames->second);
			bUseParamVec = 1;

			//G_COUT << "Using vector parameters" << std::endl;
		}
	}

	GenericModel mod(pFkt, info, pSymTab, bUseParamVec?&vecParamNames:0);

	if(!bUseParamVec)
		vecParamNames = /*convert_string_vector*/(mod.GetParamNames());
	const unsigned int iParamSize = vecParamNames.size();


	std::vector<t_real> vecHints, vecHintsErr;
	std::vector<bool> vecHintsActive, vecHintsErrActive;

	std::vector<t_real> vecLimMin, vecLimMax;
	std::vector<bool> vecLimMinActive, vecLimMaxActive;

	std::vector<t_string> vecFittingSteps;
	std::vector<t_string> vecFixedParams;

	vecLimMinActive.resize(iParamSize);
	vecLimMaxActive.resize(iParamSize);

	// parameter map
	if(bHasParamMap)
	{
		SymbolMap::t_map& mapSym = ((SymbolMap*)vecSyms[4])->GetMap();

		SymbolMap::t_map::iterator iterHints = mapSym.find(T_STR"hints");
		SymbolMap::t_map::iterator iterHintsErr = mapSym.find(T_STR"hints_errors");

		SymbolMap::t_map::iterator iterLimitsMin = mapSym.find(T_STR"lower_limits");
		SymbolMap::t_map::iterator iterLimitsMax = mapSym.find(T_STR"upper_limits");

		SymbolMap::t_map::iterator iterFixed = mapSym.find(T_STR"fixed");

		SymbolMap::t_map::iterator iterSteps = mapSym.find(T_STR"steps");

		SymbolMap::t_map::iterator iterDebug = mapSym.find(T_STR"debug");


		if(iterHints != mapSym.end())
			get_values(vecParamNames, iterHints->second, vecHints, vecHintsActive);
		if(iterHintsErr != mapSym.end())
			get_values(vecParamNames, iterHintsErr->second, vecHintsErr, vecHintsErrActive);

		if(iterLimitsMin != mapSym.end())
			get_values(vecParamNames, iterLimitsMin->second, vecLimMin, vecLimMinActive);
		if(iterLimitsMax != mapSym.end())
			get_values(vecParamNames, iterLimitsMax->second, vecLimMax, vecLimMaxActive);

		if(iterSteps != mapSym.end())
			vecFittingSteps = sym_to_vec<t_stdvec, t_string>(iterSteps->second);

		if(iterFixed != mapSym.end())
			vecFixedParams = sym_to_vec<t_stdvec, t_string>(iterFixed->second);

		if(iterDebug != mapSym.end())
			iDebug = iterDebug->second->GetValInt();

//		for(const t_string& strFixed : vecFixedParams)
//			G_COUT << "fixed params: " << strFixed << std::endl;
	}

	bool bFitterDebug = (iDebug>0);


	std::vector<t_real> vecX = sym_to_vec<t_stdvec>(vecSyms[1]);
	std::vector<t_real> vecY = sym_to_vec<t_stdvec>(vecSyms[2]);
	std::vector<t_real> vecYErr = sym_to_vec<t_stdvec>(vecSyms[3]);

	unsigned int iSize = std::min<unsigned int>(vecX.size(), vecY.size());
	iSize = std::min<unsigned int>(iSize, vecYErr.size());

	Chi2Function chi2fkt(&mod, iSize, vecX.data(), vecY.data(), vecYErr.data());


	ROOT::Minuit2::MnUserParameters params;
	for(unsigned int iParam=0; iParam<iParamSize; ++iParam)
	{
		const t_string& wstrParam = vecParamNames[iParam];
		std::string strParam = WSTR_TO_STR(wstrParam);

		t_real dHint = 0.;
		t_real dErr = 0.;

		if(iParam < vecHints.size())
			dHint = vecHints[iParam];
		if(iParam < vecHintsErr.size())
			dErr = vecHintsErr[iParam];

		//G_COUT << "hints for " << vecParamNames[iParam] << ": "
		//		<< dHint << " +- " << dErr << std::endl;
		params.Add(strParam, dHint, dErr);



		t_real dLimMin = 0.;
		t_real dLimMax = 0.;

		if(iParam < vecLimMin.size())
			dLimMin = vecLimMin[iParam];
		if(iParam < vecLimMax.size())
			dLimMax = vecLimMax[iParam];

/*		if(vecLimMinActive[iParam])
			G_COUT << "lower limit for " << vecParamNames[iParam] << ": "
						<< dLimMin << std::endl;
		if(vecLimMaxActive[iParam])
			G_COUT << "upper limit for " << vecParamNames[iParam] << ": "
						<< dLimMax << std::endl;*/

		if(vecLimMinActive[iParam] && vecLimMaxActive[iParam])
		{
			//std::cout << "Limits for " << strParam << ": " << dLimMin << ", " << dLimMax << std::endl;
			params.SetLimits(strParam, dLimMin, dLimMax);
		}
		else if(vecLimMinActive[iParam] && vecLimMaxActive[iParam]==0)
			params.SetLowerLimit(strParam, dLimMin);
		else if(vecLimMinActive[iParam]==0 && vecLimMaxActive[iParam])
			params.SetUpperLimit(strParam, dLimMax);

		if(std::find(vecFixedParams.begin(), vecFixedParams.end(), /*w*/strParam) != vecFixedParams.end())
			params.Fix(strParam);
	}


	bool bValidFit = 1;
	std::vector<ROOT::Minuit2::FunctionMinimum> minis;

	if(vecFittingSteps.size() == 0)
	{
		if(bFitterDebug)
			G_COUT << "Using default fitting steps." << std::endl;

		minis.reserve(2);

		{
			// step 1: free fit (limited)

			ROOT::Minuit2::MnMigrad migrad(chi2fkt, params, 2);
			ROOT::Minuit2::FunctionMinimum mini = migrad();
			bValidFit = mini.IsValid() && mini.HasValidParameters();

			for(const t_string& strSym : vecParamNames)
			{
				std::string _strSym = WSTR_TO_STR(strSym);

				params.SetValue(_strSym, mini.UserState().Value(_strSym));
				params.SetError(_strSym, mini.UserState().Error(_strSym));
			}

			minis.push_back(mini);
		}


		{
			// step 2: free fit (unlimited)

			for(const t_string& strSym : vecParamNames)
				params.RemoveLimits(WSTR_TO_STR(strSym));

			ROOT::Minuit2::MnMigrad migrad(chi2fkt, params, 2);
			ROOT::Minuit2::FunctionMinimum mini = migrad();
			bValidFit = mini.IsValid() && mini.HasValidParameters();

			minis.push_back(mini);
		}
	}
	else	// custom fitting steps
	{
		if(bFitterDebug)
			G_COUT << "Using " << vecFittingSteps.size()
					<< " custom fitting steps." << std::endl;

		minis.reserve(vecFittingSteps.size());

		// one strip per fitting step
		unsigned int iFitStep = 0;
		for(const t_string& strStep : vecFittingSteps)
		{
			if(bFitterDebug)
				G_COUT << "Performing fitting step " << iFitStep+1
						<< ": " << strStep << std::endl;

			// one character for each parameter
			for(unsigned int iStepParam = 0; iStepParam<strStep.length(); ++iStepParam)
			{
				t_char cOp = strStep[iStepParam];
				const t_string& strParam = vecParamNames[iStepParam];

				switch(cOp)
				{
					case 'f':			// fix param
						params.Fix(WSTR_TO_STR(strParam));
						break;
					case 'r': 			// release param
						params.Release(WSTR_TO_STR(strParam));
						break;
					case 'u': 			// remove limits
						params.RemoveLimits(WSTR_TO_STR(strParam));
						break;
					case 'x': 			// release param and remove limits
						params.Release(WSTR_TO_STR(strParam));
						params.RemoveLimits(WSTR_TO_STR(strParam));
						break;
					case 'n':			// nop
						break;
					default:
						G_CERR << "Error: Unknow fitting step operation \'"
								<< cOp << "\'." << std::endl;
						break;
				}
			} // ops


			ROOT::Minuit2::MnMigrad migrad(chi2fkt, params, 2);
			ROOT::Minuit2::FunctionMinimum mini = migrad();
			bValidFit = mini.IsValid() && mini.HasValidParameters();

			for(unsigned int iStepParam = 0; iStepParam<strStep.length(); ++iStepParam)
			{
				char cOp = strStep[iStepParam];
				const t_string& strParam = vecParamNames[iStepParam];
				std::string _strSym = WSTR_TO_STR(strParam);

				// copy fitted values only for non-fixed params
				if(cOp != 'f')
				{
					params.SetValue(_strSym, mini.UserState().Value(_strSym));
					params.SetError(_strSym, mini.UserState().Error(_strSym));
				}
			}

			minis.push_back(mini);
			++iFitStep;
		} // steps
	}




	const ROOT::Minuit2::FunctionMinimum& lastmini = *minis.rbegin();


	SymbolMap *pSymMap = new SymbolMap();
	for(const t_string& strSym : vecParamNames)
	{
		std::string _strSym = WSTR_TO_STR(strSym);

		t_real dVal = lastmini.UserState().Value(_strSym);
		t_real dErr = lastmini.UserState().Error(_strSym);
		dErr = fabs(dErr);

		SymbolArray* pArr = new SymbolArray();
		pArr->GetArr().push_back(new SymbolDouble(dVal));
		pArr->GetArr().push_back(new SymbolDouble(dErr));

		pSymMap->GetMap().insert(SymbolMap::t_map::value_type(strSym, pArr));
	}


	if(bFitterDebug)
	{
		G_CERR << "--------------------------------------------------------------------------------" << std::endl;
		unsigned int uiMini=0;
		for(const auto& mini : minis)
		{
			std::ostringstream ostrMini;
			ostrMini << mini;

			G_CERR << "result of user-defined fit step " << (++uiMini) << std::endl;
			G_CERR << "=================================" << std::endl;
			G_CERR << STR_TO_WSTR(ostrMini.str()) << std::endl;
		}
		G_CERR << "--------------------------------------------------------------------------------" << std::endl;
	}

	SymbolInt *pSymFitValid = new SymbolInt(bValidFit);
	pSymMap->GetMap()["<valid>"] = pSymFitValid;
	//if(!bValidFit)
	//	G_CERR << "Error: Fit invalid!" << std::endl;
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
		std::ostringstream ostrErr;
		ostrErr << linenr(T_STR"Error", info) << "Function needs x and y vector arguments."
			<< std::endl;
		throw Err(ostrErr.str(),0);
	}

	unsigned int ideg = 3;
	unsigned int N = 128;

	if(vecSyms.size() > 2)
		N = vecSyms[2]->GetValInt();
	if(vecSyms.size() > 3)
		ideg = vecSyms[3]->GetValInt();

	std::vector<t_real> vecX = sym_to_vec<t_stdvec>(vecSyms[0]);
	std::vector<t_real> vecY = sym_to_vec<t_stdvec>(vecSyms[1]);

	unsigned int iSize = std::min(vecX.size(), vecY.size());

	FunctionModel_param* pfkt = 0;
	if(whichfkt == FKT_SPLINE)
		pfkt = new BSpline(iSize, vecX.data(), vecY.data(), ideg);
	else if(whichfkt == FKT_BEZIER)
		pfkt = new Bezier(iSize, vecX.data(), vecY.data());
	else
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(T_STR"Error", info) << "Unknown parametric function selected."
					<< std::endl;
		throw Err(ostrErr.str(),0);
	}

	SymbolArray* pArrX = new SymbolArray();
	SymbolArray* pArrY = new SymbolArray();
	pArrX->GetArr().reserve(N);
	pArrY->GetArr().reserve(N);

	for(unsigned int i=0; i<N; ++i)
	{
		t_real t = t_real(i)/t_real(N-1);
		ublas::vector<t_real> vecSpl = (*pfkt)(t);

		if(vecSpl.size() < 2)
			continue;

		pArrX->GetArr().push_back(new SymbolDouble(vecSpl[0]));
		pArrY->GetArr().push_back(new SymbolDouble(vecSpl[1]));
	}

	delete pfkt;

	SymbolArray *pArr = new SymbolArray();
	pArr->GetArr().push_back(pArrX);
	pArr->GetArr().push_back(pArrY);
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

static Symbol* fkt_find_peaks(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() < 2)
	{
		std::ostringstream ostrErr;
		ostrErr << linenr(T_STR"Error", info) << "find_peaks needs x and y arrays."
					<< std::endl;
		throw Err(ostrErr.str(),0);
	}

	const unsigned int iOrder = 5;

	std::vector<t_real> vecX = sym_to_vec<t_stdvec>(vecSyms[0]);
	std::vector<t_real> vecY = sym_to_vec<t_stdvec>(vecSyms[1]);
	const unsigned int iLen = std::min(vecX.size(), vecY.size());

	std::vector<t_real> vecMaximaX, vecMaximaSize, vecMaximaWidth;
	::find_peaks<t_real>(iLen, vecX.data(), vecY.data(), iOrder,
						vecMaximaX, vecMaximaSize, vecMaximaWidth);

	SymbolArray* pArrX = (SymbolArray*)vec_to_sym<t_stdvec>(vecMaximaX);
	SymbolArray* pArrSizes = (SymbolArray*)vec_to_sym<t_stdvec>(vecMaximaSize);
	SymbolArray* pArrWidths = (SymbolArray*)vec_to_sym<t_stdvec>(vecMaximaWidth);

	SymbolArray *pArr = new SymbolArray();
	pArr->GetArr().push_back(pArrX);
	pArr->GetArr().push_back(pArrSizes);
	pArr->GetArr().push_back(pArrWidths);
	return pArr;
}
// --------------------------------------------------------------------------------


// ["a" : [val0, val1]]  =>  ["a" : val0]
static Symbol* fkt_map_vec_to_val(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1 || vecSyms[0]->GetType()!=SYMBOL_MAP)
	{
		G_CERR << linenr(T_STR"Error", info)
					<< "Need a map of vectors."
					<< std::endl;
		return 0;
	}


	// index parameter
	t_int iIdx = 0;
	if(vecSyms.size()>1)
	{
		iIdx = vecSyms[1]->GetValInt();
		if(iIdx < 0)
		{
			G_CERR << linenr(T_STR"Warning", info)
						<< "Ignoring negative index."
						<< std::endl;
			iIdx = 0;
		}
	}


	SymbolMap* pMapRet = new SymbolMap();

	for(const SymbolMap::t_map::value_type& pair : ((SymbolMap*)vecSyms[0])->GetMap())
	{
		const t_string& strKey = pair.first;
		const Symbol* pSym = pair.second;

		if(pSym->GetType() != SYMBOL_ARRAY)
		{
			// this also ignores the "<valid>" fit variable automatically

			//G_CERR << linenr(T_STR"Warning", info)
			//			<< "Ignoring non-vector variable " 
			//			<< "\"" << strKey << "\""
			//			<< " in map." << std::endl;
			continue;
		}

		const std::vector<Symbol*>& arr = ((SymbolArray*)pSym)->GetArr();
		if(iIdx >= arr.size())
		{
			G_CERR << linenr(T_STR"Warning", info)
						<< "Ignoring invalid index."
						<< std::endl;
			continue;
		}

		const Symbol* pElem = arr[iIdx];

		pMapRet->GetMap().insert(SymbolMap::t_map::value_type(strKey, pElem->clone()));
	}

	return pMapRet;
}

extern void init_ext_fit_calls()
{
	t_mapFkts mapFkts =
	{
		t_mapFkts::value_type(T_STR"fit", fkt_fit),
//		t_mapFkts::value_type(T_STR"fit_sin", fkt_fit_sin),
//		t_mapFkts::value_type(T_STR"fit_gauss", fkt_fit_sin),

		t_mapFkts::value_type(T_STR"map_vec_to_val", fkt_map_vec_to_val),

		t_mapFkts::value_type(T_STR"spline", fkt_spline),
		t_mapFkts::value_type(T_STR"bezier", fkt_bezier),

		t_mapFkts::value_type(T_STR"find_peaks", fkt_find_peaks),
	};

	add_ext_calls(mapFkts);
}
