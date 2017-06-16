/**
 * julia module
 * gcc -I. -I.. -I/usr/include/julia -I/usr/include/root -std=c++11 -shared -fPIC -o tlibs_jl.so ../log/log.cpp tlibs_jl.cpp -lstdc++ -lboost_system -lboost_iostreams -L/usr/lib64/root -lMinuit2
 * gcc -I. -I.. -I/usr/include/julia -I/usr/local/include -std=c++11 -shared -fPIC -o tlibs_jl.so ../log/log.cpp tlibs_jl.cpp -lstdc++ -lboost_system -lboost_iostreams -L/usr/local/lib64 -lMinuit2 -lgomp
 *
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 23-apr-2017
 * @license GPLv2 or GPLv3
 */

#define TLIBS_INC_HDR_IMPLS
#include "version.h"
#include "ext/jl.h"
#include "log/log.h"
#include "file/loadinstr.h"
#include "fit/minuit.h"


using t_real = double;


extern "C" void load_tlibs()
{
	tl::log_debug("Loaded tlibs version ", TLIBS_VERSION, ".");
}


/**
 * loads an instrument data file
 */
extern "C" jl_array_t* load_instr(const char* pcFile)
{
	// [ column names, data, keys, values ]
	jl_array_t *pArr = jl_alloc_array_1d(jl_apply_array_type(jl_any_type, 1), 4);
	jl_array_t** pArrDat = reinterpret_cast<jl_array_t**>(jl_array_data(pArr));

	tl::FileInstrBase<t_real>* pInstr = tl::FileInstrBase<t_real>::LoadInstr(pcFile);
	if(!pInstr)
		tl::log_err("In ", __func__, ": Cannot load ", pcFile, ".");

	// data column names
	pArrDat[0] = tl::make_jl_str_arr(pInstr->GetColNames());

	// data matrix
	pArrDat[1] = tl::make_jl_2darr(pInstr->GetData());

	// scan property map
	std::tie(pArrDat[2], pArrDat[3]) = tl::make_jl_strmap_arr(pInstr->GetAllParams());

	return pArr;
}


/**
 * function fitting
 */
extern "C" int fit(void *_pFkt, std::size_t iNumParams,
	const t_real *pX, const t_real *pY, const t_real *pYerr, std::size_t iArrLen,
	jl_array_t *parrParamNames, jl_array_t *parrFixed,
	const t_real *pValues, const t_real *pErrors)
{
	std::vector<tl::t_real_min> vecX, vecY, vecYerr;

	vecX.reserve(iArrLen);
	vecY.reserve(iArrLen);
	vecYerr.reserve(iArrLen);

	// copy arrays to vectors
	for(std::size_t i=0; i<iArrLen; ++i)
	{
		vecX.push_back(pX[i]);
		vecY.push_back(pY[i]);
		vecYerr.push_back(pYerr[i]);
	}


	// parameter names
	std::vector<std::string> vecParamNames =
		tl::from_jl_arr<std::vector, std::string>(parrParamNames, 1);

	// fixed parameters
	std::vector<std::string> vecFixedParams =
		tl::from_jl_arr<std::vector, std::string>(parrFixed);

	std::vector<bool> vecFixed;
	for(const std::string& strParam : vecParamNames)
	{
		bool bFixed = std::find(vecFixedParams.begin(), vecFixedParams.end(), strParam)
			!= vecFixedParams.end();
		vecFixed.push_back(bFixed);
	}

	// values & errors
	std::vector<tl::t_real_min> vecVals, vecErrs;
	for(std::size_t iIdx=0; iIdx<iNumParams; ++iIdx)
	{
		vecVals.push_back(pValues[iIdx]);
		vecErrs.push_back(pErrors[iIdx]);
	}


	// ------------------------------------------------------------------------
	// fill up missing parameters and hints
	if(vecParamNames.size() < iNumParams)
	{
		for(std::size_t iArg=vecParamNames.size(); iArg<iNumParams; ++iArg)
		{
			std::ostringstream ostrArg;
			ostrArg << "arg_" << iArg;
			vecParamNames.push_back(ostrArg.str());
		}
	}

	while(vecVals.size() < iNumParams) vecVals.push_back(0);
	while(vecErrs.size() < iNumParams) vecErrs.push_back(0);
	while(vecFixed.size() < iNumParams) vecFixed.push_back(0);
	// ------------------------------------------------------------------------


	bool bDebug = 1;
	bool bOk = 0;

	if(iNumParams == 1)
	{
		auto *pFkt = reinterpret_cast<t_real (*)(t_real, t_real)>(_pFkt);
		bOk = tl::fit<2>([pFkt](t_real x, t_real arg1) -> t_real { return pFkt(x, arg1); },
			vecX, vecY, vecYerr, vecParamNames, vecVals, vecErrs, &vecFixed, bDebug);
	}
	else if(iNumParams == 2)
	{
		auto *pFkt = reinterpret_cast<t_real (*)(t_real, t_real, t_real)>(_pFkt);
		bOk = tl::fit<3>([pFkt](t_real x, t_real arg1, t_real arg2) -> t_real { return pFkt(x, arg1, arg2); },
			vecX, vecY, vecYerr, vecParamNames, vecVals, vecErrs, &vecFixed, bDebug);
	}
	else if(iNumParams == 3)
	{
		auto *pFkt = reinterpret_cast<t_real (*)(t_real, t_real, t_real, t_real)>(_pFkt);
		bOk = tl::fit<4>([pFkt](t_real x, t_real arg1, t_real arg2, t_real arg3) -> t_real { return pFkt(x, arg1, arg2, arg3); },
			vecX, vecY, vecYerr, vecParamNames, vecVals, vecErrs, &vecFixed, bDebug);
	}
	else if(iNumParams == 4)
	{
		auto *pFkt = reinterpret_cast<t_real (*)(t_real, t_real, t_real, t_real, t_real)>(_pFkt);
		bOk = tl::fit<5>([pFkt](t_real x, t_real arg1, t_real arg2, t_real arg3, t_real arg4) -> t_real { return pFkt(x, arg1, arg2, arg3, arg4); },
			vecX, vecY, vecYerr, vecParamNames, vecVals, vecErrs, &vecFixed, bDebug);
	}
	else if(iNumParams == 5)
	{
		auto *pFkt = reinterpret_cast<t_real (*)(t_real, t_real, t_real, t_real, t_real, t_real)>(_pFkt);
		bOk = tl::fit<6>([pFkt](t_real x, t_real arg1, t_real arg2, t_real arg3, t_real arg4, t_real arg5) -> t_real { return pFkt(x, arg1, arg2, arg3, arg4, arg5); },
			vecX, vecY, vecYerr, vecParamNames, vecVals, vecErrs, &vecFixed, bDebug);
	}
	else if(iNumParams == 6)
	{
		auto *pFkt = reinterpret_cast<t_real (*)(t_real, t_real, t_real, t_real, t_real, t_real, t_real)>(_pFkt);
		bOk = tl::fit<7>([pFkt](t_real x, t_real arg1, t_real arg2, t_real arg3, t_real arg4, t_real arg5, t_real arg6) -> t_real { return pFkt(x, arg1, arg2, arg3, arg4, arg5, arg6); },
			vecX, vecY, vecYerr, vecParamNames, vecVals, vecErrs, &vecFixed, bDebug);
	}
	else if(iNumParams == 7)
	{
		auto *pFkt = reinterpret_cast<t_real (*)(t_real, t_real, t_real, t_real, t_real, t_real, t_real, t_real)>(_pFkt);
		bOk = tl::fit<8>([pFkt](t_real x, t_real arg1, t_real arg2, t_real arg3, t_real arg4, t_real arg5, t_real arg6, t_real arg7) -> t_real { return pFkt(x, arg1, arg2, arg3, arg4, arg5, arg6, arg7); },
			vecX, vecY, vecYerr, vecParamNames, vecVals, vecErrs, &vecFixed, bDebug);
	}
	else if(iNumParams == 8)
	{
		auto *pFkt = reinterpret_cast<t_real (*)(t_real, t_real, t_real, t_real, t_real, t_real, t_real, t_real, t_real)>(_pFkt);
		bOk = tl::fit<9>([pFkt](t_real x, t_real arg1, t_real arg2, t_real arg3, t_real arg4, t_real arg5, t_real arg6, t_real arg7, t_real arg8) -> t_real { return pFkt(x, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8); },
			vecX, vecY, vecYerr, vecParamNames, vecVals, vecErrs, &vecFixed, bDebug);
	}

	return int(bOk);
}
