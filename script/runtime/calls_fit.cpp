/*
 * external fit functions
 * @author tweber
 * @date jan 2014
 */

#include "calls_fit.h"
#include "calls_math.h"
#include "../calls.h"

#include "../fitter/models/interpolation.h"
namespace ublas = boost::numeric::ublas;

// --------------------------------------------------------------------------------
// fitting

static Symbol* fkt_fit(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	// TODO
	return 0;
}
// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// interpolation

template<typename T=double> using t_stdvec = std::vector<T>;

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


extern void init_ext_fit_calls()
{
	t_mapFkts mapFkts =
	{
		t_mapFkts::value_type("fit", fkt_fit),
//		t_mapFkts::value_type("fit_sin", fkt_fit_sin),
//		t_mapFkts::value_type("fit_gauss", fkt_fit_sin),

		t_mapFkts::value_type("spline", fkt_spline),
		t_mapFkts::value_type("bezier", fkt_bezier),
	};

	add_ext_calls(mapFkts);
}
