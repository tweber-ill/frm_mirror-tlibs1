/*
 * external math functions
 * @author tweber
 * @date 2013
 */

#include "calls_math.h"
#include "calls.h"
#include "helper/fourier.h"

static inline Symbol* _fkt_linlogspace(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab, bool bLog)
{
	if(vecSyms.size()<3)
	{
		std::cerr << linenr("Error", info) << "Invalid call to linspace(start, end, count)." << std::endl;
		return 0;
	}

	int iNum = vecSyms[2]->GetValInt();
	double dmin = vecSyms[0]->GetValDouble();
	double dmax = vecSyms[1]->GetValDouble();

	SymbolArray* pSymRet = new SymbolArray;
	pSymRet->m_arr.reserve(iNum);
	for(int i=0; i<iNum; ++i)
	{
		SymbolDouble *pSymD = new SymbolDouble();
		double dDiv = (iNum!=1 ? double(iNum-1) : 1);
		pSymD->m_dVal = double(i)*(dmax-dmin)/dDiv + dmin;
		if(bLog)
		{
			const double dBase = 10.;
			pSymD->m_dVal = std::pow(dBase, pSymD->m_dVal);
		}
		pSymRet->m_arr.push_back(pSymD);
	}

	pSymRet->UpdateIndices();

	return pSymRet;
}

static Symbol* fkt_linspace(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	return _fkt_linlogspace(vecSyms, info, pSymTab, false);
}

static Symbol* fkt_logspace(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	return _fkt_linlogspace(vecSyms, info, pSymTab, true);
}



// --------------------------------------------------------------------------------
// math
enum MathFkts
{
	MATH_MAX,
	MATH_MIN
};

template<MathFkts fkt, typename T>
const T& math_fkt(const T& t1, const T& t2)
{
	if(fkt == MATH_MAX)
		return std::max<T>(t1, t2);
	else if(fkt == MATH_MIN)
		return std::min<T>(t1, t2);

	static const T tErr = T(0);
	std::cerr << "Error: Invalid function selected in math_fkt." << std::endl;
	return tErr;
}

template<MathFkts fkt>
static Symbol* fkt_math_for_every(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() < 1)
	{
		std::cerr << linenr("Error", info) << "fkt_math_for_every needs at least one argument" << std::endl;
		return 0;
	}

	double dRes;
	int iRes;

	bool bHadInt = 0,
		bHadDouble = 0;

	for(Symbol* pSym : vecSyms)
	{
		Symbol *pThisSym = pSym;
		bool bCleanSym = 0;
		if(pSym->GetType() == SYMBOL_ARRAY)
		{
			pThisSym = fkt_math_for_every<fkt>(
					((SymbolArray*)pSym)->m_arr,
					info, pSymTab);

			bCleanSym = 1;
		}

		if(pThisSym->GetType() == SYMBOL_INT)
		{
			if(!bHadInt)
				iRes = ((SymbolInt*)pThisSym)->m_iVal;
			else
				iRes = math_fkt<fkt, int>(iRes, ((SymbolInt*)pThisSym)->m_iVal);

			bHadInt = 1;
		}
		else if(pThisSym->GetType() == SYMBOL_DOUBLE)
		{
			if(!bHadDouble)
				dRes = ((SymbolDouble*)pThisSym)->m_dVal;
			else
				dRes = math_fkt<fkt, double>(dRes, ((SymbolDouble*)pThisSym)->m_dVal);

			bHadDouble = 1;
		}

		if(bCleanSym)
			delete pThisSym;
	}

	if(bHadInt && !bHadDouble)
		return new SymbolInt(iRes);
	else if(bHadInt && bHadDouble)
	{
		dRes = math_fkt<fkt, double>(dRes, double(iRes));
		return new SymbolDouble(dRes);
	}
	else if(!bHadInt && bHadDouble)
		return new SymbolDouble(dRes);

	std::cerr << linenr("Error", info) << "No valid arguments given for fkt_math_for_every." << std::endl;
	return 0;
}


template<double (*FKT)(double)>
static Symbol* fkt_math_1arg(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		std::cerr << linenr("Error", info) << "fkt_math_1arg takes exactly one argument." << std::endl;
		return 0;
	}

	if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray* pArrRet = new SymbolArray();

		SymbolArray* pSymArr = (SymbolArray*)vecSyms[0];
		for(Symbol* pArrElem : pSymArr->m_arr)
		{
			std::vector<Symbol*> vecDummy;
			vecDummy.push_back(pArrElem);

			pArrRet->m_arr.push_back(fkt_math_1arg<FKT>(vecDummy, info, pSymTab));
		}

		pArrRet->UpdateIndices();
		return pArrRet;
	}
	else
	{
		double dResult = FKT(vecSyms[0]->GetValDouble());
		return new SymbolDouble(dResult);
	}

	return 0;
}

// TODO: for arrays
template<double (*FKT)(double, double)>
static Symbol* fkt_math_2args(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2)
	{
		std::cerr << linenr("Error", info) << "fkt_math_2args takes exactly two arguments." << std::endl;
		return 0;
	}

	double dResult = FKT(vecSyms[0]->GetValDouble(), vecSyms[1]->GetValDouble());
	return new SymbolDouble(dResult);
}


template<typename T> T myabs(T t)
{
	if(t < T(0))
		return -t;
	return t;
}

template<> double myabs(double t) { return ::fabs(t); }


// TODO: integrate with fkt_math_1arg
static Symbol* fkt_math_abs(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		std::cerr << linenr("Error", info) << "abs takes exactly one argument." << std::endl;
		return 0;
	}

	if(vecSyms[0]->GetType() == SYMBOL_ARRAY)
	{
		SymbolArray* pArrRet = new SymbolArray();

		SymbolArray* pSymArr = (SymbolArray*)vecSyms[0];
		for(Symbol* pArrElem : pSymArr->m_arr)
		{
			std::vector<Symbol*> vecDummy;
			vecDummy.push_back(pArrElem);

			pArrRet->m_arr.push_back(fkt_math_abs(vecDummy, info, pSymTab));
		}

		pArrRet->UpdateIndices();
		return pArrRet;
	}
	else if(vecSyms[0]->GetType() == SYMBOL_INT)
	{
		SymbolInt* pSymInt = (SymbolInt*)vecSyms[0];
		return new SymbolInt(myabs(pSymInt->m_iVal));
	}
	else if(vecSyms[0]->GetType() == SYMBOL_DOUBLE)
	{
		SymbolDouble* pSymD = (SymbolDouble*)vecSyms[0];
		return new SymbolDouble(myabs(pSymD->m_dVal));
	}


	std::cerr << linenr("Error", info) << "abs received unsupported symbol type." << std::endl;
	return 0;
}


// --------------------------------------------------------------------------------
// FFT

static Symbol* _fkt_fft(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab,
						bool bInv)
{
	bool (Fourier::*pFkt)(const double*, const double*, double*, double*)
								= (bInv ? &Fourier::ifft : &Fourier::fft);

	bool bArgsOk=1;
	std::vector<double> vecRealIn, vecImagIn;

	// real and imag part as separate arguments
	if(vecSyms.size()==2 && vecSyms[0]->GetType()==SYMBOL_ARRAY
			&& vecSyms[1]->GetType()==SYMBOL_ARRAY)
	{
		vecRealIn = ((SymbolArray*)vecSyms[0])->ToDoubleArray();
		vecImagIn = ((SymbolArray*)vecSyms[1])->ToDoubleArray();
	}
	// arrays in one array
	else if(vecSyms.size()==1 && vecSyms[0]->GetType()==SYMBOL_ARRAY)
	{
		SymbolArray* pSymArr = (SymbolArray*)vecSyms[0];
		unsigned int iSymArrSize = pSymArr->m_arr.size();

		if(iSymArrSize==0)
			bArgsOk = 0;
		if(iSymArrSize>=1 && pSymArr->m_arr[0]->GetType()==SYMBOL_ARRAY)
			vecRealIn = ((SymbolArray*)pSymArr->m_arr[0])->ToDoubleArray();
		if(iSymArrSize>=2 && pSymArr->m_arr[1]->GetType()==SYMBOL_ARRAY)
			vecImagIn = ((SymbolArray*)pSymArr->m_arr[1])->ToDoubleArray();

		// simple array containing real data
		if(pSymArr->m_arr[0]->GetType()!=SYMBOL_ARRAY)
			vecRealIn = pSymArr->ToDoubleArray();
	}

	if(!bArgsOk)
	{
		std::cerr << linenr("Error", info)
				<< "fft received invalid arguments."
				<< std::endl;
		return 0;
	}

	if(vecRealIn.size() != vecImagIn.size())
	{
		unsigned int iSize = std::max(vecRealIn.size(), vecImagIn.size());
		vecRealIn.resize(iSize);
		vecImagIn.resize(iSize);
	}

	std::vector<double> vecRealOut, vecImagOut;
	vecRealOut.resize(vecRealIn.size());
	vecImagOut.resize(vecImagIn.size());

	Fourier fourier(vecRealIn.size());
	(fourier.*pFkt)(vecRealIn.data(), vecImagIn.data(),
			vecRealOut.data(), vecImagOut.data());

	SymbolArray* pArrReal = new SymbolArray();
	SymbolArray* pArrImag = new SymbolArray();

	pArrReal->FromDoubleArray(vecRealOut);
	pArrImag->FromDoubleArray(vecImagOut);

	SymbolArray* pRet = new SymbolArray();
	pRet->m_arr.push_back(pArrReal);
	pRet->m_arr.push_back(pArrImag);

	return pRet;
}

static Symbol* fkt_fft(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{ return _fkt_fft(vecSyms, info, pSymTab, false); }
static Symbol* fkt_ifft(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{ return _fkt_fft(vecSyms, info, pSymTab, true); }

// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// linalg stuff

static Symbol* fkt_dot(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2)
	{
		std::cerr << linenr("Error", info) << "dot needs two arguments." << std::endl;
		return 0;
	}

	if(vecSyms[0]->GetType()!=SYMBOL_ARRAY || vecSyms[1]->GetType()!=SYMBOL_ARRAY)
	{
		std::cerr << linenr("Error", info) << "dot needs vector arguments." << std::endl;
		return 0;
	}

	SymbolArray* pLeft = (SymbolArray*)vecSyms[0];
	SymbolArray* pRight = (SymbolArray*)vecSyms[1];

	double dRetVal = 0.;
	unsigned int iSize = std::min(pLeft->m_arr.size(), pRight->m_arr.size());
	for(unsigned int i=0; i<iSize; ++i)
	{
		dRetVal += pLeft->m_arr[i]->GetValDouble()*pRight->m_arr[i]->GetValDouble();
	}

	return new SymbolDouble(dRetVal);
}

// --------------------------------------------------------------------------------



extern void init_ext_math_calls()
{
	t_mapFkts mapFkts =
	{
		// math stuff
		t_mapFkts::value_type("sqrt", fkt_math_1arg< ::sqrt >),
		t_mapFkts::value_type("cbrt", fkt_math_1arg< ::cbrt >),
		t_mapFkts::value_type("exp", fkt_math_1arg< ::exp >),
		t_mapFkts::value_type("exp2", fkt_math_1arg< ::exp2 >),
		t_mapFkts::value_type("expm1", fkt_math_1arg< ::expm1 >),
		t_mapFkts::value_type("log", fkt_math_1arg< ::log >),
		t_mapFkts::value_type("log1p", fkt_math_1arg< ::log1p >),
		t_mapFkts::value_type("log10", fkt_math_1arg< ::log10 >),
		t_mapFkts::value_type("log2", fkt_math_1arg< ::log2 >),
		t_mapFkts::value_type("logb", fkt_math_1arg< ::logb >),
		t_mapFkts::value_type("pow", fkt_math_2args< ::pow >),

		t_mapFkts::value_type("sin", fkt_math_1arg< ::sin >),
		t_mapFkts::value_type("cos", fkt_math_1arg< ::cos >),
		t_mapFkts::value_type("tan", fkt_math_1arg< ::tan >),
		t_mapFkts::value_type("asin", fkt_math_1arg< ::asin >),
		t_mapFkts::value_type("acos", fkt_math_1arg< ::acos >),
		t_mapFkts::value_type("atan", fkt_math_1arg< ::atan >),
		t_mapFkts::value_type("atan2", fkt_math_2args< ::atan2 >),
		t_mapFkts::value_type("hypot", fkt_math_2args< ::hypot >),

		t_mapFkts::value_type("sinh", fkt_math_1arg< ::sinh >),
		t_mapFkts::value_type("cosh", fkt_math_1arg< ::cosh >),
		t_mapFkts::value_type("tanh", fkt_math_1arg< ::tanh >),
		t_mapFkts::value_type("asinh", fkt_math_1arg< ::asinh >),
		t_mapFkts::value_type("acosh", fkt_math_1arg< ::acosh >),
		t_mapFkts::value_type("atanh", fkt_math_1arg< ::atanh >),

		t_mapFkts::value_type("erf", fkt_math_1arg< ::erf >),
		t_mapFkts::value_type("erfc", fkt_math_1arg< ::erfc >),
		t_mapFkts::value_type("tgamma", fkt_math_1arg< ::tgamma >),
		t_mapFkts::value_type("lgamma", fkt_math_1arg< ::lgamma >),

		t_mapFkts::value_type("round", fkt_math_1arg< ::round >),
		t_mapFkts::value_type("trunc", fkt_math_1arg< ::trunc >),
		t_mapFkts::value_type("rint", fkt_math_1arg< ::rint >),
		t_mapFkts::value_type("nearbyint", fkt_math_1arg< ::nearbyint >),
		t_mapFkts::value_type("fmod", fkt_math_2args< ::fmod >),
		t_mapFkts::value_type("nextafter", fkt_math_2args< ::nextafter >),
		//t_mapFkts::value_type("nexttoward", fkt_math_2args< ::nexttoward >),
		t_mapFkts::value_type("ceil", fkt_math_1arg< ::ceil >),
		t_mapFkts::value_type("floor", fkt_math_1arg< ::floor >),
		t_mapFkts::value_type("abs", fkt_math_abs),
		t_mapFkts::value_type("max", fkt_math_for_every<MATH_MAX>),
		t_mapFkts::value_type("min", fkt_math_for_every<MATH_MIN>),
		t_mapFkts::value_type("fdim", fkt_math_2args< ::fdim >),
		t_mapFkts::value_type("remainder", fkt_math_2args< ::remainder >),

		// fft
		t_mapFkts::value_type("fft", fkt_fft),
		t_mapFkts::value_type("ifft", fkt_ifft),

		// arrays
		t_mapFkts::value_type("linspace", fkt_linspace),
		t_mapFkts::value_type("logspace", fkt_logspace),

		// vector / matrix operations
		t_mapFkts::value_type("dot", fkt_dot),
	};

	add_ext_calls(mapFkts);
}
