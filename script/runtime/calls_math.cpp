/*
 * external math functions
 * @author tweber
 * @date 2013
 */

#include "../types.h"
#include "calls_math.h"
#include "../calls.h"
#include "../helper/fourier.h"
#include "../helper/linalg.h"
#include "../helper/rand.h"

static inline Symbol* _fkt_linlogspace(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab, bool bLog)
{
	if(vecSyms.size()<3)
	{
		G_CERR << linenr(T_STR"Error", info)
				<< "Invalid call to linspace(start, end, count)."
				<< std::endl;
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
	G_CERR << "Error: Invalid function selected in math_fkt." << std::endl;
	return tErr;
}

template<MathFkts fkt>
static Symbol* fkt_math_for_every(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() < 1)
	{
		G_CERR << linenr(T_STR"Error", info) << "fkt_math_for_every needs at least one argument" << std::endl;
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

	G_CERR << linenr(T_STR"Error", info) << "No valid arguments given for fkt_math_for_every." << std::endl;
	return 0;
}


template<double (*FKT)(double)>
static Symbol* fkt_math_1arg(const std::vector<Symbol*>& vecSyms,
						ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		G_CERR << linenr(T_STR"Error", info) << "fkt_math_1arg takes exactly one argument." << std::endl;
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
		G_CERR << linenr(T_STR"Error", info) << "fkt_math_2args takes exactly two arguments." << std::endl;
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
		G_CERR << linenr(T_STR"Error", info) << "abs takes exactly one argument." << std::endl;
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


	G_CERR << linenr(T_STR"Error", info) << "abs received unsupported symbol type." << std::endl;
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
		G_CERR << linenr(T_STR"Error", info)
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

template<typename T=double> using t_vec = ublas::vector<T>;
template<typename T=double> using t_mat = ublas::matrix<T>;


static Symbol* fkt_length(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 1)
	{
		G_CERR << linenr(T_STR"Error", info) << "Length needs one argument."
			<< std::endl;
		return 0;
	}

	if(!is_vec(vecSyms[0]))
	{
		G_CERR << linenr(T_STR"Error", info) << "Length needs a vector argument."
			<< std::endl;
		return 0;
	}

	t_vec<double> vec = sym_to_vec<t_vec>(vecSyms[0]);
	double dLen = std::sqrt(ublas::inner_prod(vec, vec));
	return new SymbolDouble(dLen);
}

static Symbol* fkt_cross(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2)
	{
		G_CERR << linenr(T_STR"Error", info) << "Cross product needs two arguments."
			<< std::endl;
		return 0;
	}

	if(!is_vec(vecSyms[0]) || !is_vec(vecSyms[1]))
	{
		G_CERR << linenr(T_STR"Error", info) << "Cross product needs vector arguments."
			<< std::endl;
		return 0;
	}

	t_vec<double> vecLeft = sym_to_vec<t_vec>(vecSyms[0]);
	t_vec<double> vecRight = sym_to_vec<t_vec>(vecSyms[1]);

	if(vecLeft.size()!=3 || vecRight.size()!=3)
	{
		G_CERR << linenr(T_STR"Error", info) << "Cross product needs 3-vectors."
			<< std::endl;
		return 0;
	}

	t_vec<double> vecCross = cross_3(vecLeft, vecRight);
	return vec_to_sym<t_vec>(vecCross);
}

// matrix(rows, cols)
// matrix(dim)
static Symbol* fkt_matrix(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()<1)
	{
		G_CERR << linenr(T_STR"Error", info) << "Need size of matrix." << std::endl;
		return 0;
	}

	unsigned int iRows = vecSyms[0]->GetValInt();
	unsigned int iCols = iRows;

	// cols also given
	if(vecSyms.size() >= 2)
		iCols = vecSyms[1]->GetValInt();

	t_mat<double> mat = ublas::zero_matrix<double>(iRows, iCols);
	return mat_to_sym<t_mat>(mat);
}

static Symbol* fkt_transpose(const std::vector<Symbol*>& vecSyms, 
				ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=1)
	{
		G_CERR << linenr(T_STR"Error", info) << "Transpose needs one argument" << std::endl;
		return 0;
	}

	bool bIsMat = 0;
	t_mat<double> mat = sym_to_mat<t_mat, t_vec>(vecSyms[0], &bIsMat);
	if(!bIsMat)
	{
		G_CERR << linenr(T_STR"Error", info) << "Transpose needs a matrix." << std::endl;
		return 0;
	}

	t_mat<double> mat_trans = ublas::trans(mat);
	return mat_to_sym<t_mat>(mat_trans);
}

static Symbol* fkt_inverse(const std::vector<Symbol*>& vecSyms, 
				ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=1)
	{
		G_CERR << linenr(T_STR"Error", info) << "Inverse needs one argument" << std::endl;
		return 0;
	}

	bool bIsMat = 0;
	t_mat<double> mat = sym_to_mat<t_mat, t_vec>(vecSyms[0], &bIsMat);
	if(!bIsMat)
	{
		G_CERR << linenr(T_STR"Error", info) << "Inverse needs a matrix." << std::endl;
		return 0;
	}

	t_mat<double> mat_inv;
	if(!inverse(mat, mat_inv))
		G_CERR << linenr(T_STR"Warning", info) << "Matrix inversion failed." << std::endl;

	return mat_to_sym<t_mat>(mat_inv);
}

static Symbol* fkt_determinant(const std::vector<Symbol*>& vecSyms,
                                ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=1)
	{
		G_CERR << linenr(T_STR"Error", info) << "Determinant needs one argument" << std::endl;
		return 0;
	}

	bool bIsMat = 0;
	t_mat<double> mat = sym_to_mat<t_mat, t_vec>(vecSyms[0], &bIsMat);
	if(!bIsMat || mat.size1()!=mat.size2())
	{
		G_CERR << linenr(T_STR"Error", info)
			<< "Determinant needs a square matrix."
			<< std::endl;
		return 0;
	}

	double dDet = determinant(mat);
	return new SymbolDouble(dDet);
}

static Symbol* fkt_unitmatrix(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size()!=1 || vecSyms[0]->GetType()!=SYMBOL_INT)
	{
		G_CERR << linenr(T_STR"Error", info) << "Need size of unit matrix." << std::endl;
		return 0;
	}

	int iSize = vecSyms[0]->GetValInt();
	t_mat<double> mat = unit_matrix<double>(iSize);

	return mat_to_sym<t_mat>(mat);
}

static Symbol* fkt_product(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2)
	{
		G_CERR << linenr(T_STR"Error", info) << "Product needs two arguments." << std::endl;
		return 0;
	}

	Symbol* pRet = 0;

	bool bFirstIsVec = is_vec(vecSyms[0]);
	bool bSecondIsVec = is_vec(vecSyms[1]);

	// dot product
	if(bFirstIsVec && bSecondIsVec)
	{
		t_vec<double> vec1 = sym_to_vec<t_vec>(vecSyms[0]);
		t_vec<double> vec2 = sym_to_vec<t_vec>(vecSyms[1]);

		pRet = new SymbolDouble(ublas::inner_prod(vec1, vec2));
	}
	else
	{
		unsigned int iCols1, iRows1, iCols2, iRows2;
		bool bFirstIsMat = is_mat(vecSyms[0], &iCols1, &iRows1);
		bool bSecondIsMat = is_mat(vecSyms[1], &iCols2, &iRows2);

		// matrix product
		if(bFirstIsMat && bSecondIsMat)
		{
			if(iRows1 != iCols2 || iCols1 != iRows2)
			{
				G_CERR << linenr(T_STR"Error", info)
					<< "Row and column counts of matrices do not match."
					<< std::endl;
				return 0;
			}

			t_mat<double> mat1 = sym_to_mat<t_mat, t_vec>(vecSyms[0]);
			t_mat<double> mat2 = sym_to_mat<t_mat, t_vec>(vecSyms[1]);

			t_mat<double> matProd = ublas::prod(mat1, mat2);
			pRet = mat_to_sym<t_mat>(matProd);
		}
		else if(bFirstIsMat && bSecondIsVec)
		{
			t_mat<double> mat = sym_to_mat<t_mat, t_vec>(vecSyms[0]);
			t_vec<double> vec = sym_to_vec<t_vec>(vecSyms[1]);

			t_vec<double> vecProd = ublas::prod(mat, vec);
			pRet = vec_to_sym<t_vec>(vecProd);
		}
		else if(bFirstIsVec && bSecondIsMat)
		{
			t_vec<double> vec = sym_to_vec<t_vec>(vecSyms[0]);
			t_mat<double> mat = sym_to_mat<t_mat, t_vec>(vecSyms[1]);

			t_vec<double> vecProd = ublas::prod(vec, mat);
			pRet = vec_to_sym<t_vec>(vecProd);
		}
	}

	if(!pRet)
		G_CERR << linenr(T_STR"Error", info) << "Invalid call to prod." << std::endl;
	return pRet;
}

// --------------------------------------------------------------------------------





// --------------------------------------------------------------------------------
// rand stuff


static Symbol* fkt_rand01(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{
	double dRand = rand01();
	return new SymbolDouble(dRand);
}

static Symbol* fkt_rand_real(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2)
	{
		G_CERR << linenr(T_STR"Error", info) << "rand_real needs two arguments." << std::endl;
		return 0;
	}

	double dMin = vecSyms[0]->GetValDouble();
	double dMax = vecSyms[1]->GetValDouble();

	double dRand = rand_double(dMin, dMax);
	return new SymbolDouble(dRand);
}

static Symbol* fkt_rand_int(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2)
	{
		G_CERR << linenr(T_STR"Error", info) << "rand_int needs two arguments." << std::endl;
		return 0;
	}

	int iMin = vecSyms[0]->GetValInt();
	int iMax = vecSyms[1]->GetValInt();

	int iRand = rand_int(iMin, iMax);
	return new SymbolInt(iRand);
}

static Symbol* fkt_rand_norm(const std::vector<Symbol*>& vecSyms, ParseInfo& info, SymbolTable* pSymTab)
{
	if(vecSyms.size() != 2)
	{
		G_CERR << linenr(T_STR"Error", info) << "rand_norm needs two arguments." << std::endl;
		return 0;
	}

	double dMu = 0.;
	double dSigma = 1.;

	if(vecSyms.size() >= 1)
		dMu = vecSyms[0]->GetValDouble();
	if(vecSyms.size() >= 2)
		dSigma = vecSyms[1]->GetValDouble();

	double dRand = rand_norm(dMu, dSigma);
	return new SymbolDouble(dRand);
}
// --------------------------------------------------------------------------------



extern void init_ext_math_calls()
{
	init_rand();

	t_mapFkts mapFkts =
	{
		// math stuff
		t_mapFkts::value_type(T_STR"sqrt", fkt_math_1arg< ::sqrt >),
		t_mapFkts::value_type(T_STR"cbrt", fkt_math_1arg< ::cbrt >),
		t_mapFkts::value_type(T_STR"exp", fkt_math_1arg< ::exp >),
		t_mapFkts::value_type(T_STR"exp2", fkt_math_1arg< ::exp2 >),
		t_mapFkts::value_type(T_STR"expm1", fkt_math_1arg< ::expm1 >),
		t_mapFkts::value_type(T_STR"log", fkt_math_1arg< ::log >),
		t_mapFkts::value_type(T_STR"log1p", fkt_math_1arg< ::log1p >),
		t_mapFkts::value_type(T_STR"log10", fkt_math_1arg< ::log10 >),
		t_mapFkts::value_type(T_STR"log2", fkt_math_1arg< ::log2 >),
		t_mapFkts::value_type(T_STR"logb", fkt_math_1arg< ::logb >),
		t_mapFkts::value_type(T_STR"pow", fkt_math_2args< ::pow >),

		t_mapFkts::value_type(T_STR"sin", fkt_math_1arg< ::sin >),
		t_mapFkts::value_type(T_STR"cos", fkt_math_1arg< ::cos >),
		t_mapFkts::value_type(T_STR"tan", fkt_math_1arg< ::tan >),
		t_mapFkts::value_type(T_STR"asin", fkt_math_1arg< ::asin >),
		t_mapFkts::value_type(T_STR"acos", fkt_math_1arg< ::acos >),
		t_mapFkts::value_type(T_STR"atan", fkt_math_1arg< ::atan >),
		t_mapFkts::value_type(T_STR"atan2", fkt_math_2args< ::atan2 >),
		t_mapFkts::value_type(T_STR"hypot", fkt_math_2args< ::hypot >),

		t_mapFkts::value_type(T_STR"sinh", fkt_math_1arg< ::sinh >),
		t_mapFkts::value_type(T_STR"cosh", fkt_math_1arg< ::cosh >),
		t_mapFkts::value_type(T_STR"tanh", fkt_math_1arg< ::tanh >),
		t_mapFkts::value_type(T_STR"asinh", fkt_math_1arg< ::asinh >),
		t_mapFkts::value_type(T_STR"acosh", fkt_math_1arg< ::acosh >),
		t_mapFkts::value_type(T_STR"atanh", fkt_math_1arg< ::atanh >),

		t_mapFkts::value_type(T_STR"erf", fkt_math_1arg< ::erf >),
		t_mapFkts::value_type(T_STR"erfc", fkt_math_1arg< ::erfc >),
		t_mapFkts::value_type(T_STR"tgamma", fkt_math_1arg< ::tgamma >),
		t_mapFkts::value_type(T_STR"lgamma", fkt_math_1arg< ::lgamma >),

		t_mapFkts::value_type(T_STR"round", fkt_math_1arg< ::round >),
		t_mapFkts::value_type(T_STR"trunc", fkt_math_1arg< ::trunc >),
		t_mapFkts::value_type(T_STR"rint", fkt_math_1arg< ::rint >),
		t_mapFkts::value_type(T_STR"nearbyint", fkt_math_1arg< ::nearbyint >),
		t_mapFkts::value_type(T_STR"fmod", fkt_math_2args< ::fmod >),
		t_mapFkts::value_type(T_STR"nextafter", fkt_math_2args< ::nextafter >),
		//t_mapFkts::value_type(T_STR"nexttoward", fkt_math_2args< ::nexttoward >),
		t_mapFkts::value_type(T_STR"ceil", fkt_math_1arg< ::ceil >),
		t_mapFkts::value_type(T_STR"floor", fkt_math_1arg< ::floor >),
		t_mapFkts::value_type(T_STR"abs", fkt_math_abs),
		t_mapFkts::value_type(T_STR"max", fkt_math_for_every<MATH_MAX>),
		t_mapFkts::value_type(T_STR"min", fkt_math_for_every<MATH_MIN>),
		t_mapFkts::value_type(T_STR"fdim", fkt_math_2args< ::fdim >),
		t_mapFkts::value_type(T_STR"remainder", fkt_math_2args< ::remainder >),

		// fft
		t_mapFkts::value_type(T_STR"fft", fkt_fft),
		t_mapFkts::value_type(T_STR"ifft", fkt_ifft),

		// arrays
		t_mapFkts::value_type(T_STR"linspace", fkt_linspace),
		t_mapFkts::value_type(T_STR"logspace", fkt_logspace),

		// vector operations
		//t_mapFkts::value_type(T_STR"dot", fkt_dot), -> use prod
		t_mapFkts::value_type(T_STR"cross", fkt_cross),
		t_mapFkts::value_type(T_STR"len", fkt_length),

		// matrix operations
		t_mapFkts::value_type(T_STR"mat", fkt_matrix),
		t_mapFkts::value_type(T_STR"unity", fkt_unitmatrix),
		t_mapFkts::value_type(T_STR"trans", fkt_transpose),
		t_mapFkts::value_type(T_STR"inv", fkt_inverse),
		t_mapFkts::value_type(T_STR"det", fkt_determinant),

		// matrix-vector operations
		t_mapFkts::value_type(T_STR"prod", fkt_product),


		// random numbers
		t_mapFkts::value_type(T_STR"rand01", fkt_rand01),
		t_mapFkts::value_type(T_STR"rand_real", fkt_rand_real),
		t_mapFkts::value_type(T_STR"rand_int", fkt_rand_int),
		t_mapFkts::value_type(T_STR"rand_norm", fkt_rand_norm),
	};

	add_ext_calls(mapFkts);
}
