/**
 * advanced linalg helpers
 *
 * @author: Tobias Weber
 * @date: 2013-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_LINALG2_IMPL_H__
#define __TLIBS_LINALG2_IMPL_H__

#include "linalg2.h"
#include <memory>

#ifndef NO_LAPACK
extern "C"
{
	#include <lapacke.h>
}

namespace tl {

// selects the float or double version of a lapack function
template<class T1, class T2, class F1, class F2>
struct select_func
{
	F1* m_f1 = nullptr;
	F2* m_f2 = nullptr;

	select_func(F1* f1, F2* f2) : m_f1(f1), m_f2(f2)
	{}

	template<class T>
	typename std::enable_if<std::is_same<T, T1>::value, F1*>::type 
		get_func() { return m_f1; }
	template<class T>
	typename std::enable_if<std::is_same<T, T2>::value, F2*>::type 
		get_func() { return m_f2; }
};


template<class T>
bool eigenvec(const ublas::matrix<T>& mat,
	std::vector<ublas::vector<T> >& evecs_real,
	std::vector<ublas::vector<T> >& evecs_imag,
	std::vector<T>& evals_real,
	std::vector<T>& evals_imag)
{
	select_func<float, double, decltype(LAPACKE_sgeev), decltype(LAPACKE_dgeev)> 
		sfunc(LAPACKE_sgeev, LAPACKE_dgeev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evecs_real.resize(iOrder);
	evecs_imag.resize(iOrder);
	evals_real.resize(iOrder);
	evals_imag.resize(iOrder);
	for(std::size_t i=0; i<iOrder; ++i)
	{
		evecs_real[i].resize(iOrder);
		evecs_imag[i].resize(iOrder);
	}

	/*
	// is matrix already diagonal?
	if(is_diag_matrix(mat))
	{
		for(std::size_t i=0; i<iOrder; ++i)
		{
			evals_real[i] = mat(i,i);

			evecs_real[i] = ublas::zero_vector<T>(iOrder);
			evecs_real[i][i] = 1.;


			evals_imag[i] = 0.;
			evecs_imag[i] = ublas::zero_vector<T>(iOrder);
		}

		return true;
	}
	*/

	bool bOk = true;

	std::unique_ptr<T, std::default_delete<T[]>> 
		uptrMem(new T[iOrder*iOrder + iOrder*iOrder + iOrder*iOrder + iOrder + iOrder]);
	T *pMem = uptrMem.get();

	T *pMatrix = pMem;
	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	T *p_eigenvecs = pMatrix + iOrder*iOrder;
	//T *p_eigenvecs_l = p_eigenvecs + iOrder*iOrder;
	T *peigenvals_real = p_eigenvecs/*_l*/ + iOrder*iOrder;
	T *peigenvals_imag = peigenvals_real + iOrder;

	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'N', 'V', iOrder, (T*)pMatrix, iOrder,
		(T*)peigenvals_real, (T*)peigenvals_imag,
		/*(T*)p_eigenvecs_l*/0, iOrder, (T*)p_eigenvecs, iOrder);

	if(iInfo!=0)
	{
		log_err("Could not solve eigenproblem", " (lapack error ", iInfo , ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iOrder; ++i)
	{
		bool bIsReal = 0;
		if(float_equal<T>(peigenvals_imag[i], 0.))
			bIsReal = 1;

		if(bIsReal)
		{
			for(std::size_t j=0; j<iOrder; ++j)
			{
				evecs_real[i][j] = p_eigenvecs[j*iOrder + i];
				evecs_imag[i][j] = 0.;

				evals_real[i] = peigenvals_real[i];
				evals_imag[i] = peigenvals_imag[i];
			}
		}
		else
		{
			for(std::size_t j=0; j<iOrder; ++j)
			{
				evecs_real[i][j] = p_eigenvecs[j*iOrder + i];
				evecs_imag[i][j] = p_eigenvecs[j*iOrder + i+1];

				evecs_real[i+1][j] = p_eigenvecs[j*iOrder + i];
				evecs_imag[i+1][j] = -p_eigenvecs[j*iOrder + i+1];
			}

			evals_real[i] = peigenvals_real[i];
			evals_imag[i] = peigenvals_imag[i];
			evals_real[i+1] = peigenvals_real[i+1];
			evals_imag[i+1] = peigenvals_imag[i+1];

			++i; // check: (next eigenval) == -(currrent eigenval)
		}

		//evecs_real[i] /= ublas::norm_2(evecs_real[i]);
		//evecs_imag[i] /= ublas::norm_2(evecs_imag[i]);
	}

	return bOk;
}

template<class T>
bool eigenvec_sym(const ublas::matrix<T>& mat,
	std::vector<ublas::vector<T>>& evecs,
	std::vector<T>& evals)
{
	select_func<float, double, decltype(LAPACKE_ssyev), decltype(LAPACKE_dsyev)> 
		sfunc(LAPACKE_ssyev, LAPACKE_dsyev);
	auto pfunc = sfunc.get_func<T>();

	if(mat.size1() != mat.size2())
		return false;
	if(mat.size1()==0 || mat.size1()==1)
		return false;

	const std::size_t iOrder = mat.size1();
	evecs.resize(iOrder);
	evals.resize(iOrder);
	for(std::size_t i=0; i<iOrder; ++i)
		evecs[i].resize(iOrder);


	bool bOk = true;

	std::unique_ptr<T, std::default_delete<T[]>> 
		uptrMem(new T[iOrder*iOrder + iOrder]);
	T *pMem = uptrMem.get();
	T *pMatrix = pMem;

	for(std::size_t i=0; i<iOrder; ++i)
		for(std::size_t j=0; j<iOrder; ++j)
			pMatrix[i*iOrder + j] = mat(i,j);

	T *peigenvals_real = pMatrix + iOrder*iOrder;
	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, 'V', 'U',
		iOrder, (T*)pMatrix, iOrder, (T*)peigenvals_real);

	if(iInfo!=0)
	{
		log_err("Could not solve eigenproblem", " (lapack error ", iInfo, ").");
		bOk = false;
	}

	for(std::size_t i=0; i<iOrder; ++i)
	{
		for(std::size_t j=0; j<iOrder; ++j)
			evecs[i][j] = pMatrix[j*iOrder + i];
		//evecs[i] /= ublas::norm_2(evecs[i]);
		evals[i] = peigenvals_real[i];
	}

	if(determinant<ublas::matrix<T>>(column_matrix(evecs)) < 0.)
		evecs[0] = -evecs[0];

	return bOk;
}


template<class T>
bool qr(const ublas::matrix<T>& M,
	ublas::matrix<T>& Q, ublas::matrix<T>& R)
{
	select_func<float, double, decltype(LAPACKE_sgeqrf), decltype(LAPACKE_dgeqrf)> 
		sfunc(LAPACKE_sgeqrf, LAPACKE_dgeqrf);
	auto pfunc = sfunc.get_func<T>();

	const typename ublas::matrix<T>::size_type m = M.size1();
	const typename ublas::matrix<T>::size_type n = M.size2();

	const std::size_t iTauSize = m;//std::min<std::size_t>(m,n);

	std::unique_ptr<T, std::default_delete<T[]>> 
		uptrMem(new T[n*m + iTauSize]);
	T *pMem = uptrMem.get();

	T *pMat = pMem;
	T *pTau = pMem + n*m;

	for(std::size_t i=0; i<m; ++i)
		for(std::size_t j=0; j<n; ++j)
			pMat[i*n + j] = M(i,j);

	// see: http://www.math.utah.edu/software/lapack/lapack-d/dgeqrf.html
	int iInfo = (*pfunc)(LAPACK_ROW_MAJOR, m, n, pMat, n, pTau);
	//std::cout << "dgeqrt: " << iInfo << std::endl;

	R = ublas::matrix<T>(m,n);
	for(std::size_t i=0; i<m; ++i)
		for(std::size_t j=0; j<n; ++j)
		{
			if(j>=i)
				R(i,j) = pMat[i*n + j];
			else
				R(i,j) = 0.;
		}
	//std::cout << "R = " << R << std::endl;

	ublas::vector<T> v(iTauSize);

	const ublas::matrix<T> ident = ublas::identity_matrix<T>(iTauSize);
	Q = ident;

	for(std::size_t k=1; k<=iTauSize; ++k)
	{
		T dTau = pTau[k-1];
		//std::cout << "tau " << k << " = " << dTau << std::endl;

		for(std::size_t i=1; i<=k-1; ++i)
			v[i-1] = 0.;
		v[k-1] = 1.;

		for(std::size_t i=k+1; i<=iTauSize; ++i)
			v[i-1] = pMat[(i-1)*n + (k-1)];

		ublas::matrix<T> VV = ublas::outer_prod(v, ublas::trans(v));
		ublas::matrix<T> H = ident - dTau*VV;

		Q = ublas::prod(Q, H);
	}

	//std::cout << "Q = " << Q << std::endl;
	return (iInfo==0);
}

}

#endif

#endif
