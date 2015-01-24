/*
 * advanced linalg helpers
 *
 * @author: tweber
 * @date: 30-apr-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __MIEZE_LINALG2__
#define __MIEZE_LINALG2__


#include "math.h"
#include "linalg.h"

namespace tl {

#ifndef NO_LAPACK
	#define USE_LAPACK
#endif

#ifdef NO_LAPACK

template<typename T=double>
bool qr(const ublas::matrix<T>& M, ublas::matrix<T>& Q, ublas::matrix<T>& R)
{
	return qr_decomp(M, Q, R);
}

template<typename T=double>
bool eigenvec_sym(const ublas::matrix<T>& mat, std::vector<ublas::vector<T>>& evecs, std::vector<T>& evals)
{
	return eigenvec_sym_simple(mat, evecs, evals);
}

#else

template<typename T=double>
bool qr(const ublas::matrix<T>& M, ublas::matrix<T>& Q, ublas::matrix<T>& R)
{
	log_err("No specialisation of \"qr\" available for this type.");
	return false;
}

template<>
bool qr<double>(const ublas::matrix<double>& M,
		ublas::matrix<double>& Q,
		ublas::matrix<double>& R);


template<typename T=double>
bool eigenvec(const ublas::matrix<T>& mat, 
		std::vector<ublas::vector<T> >& evecs_real, std::vector<ublas::vector<T> >& evecs_imag, 
		std::vector<T>& evals_real, std::vector<T>& evals_imag)
{
	log_err("No specialisation of \"eigenvec\" available for this type.");
	return false;
}

template<typename T=double>
bool eigenvec_sym(const ublas::matrix<T>& mat, std::vector<ublas::vector<T> >& evecs, std::vector<T>& evals)
{
	log_err("No specialisation of \"eigenvec_sym\" available for this type.");
	return false;
}

template<>
bool eigenvec<double>(const ublas::matrix<double>& mat,
			std::vector<ublas::vector<double> >& evecs_real,
			std::vector<ublas::vector<double> >& evecs_imag,
			std::vector<double>& evals_real, 
			std::vector<double>& evals_imag);
template<>
bool eigenvec_sym<double>(const ublas::matrix<double>& mat,
			std::vector<ublas::vector<double> >& evecs,
			std::vector<double>& evals);
#endif

}

#endif
