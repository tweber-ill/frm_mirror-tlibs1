/**
 * statistical functions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 2013-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_STAT_H__
#define __TLIBS_STAT_H__

#include "linalg.h"

namespace tl {

/**
 * mean value
 */
template<class vec_type>
typename vec_type::value_type mean_value(const vec_type& vec)
{
	typedef typename vec_type::value_type T;
	if(vec.size()==0) return T(0);
	else if(vec.size()==1) return vec[0];

	T tMean = vec[0];
	for(std::size_t i=1; i<vec.size(); ++i)
		tMean += vec[i];
	tMean /= vec.size();

	return tMean;
}

/**
 * mean value with given probability
 */
template<class vec_type_prob, class vec_type>
typename vec_type::value_type mean_value(const vec_type_prob& vecP, const vec_type& vec)
{
	typedef typename vec_type::value_type T;
	typedef typename vec_type_prob::value_type Tprob;
	std::size_t iSize = std::min(vecP.size(), vec.size());

	if(iSize==0) return T(0);

	T tMean = vecP[0]*vec[0];
	Tprob tProbTotal = vecP[0];
	for(std::size_t i=1; i<iSize; ++i)
	{
		tMean += vecP[i]*vec[i];
		tProbTotal += vecP[i];
	}
	tMean /= tProbTotal;

	return tMean;
}


/**
 * standard deviation of mean value, with correction factor
 * see e.g.: https://en.wikipedia.org/wiki/Bessel%27s_correction
 */
template<class vec_type>
typename vec_type::value_type std_dev(const vec_type& vec, bool bCorr=1)
{
	typedef typename vec_type::value_type T;
	if(vec.size()<=1) return T(0);

	T tProb = T(vec.size());
	if(bCorr) tProb -= T(1);

	T tMean = mean_value(vec);
	T t = T(0);
	for(const T& tval : vec)
		t += (tval-tMean) * (tval-tMean);
	t /= tProb;

	return std::sqrt(t);
}

/**
 * standard deviation with given probability
 */
template<class vec_type_prob, class vec_type>
typename vec_type::value_type std_dev(const vec_type_prob& vecP, const vec_type& vec)
{
	typedef typename vec_type::value_type T;
	std::size_t iSize = std::min(vecP.size(), vec.size());
	if(iSize<=1) return T(0);

	T tMean = mean_value<vec_type_prob, vec_type>(vecP, vec);
	T t = T(0);
	T tProbTotal = T(0);

	for(std::size_t iIdx = 0; iIdx<iSize; ++iIdx)
	{
		t += (vec[iIdx]-tMean)*(vec[iIdx]-tMean) * vecP[iIdx];
		tProbTotal += vecP[iIdx];
	}
	t /= tProbTotal;

	return std::sqrt(t);
}


// -----------------------------------------------------------------------------


/**
 * calculates the covariance and the correlation matrices
 * covariance: C_ij = cov(X_i, X_j) = < (X_i - <X_i>) * (X_j - <X_j>) >
 * correlation: K_ij = C_ij / (sigma_i sigma_j)
 * see e.g.: http://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm
 * see also e.g.: (Arfken 2013) p. 1142
 */
template<typename T=double>
std::tuple<ublas::matrix<T>, ublas::matrix<T>>
covariance(const std::vector<ublas::vector<T>>& vecVals, const std::vector<T>* pProb = 0)
{
	using t_mat = ublas::matrix<T>;
	using t_vecvec = typename std::remove_reference<decltype(vecVals)>::type;
	using t_innervec_org = decltype(vecVals[0]);
	using t_innervec = typename std::remove_const<
		typename std::remove_reference<t_innervec_org>::type>::type;

	if(vecVals.size() == 0) return std::make_tuple(t_mat(), t_mat());

	// mean vector <X_i>
	t_innervec vecMean;
	if(pProb)
		vecMean = mean_value<std::vector<T>, t_vecvec>(*pProb, vecVals);
	else
		vecMean = mean_value<t_vecvec>(vecVals);

	t_mat matCov(vecVals[0].size(), vecVals[0].size());
	T tSum = T(0);
	const std::size_t N = vecVals.size();

	for(std::size_t i=0; i<N; ++i)
	{
		T tprob = T(1);

		// X_i - <X_i>
		t_innervec vec = vecVals[i] - vecMean;

		// matrix elements, AA^t
		t_mat matOuter = ublas::outer_prod(vec, vec);

		// probabilities for final averaging, <...>
		if(pProb)
		{
			tprob = (*pProb)[i];
			matOuter *= tprob;
		}

		matCov += matOuter;
		tSum += tprob;
	}

	// average, sometimes defined as C /= (N-1)
	matCov /= tSum /*-T(1)*/;


	// --------------------------------------------------------------------------------
	// correlation matrix
	t_innervec vecVar = diag_vec(matCov);
	t_innervec vecStdDev(vecVar.size());

	std::transform(vecVar.begin(), vecVar.end(), vecStdDev.begin(),
		[](typename t_innervec::value_type d) -> typename t_innervec::value_type
		{ return std::sqrt(d); });

	t_mat matStdDev = ublas::outer_prod(vecStdDev, vecStdDev);
	t_mat matCorr = ublas::element_div(matCov, matStdDev);
	// --------------------------------------------------------------------------------
	
	return std::make_tuple(matCov, matCorr);
}
}

#endif
