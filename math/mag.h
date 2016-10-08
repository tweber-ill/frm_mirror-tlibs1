/**
 * magnetic dispersion relations
 * @author tweber
 * @date 7-jul-15
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_MAGDISP_H__
#define __TLIBS_MAGDISP_H__

#include <initializer_list>
#include <vector>
#include <cmath>
#include <complex>
#include <cassert>

#include "linalg.h"
#include "atoms.h"
#include "nn.h"
#include "rand.h"


namespace tl {
// ----------------------------------------------------------------------------

/**
 * Simple ferromagnetic dispersion (see e.g. Squires p. 161)
 * @param lstNeighbours list of distances to neighbour atoms and their coupling constants
 * @param vecq q position
 * @param tS spin
 * @return E(q)
 */
template<class t_vec = ublas::vector<double>,
	typename T = typename t_vec::value_type,
	template<class...> class t_cont = std::vector>
T ferromag(const t_cont<t_vec>& vecNeighbours, const t_cont<std::complex<T>>& vecJ,
	const ublas::vector<T>& vecq, T tS)
{
	std::complex<T> J(0., 0.), J0(0., 0.);
	J = structfact<T, std::complex<T>, t_vec, t_cont>
		(vecNeighbours, vecq, vecJ, &J0).real();
	return T(2)*tS*(J0 - J).real();
}

template<class t_vec = ublas::vector<double>,
	typename T = typename t_vec::value_type,
	typename t_cont = std::initializer_list<std::tuple<t_vec, std::complex<T>>>>
T ferromag(const t_cont& lstNeighbours, const ublas::vector<T>& vecq, T tS)
{
	return ferromag(vec_from_pairvec<0,std::vector,t_cont>()(lstNeighbours),
		vec_from_pairvec<1,std::vector,t_cont>()(lstNeighbours),
		vecq, tS);
}
// ----------------------------------------------------------------------------


/**
 * Magnetic form factors
 * see: ILL Neutron Data Booklet sec. 2.5-1 (p. 60)
 * also see: https://www.ill.eu/sites/ccsl/ffacts/
 */
template<class T=double, template<class...> class t_vec=std::initializer_list>
T j0_avg(T Q, const t_vec<T>& A, const t_vec<T>& a)
{
	assert(A.size() == a.size()+1);

	T tJ = T(0);
	for(std::size_t i=0; i<a.size(); ++i)
		tJ += A[i] * std::exp(-a[i] * Q/(T(4)*get_pi<T>())*Q/(T(4)*get_pi<T>()));
	tJ += *A.rbegin();
	return tJ;
}

template<class T=double, template<class...> class t_vec=std::initializer_list>
T j2_avg(T Q, const t_vec<T>& A, const t_vec<T>& a)
{
	return j0_avg<T, t_vec>(Q, A, a) * Q/(T(4)*get_pi<T>()) * Q/(T(4)*get_pi<T>());
}

template<class T=double, template<class...> class t_vec=std::initializer_list>
T mag_formfact(T Q, T L, T S,
	const t_vec<T>& A0, const t_vec<T>& a0,
	const t_vec<T>& A2, const t_vec<T>& a2)
{
	return (L+T(2)*S) * j0_avg<T, t_vec>(Q, A0, a0) * L * j2_avg<T, t_vec>(Q, A2, a2);
}

/**
 * @desc see: Squires, p. 139
 */
template<class T=double, template<class...> class t_vec=std::initializer_list>
T mag_formfact(T Q, T L, T S, T J,
	const t_vec<T>& A0, const t_vec<T>& a0,
	const t_vec<T>& A2, const t_vec<T>& a2)
{
	T j0 = j0_avg<T, t_vec>(Q, A0, a0);
	T j2 = j2_avg<T, t_vec>(Q, A2, a2);

	T gL = T(0.5) + (L*(L+T(1)) - S*(S+T(1))) / (T(2)*J* (J+T(1)));
	T gS = T(1) + (S*(S+T(1)) - L*(L+T(1))) / (J * (J+T(1)));

	return (gS*j0 + gL*(j0+j2)) / (gL + gS);
}
// ----------------------------------------------------------------------------


/**
 * metropolis algorithm
 */
template<class t_real, std::size_t DIM,
	template<class, std::size_t, class...> class t_arr_1d = boost::array,
	template<class, std::size_t, class...> class t_arr_nd = boost::multi_array>
t_arr_nd<bool, DIM> metrop(
	const t_arr_1d<typename t_arr_nd<bool,DIM>::index, DIM>& arrDims,
	std::size_t iNumIters, t_real dJ, t_real dk, t_real dT,
	t_real *pEtot = nullptr)
{
	using T = bool;
	using t_arr = t_arr_nd<T, DIM>;
	using t_idx = typename t_arr::index;
	using t_dim = t_arr_1d<t_idx, DIM>;

	const t_real dBeta = t_real(1)/(dk*dT);


	// get next neighbours
	auto getNN = [](const t_dim& dim, const t_dim& idx) -> std::vector<t_dim>
	{
		std::vector<t_dim> vecNN;
		for(std::size_t iDim=0; iDim<DIM; ++iDim)
		{
			if(idx[iDim] > 0)
			{
				t_dim idxNew = idx;
				--idxNew[iDim];
				vecNN.emplace_back(std::move(idxNew));
			}
			if(idx[iDim] < dim[iDim]-1)
			{
				t_dim idxNew = idx;
				++idxNew[iDim];
				vecNN.emplace_back(std::move(idxNew));
			}
		}
		return vecNN;
	};

	// calculate energy
	auto calcE = [dJ](const t_arr& arr,
		const t_dim& idxSpin, const std::vector<t_dim>& vecNN,
		bool bFlip=0) -> t_real
	{
		t_real dE = t_real(0);
		bool bSpin = arr(idxSpin);
		if(bFlip) bSpin = !bSpin;
		t_real dSpin = (bSpin ? t_real(1) : t_real(-1));

		for(const t_dim& idxNN : vecNN)
		{
			t_real dSpinNN = (arr(idxNN) ? t_real(1) : t_real(-1));
			dE += dJ*dSpin*dSpinNN;
		}
		return dE;
	};


	// TODO: search radius
	t_dim dimMin; dimMin.fill(0);
	t_dim dimMax = arrDims;

	t_arr arrSpins = rand_array<T, DIM, t_arr_1d, t_arr_nd>(arrDims);

	for(std::size_t iIter=0; iIter<iNumIters; ++iIter)
	{
		t_dim idx = rand_idx<t_idx, DIM, t_arr_1d>(dimMin, dimMax);
		std::vector<t_dim> vecNN = getNN(arrDims, idx);

		t_real dENoFlip = calcE(arrSpins, idx, vecNN, 0);
		t_real dEFlip = calcE(arrSpins, idx, vecNN, 1);

		if(dEFlip < dENoFlip)
		{
			arrSpins(idx) = !arrSpins(idx);
		}
		else
		{
			t_real dEDiff = dEFlip - dENoFlip;
			t_real dProb = std::exp(-dBeta * dEDiff);

			if(rand_prob<t_real>(dProb))
				arrSpins(idx) = !arrSpins(idx);
		}
	}

	// calculate total energy
	if(pEtot)
	{
		t_dim idxCur; idxCur.fill(0);
		while(1)
		{
			std::vector<t_dim> vecNN = getNN(arrDims, idxCur);
			*pEtot += calcE(arrSpins, idxCur, vecNN, 0);

			++idxCur[0];
			for(std::size_t iDim=0; iDim<DIM-1; ++iDim)
			{
				if(idxCur[iDim] >= arrDims[iDim])
				{
					idxCur[iDim] = 0;
					++idxCur[iDim+1];
				}
			}
			if(idxCur[DIM-1] >= arrDims[DIM-1])
				break;
		}

		for(std::size_t iDim=0; iDim<DIM; ++iDim)
			*pEtot /= t_real(arrDims[iDim]);
	}

	return arrSpins;
}


// ----------------------------------------------------------------------------

}
#endif
