/*
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

#include "linalg.h"
#include "atoms.h"
#include "nn.h"


namespace tl {
// ----------------------------------------------------------------------------

// coupling J and atom position
template<typename T=double> using t_magatompos = std::pair<std::complex<T>, ublas::vector<T>>;


/**
 * Simple ferromagnetic dispersion
 * @param lstNeighbours list of distances to neighbour atoms and their coupling constants
 * @param vecq q position
 * @param tS spin
 * @return E(q)
 */
template<typename T=double, typename t_cont=std::initializer_list<t_magatompos<T>>>
T ferromag(const t_cont& lstNeighbours, const ublas::vector<T>& vecq, T tS)
{
	std::complex<T> J(0., 0.), J0(0., 0.);

	J = structfact(vec_from_pairvec<1,std::vector,t_cont>()(lstNeighbours), vecq,
		vec_from_pairvec<0,std::vector,t_cont>()(lstNeighbours), &J0).real();

	return T(2)*tS*(J0 - J).real();
}

// ----------------------------------------------------------------------------


// Magnetic form factors
// see Neutron Data Booklet sec. 2.5-1 (p. 60)

template<class T=double>
T j0_avg(T q, T A, T a, T B, T b, T C, T c, T D)
{
	return A * std::exp(-a * q/(T(4)*get_pi<T>())*q/(T(4)*get_pi<T>())) +
		B * std::exp(-b * q/(T(4)*get_pi<T>())*q/(T(4)*get_pi<T>())) +
		C * std::exp(-c * q/(T(4)*get_pi<T>())*q/(T(4)*get_pi<T>())) + D;
}

template<class T=double>
T j2_avg(T q, T A, T a, T B, T b, T C, T c, T D)
{
	return j0_avg(q, A,a, B,b, C,c, D) * q/(T(4)*get_pi<T>())*q/(T(4)*get_pi<T>());
}

template<class T=double>
T mag_formfact(T q, T L, T S, T A, T a, T B, T b, T C, T c, T D)
{
	return (L+2.*S) * j0_avg(q, A,a, B,b, C,c, D) *
		L * j2_avg(q, A,a, B,b, C,c, D);
}
// ----------------------------------------------------------------------------

}
#endif
