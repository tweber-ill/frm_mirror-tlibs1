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

#include "linalg.h"
#include "atoms.h"
#include "nn.h"


namespace tl {
// ----------------------------------------------------------------------------

/**
 * Simple ferromagnetic dispersion
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
	J = structfact(vecNeighbours, vecq, vecJ, &J0).real();
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


// Magnetic form factors
// see: ILL Neutron Data Booklet sec. 2.5-1 (p. 60)
// also see: https://www.ill.eu/sites/ccsl/ffacts/

template<class T=double>
T j0_avg(T Q, T A, T a, T B, T b, T C, T c, T D)
{
	return A * std::exp(-a * Q/(T(4)*get_pi<T>())*Q/(T(4)*get_pi<T>())) +
		B * std::exp(-b * Q/(T(4)*get_pi<T>())*Q/(T(4)*get_pi<T>())) +
		C * std::exp(-c * Q/(T(4)*get_pi<T>())*Q/(T(4)*get_pi<T>())) + D;
}

template<class T=double>
T j2_avg(T Q, T A, T a, T B, T b, T C, T c, T D)
{
	return j0_avg(Q, A,a, B,b, C,c, D) *
		Q/(T(4)*get_pi<T>()) * Q/(T(4)*get_pi<T>());
}

template<class T=double>
T mag_formfact(T Q, T L, T S,
	T A0, T a0, T B0, T b0, T C0, T c0, T D0,
	T A2, T a2, T B2, T b2, T C2, T c2, T D2)
{
	return (L+T(2)*S) * j0_avg(Q, A0,a0, B0,b0, C0,c0, D0) *
		L * j2_avg(Q, A2,a2, B2,b2, C2,c2, D2);
}

// see: Squires, p. 139
template<class T=double>
T mag_formfact(T Q, T L, T S, T J,
	T A0, T a0, T B0, T b0, T C0, T c0, T D0,
	T A2, T a2, T B2, T b2, T C2, T c2, T D2)
{
	T j0 = j0_avg(Q, A0,a0, B0,b0, C0,c0, D0);
	T j2 = j2_avg(Q, A2,a2, B2,b2, C2,c2, D2);

	T gL = T(0.5) + (L*(L+T(1)) - S*(S+T(1))) / (T(2)*J* (J+T(1)));
	T gS = T(1) + (S*(S+T(1)) - L*(L+T(1))) / (J * (J+T(1)));

	return (gS*j0 + gL*(j0+j2)) / (gL + gS);
}
// ----------------------------------------------------------------------------

}
#endif
