/*
 * dispersion relations
 * @author tweber
 * @date 7-jul-15
 */

#ifndef __TLIBS_DISP_H__
#define __TLIBS_DISP_H__

#include <initializer_list>
#include <vector>
#include <cmath>
#include <complex>
#include "linalg.h"

namespace tl {

// ----------------------------------------------------------------------------

// structure factor
template<typename T=double>
std::complex<T> structfact(const std::initializer_list<ublas::vector<T>>& lstAtoms, 
	const ublas::vector<T>& vecG)
{
	constexpr std::complex<T> i(0., 1.);
	std::complex<T> F(0., 0.);

	for(const ublas::vector<T>& vecAtom : lstAtoms)
		F += std::exp(-i * ublas::inner_prod(vecG, vecAtom));

	return F;
}


// ----------------------------------------------------------------------------

// coupling J and atom position
template<typename T=double> using t_atompos = std::pair<T, ublas::vector<T>>;

template<typename T=double>
T ferromag(const std::initializer_list<t_atompos<T>>& lstAtoms, 
	const ublas::vector<T>& vecq, T tS)
{
	std::complex<T> J(0., 0.), J0(0., 0.);
	constexpr std::complex<T> i(0., 1.);

	for(const t_atompos<T>& atom : lstAtoms)
	{
		const T& tJ = atom.first;
		const ublas::vector<T>& vecPos = atom.second;

		J += tJ * std::exp(i * ublas::inner_prod(vecq, vecPos));
		J0 += tJ;
	}

	T tE = 2.*tS*(J0 - J).real();
	return tE;
}

// ----------------------------------------------------------------------------

}
#endif
 
