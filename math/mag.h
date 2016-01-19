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

namespace tl {


// ----------------------------------------------------------------------------

// atom positions
enum class UCType
{
	SIMPLE,
	FCC,
	BCC,
};

/**
 * Next neighbours
 * iDist == 0: nearest neighbours
 * iDist == 1: next-nearest neighbours
 */
template<typename T=double>
std::vector<ublas::vector<T>> get_neighbour_atoms(UCType crys, int iDist=0, T a=1.)
{
	std::vector<ublas::vector<T>> vecAtoms;

	if(crys == UCType::SIMPLE)
	{
		if(iDist == 0)
		{
			vecAtoms.push_back(make_vec({1., 0., 0.}));
			vecAtoms.push_back(make_vec({0., 1., 0.}));
			vecAtoms.push_back(make_vec({0., 0., 1.}));
			vecAtoms.push_back(make_vec({-1., 0., 0.}));
			vecAtoms.push_back(make_vec({0., -1., 0.}));
			vecAtoms.push_back(make_vec({0., 0., -1.}));
		}
		else if(iDist == 1)
		{
			vecAtoms.push_back(make_vec({1., 1., 0.}));
			vecAtoms.push_back(make_vec({1., -1., 0.}));
			vecAtoms.push_back(make_vec({-1., 1., 0.}));
			vecAtoms.push_back(make_vec({-1., -1., 0.}));
			vecAtoms.push_back(make_vec({1., 0., 1.}));
			vecAtoms.push_back(make_vec({1., 0., -1.}));
			vecAtoms.push_back(make_vec({-1., 0., 1.}));
			vecAtoms.push_back(make_vec({-1., 0., -1.}));
			vecAtoms.push_back(make_vec({0., 1., 1.}));
			vecAtoms.push_back(make_vec({0., 1., -1.}));
			vecAtoms.push_back(make_vec({0., -1., 1.}));
			vecAtoms.push_back(make_vec({0., -1., -1.}));
		}
	}
	else if(crys == UCType::FCC)
	{
		if(iDist == 0)
		{
			vecAtoms.push_back(make_vec({0.5, 0.5, 0.}));
			vecAtoms.push_back(make_vec({0.5, -0.5, 0.}));
			vecAtoms.push_back(make_vec({-0.5, 0.5, 0.}));
			vecAtoms.push_back(make_vec({-0.5, -0.5, 0.}));
			vecAtoms.push_back(make_vec({0.5, 0., 0.5}));
			vecAtoms.push_back(make_vec({0.5, 0., -0.5}));
			vecAtoms.push_back(make_vec({-0.5, 0., 0.5}));
			vecAtoms.push_back(make_vec({-0.5, 0., -0.5}));
			vecAtoms.push_back(make_vec({0., 0.5, 0.5}));
			vecAtoms.push_back(make_vec({0., 0.5, -0.5}));
			vecAtoms.push_back(make_vec({0., -0.5, 0.5}));
			vecAtoms.push_back(make_vec({0., -0.5, -0.5}));
		}
		else if(iDist == 1)
		{
			vecAtoms.push_back(make_vec({1., 0., 0.}));
			vecAtoms.push_back(make_vec({0., 1., 0.}));
			vecAtoms.push_back(make_vec({0., 0., 1.}));
			vecAtoms.push_back(make_vec({-1., 0., 0.}));
			vecAtoms.push_back(make_vec({0., -1., 0.}));
			vecAtoms.push_back(make_vec({0., 0., -1.}));
		}
	}
	else if(crys == UCType::BCC)
	{
		if(iDist == 0)
		{
			vecAtoms.push_back(make_vec({0.5, 0.5, 0.5}));
			vecAtoms.push_back(make_vec({0.5, 0.5, -0.5}));
			vecAtoms.push_back(make_vec({0.5, -0.5, 0.5}));
			vecAtoms.push_back(make_vec({0.5, -0.5, -0.5}));
			vecAtoms.push_back(make_vec({-0.5, 0.5, 0.5}));
			vecAtoms.push_back(make_vec({-0.5, 0.5, -0.5}));
			vecAtoms.push_back(make_vec({-0.5, -0.5, 0.5}));
			vecAtoms.push_back(make_vec({-0.5, -0.5, -0.5}));
		}
		else if(iDist == 1)
		{
			vecAtoms.push_back(make_vec({1., 0., 0.}));
			vecAtoms.push_back(make_vec({0., 1., 0.}));
			vecAtoms.push_back(make_vec({0., 0., 1.}));
			vecAtoms.push_back(make_vec({-1., 0., 0.}));
			vecAtoms.push_back(make_vec({0., -1., 0.}));
			vecAtoms.push_back(make_vec({0., 0., -1.}));
		}
	}

	if(a != 1.)
	for(ublas::vector<T>& vec : vecAtoms)
		vec *= a;

	return vecAtoms;
}



// ----------------------------------------------------------------------------

// coupling J and atom position
template<typename T=double> using t_magatompos = std::pair<T, ublas::vector<T>>;


/**
 * Simple ferromagnetic dispersion
 * @param lstAtoms list of atoms and their coupling constants
 * @param vecq q position
 * @param tS spin
 * @return E(q)
 */
template<typename T=double, typename t_cont=std::initializer_list<t_magatompos<T>>>
T ferromag(const t_cont& lstAtoms, const ublas::vector<T>& vecq, T tS)
{
	std::complex<T> J(0., 0.), J0(0., 0.);
	constexpr std::complex<T> i(0., 1.);

	for(const t_magatompos<T>& atom : lstAtoms)
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
