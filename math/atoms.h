/**
 * atoms and structural calculations
 * @author tweber
 * @date nov-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __ATOMS_H__
#define __ATOMS_H__


#include "linalg.h"
#include "linalg_ops.h"


namespace tl{

/**
 * Maps atom position back to units cell
 */
template<class t_vec>
void restrict_to_uc(t_vec& vec,
	typename t_vec::value_type tMin=0, typename t_vec::value_type tMax=1)
{
	using T = typename t_vec::value_type;

	for(std::size_t i=0; i<vec.size(); ++i)
	{
		vec[i] = std::fmod(vec[i], T(1));

		while(vec[i] < tMin) vec[i] += T(1);
		while(vec[i] >= tMax) vec[i] -= T(1);
	}
}

/**
 * Generates atom positions using trafo matrices
 */
template<class t_mat, class t_vec, template<class ...Args> class t_cont>
t_cont<t_vec> generate_atoms(const t_cont<t_mat>& trafos, const t_vec& vecAtom,
	typename t_vec::value_type tUCMin=0, typename t_vec::value_type tUCMax=1)
{
	//typedef typename t_vec::value_type t_real;
	t_cont<t_vec> vecvecRes;

	for(const t_mat& mat : trafos)
	{
		t_vec vecRes = mat * vecAtom;
		restrict_to_uc<t_vec>(vecRes, tUCMin, tUCMax);

		bool bPushBack = 1;
		// already have pos?
		for(const t_vec& vecOld : vecvecRes)
		{
			if(vec_equal(vecOld, vecRes))
			{
				bPushBack = 0;
				break;
			}
		}

		if(bPushBack)
			vecvecRes.push_back(std::move(vecRes));
	}

	return vecvecRes;
}



// ----------------------------------------------------------------------------


/**
 * calculates atomic form factors
 * @param G Length of lattice vector
 * @param vecA "a" coefficients
 * @param vecB "b" coefficients
 * @param c "c" coefficient
 * @return form factor
 * see: Waasmaier and Kirfel, Acta Cryst. A51, 416-431 (1995)
 */
template<class T=double, template<class...> class t_cont>
T formfact(T G, const t_cont<T>& vecA, const t_cont<T>& vecB, T c)
{
	T ff = T(0);
	T s = G / T(4.*M_PI);

	typename t_cont<T>::const_iterator iterA = vecA.begin();
	typename t_cont<T>::const_iterator iterB = vecB.begin();

	for(; iterA!=vecA.end() && iterB!=vecB.end(); ++iterA, ++iterB)
		ff += (*iterA)*std::exp(-(*iterB)*s*s);
	ff += c;

	return ff;
}



/**
 * calculates the structure factor F
 * @param lstAtoms List of atom positions
 * @param vecG Lattice vector
 * @param lstf G-dependent Atomic form factors (x-rays) or coherent scattering length (neutrons)
 * @return structure factor
 */
template<typename T=double, typename t_ff = std::complex<T>,
	class t_vec=ublas::vector<T>,
	template<class ...> class t_cont=std::initializer_list>
std::complex<T> structfact(const t_cont<t_vec>& lstAtoms, const t_vec& vecG, 
	const t_cont<t_ff>& lstf = t_cont<t_ff>())
{
	constexpr std::complex<T> i(0., 1.);
	std::complex<T> F(0., 0.);

	using t_iter_atoms = typename t_cont<t_vec>::const_iterator;
	using t_iter_ffact = typename t_cont<t_ff>::const_iterator;

	t_iter_atoms iterAtom = lstAtoms.begin();
	t_iter_ffact iterFFact = lstf.begin();

	for(; iterAtom!=lstAtoms.end(); ++iterAtom)
	{
		// only use form factors or scattering lengths when available
		t_ff tFF = T(1);
		if(iterFFact != lstf.end())
			tFF = *iterFFact;

		F += tFF * std::exp(i * ublas::inner_prod(vecG, *iterAtom));

		if(iterFFact != lstf.end())
			++iterFFact;
	}

	return F;
}


/**
 * Lorentz factor
 * @param twotheta Scattering angle in rad
 */
template<typename T=double>
T lorentz_factor(T twotheta)
{
	T theta = 0.5*twotheta;
	return T(0.25) / (std::sin(theta)*std::sin(theta) * std::cos(theta));
}

/**
 * Lorentz polarisation factor (only for x-rays)
 * @param twotheta Scattering angle in rad
 */
template<typename T=double>
T lorentz_pol_factor(T twotheta)
{
	return T(0.5) + T(0.5)*std::cos(twotheta)*std::cos(twotheta);
}

}
#endif
