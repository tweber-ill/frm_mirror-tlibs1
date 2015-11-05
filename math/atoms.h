/**
 * atoms
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
 * Generates atom positions using trafo matrices
 */
template<class t_mat, class t_vec, template<class ...Args> class t_cont>
t_cont<t_vec> generate_atoms(const t_cont<t_mat>& trafos, const t_vec& vecAtom)
{
	typedef typename t_mat::value_type t_real;
	t_cont<t_vec> vecvecRes;

	for(const t_mat& mat : trafos)
	{
		t_vec vecRes = mat * vecAtom;

		bool bPushBack = 1;
		// already have pos?
		for(const t_vec& vecOld : vecvecRes)
			if(vec_equal(vecOld, vecRes))
			{
				bPushBack = 0;
				break;
			}

		if(bPushBack)
			vecvecRes.push_back(vecRes);
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
 * calculates the structure factor
 * @param lstAtoms: List of atom positions
 * @param lstf: Atomic form factors
 * @param vecG: Lattice vector
 * @return structure factor
 */
template<typename T=double,
	class t_vec=ublas::vector<T>,
	template<class ...> class t_cont=std::initializer_list>
std::complex<T> structfact(const t_cont<t_vec>& lstAtoms, const t_cont<T>& lstf, const t_vec& vecG)
{
	constexpr std::complex<T> i(0., 1.);
	std::complex<T> F(0., 0.);

	using t_iter_atoms = typename t_cont<t_vec>::const_iterator;
	using t_iter_ffact = typename t_cont<T>::const_iterator;

	t_iter_atoms iterAtom=lstAtoms.begin();
	t_iter_ffact iterFFact=lstf.begin();

	for(; iterAtom!=lstAtoms.end(); ++iterAtom, ++iterFFact)
		F += *iterFFact * std::exp(-i * ublas::inner_prod(vecG, *iterAtom));

	return F;
}


}
#endif
