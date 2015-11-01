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


}
#endif
