/**
 * Swarming algorithms
 *
 * @author Tobias Weber
 * @date Feb-17
 * @license GPLv2 or GPLv3
 */

#ifndef __SWARM_H__
#define __SWARM_H__

#include <vector>
#include <functional>

#include "../math/linalg.h"
#include "../math/rand.h"
#include "../log/log.h"


namespace tl {

template<class t_real, template<class...> class t_vec>
struct Raven
{
	t_vec<t_real> vecPos, vecBestPos;
	t_vec<t_real> vecVel;
};


/**
 * swarm minimisation
 * algorithm: https://en.wikipedia.org/wiki/Particle_swarm_optimization
 */
template<class t_real, template<class...> class t_vec>
class Unkindness
{
protected:
	std::vector<Raven<t_real, t_vec>> m_vecRavens;
	t_vec<t_real> m_vecBestPos;

	// function to minimise
	std::function<t_real(t_vec<t_real>)> m_func;

public:
	void Init(std::size_t iNumRavens,
		const t_vec<t_real>& vecMin, const t_vec<t_real>& vecMax)
	{
		m_vecRavens.clear();
		m_vecRavens.reserve(iNumRavens);

		// random initial best position
		m_vecBestPos = convert_vec_full<t_real, t_real, std::vector, t_vec>(
			rand_minmax_nd<t_real, std::vector>(
			convert_vec_full<t_real, t_real, t_vec, std::vector>(vecMin),
			convert_vec_full<t_real, t_real, t_vec, std::vector>(vecMax)));

		t_vec<t_real> vecVelMin = vecMin-vecMax;
		t_vec<t_real> vecVelMax = vecMax-vecMin;

		for(std::size_t i=0; i<iNumRavens; ++i)
		{
			Raven<t_real, t_vec> raven;

			// random positions
			raven.vecPos = raven.vecBestPos =
				convert_vec_full<t_real, t_real, std::vector, t_vec>(
					rand_minmax_nd<t_real, std::vector>(
					convert_vec_full<t_real, t_real, t_vec, std::vector>(vecMin),
					convert_vec_full<t_real, t_real, t_vec, std::vector>(vecMax)));

			// random velocities
			raven.vecVel =
				convert_vec_full<t_real, t_real, std::vector, t_vec>(
					rand_minmax_nd<t_real, std::vector>(
					convert_vec_full<t_real, t_real, t_vec, std::vector>(vecVelMin),
					convert_vec_full<t_real, t_real, t_vec, std::vector>(vecVelMax)));

			if(m_func(raven.vecPos) < m_func(m_vecBestPos))
				m_vecBestPos = raven.vecPos;

			m_vecRavens.emplace_back(std::move(raven));
		}
	}


	void Run()
	{
		// TODO
	}
};

}

#endif
