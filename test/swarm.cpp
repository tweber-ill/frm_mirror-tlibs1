/**
 * Swarming algorithms
 *
 * @author Tobias Weber
 * @date Feb-17
 * @license GPLv2 or GPLv3
 * 
 * gcc -o swarm swarm.cpp ../math/rand.cpp ../log/log.cpp -lstdc++ -lm -lpthread
 */

#include "../fit/swarm.h"
namespace ublas = tl::ublas;

using t_real = double;
template<class T> using t_vec = ublas::vector<T>;

int main()
{
	tl::Unkindness<t_real, t_vec> unk;
	unk.Init(1000, tl::make_vec({0.,0.,0.}), tl::make_vec({1.,1.,1.}));

	
    return 0;
}
