#include "../math/linalg.h"

namespace ublas = boost::numeric::ublas;

int main()
{
	ublas::vector<double> vec1 = tl::make_vec({1.,2.,3.});
	ublas::vector<double> vec2 = tl::make_vec({2.,-3.,4.});
	ublas::vector<double> vec3 = tl::make_vec({3.,4.,5.});

	std::vector<ublas::vector<double>> vecsIn = {vec1, vec2, vec3};
	std::vector<ublas::vector<double>> vecsOut = tl::gram_schmidt(vecsIn);

	for(const ublas::vector<double>& vec : vecsOut)
		std::cout << vec << std::endl;

	std::cout << "v0 * v1 = " << ublas::inner_prod(vecsOut[0], vecsOut[1]) << std::endl;
	std::cout << "v0 * v2 = " << ublas::inner_prod(vecsOut[0], vecsOut[2]) << std::endl;
	std::cout << "v1 * v2 = " << ublas::inner_prod(vecsOut[1], vecsOut[2]) << std::endl;

	return 0;
}
