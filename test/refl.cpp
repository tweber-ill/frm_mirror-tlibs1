#include "../math/linalg.h"

namespace ublas = boost::numeric::ublas;

int main()
{
	ublas::vector<double> vec = tl::make_vec({1.,2.,3.});
	ublas::vector<double> vecN = tl::make_vec({1.,2.,0.});
	ublas::vector<double> vecRefl = tl::reflection(vec, vecN);

	std::cout << "vec: " << vec << std::endl;
	std::cout << "norm: " << vecN << std::endl;
	std::cout << "refl: " << vecRefl << std::endl;

	return 0;
}
