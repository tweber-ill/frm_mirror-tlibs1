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


	ublas::matrix<double> M = tl::make_mat({{1.1,-2.2,3.3},{-4.4,5.5,6.6}, {7.7,8.8,-9.9}});
	std::cout << M << std::endl;
	std::cout << tl::insert_unity(M, 2) << std::endl;


	ublas::matrix<double> Q, R;
	if(tl::qr_decomp(M, Q, R))
	{
		std::cout << "M = " << M << std::endl;
		std::cout << "Q = " << Q << std::endl;
		std::cout << "R = " << R << std::endl;
		std::cout << "QR = " << ublas::prod(Q,R) << std::endl;
	}

	return 0;
}
