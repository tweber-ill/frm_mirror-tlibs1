#include "../math/math.h"
#include "../math/linalg.h"
#include <complex>
#include <iostream>

int main()
{
	std::cout << tl::get_epsilon<double>() << std::endl;
	std::cout << tl::get_epsilon<std::complex<double>>() << std::endl;
	std::cout << tl::get_epsilon<tl::ublas::matrix<std::complex<double>>>() << std::endl;

	std::cout << tl::float_equal<double>(0., 1.) << std::endl;
	std::cout << tl::float_equal<std::complex<double>>(
		std::complex<double>(0.,1.), std::complex<double>(1.,1.))
		<< std::endl;
	return 0;
}
