// clang -o struct test/struct.cpp -lstdc++ -lm -std=c++11

#include "../math/disp.h"
#include <iostream>

int main()
{
	std::complex<double> F = tl::structfact({
		tl::make_vec({1., 0., 0.}),
		tl::make_vec({-1., 0., 0.}),
		tl::make_vec({0., 1., 0.}),
		tl::make_vec({0., -1., 0.}),
		tl::make_vec({0., 0., 1.}),
		tl::make_vec({0., 0., -1.}),
			},
		tl::make_vec({1, 0., 0.}));

	std::cout << F << std::endl;
	return 0;
}
