// clang -o struct test/struct.cpp -lstdc++ -lm -std=c++11

#include "../math/atoms.h"
#include "../math/disp.h"
#include <iostream>

int main()
{
	double a = 5.;
	double h = 0., k = 0., l = 0.;

	while(1)
	{
		std::cout << "Enter hkl: ";
		std::cin >> h >> k >> l;

		std::complex<double> F = tl::structfact(
			{
				tl::make_vec({0., 0., 0.}),
				tl::make_vec({a*0.5, a*0.5, a*0.5}),
			},
			tl::make_vec({h*2.*M_PI/a, k*2.*M_PI/a, l*2.*M_PI/a}),
			{1., 1.});

		std::cout << "F = " << F << std::endl;
	}
	return 0;
}
