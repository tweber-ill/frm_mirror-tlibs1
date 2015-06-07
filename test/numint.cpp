#include <iostream>
#include <cmath>
#include "../math/numint.hpp"


double fkt(double x)
{
	return x*x*x;
}

int main()
{
	std::function<double(double)> f(fkt);
	std::cout << "rect: " << tl::numint_rect(f, 3., 5., 128) << std::endl;
	std::cout << "trap: " << tl::numint_trap(f, 3., 5.) << std::endl;
	std::cout << "trapN: " << tl::numint_trapN(f, 3., 5., 128) << std::endl;
	std::cout << "simp: " << tl::numint_simp(f, 3., 5.) << std::endl;
	std::cout << "simpN: " << tl::numint_simpN(f, 3., 5., 128) << std::endl;

	std::cout << "calc: " << 0.25*5.*5.*5.*5. - 0.25*3.*3.*3.*3. << std::endl;

	return 0;
}
