// gcc -o units test/units.cpp -std=c++14 -lstdc++

#include <iostream>
#include "../math/units.h"

int main()
{
	tl::t_length_si<double> len = 1.*tl::t_meters<double>;
	std::cout << len << std::endl;

	tl::t_energy_si<double> E = 1.*tl::t_meV<double>;
	std::cout << E << std::endl;

	tl::t_energy_si<float> fE = 1.f*tl::t_meV<float>;
	std::cout << fE << std::endl;

	return 0;
}
