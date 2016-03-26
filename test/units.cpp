// gcc -o units test/units.cpp -std=c++14 -lstdc++

#include <iostream>
#include "../math/units.h"
#include "../math/neutrons.h"
#include "../log/debug.h"

int main()
{
	tl::t_length_si<double> len = 1.*tl::t_meters<double>;
	std::cout << len << std::endl;

	tl::t_energy_si<double> E = 1.*tl::t_meV<double>;
	std::cout << E << std::endl;

	tl::t_energy_si<float> fE = 1.f*tl::t_meV<float>;
	std::cout << fE << std::endl;

	std::cout << "d: " << tl::t_KSQ2E<double> << std::endl;
	std::cout << "ld: " << tl::t_KSQ2E<long double> << std::endl;

	std::cout << tl::get_typename<decltype(tl::co::hbar)>() << std::endl;
	std::cout << tl::get_typename<tl::co::hbar_t>() << std::endl;

	std::cout << tl::get_m_n<double>() << std::endl;
	std::cout << tl::get_m_n<float>() << std::endl;

	std::cout << tl::get_kB<double>() << std::endl;
	std::cout << tl::get_kB<float>() << std::endl;
	return 0;
}
