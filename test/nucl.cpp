#include "../math/nuclear.h"
#include <fstream>

int main()
{
	std::ofstream ofstr("collisions.dat");

	for(double A=1.001; A<6.; A+=0.001)
	{
		tl::t_energy_si<double> E_from = 2. * tl::get_one_MeV<double>();
		tl::t_energy_si<double> E_to = 25. * tl::get_one_meV<double>();

		ofstr << std::setw(25) << std::left << A << " "
			<< tl::mean_collisions(A, E_from, E_to) << std::endl;
	}

	return 0;
}
