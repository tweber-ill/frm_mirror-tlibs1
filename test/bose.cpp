// gcc -I. -o bose test/bose.cpp -std=c++11 -lstdc++ -lm

#include <iostream>
#include <fstream>
#include "math/neutrons.hpp"

int main()
{
	std::ofstream ofstr("bosetst.dat");
	double dT = 80.;

	for(double dE=1.; dE<10.; dE+=0.1)
	{
		double dBose = tl::bose<double>(dE, dT);
		ofstr << dE << "\t" << dBose << "\n";
	}

	return 0;
}
