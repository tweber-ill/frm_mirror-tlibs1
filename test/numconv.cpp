#include <iostream>
#include <vector>
#include "../math/numint.hpp"


int main()
{
	std::vector<double> vecA = {1.,2.,3.,4.,5.};
	std::vector<double> vecB = {3.,2.,1.};

	std::vector<double> vecC = tl::convolute_discrete(vecA, vecB);

	for(double d : vecC)
		std::cout << d << ", ";
	std::cout << std::endl;

	return 0;
}
