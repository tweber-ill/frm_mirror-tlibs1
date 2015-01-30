// gcc -o dft2 dft2.cpp ../math/fourier.cpp ../helper/log.cpp -lstdc++ -lm -std=c++11

#include "../math/fourier.h"
#include <iostream>

int main()
{
	typedef std::complex<double> t_c;
	std::vector<t_c> vecIn = {t_c(1,4), t_c(2,3), t_c(3,2), t_c(4,1)};

	std::cout << "in: ";
	for(const t_c& c : vecIn) std::cout << c << ", ";
	std::cout << std::endl;


	std::vector<t_c> vecOut = tl::dft_direct(vecIn, 0, 0);
	std::vector<t_c> vecOut2 = tl::fft_direct(vecIn);

	std::cout << "dft: ";
	for(const t_c& c : vecOut) std::cout << c << ", ";
	std::cout << std::endl;

	std::cout << "fft: ";
	for(const t_c& c : vecOut) std::cout << c << ", ";
	std::cout << std::endl;


	vecOut = tl::dft_double(vecOut);
	vecIn = tl::dft_direct(vecOut, 1, 1);

	std::cout << "idft: ";
	for(const t_c& c : vecIn) std::cout << c << ", ";
	std::cout << std::endl;

	return 0;
}
