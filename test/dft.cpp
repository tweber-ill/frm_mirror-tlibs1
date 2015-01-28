// gcc -DUSE_FFTW -o dft dft.cpp ../math/fourier.cpp ../helper/log.cpp -lstdc++ -lm -std=c++11 -lfftw3
// gcc -o dft dft.cpp ../math/fourier.cpp ../helper/log.cpp -lstdc++ -lm -std=c++11

#include "../math/fourier.h"
#include <iostream>

int main()
{
	double dInR[8] = {1., 2., 3., 4., 5., 6., 7., 8.};
	double dInI[8] = {8., 7., 6., 5., 4., 3., 2., 1.};
	double dOutR[8], dOutI[8];

	tl::Fourier fourier(8);
	fourier.fft(dInR, dInI, dOutR, dOutI);

	std::cout << "fft: ";
	for(int i=0; i<8; ++i)
	{
		std::cout << dOutR[i];
		if(!tl::float_equal(dOutI[i], 0.))
			std::cout << " + " << dOutI[i] << "*i";
		std::cout << ", ";
	}
	std::cout << std::endl;



	fourier.ifft(dOutR, dOutI, dInR, dInI);

	std::cout << "ifft: ";
	for(int i=0; i<8; ++i)
	{
		std::cout << dInR[i]/8.;
		if(!tl::float_equal(dInI[i], 0., 1e-7))
			std::cout << " + " << dInI[i]/8. << "*i";
		std::cout << ", ";
	}
	std::cout << std::endl;

	return 0;
}
