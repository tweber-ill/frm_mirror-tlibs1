/*
 * Wrapper to choose between FFT and FFTw
 *
 * @author tweber
 * @date jan-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __FOURIER_CHOOSER__
#define __FOURIER_CHOOSER__

#ifdef USE_FFTW
	#include "fftw.h"
#else
	#include "dft.h"
#endif


namespace tl {

class Fourier
{
	protected:
		unsigned int m_iSize;

#ifdef USE_FFTW
		FFTw m_dft;
#else
		FFT<double> m_dft;
#endif


	public:
		Fourier(unsigned int iSize);
		virtual ~Fourier();

		void fft(const double* pRealIn, const double *pImagIn,
				double *pRealOut, double *pImagOut);
		void ifft(const double* pRealIn, const double *pImagIn,
				double *pRealOut, double *pImagOut);
};
//------------------------------------------------------------------------------

}

#endif
