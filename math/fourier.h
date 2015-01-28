/*
 * FFT & DFT routines
 *
 * @author Tobias Weber <tweber@frm2.tum.de>
 * @date August 2012
 * @license GPLv2 or GPLv3
 */

#ifndef __FOURIER__
#define __FOURIER__

#include "dft.h"

namespace tl {


//------------------------------------------------------------------------------
// algorithms
// perform a zero-order phase correction;
// for a description (albeit in the context of NMR) see e.g. here:
// http://www-keeler.ch.cam.ac.uk/lectures/Irvine/chapter4.pdf
template<typename T>
std::complex<T> phase_correction_0(const std::complex<T>& c, T dPhase)
{
	return c * std::complex<T>(cos(-dPhase), sin(-dPhase));
}

// perform a first-order phase correction:
// dPhase = dPhaseOffs + x*dPhaseSlope
template<typename T>
std::complex<T> phase_correction_1(const std::complex<T>& c,
				T dPhaseOffs, T dPhaseSlope, T x)
{
	return phase_correction_0<T>(c, dPhaseOffs + x*dPhaseSlope);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// fft using fftw
class Fourier
{
	protected:
		unsigned int m_iSize;
		
#ifdef USE_FFTW
		void *m_pIn, *m_pOut;
		void *m_pPlan, *m_pPlan_inv;
#else
		DFT<double> m_dft;
#endif


	public:
		Fourier(unsigned int iSize);
		virtual ~Fourier();

		bool fft(const double* pRealIn, const double *pImagIn,
									double *pRealOut, double *pImagOut);
		bool ifft(const double* pRealIn, const double *pImagIn,
									double *pRealOut, double *pImagOut);

		// shift a sine given in pDatIn by dPhase
		// we cannot phase shift a sine directly in the time domain due to
		// binning constraints; but in the frequency domain the phase is
		// a continuous variable which we can arbitrarily change and then
		// transform the data back into the time domain to get a shifted rebinning
		bool shift_sin(double dNumOsc, const double* pDatIn,
						double *pDataOut, double dPhase);

		bool phase_correction_0(const double* pDatIn, double *pDataOut,
								double dPhase);

		bool phase_correction_1(const double* pDatIn, double *pDataOut,
								double dPhaseOffs, double dPhaseSlope);

		bool get_contrast(double dNumOsc, const double* pDatIn,
						  double& dC, double& dPh);
};
//------------------------------------------------------------------------------

}

#endif
