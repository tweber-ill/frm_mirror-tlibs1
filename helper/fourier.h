/*
 * FFT & DFT routines
 *
 * Author: Tobias Weber <tweber@frm2.tum.de>
 * Date: August 2012
 */

#ifndef __FOURIER__
#define __FOURIER__

#include <complex>
#define USE_FFTW


//------------------------------------------------------------------------------
// standard dft
// dft formulas from here:
// http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029
template<typename T>
std::complex<T> dft_coeff(int k,
					const T *pReal, const T *pImag,
					unsigned int n)
{
	std::complex<T> imag(0., 1.);

	std::complex<T> f(0.,0.);
	for(unsigned int j=0; j<n; ++j)
	{
		std::complex<T> t(pReal?pReal[j]:T(0), pImag?pImag[j]:T(0));

		T dv = -2.*M_PI*T(j)*T(k)/T(n);
		f += t * (cos(dv) + imag*sin(dv));
	}

	return f;
}

template<typename T>
void dft(const T *pRealIn, const T *pImagIn,
			   T *pRealOut, T *pImagOut, unsigned int n)
{
	for(unsigned int k=0; k<n; ++k)
	{
		std::complex<T> f = dft_coeff<T>(k, pRealIn, pImagIn, n);
		pRealOut[k] = f.real();
		pImagOut[k] = f.imag();
	}
}

template<typename T>
std::complex<T> idft_coeff(int k,
					const T *pReal, const T *pImag,
					unsigned int n)
{
	std::complex<T> imag(0., 1.);

	std::complex<T> t(0.,0.);
	for(unsigned int j=0; j<n; ++j)
	{
		std::complex<T> f(pReal?pReal[j]:T(0), pImag?pImag[j]:T(0));

		T dv = 2.*M_PI*T(j)*T(k)/T(n);
		t += f * (cos(dv) + imag*sin(dv));
	}

	return t;
}

template<typename T>
void idft(const T *pRealIn, const T *pImagIn,
				T *pRealOut, T *pImagOut, unsigned int n)
{
	for(unsigned int k=0; k<n; ++k)
	{
		std::complex<T> t = idft_coeff<T>(k, pRealIn, pImagIn, n);
		pRealOut[k] = t.real();
		pImagOut[k] = t.imag();
	}
}
//------------------------------------------------------------------------------


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
		void *m_pIn, *m_pOut;
		void *m_pPlan, *m_pPlan_inv;

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

#endif
