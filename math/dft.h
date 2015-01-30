/*
 * dft stuff
 * @author tweber
 * @date 2012, jan-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_DFT_H__
#define __TLIBS_DFT_H__

#include <complex>
#include <vector>
#include <memory>

#include "math.h"
#include "linalg.h"

namespace tl
{
//------------------------------------------------------------------------------
// standard dft
// dft formulas from here:
// http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029
template<typename T=double>
std::complex<T> dft_coeff(int k,
			const T *pReal, const T *pImag, std::size_t n,
			bool bInv=0)
{
	std::complex<T> imag(0., 1.);

	std::complex<T> f(0.,0.);
	for(std::size_t j=0; j<n; ++j)
	{
		std::complex<T> t(pReal?pReal[j]:T(0), pImag?pImag[j]:T(0));

		T dv = -2.*M_PI*T(j)*T(k)/T(n);
		if(bInv) dv = -dv;
		f += t * (cos(dv) + imag*sin(dv));
	}

	return f;
}

template<typename T=double>
void dft_direct(const T *pRealIn, const T *pImagIn,
		T *pRealOut, T *pImagOut, std::size_t n,
		bool bInv=0, bool bNorm=0)
{
	for(std::size_t k=0; k<n; ++k)
	{
		std::complex<T> f = dft_coeff<T>(k, pRealIn, pImagIn, n, bInv);
		pRealOut[k] = f.real();
		pImagOut[k] = f.imag();

		if(bNorm && bInv)
		{
			pRealOut[k] /= n;
			pImagOut[k] /= n;
		}
	}
}

template<typename T=double>
std::vector<std::complex<T>> dft_direct(const std::vector<std::complex<T>>& vecIn,
					bool bInv=0, bool bNorm=0)
{
	const std::size_t n = vecIn.size();
	std::vector<std::complex<T>> vecOut;
	vecOut.reserve(n);

	std::unique_ptr<T[]> in_real(new T[n]);
	std::unique_ptr<T[]> in_imag(new T[n]);
	std::unique_ptr<T[]> out_real(new T[n]);
	std::unique_ptr<T[]> out_imag(new T[n]);

	T* pInR = in_real.get();
	T* pInI = in_imag.get();
	T* pOutR = out_real.get();
	T* pOutI = out_imag.get();

	for(std::size_t i=0; i<n; ++i)
	{
		pInR[i] = vecIn[i].real();
		pInI[i] = vecIn[i].imag();
	}

	dft_direct(pInR, pInI, pOutR, pOutI, n, bInv, bNorm);

	for(std::size_t i=0; i<n; ++i)
		vecOut.push_back(std::complex<T>(pOutR[i], pOutI[i]));

	return vecOut;
}


//------------------------------------------------------------------------------

template<typename T=double>
std::vector<std::complex<T>> arrs_to_cvec(const T* pReal, const T* pImag, std::size_t N)
{
	std::vector<std::complex<T>> vec;
	vec.reserve(N);

	for(std::size_t n=0; n<N; ++n)
		vec.push_back(std::complex<T>(pReal?pReal[n]:0., pImag?pImag[n]:0.));

	return vec;
}

template<typename T=double>
void cvec_to_arrs(const std::vector<std::complex<T>>& vec, T* pReal, T* pImag)
{
	for(std::size_t n=0; n<vec.size(); ++n)
	{
		if(pReal) pReal[n] = vec[n].real();
		if(pImag) pImag[n] = vec[n].imag();
	}
}


//------------------------------------------------------------------------------


// tShift=1: shift 1 sample to the right
// tShift=-11: shift 1 sample to the left
template<typename T=double>
std::vector<std::complex<T>> dft_shift(const std::vector<std::complex<T>>& vecIn, T tShift=1.)
{
	const std::size_t N = vecIn.size();
	std::vector<std::complex<T>> vecOut;
	vecOut.reserve(N);

	for(std::size_t i=0; i<N; ++i)
	{
		const T c = std::cos(2.*M_PI*i*tShift / N);
		const T s = std::sin(2.*M_PI*i*tShift / N);
		std::complex<T> ph(c, -s);

		vecOut.push_back(vecIn[i]*ph);
	}

	return vecOut;
}

template<typename T=double>
std::vector<std::complex<T>> dft_double(const std::vector<std::complex<T>>& vecIn)
{
	const std::size_t N = vecIn.size();
	std::vector<std::complex<T>> vecOut;
	vecOut.reserve(N*2);

	for(std::size_t i=0; i<2*N; ++i)
		vecOut.push_back(vecIn[i%N]);

	return vecOut;
}


//------------------------------------------------------------------------------

template<typename T=unsigned int>
T count_bits(T imax)
{
	T inum = 0;
	for(; imax!=0; imax>>=1) ++inum;
	return inum;
}

template<typename T=unsigned int>
T bit_reverse(T imax, T inum)
{
	if(imax<2) return inum;

	T irev = 0;
	T ibitcnt = count_bits(imax)-2;

	for(T i=1; i<imax; i<<=1)
	{
		if(inum & i)
			irev |= (1 << ibitcnt);
		--ibitcnt;
	}

	return irev;
}

template<typename T=unsigned int>
std::vector<T> bit_reverse_indices(T imax)
{
	std::vector<T> vec;
	vec.reserve(imax);

	for(T i=0; i<imax; ++i)
		vec.push_back(bit_reverse(imax,i));

	return vec;
}

template<typename T=double>
std::complex<T> fft_factor(T N, T k)
{
	T c = std::cos(2.*M_PI*k/N);
	T s = std::sin(2.*M_PI*k/N);
	return std::complex<T>(c, -s);
}

template<typename T=double>
std::vector<std::complex<T>> fft_reorder(const std::vector<std::complex<T>>& vecIn)
{
	std::vector<std::size_t> vecIdx = bit_reverse_indices(vecIn.size());

	std::vector<std::complex<T>> vecInRev;
	vecInRev.reserve(vecIn.size());

	for(std::size_t i=0; i<vecIn.size(); ++i)
		vecInRev.push_back(vecIn[vecIdx[i]]);

	return vecInRev;
}

template<typename T=double>
std::vector<std::complex<T>> fft_twopoint(const std::vector<std::complex<T>>& vecIn)
{
	std::vector<std::complex<T>> vecOut;
	vecOut.reserve(vecIn.size());

	for(std::size_t i=0; i<vecIn.size(); i+=2)
	{
		vecOut.push_back(vecIn[i] + vecIn[i+1]);
		vecOut.push_back(vecIn[i] - vecIn[i+1]);
	}

	return vecOut;
}

template<typename T=double>
std::vector<std::complex<T>> fft_merge(const std::vector<std::complex<T>>& vecIn)
{
	const std::size_t N = vecIn.size();
	const std::size_t N2 = N/2;

	std::vector<std::complex<T>> vecOut;
	vecOut.resize(N);

	vecOut[0] = vecIn[0] + vecIn[N2+0]*fft_factor<T>(N,0);
	vecOut[1] = vecIn[1] + vecIn[N2+1]*fft_factor<T>(N,1);
	vecOut[N2+0] = vecIn[0] + vecIn[N2+0]*fft_factor<T>(N,N2+0);
	vecOut[N2+1] = vecIn[1] + vecIn[N2+1]*fft_factor<T>(N,N2+1);

	// TODO

	return vecOut;
}

template<typename T=double>
std::vector<std::complex<T>> fft_direct(const std::vector<std::complex<T>>& vecIn)
{
	const std::size_t n = vecIn.size();
	std::vector<std::complex<T>> vecOut;
	vecOut.resize(n);

	if(n==0);
	else if(n==1)
	{
		vecOut[0] = vecIn[0];
	}
	else if(n==2)
	{
		vecOut[0] = vecIn[0] + vecIn[1];
		vecOut[1] = vecIn[0] - vecIn[1];
	}
	else
	{
		vecOut = fft_reorder(vecIn);
		vecOut = fft_twopoint(vecOut);
		vecOut = fft_merge(vecOut);
	}

	return vecOut;
}


//------------------------------------------------------------------------------


// dft with pre-calculated coefficients
template<class T=double>
class DFT
{
	protected:
		ublas::matrix<std::complex<T>> m_matCoeff;
		ublas::matrix<std::complex<T>> m_matCoeffInv;

	protected:
		void InitCoeffMatrices(std::size_t n)
		{
			m_matCoeff.resize(n,n,0);

			for(std::size_t i=0; i<n; ++i)
				for(std::size_t j=0; j<n; ++j)
				{
					const T c = std::cos(2.*M_PI*i*j / n);
					const T s = std::sin(2.*M_PI*i*j / n);
					m_matCoeff(i,j) = std::complex<T>(c, -s);
				}

				m_matCoeffInv = ublas::conj(m_matCoeff);
		}

	public:
		DFT(std::size_t n)
		{
			InitCoeffMatrices(n);
		}

		virtual ~DFT() = default;

		ublas::vector<std::complex<T>> trafo(const ublas::vector<std::complex<T>>& vec, bool bInv=0)
		{
			if(!bInv)
				return ublas::prod(m_matCoeff, vec);
			else
				return ublas::prod(m_matCoeffInv, vec);
		}

		void trafo(const double* pInR, const double *pInI, 
					double *pOutR, double *pOutI, 
					bool bInv=0)
		{
			const std::size_t iSize = m_matCoeff.size1();

			typedef ublas::vector<std::complex<T>> t_vec;
			t_vec vecIn(iSize), vecOut;

			for(std::size_t i=0; i<iSize; ++i)
			{
				T dR = pInR ? pInR[i] : 0.;
				T dI = pInI ? pInI[i] : 0.;
				vecIn[i] = std::complex<T>(dR, dI);
			}

			vecOut = trafo(vecIn, bInv);

			for(std::size_t i=0; i<iSize; ++i)
			{
				if(pOutR) pOutR[i] = vecOut[i].real();
				if(pOutI) pOutI[i] = vecOut[i].imag();
			}
		}
};

}

#endif
