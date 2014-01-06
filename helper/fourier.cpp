/*
 * FFT & DFT routines
 * 
 * Author: Tobias Weber <tweber@frm2.tum.de>
 * Date: August 2012
 */

#include "fourier.h"
#include "misc.h"
#include "math.h"
#include <iostream>
#include <string.h>

#ifdef USE_FFTW
	#include <fftw3.h>
#endif

#ifdef USE_CUDA
	#undef USE_FFTW
	#include <cuda.h>
	#include <cuda_runtime.h>
	#include <cufft.h>
#endif

//------------------------------------------------------------------------------
// fft using fftw
#ifdef USE_FFTW

Fourier::Fourier(unsigned int iSize) : m_iSize(iSize)
{
	m_pIn = fftw_malloc(m_iSize * sizeof(fftw_complex));
	m_pOut = fftw_malloc(m_iSize * sizeof(fftw_complex));

	m_pPlan = malloc(sizeof(fftw_plan));
	m_pPlan_inv = malloc(sizeof(fftw_plan));

	fftw_plan* pPlan = (fftw_plan*) m_pPlan;
	fftw_plan* pPlan_inv = (fftw_plan*) m_pPlan_inv;

	*pPlan = fftw_plan_dft_1d(iSize,
								(fftw_complex*)m_pIn, (fftw_complex*)m_pOut,
								FFTW_FORWARD, FFTW_MEASURE);

	*pPlan_inv = fftw_plan_dft_1d(iSize,
								(fftw_complex*)m_pIn, (fftw_complex*)m_pOut,
								FFTW_BACKWARD, FFTW_MEASURE);

	if(!*pPlan || !*pPlan_inv)
		std::cerr << "Fourier: Could not create plan." << std::endl;
}

Fourier::~Fourier()
{
	fftw_destroy_plan(*(fftw_plan*)m_pPlan);
	fftw_destroy_plan(*(fftw_plan*)m_pPlan_inv);

	free(m_pPlan);
	free(m_pPlan_inv);

	fftw_free(m_pIn);
	fftw_free(m_pOut);
}

bool Fourier::fft(const double* pRealIn, const double *pImagIn,
					double *pRealOut, double *pImagOut)
{
	fftw_plan* pPlan = (fftw_plan*) m_pPlan;

	fftw_complex* pIn = (fftw_complex*)m_pIn;
	fftw_complex* pOut = (fftw_complex*)m_pOut;

	for(unsigned int i=0; i<m_iSize; ++i)
	{
		pIn[i][0] = pRealIn ? pRealIn[i] : 0.;
		pIn[i][1] = pImagIn ? pImagIn[i] : 0.;

		pOut[i][0] = 0.;
		pOut[i][1] = 0.;
	}

	fftw_execute(*pPlan);

	for(unsigned int i=0; i<m_iSize; ++i)
	{
		if(pRealOut) pRealOut[i] = pOut[i][0];
		if(pImagOut) pImagOut[i] = pOut[i][1];
	}

	return true;
}

bool Fourier::ifft(const double* pRealIn, const double *pImagIn,
					double *pRealOut, double *pImagOut)
{
	fftw_plan* pPlan_inv = (fftw_plan*) m_pPlan_inv;

	fftw_complex* pIn = (fftw_complex*)m_pIn;
	fftw_complex* pOut = (fftw_complex*)m_pOut;

	for(unsigned int i=0; i<m_iSize; ++i)
	{
		pIn[i][0] = pRealIn ? pRealIn[i] : 0.;
		pIn[i][1] = pImagIn ? pImagIn[i] : 0.;

		pOut[i][0] = 0.;
		pOut[i][1] = 0.;
	}

	fftw_execute(*pPlan_inv);

	for(unsigned int i=0; i<m_iSize; ++i)
	{
		if(pRealOut) pRealOut[i] = pOut[i][0];
		if(pImagOut) pImagOut[i] = pOut[i][1];
	}

	return true;
}

#else
#ifdef USE_CUDA


Fourier::Fourier(unsigned int iSize) : m_iSize(iSize)
{
	cudaError_t errIn = cudaMalloc((void**)&m_pIn, m_iSize*sizeof(cufftComplex));
	cudaError_t errOut = cudaMalloc((void**)&m_pOut, m_iSize*sizeof(cufftComplex));

	if(errIn==cudaErrorMemoryAllocation || errOut==cudaErrorMemoryAllocation)
	{
		std::cerr << "Error: Cuda could not allocate memory." << std::endl;
		return;
	}

	m_pTmp = malloc(m_iSize*sizeof(cufftComplex));

	m_pPlan = malloc(sizeof(cufftHandle));
	m_pPlan_inv = malloc(sizeof(cufftHandle));

	cufftPlan1d((cufftHandle*)m_pPlan, m_iSize, CUFFT_C2C, 1);
	cufftPlan1d((cufftHandle*)m_pPlan_inv, m_iSize, CUFFT_C2C, 1);

	if(!*(cufftHandle*)m_pPlan || !*(cufftHandle*)m_pPlan_inv)
		std::cerr << "Fourier: Could not create plan." << std::endl;
}

Fourier::~Fourier()
{
	cufftDestroy(*(cufftHandle*)m_pPlan);
	cufftDestroy(*(cufftHandle*)m_pPlan_inv);

	free(m_pTmp);

	free(m_pPlan);
	free(m_pPlan_inv);

	cudaFree(m_pIn);
	cudaFree(m_pOut);
}

bool Fourier::fft(const double* pRealIn, const double *pImagIn,
					double *pRealOut, double *pImagOut)
{
	cufftHandle* pPlan = (cufftHandle*)m_pPlan;
	cufftComplex* pTmp = (cufftComplex*)m_pTmp;

	for(unsigned int i=0; i<m_iSize; ++i)
	{
		pTmp[i].x = pRealIn ? pRealIn[i] : 0.;
		pTmp[i].y = pImagIn ? pImagIn[i] : 0.;
	}

	cudaMemcpy(m_pIn, pTmp, m_iSize*sizeof(cufftComplex), cudaMemcpyHostToDevice);
	cufftExecC2C(*pPlan, (cufftComplex*)m_pIn, (cufftComplex*)m_pOut, CUFFT_FORWARD);
	cudaMemcpy(pTmp, m_pOut, m_iSize*sizeof(cufftComplex), cudaMemcpyDeviceToHost);

	for(unsigned int i=0; i<m_iSize; ++i)
	{
		if(pRealOut) pRealOut[i] = pTmp[i].x;
		if(pImagOut) pImagOut[i] = pTmp[i].y;
	}

	return true;
}

bool Fourier::ifft(const double* pRealIn, const double *pImagIn,
					double *pRealOut, double *pImagOut)
{
	cufftHandle* pPlan = (cufftHandle*)m_pPlan;
	cufftComplex* pTmp = (cufftComplex*)m_pTmp;

	for(unsigned int i=0; i<m_iSize; ++i)
	{
		pTmp[i].x = pRealIn ? pRealIn[i] : 0.;
		pTmp[i].y = pImagIn ? pImagIn[i] : 0.;
	}

	cudaMemcpy(m_pIn, pTmp, m_iSize*sizeof(cufftComplex), cudaMemcpyHostToDevice);
	cufftExecC2C(*pPlan, (cufftComplex*)m_pIn, (cufftComplex*)m_pOut, CUFFT_INVERSE);
	cudaMemcpy(pTmp, m_pOut, m_iSize*sizeof(cufftComplex), cudaMemcpyDeviceToHost);

	for(unsigned int i=0; i<m_iSize; ++i)
	{
		if(pRealOut) pRealOut[i] = pTmp[i].x;
		if(pImagOut) pImagOut[i] = pTmp[i].y;
	}

	return true;
}


#else

Fourier::Fourier(unsigned int iSize) : m_iSize(iSize)
{}

Fourier::~Fourier()
{}

bool Fourier::fft(const double *pRealIn, const double *pImagIn,
					double *pRealOut, double *pImagOut)
{
	// if we don't have the fftw available, use dft instead
	dft<double>(pRealIn, pImagIn, pRealOut, pImagOut, m_iSize);

	return true;
}

bool Fourier::ifft(const double *pRealIn, const double *pImagIn,
					double *pRealOut, double *pImagOut)
{
	// if we don't have the fftw available, use dft instead
	idft<double>(pRealIn, pImagIn, pRealOut, pImagOut, m_iSize);

	return true;	
}

#endif	// USE_CUDA
#endif	// USE_FFTW


bool Fourier::shift_sin(double dNumOsc, const double* pDatIn,
				double *pDataOut, double dPhase)
{
	unsigned int iSize = m_iSize;
	const int iNumOsc = int(dNumOsc);
	dNumOsc = double(iNumOsc);			// consider only full oscillations

	//double dShiftSamples = dPhase/(2.*M_PI) * dSize;

	double *pZero = new double[iSize];
	memset(pZero, 0, sizeof(double)*iSize);

	double *pDatFFT_real = new double[iSize];
	double *pDatFFT_imag = new double[iSize];

	autodeleter<double> _a0(pZero, true);
	autodeleter<double> _a1(pDatFFT_real, true);
	autodeleter<double> _a2(pDatFFT_imag, true);

	if(!fft(pDatIn, pZero, pDatFFT_real, pDatFFT_imag))
		return false;

	// filter out everything not concerning the sine with iNumOsc oscillations
	for(unsigned int i=0; i<iSize; ++i)
	{
		if(int(i)==iNumOsc || i==0)
			continue;

		pDatFFT_real[i] = 0.;
		pDatFFT_imag[i] = 0.;
	}

	// amp & phase
	std::complex<double> c(pDatFFT_real[iNumOsc], pDatFFT_imag[iNumOsc]);

	// since the signal is real we can take the first half of the fft data
	// and multiply by two.
	c *= 2.;
	c = ::phase_correction_0<double>(c, dPhase);

	pDatFFT_real[iNumOsc] = c.real();
	pDatFFT_imag[iNumOsc] = c.imag();

	// offset
	pDatFFT_imag[0] = 0.;


	if(!ifft(pDatFFT_real, pDatFFT_imag, pDataOut, pZero))
		return false;


	// normalization
	for(unsigned int i=0; i<iSize; ++i)
		pDataOut[i] /= double(iSize);

	//save_dat("in.dat", pDatIn, iSize);
	//save_dat("out.dat", pDataOut, iSize);
	return true;
}

bool Fourier::phase_correction_0(const double* pDatIn, double *pDataOut,
								double dPhase)
{
	unsigned int iSize = m_iSize;

	double *pZero = new double[iSize];
	memset(pZero, 0, sizeof(double)*iSize);

	double *pDatFFT_real = new double[iSize];
	double *pDatFFT_imag = new double[iSize];

	autodeleter<double> _a0(pZero, true);
	autodeleter<double> _a1(pDatFFT_real, true);
	autodeleter<double> _a2(pDatFFT_imag, true);

	if(!fft(pDatIn, pZero, pDatFFT_real, pDatFFT_imag))
		return false;

	for(unsigned int i=1; i<iSize; ++i)
	{
		std::complex<double> c(pDatFFT_real[i], pDatFFT_imag[i]);
		if(i<iSize/2)
		{
			c *= 2.;
			c = ::phase_correction_0<double>(c, dPhase*double(i));
		}
		else
		{
			// not needed in real input data
			c = std::complex<double>(0., 0.);
		}

		pDatFFT_real[i] = c.real();
		pDatFFT_imag[i] = c.imag();
	}

	if(!ifft(pDatFFT_real, pDatFFT_imag, pDataOut, pZero))
		return false;

	// normalization
	for(unsigned int i=0; i<iSize; ++i)
		pDataOut[i] /= double(iSize);

	//save_dat("in.dat", pDatIn, iSize);
	//save_dat("out.dat", pDataOut, iSize);
	return true;
}

bool Fourier::phase_correction_1(const double* pDatIn,
				double *pDataOut, double dPhaseOffs, double dPhaseSlope)
{
	unsigned int iSize = m_iSize;

	double *pZero = new double[iSize];
	memset(pZero, 0, sizeof(double)*iSize);

	double *pDatFFT_real = new double[iSize];
	double *pDatFFT_imag = new double[iSize];

	autodeleter<double> _a0(pZero, true);
	autodeleter<double> _a1(pDatFFT_real, true);
	autodeleter<double> _a2(pDatFFT_imag, true);

	if(!fft(pDatIn, pZero, pDatFFT_real, pDatFFT_imag))
		return false;

	for(unsigned int i=1; i<iSize; ++i)
	{	
		std::complex<double> c(pDatFFT_real[i], pDatFFT_imag[i]);
		if(i<iSize/2)
		{
			double dX = double(i)/double(iSize);
			
			c *= 2.;
			c = ::phase_correction_1<double>(c,
							dPhaseOffs*double(i), dPhaseSlope*double(i),
							dX);
		}
		else
		{
			// not needed in real input data
			c = std::complex<double>(0., 0.);
		}

		pDatFFT_real[i] = c.real();
		pDatFFT_imag[i] = c.imag();
	}

	if(!ifft(pDatFFT_real, pDatFFT_imag, pDataOut, pZero))
		return false;

	// normalization
	for(unsigned int i=0; i<iSize; ++i)
		pDataOut[i] /= double(iSize);

	//save_dat("in.dat", pDatIn, iSize);
	//save_dat("out.dat", pDataOut, iSize);
	return true;
}


bool Fourier::get_contrast(double dNumOsc, const double* pDatIn,
						   double& dC, double& dPh)
{
	unsigned int iSize = m_iSize;
	const int iNumOsc = int(dNumOsc);
	dNumOsc = double(iNumOsc);			// consider only full oscillations

	double *pZero = new double[iSize];
	memset(pZero, 0, sizeof(double)*iSize);

	double *pDatFFT_real = new double[iSize];
	double *pDatFFT_imag = new double[iSize];

	autodeleter<double> _a0(pZero, true);
	autodeleter<double> _a1(pDatFFT_real, true);
	autodeleter<double> _a2(pDatFFT_imag, true);

	if(!fft(pDatIn, pZero, pDatFFT_real, pDatFFT_imag))
		return false;

	double dReal = 2.*pDatFFT_real[iNumOsc]/double(m_iSize);
	double dImag = 2.*pDatFFT_imag[iNumOsc]/double(m_iSize);

	double dAmp = sqrt(dReal*dReal + dImag*dImag);
	double dOffs = pDatFFT_real[0] / double(m_iSize);

	//std::cout << "amp = " << dAmp << std::endl;
	//std::cout << "offs = " << dOffs << std::endl;

	dC = dAmp/dOffs;
	dPh = atan2(dImag, dReal) + M_PI/2.;

	// half a bin
	dPh -= 2.*M_PI*dNumOsc*0.5/double(m_iSize);

	// half of first bin
	//dPh += - 0.5/dNumOsc * 2.*M_PI;

	if(dPh<0.)
		dPh += 2.*M_PI;

	return true;
}
