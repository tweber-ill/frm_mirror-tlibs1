/*
 * FFT & DFT routines
 *
 * @author Tobias Weber <tweber@frm2.tum.de>
 * @date August 2012
 * @license GPLv2 or GPLv3
 */

#include "fourier.h"
#include "math.h"
#include "../helper/log.h"

#include <iostream>
#include <cstring>
#include <memory>

namespace tl {

Fourier::Fourier(unsigned int iSize) : m_iSize(iSize), m_dft(m_iSize)
{}

Fourier::~Fourier()
{}

void Fourier::fft(const double *pRealIn, const double *pImagIn,
			double *pRealOut, double *pImagOut)
{
	m_dft.trafo(pRealIn, pImagIn, pRealOut, pImagOut, 0);
}

void Fourier::ifft(const double *pRealIn, const double *pImagIn,
			double *pRealOut, double *pImagOut)
{
	m_dft.trafo(pRealIn, pImagIn, pRealOut, pImagOut, 1);
}

}
