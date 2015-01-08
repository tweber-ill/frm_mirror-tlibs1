/*
 * data point interpolation
 *
 * @author: Tobias Weber
 * @date: 25-04-2013
 */

#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

#include <cmath>
#include <boost/math/special_functions/binomial.hpp>
#include <vector>
#include <algorithm>

#include "math.h"
#include "geo.h"
#include "misc.h"
#include "log.h"
#include "funcmod.h"


// see:
// http://mathworld.wolfram.com/BernsteinPolynomial.html
template<typename T> T bernstein(int i, int n, T t)
{
	T bino = boost::math::binomial_coefficient<T>(n, i);
	return bino * pow(t, i) * pow(1-t, n-i);
}

// see:
// http://mathworld.wolfram.com/BezierCurve.html
template<typename T>
boost::numeric::ublas::vector<T> bezier(const boost::numeric::ublas::vector<T>* P, unsigned int N, T t)
{
	if(N==0) return boost::numeric::ublas::vector<T>(0);
	const int n = N-1;

	boost::numeric::ublas::vector<T> vec(P[0].size());
	for(unsigned int i=0; i<vec.size(); ++i) vec[i] = T(0);

	for(int i=0; i<=n; ++i)
		vec += P[i]*bernstein(i, n, t);

	return vec;
}


// see:
// http://mathworld.wolfram.com/B-Spline.html
template<typename T> T bspline_base(int i, int j, T t, const std::vector<T>& knots)
{
	if(j==0)
	{
		if((knots[i] <= t) && (t < knots[i+1]) && (knots[i]<knots[i+1]))
			return 1.;
		return 0.;
	}

	T val11 = (t - knots[i]) / (knots[i+j]-knots[i]);
	T val12 = bspline_base(i, j-1, t, knots);
	T val1 = val11 * val12;

	T val21 = (knots[i+j+1]-t) / (knots[i+j+1]-knots[i+1]);
	T val22 = bspline_base(i+1, j-1, t, knots);
	T val2 = val21 * val22;

	T val = val1 + val2;
	return val;
}


// see:
// http://mathworld.wolfram.com/B-Spline.html
template<typename T>
boost::numeric::ublas::vector<T> bspline(const boost::numeric::ublas::vector<T>* P, unsigned int N, T t, const std::vector<T>& knots)
{
	if(N==0) return boost::numeric::ublas::vector<T>(0);
	const int n = N-1;
	const int m = knots.size()-1;
	const int degree = m-n-1;

	boost::numeric::ublas::vector<T> vec(P[0].size());
	for(unsigned int i=0; i<vec.size(); ++i) vec[i] = T(0);

	for(int i=0; i<=n; ++i)
		vec += P[i]*bspline_base(i, degree, t, knots);

	return vec;
}


class Bezier : public FunctionModel_param
{
	protected:
		boost::numeric::ublas::vector<double> *m_pvecs;
		unsigned int m_iN;

	public:
		Bezier(unsigned int N, const double *px, const double *py);
		virtual ~Bezier();

		virtual boost::numeric::ublas::vector<double> operator()(double t) const;

		virtual const char* GetModelName() const { return "bezier"; };
};


class BSpline : public FunctionModel_param
{
	protected:
		boost::numeric::ublas::vector<double> *m_pvecs;
		unsigned int m_iN, m_iDegree;
		std::vector<double> m_vecKnots;

	public:
		BSpline(unsigned int N, const double *px, const double *py, unsigned int iDegree=3);
		virtual ~BSpline();

		virtual boost::numeric::ublas::vector<double> operator()(double t) const;

		virtual const char* GetModelName() const { return "bspline"; };
};


template<typename T>
void find_peaks(unsigned int iLen, const T* px, const T* py, unsigned int iOrder,
						std::vector<T>& vecMaximaX, std::vector<T>& vecMaximaSize, std::vector<T>& vecMaximaWidth)
{
	BSpline spline(iLen, px, py, iOrder);
	const unsigned int iNumSpline = 512;    // TODO: in config

	T *pSplineX = new T[iNumSpline];
	T *pSplineY = new T[iNumSpline];
	T *pSplineDiff = new T[iNumSpline];
	T *pSplineDiff2 = new T[iNumSpline];

	const double* pdyMin = std::min_element(py, py+iLen);

	for(unsigned int iSpline=0; iSpline<iNumSpline; ++iSpline)
	{
		const T dT = T(iSpline) / T(iNumSpline-1);
		boost::numeric::ublas::vector<T> vec = spline(dT);

		pSplineX[iSpline] = vec[0];
		pSplineY[iSpline] = vec[1];
	}

	::diff(iNumSpline, pSplineX, pSplineY, pSplineDiff);
	::diff(iNumSpline, pSplineX, pSplineDiff, pSplineDiff2);
	std::vector<unsigned int> vecZeroes = ::find_zeroes<T>(iNumSpline, pSplineDiff);


	for(unsigned int iZeroIdx = 0; iZeroIdx<vecZeroes.size(); ++iZeroIdx)
	{
		const unsigned int iZero = vecZeroes[iZeroIdx];

		// minima / saddle points
		if(pSplineDiff2[iZero] >= 0.)
			continue;

		vecMaximaX.push_back(pSplineX[iZero]);

		int iMinIdxLeft = -1;
		int iMinIdxRight = -1;
		if(iZeroIdx > 0)
			iMinIdxLeft = vecZeroes[iZeroIdx-1];
		if(iZeroIdx+1 < vecZeroes.size())
			iMinIdxRight = vecZeroes[iZeroIdx+1];

		T dHeight = 0.;
		T dWidth = 0.;
		T dDiv = 0.;

		// minimum left of the peak
		if(iMinIdxLeft>=0)
		{
			dHeight += (pSplineY[iZero]-pSplineY[iMinIdxLeft]);
			dWidth += fabs((pSplineX[iZero]-pSplineX[iMinIdxLeft]));
			dDiv += 1.;
		}

		// minimum right of the peak
		if(iMinIdxRight>=0)
		{
			dHeight += (pSplineY[iZero]-pSplineY[iMinIdxRight]);
			dWidth += fabs((pSplineX[iZero]-pSplineX[iMinIdxRight]));
			dDiv += 1.;
		}

		// no adjacent minima...
		if(iMinIdxLeft<0 && iMinIdxRight<0)
		{
			dHeight = pSplineY[iZero]- *pdyMin;
			dWidth = (px[iLen-1] - px[0]) / 10.;	// guess something...
			dDiv = 1.;
		}

		if(dDiv != 0.)
		{
			dHeight /= dDiv;
			dWidth /= dDiv;
		}

		vecMaximaSize.push_back(dHeight);
		vecMaximaWidth.push_back(dWidth);
	}

	::sort_3<std::vector<double>::iterator>(vecMaximaSize.begin(), vecMaximaSize.end(), vecMaximaWidth.begin(), vecMaximaX.begin());
	std::reverse(vecMaximaSize.begin(), vecMaximaSize.end());
	std::reverse(vecMaximaWidth.begin(), vecMaximaWidth.end());
	std::reverse(vecMaximaX.begin(), vecMaximaX.end());



	std::ostringstream ostrDbg;
	ostrDbg << "Prefitter found peaks at: ";
	for(double dValX : vecMaximaX)
		ostrDbg << dValX << ", ";
	log_debug(ostrDbg.str());


	delete[] pSplineX;
	delete[] pSplineY;
	delete[] pSplineDiff;
	delete[] pSplineDiff2;
}

#endif
