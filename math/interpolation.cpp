/*
 * data point interpolation
 *
 * @author: Tobias Weber
 * @date: 25-04-2013
 */

#include "interpolation.h"
//#include <boost/algorithm/minmax_element.hpp>
#include <limits>

namespace ublas = boost::numeric::ublas;


Bezier::Bezier(unsigned int N, const double *px, const double *py)
	: m_pvecs(0), m_iN(N)
{
	m_pvecs = new ublas::vector<double>[m_iN];

	for(unsigned int i=0; i<N; ++i)
	{
		m_pvecs[i].resize(2);
		m_pvecs[i][0] = px[i];
		m_pvecs[i][1] = py[i];
	}

	//auto MinMax = boost::minmax_element(px, px+N);
	//m_dMin = *MinMax.first;
	//m_dMax = *MinMax.second;
}

Bezier::~Bezier()
{
	if(m_pvecs) delete[] m_pvecs;
}

ublas::vector<double> Bezier::operator()(double t) const
{
	return ::bezier<double>(m_pvecs, m_iN, t);
}




BSpline::BSpline(unsigned int N, const double *px, const double *py, unsigned int iDegree)
	: m_pvecs(0), m_iN(N), m_iDegree(iDegree)
{
	m_pvecs = new ublas::vector<double>[m_iN];

	for(unsigned int i=0; i<m_iN; ++i)
	{
		m_pvecs[i].resize(2);
		m_pvecs[i][0] = px[i];
		m_pvecs[i][1] = py[i];
	}

	unsigned int iM = m_iDegree + m_iN + 1;
	m_vecKnots.resize(iM);

	const double eps = std::numeric_limits<double>::epsilon();

	// set knots to uniform, nonperiodic B-Spline
	for(unsigned int i=0; i<m_iDegree+1; ++i)
		m_vecKnots[i] = 0.+i*eps;
	for(unsigned int i=iM-m_iDegree-1; i<iM; ++i)
		m_vecKnots[i] = 1.-i*eps;
	for(unsigned int i=m_iDegree+1; i<iM-m_iDegree-1; ++i)
		m_vecKnots[i] = double(i+1-m_iDegree-1) / double(iM-2*m_iDegree-2 + 1);

	//for(unsigned int i=0; i<iM; ++i)
	//	m_vecKnots[i] = double(i) / double(iM-1);

	/*std::cout << "knots: ";
	for(double d: m_vecKnots)
		std::cout << d << " ";
	std::cout << std::endl;*/
}

BSpline::~BSpline()
{
	if(m_pvecs) delete[] m_pvecs;
}

boost::numeric::ublas::vector<double> BSpline::operator()(double t) const
{
	if(m_iN==0)
	{
		boost::numeric::ublas::vector<double> vecNull(2);
		vecNull[0] = vecNull[1] = 0.;
		return vecNull;
	}

	boost::numeric::ublas::vector<double> vec =  ::bspline<double>(m_pvecs, m_iN, t, m_vecKnots);

	// remove epsilon dependence
	if(t<=0.) vec = m_pvecs[0];
	if(t>=1.) vec = m_pvecs[m_iN-1];

	return vec;
}
