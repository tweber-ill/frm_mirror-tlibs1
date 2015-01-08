/*
 * Minuit interface
 *
 * Author: Tobias Weber
 * Date: April 2012
 * 
 * general fitter structure (i.e. function => chi^2 calculation => calling
 * minuit) originally based on the examples in the Minuit user's guide:
 * http://seal.cern.ch/documents/minuit/mnusersguide.pdf
 *
 */

#include "minuit.h"

#include <limits>
#include <cmath>
#include <algorithm>
#include <sstream>

#include "../helper/misc.h"

std::ostream& operator<<(std::ostream& ostr, const MinuitFuncModel& fkt)
{
	ostr << fkt.print();
	return ostr;
}

std::ostream& operator<<(std::ostream& ostr, const MinuitFuncModel_nd& fkt)
{
	ostr << fkt.print();
	return ostr;
}

// chi^2 calculation
// based on the example in the Minuit user's guide:
// http://seal.cern.ch/documents/minuit/mnusersguide.pdf
double Chi2Function::chi2(const std::vector<double>& vecParams) const
{
	// cannot operate on m_pfkt directly, because Minuit
	// uses more than one thread!
	MinuitFuncModel* pfkt = m_pfkt->copy();
	autodeleter<MinuitFuncModel> a(pfkt);

	bool bParamsOk = pfkt->SetParams(vecParams);
	//if(!bParamsOk)
	//	return std::numeric_limits<double>::max();

	double dChi2 = 0.;
	for(unsigned int i=0; i<m_uiLen; ++i)
	{
		double d = (*pfkt)(m_px[i]) - m_py[i];
		double dy = m_pdy ? m_pdy[i] : 0.1*d;	// assume 10% error if none given
		if(fabs(dy) < std::numeric_limits<double>::min())
			dy = std::numeric_limits<double>::min();

		d /= dy;
		dChi2 += d*d;
	}
	return dChi2;
}

double Chi2Function_nd::chi2(const std::vector<double>& vecParams) const
{
	MinuitFuncModel_nd* pfkt = m_pfkt->copy();
	autodeleter<MinuitFuncModel_nd> _a0(pfkt);

	bool bParamsOk = pfkt->SetParams(vecParams);

	double *px = new double[m_uiDim];
	autodeleter<double> _a1(px, 1);

	double dChi2 = 0.;
	for(unsigned int i=0; i<m_uiLen; ++i)
	{
		for(unsigned int iX=0; iX<m_uiDim; ++iX)
			px[iX] = m_vecpx[iX][i];

		double d = (*pfkt)(px) - m_py[i];
		double dy = m_pdy ? m_pdy[i] : 0.1*d;	// assume 10% error if none given
		if(fabs(dy) < std::numeric_limits<double>::min())
			dy = std::numeric_limits<double>::min();

		d /= dy;
		dChi2 += d*d;
	}
	return dChi2;
}
