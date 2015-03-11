/*
 * Minuit interface
 *
 * @author Tobias Weber
 * @date April 2012
 * @license GPLv2 or GPLv3
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

namespace tl {

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
	std::unique_ptr<MinuitFuncModel> uptrFkt(m_pfkt->copy());
	MinuitFuncModel* pfkt = uptrFkt.get();

	/*bool bParamsOk = */pfkt->SetParams(vecParams);
	//if(!bParamsOk)
	//	return std::numeric_limits<double>::max();

	return tl::chi2<double, decltype(*pfkt)>(*pfkt, m_uiLen, m_px, m_py, m_pdy);
}

double Chi2Function_nd::chi2(const std::vector<double>& vecParams) const
{
	std::unique_ptr<MinuitFuncModel_nd> uptrFkt(m_pfkt->copy());
	MinuitFuncModel_nd* pfkt = uptrFkt.get();

	/*bool bParamsOk = */pfkt->SetParams(vecParams);

	std::unique_ptr<double[]> uptrX(new double[m_uiDim]);
	double *px = uptrX.get();

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

}
