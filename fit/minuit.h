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
#ifndef __MINUIT_IFACE_H__
#define __MINUIT_IFACE_H__

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnFcn.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <type_traits>
#include <limits>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>

#include "funcmod.h"
#include "../helper/misc.h"
#include "../log/log.h"


namespace tl {
//using t_real_min = double;
using t_real_min = typename std::result_of<
	decltype(&ROOT::Minuit2::MnFcn::Up)(ROOT::Minuit2::MnFcn) >::type;

class MinuitFuncModel : public FunctionModel<t_real_min>
{
public:
	virtual ~MinuitFuncModel() {}

	virtual bool SetParams(const std::vector<t_real_min>& vecParams) = 0;
	virtual t_real_min operator()(t_real_min x) const = 0;

	virtual MinuitFuncModel* copy() const = 0;
	virtual std::string print(bool bFillInSyms=true) const = 0;

	virtual const char* GetModelName() const = 0;
	virtual std::vector<std::string> GetParamNames() const = 0;
	virtual std::vector<t_real_min> GetParamValues() const = 0;
	virtual std::vector<t_real_min> GetParamErrors() const = 0;


	friend std::ostream& operator<<(std::ostream& ostr, const MinuitFuncModel& fkt)
	{
		ostr << fkt.print();
		return ostr;
	}
};

class MinuitFuncModel_nd : public FunctionModel_nd<t_real_min>
{
public:
	virtual ~MinuitFuncModel_nd() {}

	virtual std::size_t GetDim() const = 0;

	virtual bool SetParams(const std::vector<t_real_min>& vecParams) = 0;
	virtual t_real_min operator()(const t_real_min* px) const = 0;

	virtual MinuitFuncModel_nd* copy() const = 0;
	virtual std::string print(bool bFillInSyms=true) const = 0;

	virtual const char* GetModelName() const = 0;


	friend std::ostream& operator<<(std::ostream& ostr, const MinuitFuncModel_nd& fkt)
	{
		ostr << fkt.print();
		return ostr;
	}
};


// generic chi^2 calculation
class Chi2Function : public ROOT::Minuit2::FCNBase
{
protected:
	const MinuitFuncModel *m_pfkt;

	std::size_t m_uiLen;
	const t_real_min* m_px;
	const t_real_min* m_py;
	const t_real_min* m_pdy;

	t_real_min m_dSigma = 1.;

	bool m_bDebug = 0;

public:
	Chi2Function(const MinuitFuncModel* fkt=0,
		std::size_t uiLen=0, const t_real_min *px=0,
		const t_real_min *py=0, const t_real_min *pdy=0)
		: m_pfkt(fkt), m_uiLen(uiLen), m_px(px), m_py(py), m_pdy(pdy)
	{}

	virtual ~Chi2Function() {}

	/* 
	 * chi^2 calculation
	 * based on the example in the Minuit user's guide:
	 * http://seal.cern.ch/documents/minuit/mnusersguide.pdf
	 */
	t_real_min chi2(const std::vector<t_real_min>& vecParams) const
	{
		// cannot operate on m_pfkt directly, because Minuit
		// uses more than one thread!
		std::unique_ptr<MinuitFuncModel> uptrFkt(m_pfkt->copy());
		MinuitFuncModel* pfkt = uptrFkt.get();

		/*bool bParamsOk = */pfkt->SetParams(vecParams);
		//if(!bParamsOk)
		//	return std::numeric_limits<t_real_min>::max();

		return tl::chi2<t_real_min, decltype(*pfkt)>(*pfkt, m_uiLen, m_px, m_py, m_pdy);
	}

	virtual t_real_min Up() const override { return m_dSigma*m_dSigma; }

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		t_real_min dChi2 = chi2(vecParams);
		if(m_bDebug) tl::log_debug("Chi2 = ", dChi2);
		return dChi2;
	}

	void SetSigma(t_real_min dSig) { m_dSigma = dSig; }
	t_real_min GetSigma() const { return m_dSigma; }

	void SetDebug(bool b) { m_bDebug = b; }
};


// in n dimensions
class Chi2Function_nd : public ROOT::Minuit2::FCNBase
{
protected:
	const MinuitFuncModel_nd *m_pfkt;
	std::size_t m_uiDim;

	std::size_t m_uiLen;
	std::vector<const t_real_min*> m_vecpx;

	const t_real_min* m_py;
	const t_real_min* m_pdy;

	bool m_bDebug = 0;

public:
	Chi2Function_nd(const MinuitFuncModel_nd* fkt=0,
		std::size_t uiLen=0, const t_real_min **ppx=0,
		const t_real_min *py=0, const t_real_min *pdy=0)
		: m_pfkt(fkt), m_uiDim(fkt->GetDim()), m_uiLen(uiLen), m_py(py), m_pdy(pdy)
	{
		m_vecpx.resize(m_uiDim);

		for(std::size_t i=0; i<m_uiDim; ++i)
			m_vecpx[i] = ppx[i];
	}

	virtual ~Chi2Function_nd() {}
	t_real_min chi2(const std::vector<t_real_min>& vecParams) const
	{
		std::unique_ptr<MinuitFuncModel_nd> uptrFkt(m_pfkt->copy());
		MinuitFuncModel_nd* pfkt = uptrFkt.get();

		/*bool bParamsOk = */pfkt->SetParams(vecParams);

		std::unique_ptr<t_real_min[]> uptrX(new t_real_min[m_uiDim]);

		t_real_min dChi2 = 0.;
		for(std::size_t i=0; i<m_uiLen; ++i)
		{
			for(std::size_t iX=0; iX<m_uiDim; ++iX)
				uptrX[iX] = m_vecpx[iX][i];

			t_real_min d = (*pfkt)(uptrX.get()) - m_py[i];
			t_real_min dy = m_pdy ? m_pdy[i] : 0.1*d;	// assume 10% error if none given
			if(fabs(dy) < std::numeric_limits<t_real_min>::min())
				dy = std::numeric_limits<t_real_min>::min();

			d /= dy;
			dChi2 += d*d;
		}
		return dChi2;
	}


	virtual t_real_min Up() const override
	{
		// 1. for chi^2
		return 1.;
	}

	virtual t_real_min operator()(const std::vector<t_real_min>& vecParams) const override
	{
		t_real_min dChi2 = chi2(vecParams);
		if(m_bDebug) tl::log_debug("Chi2 = ", dChi2);
		return dChi2;
	}

	void SetDebug(bool b) { m_bDebug = b; }
};
}

#endif
