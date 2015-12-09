/*
 * Minuit interface
 *
 * @author Tobias Weber
 * @date April 2012
 * @license GPLv2 or GPLv3
 */

#ifndef __MINUIT_IFACE_H__
#define __MINUIT_IFACE_H__

#include <Minuit2/FCNBase.h>
#include <vector>
#include <iostream>
#include <string>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>

#include "funcmod.h"


namespace tl {

class MinuitFuncModel : public FunctionModel
{
	public:
		virtual ~MinuitFuncModel() {}

		virtual bool SetParams(const std::vector<double>& vecParams) = 0;
		virtual double operator()(double x) const = 0;

		virtual MinuitFuncModel* copy() const = 0;
		virtual std::string print(bool bFillInSyms=true) const = 0;

		virtual const char* GetModelName() const = 0;
		virtual std::vector<std::string> GetParamNames() const = 0;
		virtual std::vector<double> GetParamValues() const = 0;
		virtual std::vector<double> GetParamErrors() const = 0;
};

class MinuitFuncModel_nd
{
	public:
		virtual ~MinuitFuncModel_nd() {}

		virtual unsigned int GetDim() const = 0;

		virtual bool SetParams(const std::vector<double>& vecParams) = 0;
		virtual double operator()(const double* px) const = 0;

		virtual MinuitFuncModel_nd* copy() const = 0;
		virtual std::string print(bool bFillInSyms=true) const = 0;

		virtual const char* GetModelName() const = 0;
};

std::ostream& operator<<(std::ostream& ostr, const MinuitFuncModel& fkt);
std::ostream& operator<<(std::ostream& ostr, const MinuitFuncModel_nd& fkt);


// generic chi^2 calculation
class Chi2Function : public ROOT::Minuit2::FCNBase
{
	protected:
		const MinuitFuncModel *m_pfkt;

		unsigned int m_uiLen;
		const double* m_px;
		const double* m_py;
		const double* m_pdy;

		double m_dSigma = 1.;

	public:
		Chi2Function(const MinuitFuncModel* fkt=0,
			unsigned int uiLen=0, const double* px=0,
			const double *py=0, const double *pdy=0)
			: m_pfkt(fkt), m_uiLen(uiLen), m_px(px), m_py(py), m_pdy(pdy)
		{}

		virtual ~Chi2Function() {}

		double chi2(const std::vector<double>& vecParams) const;
		virtual double Up() const override { return m_dSigma*m_dSigma; }

		virtual double operator()(const std::vector<double>& vecParams) const override
		{
			return chi2(vecParams);
		}

		void SetSigma(double dSig) { m_dSigma = dSig; }
		double GetSigma() const { return m_dSigma; }
};


// in n dimensions
class Chi2Function_nd : public ROOT::Minuit2::FCNBase
{
	protected:
		const MinuitFuncModel_nd *m_pfkt;
		unsigned int m_uiDim;

		unsigned int m_uiLen;
		std::vector<const double*> m_vecpx;

		const double* m_py;
		const double* m_pdy;

	public:
		Chi2Function_nd(const MinuitFuncModel_nd* fkt=0,
			unsigned int uiLen=0, const double** ppx=0,
			const double *py=0, const double *pdy=0)
			: m_pfkt(fkt), m_uiDim(fkt->GetDim()), m_uiLen(uiLen), m_py(py), m_pdy(pdy)
		{
			m_vecpx.resize(m_uiDim);

			for(unsigned int i=0; i<m_uiDim; ++i)
				m_vecpx[i] = ppx[i];
		}

		virtual ~Chi2Function_nd() {}
		double chi2(const std::vector<double>& vecParams) const;

		virtual double Up() const override
		{
			// 1. for chi^2
			return 1.;
		}

		virtual double operator()(const std::vector<double>& vecParams) const override
		{
			return chi2(vecParams);
		}
};
}

#endif
