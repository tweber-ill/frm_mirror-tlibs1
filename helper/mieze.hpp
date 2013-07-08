/*
 * MIEZE formulas
 * Author: Tobias Weber
 * Date: May 2012, 29-may-2013
 */

#ifndef __MIEZE_FORMULAS__
#define __MIEZE_FORMULAS__

#include <cmath>

#include <boost/units/unit.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/dimensionless_quantity.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/physical_dimensions.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/systems/angle/degrees.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/neutron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
namespace units = boost::units;
namespace co = boost::units::si::constants::codata;

//static const units::quantity<units::si::length> angstrom = 1e-10 * units::si::meter;


#include <boost/numeric/ublas/matrix.hpp>
namespace ublas = boost::numeric::ublas;



template<class Sys, class Y>
units::quantity<units::unit<units::momentum_dimension, Sys>, Y>
lam2p(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam)
{
	return co::h / lam;
}

template<class Sys, class Y>
units::quantity<units::unit<units::length_dimension, Sys>, Y>
p2lam(const units::quantity<units::unit<units::momentum_dimension, Sys>, Y>& p)
{
        return co::h / p;
}


//------------------------------------------------------------------------------
// MIEZE time (eq. 117 from [Keller, Golub, Gähler, 2000])
// tau = hbar * omega * Ls / (m*v^3)
template<class Sys, class Y>
units::quantity<units::unit<units::time_dimension, Sys>, Y>
mieze_tau(const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& fm,
				const units::quantity<units::unit<units::length_dimension, Sys>, Y>& Ls,
				const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam)
{
	units::quantity<units::unit<units::velocity_dimension, Sys>, Y> v = lam2p(lam) / co::m_n;
	return 2.*M_PI * fm * Ls * co::hbar / (co::m_n * v*v*v);
}

template<class Sys, class Y>
units::quantity<units::unit<units::frequency_dimension, Sys>, Y>
mieze_tau_fm(const units::quantity<units::unit<units::time_dimension, Sys>, Y>& tau,
						const units::quantity<units::unit<units::length_dimension, Sys>, Y>& Ls,
						const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam)
{
	units::quantity<units::unit<units::velocity_dimension, Sys>, Y> v = lam2p(lam) / co::m_n;
	return tau / (2.*M_PI * Ls * co::hbar) * (co::m_n * v*v*v);
}

template<class Sys, class Y>
units::quantity<units::unit<units::length_dimension, Sys>, Y>
mieze_tau_Ls(const units::quantity<units::unit<units::time_dimension, Sys>, Y>& tau,
						const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& fm,
						const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam)
{
	units::quantity<units::unit<units::velocity_dimension, Sys>, Y> v = lam2p(lam) / co::m_n;
	return tau / (2.*M_PI * fm * co::hbar) * (co::m_n * v*v*v);
}

template<class Sys, class Y>
units::quantity<units::unit<units::length_dimension, Sys>, Y>
mieze_tau_lam(const units::quantity<units::unit<units::time_dimension, Sys>, Y>& tau,
						const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& fm,
						const units::quantity<units::unit<units::length_dimension, Sys>, Y>& Ls)
{
	units::quantity<units::unit<units::velocity_dimension, Sys>, Y> v;
	v = boost::units::root<3>(2.*M_PI * fm * Ls * co::hbar / (tau * co::m_n));

	units::quantity<units::unit<units::momentum_dimension, Sys>, Y> p = v * co::m_n;
	return p2lam(p);
}
//------------------------------------------------------------------------------



template<class Sys, class Y>
double
mieze_reduction_det(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lx,
						const units::quantity<units::unit<units::length_dimension, Sys>, Y>& ly,
						const units::quantity<units::unit<units::length_dimension, Sys>, Y>& xpos,
						const units::quantity<units::unit<units::length_dimension, Sys>, Y>& ypos,
						const units::quantity<units::unit<units::length_dimension, Sys>, Y>& Ls,
						const units::quantity<units::unit<units::time_dimension, Sys>, Y>& tau,
						const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam,
						const units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>& central_phase,
						unsigned int iXPixels=128, unsigned int iYPixels=128,
						ublas::matrix<double>* pPhaseMatrix=0)
{
	using namespace units;
	using namespace co;

	const quantity<unit<mass_dimension, Sys>, Y> mn = co::m_n;
	quantity<unit<velocity_dimension, Sys> > v0 = lam2p(lam)/mn;

	const units::quantity<units::unit<units::frequency_dimension, Sys>, Y> fM = mieze_tau_fm(tau, Ls, lam);
	const quantity<unit<frequency_dimension, Sys>, Y> omegaM = 2.*M_PI*fM;

	quantity<unit<length_dimension, Sys>, Y> lx_inc = lx / double(iXPixels);   // pixel size
	quantity<unit<length_dimension, Sys>, Y> ly_inc = ly / double(iYPixels);   // pixel size

	if(pPhaseMatrix)
		pPhaseMatrix->resize(iXPixels, iYPixels, false);

	quantity<unit<length_dimension, Sys>, Y> dx, dy;
	quantity<unit<area_dimension, Sys>, Y> int_red = 0.*si::meter*si::meter;


	unsigned int iX, iY;
	// integrate over detector
	for(dx=-lx/2., iX=0; dx<lx/2. && iX<iXPixels; dx+=lx_inc, ++iX)
	{
		for(dy=-ly/2., iY=0; dy<ly/2. && iY<iYPixels; dy+=ly_inc, ++iY)
		{
			quantity<unit<length_dimension, Sys>, Y> dx_new = dx-xpos;
			quantity<unit<length_dimension, Sys>, Y> dy_new = dy-ypos;

			quantity<unit<length_dimension, Sys>, Y> path_diff = sqrt(dx_new*dx_new + dy_new*dy_new + Ls*Ls) - Ls;

			// additional time needed for the neutron
			quantity<unit<time_dimension, Sys>, Y> dt = path_diff / v0;

			// additional phase
			double phase = -omegaM * dt + central_phase/units::si::radians;
			phase = fmod(phase, 2.*M_PI);

			int_red += cos(phase/2.)*lx_inc*ly_inc;

			if(pPhaseMatrix)
			{
				phase += 2.*M_PI;
				phase = fmod(phase, 2.*M_PI);
				(*pPhaseMatrix)(iX, iY) = phase;
			}
		}
	}

	double dreduction = int_red / (lx*ly);
	dreduction = fabs(dreduction);

	return dreduction;
}


template<typename T=double>
T get_mieze_freq(const T* px, unsigned int iLen, T dNumOsc=2.)
{
	if(iLen==0)
		return -1.;
	double dTLen = (px[iLen-1]-px[0])/double(iLen-1)*double(iLen);
	return dNumOsc * 2.*M_PI/dTLen;
}


//------------------------------------------------------------------------------

template<class Sys, class Y>
units::quantity<units::unit<units::length_dimension, Sys>, Y>
mieze_condition_L2(const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& f1,
					const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& f2,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& L1)
{
        return L1 / (f2/f1 - 1.);
}

template<class Sys, class Y>
units::quantity<units::unit<units::length_dimension, Sys>, Y>
mieze_condition_L1(const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& f1,
					const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& f2,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& L2)
{
        return L2 * (f2/f1 - 1.);
}

template<class Sys, class Y>
units::quantity<units::unit<units::frequency_dimension, Sys>, Y>
mieze_condition_f2(const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& f1,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& L1,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& L2)
{
        return (L1/L2 + 1.)*f1;
}

template<class Sys, class Y>
units::quantity<units::unit<units::frequency_dimension, Sys>, Y>
mieze_condition_f1(const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& f2,
				const units::quantity<units::unit<units::length_dimension, Sys>, Y>& L1,
				const units::quantity<units::unit<units::length_dimension, Sys>, Y>& L2)
{
        return f2 / (L1/L2 + 1.);
}



template<class Sys, class Y>
units::quantity<units::unit<units::length_dimension, Sys>, Y>
mieze_condition_inel_Ls(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& Ls0,
					const units::quantity<units::unit<units::energy_dimension, Sys>, Y>& dE,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam)
{
	using namespace units;

	quantity<unit<velocity_dimension, Sys> > v0 = lam2p(lam)/co::m_n;

	quantity<unit<velocity_dimension, Sys> > dv = sqrt(v0*v0 + 2.*dE/co::m_n) - v0;
	quantity<unit<velocity_dimension, Sys> > v1 = v0 + dv;
	return Ls0 * (v1*v1*v1/(v0*v0*v0));
}


//------------------------------------------------------------------------------

template<class Sys, class Y>
units::quantity<units::unit<units::frequency_dimension, Sys>, Y>
mieze_det_misaligned_df1(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& L1,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& L2,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& dL,
					const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& f1,
					const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& f2)
{
	units::quantity<units::unit<units::frequency_dimension, Sys>, Y> df;
	df = (f2-f1) - (L1*f2)/(L1+L2+dL);
	return df;
}

template<class Sys, class Y>
units::quantity<units::unit<units::frequency_dimension, Sys>, Y>
mieze_det_misaligned_df2(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& L1,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& L2,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& dL,
					const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& f1,
					const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& f2)
{
	units::quantity<units::unit<units::frequency_dimension, Sys>, Y> df;
	df = (L1 / (L2+dL) + 1.) * f1 - f2;
	return df;
}

//------------------------------------------------------------------------------

#endif
