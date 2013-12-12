/*
 * neutron formulas
 * Author: Tobias Weber
 * Date: May 2012, 11-jul-2013
 */

#ifndef __NEUTRON_FORMULAS__
#define __NEUTRON_FORMULAS__

#include "math.h"

#include <boost/units/unit.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/dimensionless_quantity.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/physical_dimensions.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/systems/angle/degrees.hpp>
//#include <boost/units/systems/temperature/celsius.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/neutron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
namespace units = boost::units;
namespace co = boost::units::si::constants::codata;

#include <boost/units/base_units/metric/angstrom.hpp>
#include <boost/units/cmath.hpp>

#ifndef one_meV
	#define one_meV (1e-3 * co::e * units::si::volts)
#endif

#ifndef one_eV
	#define one_eV (co::e * units::si::volts)
#endif

static const double SIGMA2FWHM = 2.*sqrt(2.*log(2.));
static const double SIGMA2HWHM = SIGMA2FWHM/2.;

static const units::quantity<units::si::length> angstrom = 1e-10 * units::si::meter;

static const double KSQ2E = (co::hbar*co::hbar / (2.*co::m_n)) / one_meV / (angstrom*angstrom);
static const double E2KSQ = 1./KSQ2E;

//static const units::quantity<units::si::length> angstrom = 1e-10 * units::si::meter;


// --------------------------------------------------------------------------------
// de Broglie stuff
// lam = h/p
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


// lam = 2pi/k
template<class Sys, class Y>
units::quantity<units::unit<units::length_dimension, Sys>, Y>
k2lam(const units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>& k)
{
	return 2.*M_PI / k;
}

template<class Sys, class Y>
units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>
lam2k(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam)
{
	return 2.*M_PI / lam;
}

template<class Sys, class Y>
units::quantity<units::unit<units::momentum_dimension, Sys>, Y>
k2p(const units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>& k)
{
	return co::hbar*k;
}

template<class Sys, class Y>
units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>
p2k(const units::quantity<units::unit<units::momentum_dimension, Sys>, Y>& p)
{
	return p/co::hbar;
}
// --------------------------------------------------------------------------------

// --------------------------------------------------------------------------------
// E = hbar*omega
template<class Sys, class Y>
units::quantity<units::unit<units::energy_dimension, Sys>, Y>
omega2E(const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& omega)
{
	return co::hbar * omega;
}

template<class Sys, class Y>
units::quantity<units::unit<units::frequency_dimension, Sys>, Y>
E2omega(const units::quantity<units::unit<units::energy_dimension, Sys>, Y>& en)
{
	return en / co::hbar;
}


template<class Sys, class Y>
units::quantity<units::unit<units::energy_dimension, Sys>, Y>
k2E(const units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>& k)
{
	units::quantity<units::unit<units::momentum_dimension, Sys>, Y>
		p = co::hbar*k;
	units::quantity<units::unit<units::energy_dimension, Sys>, Y>
		E = p*p / (2.*co::m_n);
	return E;
}

template<class Sys, class Y>
units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>
E2k(const units::quantity<units::unit<units::energy_dimension, Sys>, Y>& E, bool &bImag)
{
	if(E < 0.*one_meV)
		bImag = 1;
	else
		bImag = 0;

	units::quantity<units::unit<units::momentum_dimension, Sys>, Y>
		p = units::sqrt(2.*co::m_n*units::abs(E));
	units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>
		k = p / co::hbar;
	return k;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// Bragg
// real: n * lam = 2d * sin(twotheta/2)
template<class Sys, class Y>
units::quantity<units::unit<units::length_dimension, Sys>, Y>
bragg_real_lam(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& d,
				const units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>& twotheta,
				double n)
{
	return 2.*d/n * sin(twotheta/2.);
}

template<class Sys, class Y>
units::quantity<units::unit<units::length_dimension, Sys>, Y>
bragg_real_d(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam,
			const units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>& twotheta,
			double n)
{
	return n * lam / (2.* sin(twotheta/2.));
}

template<class Sys, class Y>
units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>
bragg_real_twotheta(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& d,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam,
					double n)
{
	return asin(n*lam/(2.*d)) * 2.;
}

// reciprocal: Q * lam = 4pi * sin(twotheta/2)
template<class Sys, class Y>
units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>
bragg_recip_twotheta(const units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>& Q,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam,
					double n)
{
	return asin(Q*n*lam/(4.*M_PI)) * 2.;
}

template<class Sys, class Y>
units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>
bragg_recip_Q(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam,
			const units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>& twotheta,
			double n)
{
	return 4.*M_PI / (n*lam) * sin(twotheta/2.);
}

template<class Sys, class Y>
units::quantity<units::unit<units::length_dimension, Sys>, Y>
bragg_recip_lam(const units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>& Q,
				const units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>& twotheta,
				double n)
{
	return 4.*M_PI / Q * sin(twotheta/2.) / n;
}
// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
template<class Sys, class Y>
units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>
kinematic_plane(bool bFixedKi,
				const units::quantity<units::unit<units::energy_dimension, Sys>, Y>& EiEf,
				const units::quantity<units::unit<units::energy_dimension, Sys>, Y>& DeltaE,
				const units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>& twotheta)
{
	const units::quantity<units::unit<units::energy_dimension, Sys>, Y> dE = DeltaE;
	if(bFixedKi)
		dE = -dE;

	units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y> Q =
			units::sqrt(2.*co::m_n / co::hbar) *
			(2*EiEf + dE - 2.*units::cos(twotheta)*units::sqrt(EiEf*(EiEf + dE)));

	return Q;
}

template<class Sys, class Y>
units::quantity<units::unit<units::energy_dimension, Sys>, Y>
kinematic_plane(bool bFixedKi, bool bBranch,
				const units::quantity<units::unit<units::energy_dimension, Sys>, Y>& EiEf,
				const units::quantity<units::unit<units::wavenumber_dimension, Sys>, Y>& Q,
				const units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>& twotheta)
{
	auto c = 2.*co::m_n / (co::hbar*co::hbar);
	double ctt = units::cos(twotheta);
	double c2tt = units::cos(2.*twotheta);

	double dSign = -1.;
	if(bBranch)
		dSign = 1.;

	double dSignFixedKf = 1.;
	if(bFixedKi)
		dSignFixedKf = -1.;

	auto rt = c*c*c*c * (-EiEf*EiEf)*ctt*ctt
			+ c*c*c*c*EiEf*EiEf*ctt*ctt*c2tt + 2.*c*c*c*EiEf*Q*Q*ctt*ctt;

	units::quantity<units::unit<units::energy_dimension, Sys>, Y> E =
			1./(c*c)*(dSignFixedKf*2.*c*c*EiEf*ctt*ctt
			- dSignFixedKf*2.*c*c*EiEf
			+ dSign*std::sqrt(2.) * units::sqrt(rt)
			+ dSignFixedKf*c*Q*Q);

	return E;
}
// --------------------------------------------------------------------------------


// --------------------------------------------------------------------
/* formulas.cpp */
extern void init_formulas();
extern bool get_val(const std::string& str, double& dVal, std::string& strUnit);
extern units::quantity<units::si::frequency> get_freq(const std::string& strVar);
extern units::quantity<units::si::energy> get_energy(const std::string& strVar);
extern units::quantity<units::si::length> get_length(const std::string& strUnit);
extern units::quantity<units::si::wavenumber> get_wavenumber(const std::string& strVar);
extern units::quantity<units::si::time> get_time(const std::string& strVar);
extern units::quantity<units::si::plane_angle> get_angle(const std::string& strVar);
extern units::quantity<units::si::temperature> get_temperature(const std::string& strVar);
extern double get_scalar(const std::string& strVar);
// --------------------------------------------------------------------

#endif
