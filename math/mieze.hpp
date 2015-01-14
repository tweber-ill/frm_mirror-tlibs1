/*
 * MIEZE formulas
 * @author Tobias Weber
 * @date May 2012, 29-may-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __MIEZE_FORMULAS__
#define __MIEZE_FORMULAS__

#include "neutrons.hpp"
#include "linalg.h"

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


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
namespace ublas = boost::numeric::ublas;

namespace tl {

//------------------------------------------------------------------------------
// MIEZE time (eq. 117 from [Keller, Golub, G��hler, 2000])
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


//------------------------------------------------------------------------------
// MIEZE contrast reduction due to detector geometry
template<class Sys, class Y>
Y mieze_reduction_det(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lx,
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
			Y phase = -omegaM * dt + central_phase/units::si::radians;
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

// reduction factor due to detector thickness
template<class Sys, class Y>
Y mieze_reduction_det_d(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& d,
						const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& fM,
						const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam)
{
	using namespace units;
	using namespace co;

	const quantity<unit<mass_dimension, Sys>, Y> mn = co::m_n;
	quantity<unit<velocity_dimension, Sys> > v0 = lam2p(lam)/mn;

	const quantity<unit<frequency_dimension, Sys>, Y> omegaM = 2.*M_PI*fM;

	quantity<unit<length_dimension, Sys>, Y> int_red = 0.*si::meter;
	int_red = std::sin(-0.5*omegaM/v0 * d)/(-0.5*omegaM/v0);

	Y dreduction = int_red / d;
	dreduction = fabs(dreduction);

	return dreduction;
}

//------------------------------------------------------------------------------



template<typename T=double>
T get_mieze_freq(const T* px, unsigned int iLen, T dNumOsc=2.)
{
	if(iLen==0)
		return -1.;
	double dTLen = (px[iLen-1]-px[0])/double(iLen-1)*double(iLen);
	return dNumOsc * 2.*M_PI/dTLen;
}



//------------------------------------------------------------------------------
// MIEZE contrast reduction due to sample geometry

// numerical approximation to the R_sample integral of formula (9) in [Brandl 11]
template<class Sys, class Y>
Y mieze_reduction_sample_cuboid(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& len_x,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& len_y,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& len_z,
					const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& fM,
					const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam,
					const units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>& twotheta,
					const units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>& theta_s,
					unsigned int ITERS=100)
{
	using namespace units;
	using namespace co;

	const quantity<unit<frequency_dimension, Sys>, Y> omegaM = 2.*M_PI*fM;
	quantity<unit<velocity_dimension, Sys> > v = lam2p(lam)/co::m_n;

	ublas::vector<Y> ki(3);
	ki[0] = 0.;
	ki[1] = 0.;
	ki[2] = 1.;

	ublas::vector<Y> kf(3);
	kf[0] = sin(twotheta);
	kf[1] = 0.;
	kf[2] = cos(twotheta);

	ublas::vector<Y> q_dir = ki-kf;

	quantity<unit<length_dimension, Sys>, Y> dX = len_x / Y(ITERS);
	quantity<unit<length_dimension, Sys>, Y> dY = len_y / Y(ITERS);
	quantity<unit<length_dimension, Sys>, Y> dZ = len_z / Y(ITERS);

	quantity<unit<length_dimension, Sys>, Y> x, y, z;

	quantity<unit<volume_dimension, Sys>, Y> integral = 0.*si::meter*si::meter*si::meter;
	quantity<unit<volume_dimension, Sys>, Y> vol = 0.*si::meter*si::meter*si::meter;

	const Y stheta_s = sin(theta_s);
	const Y ctheta_s = cos(theta_s);

	for(x=-len_x/2.; x<len_x/2.; x+=dX)
		for(y=-len_y/2.; y<len_y/2.; y+=dY)
			for(z=-len_z/2.; z<len_z/2.; z+=dZ)
			{
				quantity<unit<length_dimension, Sys>, Y> pos[3];
				// rotate sample
				pos[0] = ctheta_s*x + stheta_s*z;
				pos[1] = y;
				pos[2] = -stheta_s*x + ctheta_s*z;

				quantity<unit<length_dimension, Sys>, Y> path_diff = q_dir[0]*pos[0] + q_dir[1]*pos[1] + q_dir[2]*pos[2];
				Y phase = omegaM * path_diff / v;

				quantity<unit<volume_dimension, Sys>, Y> func_det = dX*dY*dZ;

				vol += func_det;
				integral += func_det * cos(phase);
			}

	//std::cout << "\ndiff_vol: " << (vol - (len_x*len_y*len_z)) << std::endl;
	//std::cout << "dxdydz = " << dX*dY*dZ << std::endl;

	vol = len_x*len_y*len_z;
	return integral / vol;
}

// Scattering with extinction
template<class Sys, class Y>
Y mieze_reduction_sample_cuboid_extinction(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& len_x,
			const units::quantity<units::unit<units::length_dimension, Sys>, Y>& len_y,
			const units::quantity<units::unit<units::length_dimension, Sys>, Y>& len_z,
			const units::quantity<units::unit<units::derived_dimension<units::length_base_dimension, -1>::type, Sys>, Y>& mu,
			const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& fM,
			const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam,
			const units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>& twotheta,
			unsigned int ITERS=100)
{
	const Y SUBDIVS = Y(ITERS);

	using namespace units;
	using namespace co;

	typedef quantity<unit<length_dimension, Sys>, Y> length;
	typedef quantity<unit<volume_dimension, Sys>, Y> volume;
	typedef const quantity<unit<frequency_dimension, Sys>, Y> frequency;
	typedef quantity<unit<velocity_dimension, Sys> > velocity;
	typedef quantity<unit<plane_angle_dimension, Sys>, Y> angle;

	frequency omegaM = 2.*M_PI*fM;
	velocity v = lam2p(lam)/co::m_n;

	ublas::vector<Y> ki(3);
	ki[0] = 0.;
	ki[1] = 0.;
	ki[2] = 1.;

	ublas::vector<Y> kf(3);
	kf[0] = sin(twotheta);
	kf[1] = 0.;
	kf[2] = cos(twotheta);

	ublas::vector<Y> q_dir = ki-kf;


	length x, y, z;
	volume integral = 0.*si::meter*si::meter*si::meter;
	volume vol = 0.*si::meter*si::meter*si::meter;


	angle theta_s = twotheta/2. - M_PI/2.*si::radians;
	const Y stheta_s = sin(theta_s);
	const Y ctheta_s = cos(theta_s);


	length zpath = 1./(mu*2.);	// reflexive: path taken twice
	//zpath /= ctheta_s;
	//if(zpath > len_z)
		zpath = len_z;

	length dX = len_x / SUBDIVS;
	length dY = len_y / SUBDIVS;
	length dZ = zpath / SUBDIVS;

	const volume func_det = dX*dY*dZ;

	for(x=-len_x/2.; x<len_x/2.; x+=dX)
		for(y=-len_y/2.; y<len_y/2.; y+=dY)
			for(z=0.*si::meter; z<zpath; z+=dZ)
			{
				//z = len_z/2.;
				length pos[3];
				// rotate sample
				pos[0] = ctheta_s*x + stheta_s*z;
				pos[1] = y;
				pos[2] = -stheta_s*x + ctheta_s*z;

				length path_diff = q_dir[0]*pos[0] + q_dir[1]*pos[1] + q_dir[2]*pos[2];
				Y phase = omegaM * path_diff / v;


				length dist = z/*+len_z/2.*/;
				dist /= ctheta_s;

				Y extinction_factor = exp(-mu * dist * 2.); // reflexive: path taken twice
				//std::cout << "extinction: " << extinction_factor << ", phase: " << phase << std::endl;
				//extinction_factor = 1.;

				//if(extinction_factor > 1 / std::exp(1))
				{
					vol += extinction_factor * func_det;
					integral += extinction_factor * func_det * cos(phase);
				}
			}

	return integral / vol;
}


// Bragg scattering with extinction
template<class Sys, class Y>
Y mieze_reduction_sample_cuboid_bragg(const units::quantity<units::unit<units::length_dimension, Sys>, Y>& len_x,
			const units::quantity<units::unit<units::length_dimension, Sys>, Y>& len_y,
			const units::quantity<units::unit<units::length_dimension, Sys>, Y>& len_z,
			const units::quantity<units::unit<units::derived_dimension<units::length_base_dimension, -1>::type, Sys>, Y>& mu,
			const units::quantity<units::unit<units::frequency_dimension, Sys>, Y>& fM,
			const units::quantity<units::unit<units::length_dimension, Sys>, Y>& lam,
			const units::quantity<units::unit<units::length_dimension, Sys>, Y>& d_ana,
			unsigned int ITERS=100)
{
	const units::quantity<units::unit<units::plane_angle_dimension, Sys>, Y>
		twotheta = bragg_real_twotheta(d_ana, lam, 1.);

	return mieze_reduction_sample_cuboid_extinction(len_x, len_y, len_z, mu, fM, lam, twotheta, ITERS);
}

//------------------------------------------------------------------------------





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
}

#endif
