/*
 * neutron formulas
 * @author Tobias Weber
 * @date May 2012, 11-jul-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_NEUTRONS__
#define __TLIBS_NEUTRONS__

#include "math.h"
#include "linalg.h"
#include "units.h"
#include "../helper/exception.h"

#include <cmath>


namespace tl {
/*
template<class T=double> T KSQ2E = T((co::hbar*co::hbar / (T(2)*co::m_n)) / one_meV / (angstrom*angstrom));
template<class T=double> T E2KSQ = T(1)/KSQ2E<T>;
*/

static const double KSQ2E = (co::hbar*co::hbar / (2.*co::m_n)) / one_meV / (angstrom*angstrom);
static const double E2KSQ = 1./KSQ2E;


// --------------------------------------------------------------------------------
// de Broglie stuff
// lam = h/p
template<class Sys, class Y>
t_momentum<Sys,Y> lam2p(const t_length<Sys,Y>& lam)
{
	return co::h / lam;
}

template<class Sys, class Y>
t_length<Sys,Y> p2lam(const t_momentum<Sys,Y>& p)
{
	return co::h / p;
}


// lam = 2pi/k
template<class Sys, class Y>
t_length<Sys,Y> k2lam(const t_wavenumber<Sys,Y>& k)
{
	return 2.*M_PI / k;
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> lam2k(const t_length<Sys,Y>& lam)
{
	return 2.*M_PI / lam;
}

template<class Sys, class Y>
t_momentum<Sys,Y> k2p(const t_wavenumber<Sys,Y>& k)
{
	return co::hbar*k;
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> p2k(const t_momentum<Sys,Y>& p)
{
	return p/co::hbar;
}

template<class Sys, class Y>
t_velocity<Sys,Y> k2v(const t_wavenumber<Sys,Y>& k)
{
	return k2p(k) / co::m_n;
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> v2k(const t_velocity<Sys,Y>& v)
{
	return co::m_n*v/co::hbar;
}
// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// E = hbar*omega
template<class Sys, class Y>
t_energy<Sys,Y> omega2E(const t_freq<Sys,Y>& omega)
{
	return co::hbar * omega;
}

template<class Sys, class Y>
t_freq<Sys,Y> E2omega(const t_energy<Sys,Y>& en)
{
	return en / co::hbar;
}


template<class Sys, class Y>
t_energy<Sys,Y> k2E(const t_wavenumber<Sys,Y>& k)
{
	t_momentum<Sys,Y> p = co::hbar*k;
	t_energy<Sys,Y> E = p*p / (Y(2.)*co::m_n);
	return E;
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> E2k(const t_energy<Sys,Y>& _E, bool &bImag)
{
	bImag = (_E < Y(0.)*get_one_meV<Y>());
	t_energy<Sys,Y> E = bImag ? -_E : _E;

	// hack
	const t_mass<Sys,Y> m_n = Y(co::m_n/units::si::kilograms)*units::si::kilograms;
	const t_action<Sys,Y> hbar = Y(co::hbar/units::si::joules/units::si::seconds)*
		units::si::joules*units::si::seconds;

	auto pp = Y(2.) * m_n * E;
	//t_momentum<Sys,Y> p = units::sqrt<typename decltype(pp)::unit_type, Y>(pp);
	t_momentum<Sys,Y> p = my_units_sqrt<t_momentum<Sys,Y>, decltype(pp)>(pp);
	t_wavenumber<Sys,Y> k = p / hbar;
	return k;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// Bragg
// real: n * lam = 2d * sin(twotheta/2)
template<class Sys, class Y>
t_length<Sys,Y> bragg_real_lam(const t_length<Sys,Y>& d,
	const t_angle<Sys,Y>& twotheta, Y n)
{
	return 2.*d/n * sin(twotheta/2.);
}

template<class Sys, class Y>
t_length<Sys,Y> bragg_real_d(const t_length<Sys,Y>& lam,
			const t_angle<Sys,Y>& twotheta,
			Y n)
{
	return n * lam / (2.* sin(twotheta/2.));
}

template<class Sys, class Y>
t_angle<Sys,Y> bragg_real_twotheta(const t_length<Sys,Y>& d,
	const t_length<Sys,Y>& lam, Y n)
{
	return asin(n*lam/(2.*d)) * 2.;
}

// reciprocal: Q * lam = 4pi * sin(twotheta/2)
template<class Sys, class Y>
t_angle<Sys,Y> bragg_recip_twotheta(const t_wavenumber<Sys,Y>& Q,
	const t_length<Sys,Y>& lam, Y n)
{
	return asin(Q*n*lam/(4.*M_PI)) * 2.;
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> bragg_recip_Q(const t_length<Sys,Y>& lam,
	const t_angle<Sys,Y>& twotheta, Y n)
{
	return 4.*M_PI / (n*lam) * sin(twotheta/2.);
}

template<class Sys, class Y>
t_length<Sys,Y> bragg_recip_lam(const t_wavenumber<Sys,Y>& Q,
	const t_angle<Sys,Y>& twotheta, Y n)
{
	return 4.*M_PI / Q * sin(twotheta/2.) / n;
}
// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------

// see e.g. ILL blue book sec. 2.6-2
template<class Sys, class Y>
t_wavenumber<Sys,Y> kinematic_plane(bool bFixedKi,
	const t_energy<Sys,Y>& EiEf, const t_energy<Sys,Y>& DeltaE,
	const t_angle<Sys,Y>& twotheta)
{
	const t_energy<Sys,Y> dE = DeltaE;
	if(bFixedKi)
		dE = -dE;

	t_wavenumber<Sys,Y> Q =
		units::sqrt(2.*co::m_n / co::hbar) *
		(2*EiEf + dE - 2.*units::cos(twotheta)*units::sqrt(EiEf*(EiEf + dE)));

	return Q;
}

template<class Sys, class Y>
t_energy<Sys,Y> kinematic_plane(bool bFixedKi, bool bBranch,
	const t_energy<Sys,Y>& EiEf, const t_wavenumber<Sys,Y>& Q,
	const t_angle<Sys,Y>& twotheta)
{
	auto c = 2.*co::m_n / (co::hbar*co::hbar);
	Y ctt = units::cos(twotheta);
	Y c2tt = units::cos(2.*twotheta);

	Y dSign = -1.;
	if(bBranch)
		dSign = 1.;

	Y dSignFixedKf = 1.;
	if(bFixedKi)
		dSignFixedKf = -1.;

	auto rt = c*c*c*c * (-EiEf*EiEf)*ctt*ctt
		+ c*c*c*c*EiEf*EiEf*ctt*ctt*c2tt
		+ 2.*c*c*c*EiEf*Q*Q*ctt*ctt;

	t_energy<Sys,Y> E =
		1./(c*c)*(dSignFixedKf*2.*c*c*EiEf*ctt*ctt
		- dSignFixedKf*2.*c*c*EiEf
		+ dSign*std::sqrt(2.) * units::sqrt(rt)
		+ dSignFixedKf*c*Q*Q);

	return E;
}
// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// Debye-Waller factor

template<class Sys, class Y>
Y debye_waller_high_T(const t_temperature<Sys,Y>& T_D,
	const t_temperature<Sys,Y>& T, const t_mass<Sys,Y>& M,
	const t_wavenumber<Sys,Y>& Q, t_length_square<Sys,Y>* pZeta_sq=0)
{
	t_length_square<Sys,Y> zeta_sq;
	zeta_sq = 9.*co::hbar*co::hbar / (co::k_B * T_D * M) * T/T_D;
	Y dwf = units::exp(-1./3. * Q*Q * zeta_sq);

	if(pZeta_sq) *pZeta_sq = zeta_sq;
	return dwf;
}


template<class Sys, class Y>
Y debye_waller_low_T(const t_temperature<Sys,Y>& T_D,
	const t_temperature<Sys,Y>& T, const t_mass<Sys,Y>& M, 
	const t_wavenumber<Sys,Y>& Q, t_length_square<Sys,Y>* pZeta_sq=0)
{
	t_length_square<Sys,Y> zeta_sq;
	zeta_sq = 9.*co::hbar*co::hbar / (4.*co::k_B*T_D*M) * (1. + 2./3. * M_PI*M_PI * (T/T_D)*(T/T_D));
	Y dwf = units::exp(-1./3. * Q*Q * zeta_sq);

	if(pZeta_sq) *pZeta_sq = zeta_sq;
	return dwf;
}

// --------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// scattering triangle / TAS stuff

// Q_vec = ki_vec - kf_vec
// kf_vec = ki_vec - Q_vec
// kf^2 = ki^2 + Q^2 - 2ki Q cos th
// cos th = (-kf^2 + ki^2 + Q^2) / (2kiQ)
template<class Sys, class Y>
t_angle<Sys,Y> get_angle_ki_Q(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf,
	const t_wavenumber<Sys,Y>& Q,
	bool bPosSense=1, bool bAngleOutsideTriag=0)
{
	t_angle<Sys,Y> angle;

	if(Q*(1e-10 * meters) == 0.)
		angle = M_PI/2. * radians;
	else
	{
		auto c = (ki*ki - kf*kf + Q*Q)/(2.*ki*Q);
		if(units::abs(c) > 1.)
			throw Err("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(bAngleOutsideTriag) angle = M_PI*radians - angle;
	if(!bPosSense) angle = -angle;

	return angle;
}

// Q_vec = ki_vec - kf_vec
// ki_vec = Q_vec + kf_vec
// ki^2 = Q^2 + kf^2 + 2Q kf cos th
// cos th = (ki^2 - Q^2 - kf^2) / (2Q kf)
template<class Sys, class Y>
t_angle<Sys,Y> get_angle_kf_Q(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf,
	const t_wavenumber<Sys,Y>& Q,
	bool bPosSense=1, bool bAngleOutsideTriag=1)
{
	t_angle<Sys,Y> angle;

	if(Q*(1e-10 * meters) == 0.)
		angle = M_PI/2. * radians;
	else
	{
		auto c = (-kf*kf + ki*ki - Q*Q)/(2.*kf*Q);
		if(units::abs(c) > 1.)
			throw Err("Scattering triangle not closed.");

		angle = units::acos(c);
	}

	if(!bAngleOutsideTriag) angle = M_PI*radians - angle;
	if(!bPosSense) angle = -angle;

	return angle;
}


template<class Sys, class Y>
t_angle<Sys,Y> get_mono_twotheta(const t_wavenumber<Sys,Y>& k,
	const t_length<Sys,Y>& d, bool bPosSense=1)
{
	const Y dOrder = 1.;
	t_angle<Sys,Y> tt = bragg_real_twotheta(d, k2lam(k), dOrder);
	if(!bPosSense)
		tt = -tt;
	return tt;
}

template<class Sys, class Y>
t_wavenumber<Sys,Y> get_mono_k(const t_angle<Sys,Y>& _theta,
	const t_length<Sys,Y>& d, bool bPosSense=1)
{
	t_angle<Sys,Y> theta = _theta;
	if(!bPosSense)
		theta = -theta;

	const Y dOrder = 1.;
	return lam2k(bragg_real_lam(d, 2.*theta, dOrder));
}


// Q_vec = ki_vec - kf_vec
// Q^2 = ki^2 + kf^2 - 2ki kf cos 2th
// cos 2th = (-Q^2 + ki^2 + kf^2) / (2ki kf)
template<class Sys, class Y>
t_angle<Sys,Y> get_sample_twotheta(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf, const t_wavenumber<Sys,Y>& Q,
	bool bPosSense=1)
{
	t_dimensionless<Sys,Y> ttCos = (ki*ki + kf*kf - Q*Q)/(2.*ki*kf);
	if(units::abs(ttCos) > 1.)
		throw Err("Scattering triangle not closed.");

	t_angle<Sys,Y> tt;
	tt = units::acos(ttCos);

	if(!bPosSense) tt = -tt;
	return tt;
}


// again cos theorem:
// Q_vec = ki_vec - kf_vec
// Q^2 = ki^2 + kf^2 - 2ki kf cos 2th
// Q = sqrt(ki^2 + kf^2 - 2ki kf cos 2th)
template<class Sys, class Y>
const t_wavenumber<Sys,Y>
get_sample_Q(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf, const t_angle<Sys,Y>& tt)
{
	t_dimensionless<Sys,Y> ctt = units::cos(tt);
	decltype(ki*ki) Qsq = ki*ki + kf*kf - 2.*ki*kf*ctt;
	if(Y(Qsq*get_one_angstrom()*get_one_angstrom()) < 0.)
	{
		// TODO

		Qsq = -Qsq;
	}

	t_wavenumber<Sys,Y> Q = units::sqrt(Qsq);
	return Q;
}



template<class Sys, class Y>
t_energy<Sys,Y> get_energy_transfer(const t_wavenumber<Sys,Y>& ki,
	const t_wavenumber<Sys,Y>& kf)
{
	return k2E<Sys,Y>(ki) - k2E<Sys,Y>(kf);
}


// (hbar*ki)^2 / (2*mn)  -  (hbar*kf)^2 / (2mn)  =  E
// 1) ki^2  =  +E * 2*mn / hbar^2  +  kf^2
// 2) kf^2  =  -E * 2*mn / hbar^2  +  ki^2
template<class Sys, class Y>
t_wavenumber<Sys,Y> get_other_k(const t_energy<Sys,Y>& E,
	const t_wavenumber<Sys,Y>& kfix, bool bFixedKi)
{
	// hack
	const t_mass<Sys,Y> m_n = Y(co::m_n/units::si::kilograms)*units::si::kilograms;
	const t_action<Sys,Y> hbar = Y(co::hbar/units::si::joules/units::si::seconds)*
		units::si::joules*units::si::seconds;

	auto kE_sq = E*Y(2.)*m_n/hbar/hbar;
	if(bFixedKi) kE_sq = -kE_sq;

	auto k_sq = kE_sq + kfix*kfix;
	if(k_sq*get_one_angstrom()*get_one_angstrom() < 0.)
		throw Err("Scattering triangle not closed.");

	return my_units_sqrt<t_wavenumber<Sys,Y>, decltype(k_sq)>(k_sq);
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// kf^3 factor, see e.g. Shirane p. 125

template<class Sys, class Y>
Y ana_effic_factor(const t_wavenumber<Sys, Y>& kf, const t_length<Sys, Y>& d)
{
	t_angle<Sys, Y> theta = 0.5*units::abs(get_mono_twotheta<Sys, Y>(kf, d, true));
	return kf*kf*kf / units::tan(theta) *
		get_one_angstrom()*get_one_angstrom()*get_one_angstrom();
}

// --------------------------------------------------------------------------------



// --------------------------------------------------------------------------------
// spurions

// Bragg tail -> see Shirane p. 152
template<class Sys, class Y>
t_energy<Sys,Y> get_bragg_tail(t_wavenumber<Sys,Y> k,
	t_wavenumber<Sys,Y> q, bool bConstEi=0)
{
	Y t = q / (2.*k);
	if(!bConstEi)
		t = -t;

	t_energy<Sys,Y> E = co::hbar*co::hbar / co::m_n * k*q*(1.+t);
	return E;
}


// higher-order inelastic spurions -> Shirane pp. 146-148
template<class Sys, class Y>
t_energy<Sys,Y> get_inelastic_spurion(bool bConstEi, t_energy<Sys,Y> E, 
	unsigned int iOrderMono, unsigned int iOrderAna)
{
	const Y dOrderMonoSq = Y(iOrderMono)*Y(iOrderMono);
	const Y dOrderAnaSq = Y(iOrderAna)*Y(iOrderAna);

	t_energy<Sys,Y> E_sp;

	// formulas from Shirane, p. 147
	if(bConstEi)
		E_sp = (1. - dOrderMonoSq/dOrderAnaSq) * E;
	else
		E_sp = (dOrderAnaSq/dOrderMonoSq - 1.) * E;

	return E_sp;
}

template<class Y=double>
struct InelasticSpurion
{
	Y dE_meV = 0.;
	unsigned int iOrderMono = 1;
	unsigned int iOrderAna = 1;
};

template<class Sys, class Y>
std::vector<InelasticSpurion<Y>> check_inelastic_spurions(bool bConstEi,
	t_energy<Sys,Y> Ei, t_energy<Sys,Y> Ef,
	t_energy<Sys,Y> E, unsigned int iMaxOrder=5)
{
	const Y dESensitivity = 0.25;	// meV

	std::vector<InelasticSpurion<Y>> vecSpuris;

	for(unsigned int iOrder=1; iOrder<=iMaxOrder; ++iOrder)
	{
		InelasticSpurion<Y> spuri;
		t_energy<Sys,Y> EiEf;

		if(bConstEi)
		{
			spuri.iOrderAna = iOrder;
			EiEf = Ei;
		}
		else
		{
			spuri.iOrderMono = iOrder;
			EiEf = Ef;
		}

		spuri.dE_meV = get_inelastic_spurion(bConstEi, EiEf,
			spuri.iOrderMono, spuri.iOrderAna) / get_one_meV<Y>();

		//std::cout << spuri.dE_meV << " *** " << Y(E/get_one_meV<Y>()) << std::endl;
		if(spuri.dE_meV!=0. && float_equal(spuri.dE_meV, Y(E/get_one_meV<Y>()), dESensitivity))
			vecSpuris.push_back(spuri);
	}

	return vecSpuris;
}

struct ElasticSpurion
{
	bool bAType = 0;
	bool bMType = 0;

	bool bAKfSmallerKi = 0;
	bool bMKfSmallerKi = 0;
};

// accidental elastic (currat-axe) spurions -> Shirane pp. 150-155 (esp. fig. 6.2)
template<typename T=double>
ElasticSpurion check_elastic_spurion(const ublas::vector<T>& ki,
	const ublas::vector<T>& kf, const ublas::vector<T>& q)
{
	const T dKi = ublas::norm_2(ki);
	const T dKf = ublas::norm_2(kf);
	const T dq = ublas::norm_2(q);

	const T dAngleSensitivity = 2.;
	const T dQSensitivity = std::max(dKi, dKf) / 50.;


	ElasticSpurion result;

	ublas::vector<T> ki_norm = ki;	ki_norm /= dKi;
	ublas::vector<T> kf_norm = kf;	kf_norm /= dKf;

	// Q, q and G point in the opposite direction in Shirane!
	// Shirane: Q = kf - ki, E = Ei - Ef
	// here: Q = ki - kf, E = Ei - Ef
	ublas::vector<T> q_norm = -q;	q_norm /= dq;

	T dAngleKfq = std::acos(ublas::inner_prod(kf_norm, q_norm));
	T dAngleKiq = std::acos(ublas::inner_prod(ki_norm, q_norm));

	//std::cout << "angle ki q: " << dAngleKiq/M_PI*180. << std::endl;
	//std::cout << "angle kf q: " << dAngleKfq/M_PI*180. << std::endl;

	bool bKiqParallel = 0, bkiqAntiParallel = 0;
	bool bKfqParallel = 0, bKfqAntiParallel = 0;

	if(float_equal(dAngleKiq, 0., tl::d2r(dAngleSensitivity)))
		bKiqParallel = 1;
	else if(float_equal(dAngleKiq, M_PI, tl::d2r(dAngleSensitivity)))
		bkiqAntiParallel = 1;
	if(float_equal(dAngleKfq, 0., tl::d2r(dAngleSensitivity)))
		bKfqParallel = 1;
	else if(float_equal(dAngleKfq, M_PI, tl::d2r(dAngleSensitivity)))
		bKfqAntiParallel = 1;

	// type A: q || kf, kf > ki
	if(bKfqParallel)
	{
		T dApparentKf = dKf - dq;

		if(float_equal(dApparentKf, dKi, dQSensitivity))
		{
			result.bAType = 1;
			result.bAKfSmallerKi = 0;
		}
	}
	// type A: q || kf, kf < ki
	else if(bKfqAntiParallel)
	{
		T dApparentKf = dKf + dq;

		if(float_equal(dApparentKf, dKi, dQSensitivity))
		{
			result.bAType = 1;
			result.bAKfSmallerKi = 1;
		}
	}

	// type M: q || ki, kf > ki
	if(bKiqParallel)
	{
		T dApparentKi = dKi + dq;

		if(float_equal(dApparentKi, dKf, dQSensitivity))
		{
			result.bMType = 1;
			result.bMKfSmallerKi = 0;
		}
	}
	// type M: q || ki, kf < ki
	else if(bkiqAntiParallel)
	{
		T dApparentKi = dKi - dq;

		if(float_equal(dApparentKi, dKf, dQSensitivity))
		{
			result.bMType = 1;
			result.bMKfSmallerKi = 1;
		}
	}

	return result;
}


// --------------------------------------------------------------------------------

template<class t_real=double>
t_real bose(t_real E, t_real T)
{
	t_real kB = co::k_B * units::si::kelvin/tl::get_one_meV();

	if(E >= 0.)
		return 1./(std::exp(std::abs(E)/(kB*T)) - 1.) + 1.;
	else
		return 1./(std::exp(std::abs(E)/(kB*T)) - 1.);
}

template<class Sys, class Y>
Y bose(const t_energy<Sys,Y>& E, const t_temperature<Sys,Y>& T)
{
	return bose<Y>(Y(E/tl::get_one_meV()), Y(T/tl::kelvin));
}

// see: B. Fak, B. Dorner, Physica B 234-236 (1997) pp. 1107-1108
template<class t_real=double>
t_real DHO_model(t_real E, t_real T, t_real E0, t_real hwhm, t_real amp, t_real offs)
{
	if(E0*E0 - hwhm*hwhm < 0.) return 0.;
	return std::abs(bose<t_real>(E, T)*amp/(E0*M_PI) *
		(hwhm/((E-E0)*(E-E0) + hwhm*hwhm) - hwhm/((E+E0)*(E+E0) + hwhm*hwhm)))
		+ offs;
}


// --------------------------------------------------------------------------------


// get macroscopic from microscopic cross-section
template<class Sys, class Y=double>
t_length_inverse<Sys, Y> macro_xsect(const t_area<Sys, Y>& xsect,
	unsigned int iNumAtoms, const t_volume<Sys, Y>& volUC)
{
	return xsect * Y(iNumAtoms) / volUC;
}



// --------------------------------------------------------------------------------

// thin lens equation: 1/f = 1/lenB + 1/lenA
template<class Sys, class Y=double>
t_length<Sys, Y> focal_len(const t_length<Sys, Y>& lenBefore, const t_length<Sys, Y>& lenAfter)
{
	const t_length_inverse<Sys, Y> f_inv = Y(1)/lenBefore + Y(1)/lenAfter;
	return Y(1) / f_inv;
}

// optimal mono/ana curvature, see e.g. Monochromator_curved.comp in McStas or Shirane p. 66
template<class Sys, class Y=double>
t_length<Sys, Y> foc_curv(const t_length<Sys, Y>& lenBefore, const t_length<Sys, Y>& lenAfter,
	const t_angle<Sys, Y>& tt, bool bVert)
{
	const t_length<Sys, Y> f = focal_len<Sys, Y>(lenBefore, lenAfter);
	const Y s = Y(units::abs(units::sin(Y(0.5)*tt)));

	const t_length<Sys, Y> curv = bVert ? Y(2)*f*s : Y(2)*f/s;

	//std::cout << "f = " << f << std::endl;
	//std::cout << "s = " << s << std::endl;
	//std::cout << "curv = " <<  curv << std::endl;
	return curv;
}


}

#endif
