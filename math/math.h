/**
 * math helpers
 *
 * @author: tweber
 * @date: 23-apr-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_MATH__
#define __TLIBS_MATH__

#include <cmath>
#include <complex>
#include <vector>
#include <limits>
#include <tuple>
#include "../helper/traits.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/constants/constants.hpp>

#ifndef M_PI
	#define M_PI (boost::math::constants::pi<long double>())
#endif

namespace tl {

//template<typename T=double> constexpr T get_pi() { return boost::math::constants::pi<T>(); }
template<typename T=double> constexpr typename get_scalar_type<T>::value_type get_pi()
{ return boost::math::constants::pi<typename get_scalar_type<T>::value_type>(); }

#if __cplusplus >= 201402L
	template<typename T=double> static constexpr T g_pi = get_pi<T>();
#endif

template<typename INT=int> bool is_even(INT i) { return (i%2 == 0); }
template<typename INT=int> bool is_odd(INT i) { return !is_even<INT>(i); }

template<class T=double> constexpr T r2d(T rad) { return rad/get_pi<T>()*T(180); }	// rad -> deg
template<class T=double> constexpr T d2r(T deg) { return deg/T(180)*get_pi<T>(); }	// deg -> rad
template<class T=double> constexpr T r2m(T rad) { return rad/get_pi<T>()*T(180*60); }	// rad -> min
template<class T=double> constexpr T m2r(T min) { return min/T(180*60)*get_pi<T>(); }	// min -> rad

template<typename T>
T sign(T t)
{
	if(t<0.) return -T(1);
	return T(1);
}

template<typename T> T cot(T t)
{
	//return T(1) / std::tan(t);
	return std::tan(T(0.5)*get_pi<T>() - t);
}

template<typename T> T coth(T t)
{
	return T(1) / std::tanh(t);
}

// -----------------------------------------------------------------------------


template<class vec_type>
typename vec_type::value_type mean_value(const vec_type& vec)
{
	typedef typename vec_type::value_type T;
	if(vec.size()==0) return T(0);
	else if(vec.size()==1) return vec[0];

	T tMean = vec[0];
	for(std::size_t i=1; i<vec.size(); ++i)
		tMean += vec[i];
	tMean /= vec.size();

	return tMean;
}

// standard deviation of mean value
template<class vec_type>
typename vec_type::value_type std_dev(const vec_type& vec)
{
	typedef typename vec_type::value_type T;
	if(vec.size()<=1) return T(0);

	T tMean = mean_value(vec);
	T t = T(0);
	for(const T& tval : vec)
		t += (tval-tMean) * (tval-tMean);

	T tN = T(vec.size());
	t /= tN-T(1);

	//std::cout << t << " " << std::sqrt(t) << std::endl;
	return std::sqrt(t);
}


// -----------------------------------------------------------------------------


template<typename T=double>
void diff(unsigned int N, const T* pXIn, const T* pYIn, T* pYOut)
{
	for(unsigned int i=0; i<N-1; ++i)
		pYOut[i] = (pYIn[i+1]-pYIn[i]) / (pXIn[i+1]-pXIn[i]);

	// copy last value
	pYOut[N-1] = pYOut[N-2];
}

// -----------------------------------------------------------------------------


template<typename T> T t_abs(const T& t)
{
	if(t < T(0))
		return -t;
	return t;
}



template<class T=double>
struct _get_epsilon_impl
{
	using t_eps = underlying_value_type_t<T>;

	t_eps operator()() const
	{
		return std::numeric_limits<t_eps>::epsilon();
	}
};

template<class T>
typename _get_epsilon_impl<T>::t_eps get_epsilon()
{
	return _get_epsilon_impl<T>()();
}


template<class T = double, class t_eps = typename _get_epsilon_impl<T>::t_eps,
	LinalgType ty = get_linalg_type<T>::value>
struct _float_equal_impl
{
	bool operator()(T t1, T t2, t_eps eps = get_epsilon<T>()) const
	{
		return t_abs<T>(t1-t2) < eps;
	}
};

template<class T, class t_eps>
struct _float_equal_impl<T, t_eps, LinalgType::COMPLEX>
{
	bool operator()(const T& t1, const T& t2, t_eps eps = get_epsilon<T>()) const
	{
		return t_abs<t_eps>(t1.real()-t2.real()) < eps &&
			t_abs<t_eps>(t1.imag()-t2.imag()) < eps;
	}
};

template<typename T = double>
bool float_equal(T t1, T t2, typename _get_epsilon_impl<T>::t_eps eps = get_epsilon<T>())
{
	return _float_equal_impl<T>()(t1, t2, eps);
}


// -----------------------------------------------------------------------------


template<class T, typename REAL=double>
T lerp(const T& a, const T& b, REAL val)
{
	return a + T((b-a)*val);
}


// x=0..1
template<typename T=double>
T linear_interp(T x0, T x1, T x)
{
	return lerp<T,T>(x0, x1, x);
}


// x=0..1, y=0..1
template<typename T=double>
T bilinear_interp(T x0y0, T x1y0, T x0y1, T x1y1, T x, T y)
{
	T top = linear_interp<T>(x0y1, x1y1, x);
	T bottom = linear_interp<T>(x0y0, x1y0, x);

	return linear_interp<T>(bottom, top, y);
}


template<typename T=double, typename REAL=double,
	template<class...> class t_vec=std::vector>
t_vec<T> linspace(const T& tmin, const T& tmax, std::size_t iNum)
{
	t_vec<T> vec;
	vec.reserve(iNum);

	for(std::size_t i=0; i<iNum; ++i)
		vec.push_back(lerp<T,REAL>(tmin, tmax, REAL(i)/REAL(iNum-1)));
	return vec;
}

template<typename T=double, typename REAL=double,
	template<class...> class t_vec=std::vector>
t_vec<T> logspace(const T& tmin, const T& tmax, std::size_t iNum, T tBase=T(10))
{
	t_vec<T> vec = linspace<T, REAL>(tmin, tmax, iNum);

	for(T& t : vec)
		t = std::pow(tBase, t);
	return vec;
}

template<typename T>
T clamp(T t, T min, T max)
{
	if(t < min) t = min;
	if(t > max) t = max;

	return t;
}


template<class T>
bool is_in_range(T val, T centre, T pm)
{
	pm = std::abs(pm);

	if(val < centre-pm) return false;
	if(val > centre+pm) return false;
	return true;
}

// -----------------------------------------------------------------------------


// solve a*x^2 + b*x + c for x
template<class T=double>
std::vector<T> quadratic_solve(T a, T b, T c)
{
	std::vector<T> vec;

	if(float_equal(a, 0.))
	{
		// b*x + c = 0
		T t = -c/b;
		if(!std::isnan(t) && !std::isinf(t))
			vec.push_back(t);
	}
	else
	{
		T D = b*b - 4.*a*c;
		if(float_equal(D, 0.))
		{
			T t = -b/(2.*a);
			if(!std::isnan(t) && !std::isinf(t))
				vec.push_back(t);
		}
		else if(D > 0.)
		{
			T r = std::sqrt(D);
			T t0 = (-b + r) / (2.*a);
			T t1 = (-b - r) / (2.*a);
			if(!std::isnan(t0) && !std::isinf(t0)) vec.push_back(t0);
			if(!std::isnan(t1) && !std::isinf(t1)) vec.push_back(t1);
		}
	}

	return vec;
}


// -----------------------------------------------------------------------------


/**
 * calculates chi^2 distance of a function model to data points
 */
template<class T, class t_func, class t_iter_dat=T*>
T chi2(const t_func& func, std::size_t N,
	const t_iter_dat x, const t_iter_dat y, const t_iter_dat dy)
{
	using t_dat = typename std::remove_pointer<t_iter_dat>::type;
	T tchi2 = T(0);

	for(std::size_t i=0; i<N; ++i)
	{
		T td = T(y[i]) - func(T(x[i]));
		T tdy = dy ? T(dy[i]) : T(0.1*td);	// 10% error if none given

		if(std::abs(tdy) < std::numeric_limits<t_dat>::min())
			tdy = std::numeric_limits<t_dat>::min();

		T tchi = T(td) / T(tdy);
		tchi2 += tchi*tchi;
	}

	return tchi2;
}

template<class t_vec, class t_func>
typename t_vec::value_type chi2(const t_func& func,
	const t_vec& x, const t_vec& y, const t_vec& dy)
{
	using T = typename t_vec::value_type;
	return chi2<T, t_func, T*>(func, x.size(), x.data(), y.data(),
		dy.size() ? dy.data() : nullptr);
}


// -----------------------------------------------------------------------------


template<typename T=double>
T log(T tbase, T tval)
{
	return T(std::log(tval)/std::log(tbase));
}

template<typename T=double>
T nextpow(T tbase, T tval)
{
	return T(std::pow(tbase, std::ceil(log(tbase, tval))));
}


// -----------------------------------------------------------------------------

/*
template<class T=double> static const T SIGMA2FWHM = T(2)*std::sqrt(T(2)*std::log(T(2)));
template<class T=double> static const T SIGMA2HWHM = std::sqrt(T(2)*std::log(T(2)));
template<class T=double> static const T FWHM2SIGMA = T(1)/SIGMA2FWHM<T>;
template<class T=double> static const T HWHM2SIGMA = T(1)/SIGMA2HWHM<T>;
*/

template<class T=double> T get_SIGMA2FWHM()
{
	return T(2)*std::sqrt(T(2)*std::log(T(2)));
}
template<class T=double> T get_SIGMA2HWHM()
{
	return std::sqrt(T(2)*std::log(T(2)));
}
template<class T=double> T get_FWHM2SIGMA()
{
	return T(1)/get_SIGMA2FWHM<T>();
}
template<class T=double> T get_HWHM2SIGMA()
{
	return T(1)/get_SIGMA2HWHM<T>();
}

static const double SIGMA2FWHM = 2.*std::sqrt(2.*std::log(2.));
static const double SIGMA2HWHM = SIGMA2FWHM/2.;
static const double HWHM2SIGMA = 1./SIGMA2HWHM;
static const double FWHM2SIGMA = 1./SIGMA2FWHM;

template<class T=double>
T gauss_model(T x, T x0, T sigma, T amp, T offs)
{
	T norm = T(1)/(std::sqrt(T(2)*get_pi<T>()) * sigma);
	return amp * norm * std::exp(-0.5 * ((x-x0)/sigma)*((x-x0)/sigma)) + offs;
}

template<class T=double>
T gauss_model_amp(T x, T x0, T sigma, T amp, T offs)
{
	return amp * std::exp(-0.5 * ((x-x0)/sigma)*((x-x0)/sigma)) + offs;
}

template<class T=double>
T lorentz_model_amp(T x, T x0, T hwhm, T amp, T offs)
{
	return amp*hwhm*hwhm / ((x-x0)*(x-x0) + hwhm*hwhm) + offs;
}


// -----------------------------------------------------------------------------
template<class t_real_to, class t_real_from,
	bool bIsEqu = std::is_same<t_real_from, t_real_to>::value>
struct complex_cast
{
	const std::complex<t_real_to>& operator()(const std::complex<t_real_from>& c) const
	{ return c; }
};

template<class t_real_to, class t_real_from>
struct complex_cast<t_real_to, t_real_from, 0>
{
	std::complex<t_real_to> operator()(const std::complex<t_real_from>& c) const
	{ return std::complex<t_real_to>(t_real_to(c.real()), t_real_to(c.imag())); }
};
// -----------------------------------------------------------------------------


#ifdef HAS_COMPLEX_ERF
}

#include <Faddeeva.hh>
using t_real_fadd = double;

namespace tl{

/**
* Complex error function
*/
template<class T=double>
std::complex<T> erf(const std::complex<T>& z)
{
	complex_cast<t_real_fadd, T> cst;
	complex_cast<T, t_real_fadd> inv_cst;
	return inv_cst(::Faddeeva::erf(cst(z)));
}

/**
* Complex complementary error function
*/
template<class T=double>
std::complex<T> erfc(const std::complex<T>& z)
{
	complex_cast<t_real_fadd, T> cst;
	complex_cast<T, t_real_fadd> inv_cst;
	return inv_cst(::Faddeeva::erfc(cst(z)));
}

/**
* Faddeeva function
*/
template<class T=double>
std::complex<T> faddeeva(const std::complex<T>& z)
{
	std::complex<T> i(0, 1.);
	return std::exp(-z*z) * erfc(-i*z);
}

/**
* Voigt profile
* see e.g.: https://en.wikipedia.org/wiki/Voigt_profile
*/
template<class T=double>
T voigt_model(T x, T x0, T sigma, T gamma, T amp, T offs)
{
	T norm = T(1)/(std::sqrt(T(2)*get_pi<T>()) * sigma);
	std::complex<T> z = std::complex<T>(x-x0, gamma) / (sigma * std::sqrt(T(2)));

	return amp*norm * faddeeva<T>(z).real() + offs;
}

template<class T=double>
T voigt_model_amp(T x, T x0, T sigma, T gamma, T amp, T offs)
{
	std::complex<T> z = std::complex<T>(x-x0, gamma) / (sigma * std::sqrt(T(2)));
	return amp * faddeeva<T>(z).real() + offs;
}

#endif


// -----------------------------------------------------------------------------


// wrapper for boost's Y function
template<class T=double>
std::complex<T> Ylm(int l /*0..i*/, int m /*-l..l*/, T th /*0..pi*/, T ph /*0..2pi*/)
{
	return boost::math::spherical_harmonic<T,T>(l,m, th, ph);
}


// -----------------------------------------------------------------------------
// coordinate trafos

template<class T = double>
std::tuple<T,T,T> cart_to_sph(T x, T y, T z)
{
	T rho = std::sqrt(x*x + y*y + z*z);
	T phi = std::atan2(y, x);
	T theta = std::acos(z/rho);

	return std::make_tuple(rho, phi, theta);
}

template<class T = double>
std::tuple<T,T,T> sph_to_cart(T rho, T phi, T theta)
{
	T x = rho * std::cos(phi)*std::sin(theta);
	T y = rho * std::sin(phi)*std::sin(theta);
	T z = rho * std::cos(theta);

	return std::make_tuple(x, y, z);
}

template<class T = double>
std::tuple<T,T,T> cyl_to_sph(T rho_cyl, T phi_cyl, T z_cyl)
{
	T rho = std::sqrt(rho_cyl*rho_cyl + z_cyl*z_cyl);
	T theta = std::acos(z_cyl/rho);

	return std::make_tuple(rho, phi_cyl, theta);
}

template<class T = double>
std::tuple<T,T,T> sph_to_cyl(T rho_sph, T phi_sph, T theta_sph)
{
	T rho = rho_sph * std::sin(theta_sph);
	T z = rho_sph * std::cos(theta_sph);

	return std::make_tuple(rho, phi_sph, z);
}

template<class T = double>
std::tuple<T,T,T> cyl_to_cart(T rho, T phi, T z)
{
	T x = rho * std::cos(phi);
	T y = rho * std::sin(phi);

	return std::make_tuple(x, y, z);
}

template<class T = double>
std::tuple<T,T,T> cart_to_cyl(T x, T y, T z)
{
	T rho = std::sqrt(x*x + y*y);
	T phi = std::atan2(y, x);

	return std::make_tuple(rho, phi, z);
}



template<class T = double>
std::tuple<T,T> crys_to_sph(T twophi_crys, T twotheta_crys)
{
	// converts the out-of-plane scattering angle '2theta' to the spherical theta
	T theta_sph = get_pi<T>()/T(2) - twotheta_crys;
	// converts in-plane scattering angle '2phi' to the spherical phi
	T phi_sph = twophi_crys - get_pi<T>()/T(2);

	return std::make_tuple(phi_sph, theta_sph);
}

template<class T = double>
std::tuple<T,T> sph_to_crys(T phi, T theta)
{
	return crys_to_sph<T>(phi, theta);
}


/**
 * gnomonic projection (similar to perspective projection with fov=90°)
 * @return [x,y]
 * @desc: see http://mathworld.wolfram.com/GnomonicProjection.html
 */
template<class T = double>
std::tuple<T,T> gnomonic_proj(T twophi_crys, T twotheta_crys)
{
	T x = -std::tan(twophi_crys);
	T y = std::tan(twotheta_crys) / std::cos(twophi_crys);

	return std::make_tuple(x, y);
}

/**
 * stereographic projection
 * @return [x,y]
 * @desc: see http://mathworld.wolfram.com/StereographicProjection.html
 */
template<class T = double>
std::tuple<T,T> stereographic_proj(T twophi_crys, T twotheta_crys, T rad)
{
	const T sth = std::sin(twotheta_crys);
	const T cth = std::cos(twotheta_crys);
	const T sph = std::sin(twophi_crys);
	const T cph = std::cos(twophi_crys);

	T x = -T(2) * rad * sph * cth / (T(1) + cth*cph);
	T y = T(2) * rad * sth / (T(1) + cth*cph);

	return std::make_tuple(x, y);
}

}
#endif
