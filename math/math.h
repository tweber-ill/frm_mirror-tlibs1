/*
 * math helpers
 *
 * @author: tweber
 * @date: 23-apr-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_MATH__
#define __TLIBS_MATH__

#include <cmath>
//#include <ctgmath>
#include <complex>
#include <vector>
#include <limits>

#ifndef M_PI
	#define M_PI (3.141592653589793238462643383279502884197169)
#endif

namespace tl {
//template<typename T=double> static constexpr T g_pi = static_cast<T>(M_PI);

template<typename INT=int> bool is_even(INT i) { return (i%2 == 0); }
template<typename INT=int> bool is_odd(INT i) { return !is_even<INT>(i); }

template<class T=double> T r2d(T rad) { return rad/T(M_PI)*T(180); }
template<class T=double> T d2r(T deg) { return deg/T(180)*T(M_PI); }

template<typename T>
T sign(T t)
{
	if(t<0.) return -T(1);
	return T(1);
}

template<typename T> T cot(T t) { return T(1)/std::tan(t); }


// -----------------------------------------------------------------------------


template<class vec_type>
typename vec_type::value_type mean_value(const vec_type& vec)
{
	typedef typename vec_type::value_type T;
	if(vec.size()==0) return T();

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

template<typename T=double>
bool float_equal(T t1, T t2, T eps = std::numeric_limits<T>::epsilon())
{
	return t_abs<T>(t1-t2) < eps;
}


// -----------------------------------------------------------------------------


// x=0..1
template<typename T=double>
T linear_interp(T x0, T x1, T x)
{
	return x0 + (x1 - x0)*x;
}


// x=0..1, y=0..1
template<typename T=double>
T bilinear_interp(T x0y0, T x1y0, T x0y1, T x1y1, T x, T y)
{
	T top = linear_interp<T>(x0y1, x1y1, x);
	T bottom = linear_interp<T>(x0y0, x1y0, x);

	return linear_interp<T>(bottom, top, y);
}

template<class T, typename REAL=double>
T lerp(const T& a, const T& b, REAL val)
{
	return a + T((b-a)*val);
}


template<typename T=double, typename REAL=double>
std::vector<T> linspace(const T& tmin, const T& tmax, unsigned int iNum)
{
	std::vector<T> vec;
	vec.reserve(iNum);

	for(unsigned int i=0; i<iNum; ++i)
		vec.push_back(REAL(i)*(tmax-tmin)/REAL(iNum-1) + tmin);

	return vec;
}

template<typename T=double, typename REAL=double>
std::vector<T> logspace(const T& tmin, const T& tmax, unsigned int iNum, T tBase=T(10))
{
	std::vector<T> vec = linspace<T, REAL>(tmin, tmax, iNum);
	for(T& t : vec)
		t = std::pow(tBase, t);
	return vec;
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


template<class T, class t_func>
T chi2(const t_func& func, std::size_t N,
			const T* x, const T* y, const T* dy)
{
	T tchi2 = T(0);

	for(std::size_t i=0; i<N; ++i)
	{
		T td = y[i] - func(x[i]);
		T tdy = dy ? dy[i] : 0.1*td;	// 10% error if none given

		if(std::fabs(tdy) < std::numeric_limits<T>::min())
			tdy = std::numeric_limits<T>::min();

		T tchi = td / tdy;
		tchi2 += tchi*tchi;
	}

	return tchi2;
}

template<class t_vec, class t_func>
typename t_vec::value_type chi2(const t_func& func,
			const t_vec& x, const t_vec& y, const t_vec& dy)
{
	using T = typename t_vec::value_type;
	return chi2<T, t_func>(func, x.size(), x.data(), y.data(), dy.size()?dy.data():nullptr);
}


// -----------------------------------------------------------------------------

/*
template<class T=double> static const T SIGMA2FWHM = T(2)*std::sqrt(T(2)*std::log(T(2)));
template<class T=double> static const T SIGMA2HWHM = std::sqrt(T(2)*std::log(T(2)));
template<class T=double> static const T FWHM2SIGMA = T(1)/SIGMA2FWHM<T>;
template<class T=double> static const T HWHM2SIGMA = T(1)/SIGMA2HWHM<T>;
*/

static const double SIGMA2FWHM = 2.*std::sqrt(2.*std::log(2.));
static const double SIGMA2HWHM = SIGMA2FWHM/2.;
static const double HWHM2SIGMA = 1./SIGMA2HWHM;
static const double FWHM2SIGMA = 1./SIGMA2FWHM;

template<class T=double>
T gauss_model(T x, T x0, T sigma, T amp, T offs)
{
	return amp * std::exp(-0.5 * ((x-x0)/sigma)*((x-x0)/sigma)) + offs;
}



#ifdef HAS_COMPLEX_ERF

}
#include <Faddeeva.hh>
namespace tl{

template<class T=double>
std::complex<T> erf(const std::complex<T>& z)
{
	return ::Faddeeva::erf(z);
}

template<class T=double>
std::complex<T> erfc(const std::complex<T>& z)
{
	return ::Faddeeva::erfc(z);
}

template<class T=double>
std::complex<T> faddeeva(const std::complex<T>& z)
{
	std::complex<T> i(0, 1.);
	return std::exp(-z*z) * erfc(-i*z);
}

template<class T=double>
T voigt_model(T x, T x0, T sigma, T gamma, T amp, T offs)
{
	std::complex<T> z = std::complex<T>(x-x0, gamma) / (sigma * std::sqrt(2.));
	return amp * (faddeeva<T>(z)).real() / (sigma * std::sqrt(2.*M_PI)) + offs;
}

#endif


}

#endif
