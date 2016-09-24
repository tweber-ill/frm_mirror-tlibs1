/**
 * distributions
 * @author tweber
 * @date sep-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_DISTR_H__
#define __TLIBS_DISTR_H__

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>

namespace tl {
namespace m = boost::math;


enum class DistrType
{
	NORMAL,
	CAUCHY,
	POISSON,
	BINOMIAL,
	HYPERGEOMETRIC,
	CHI2,
	STUDENT,

	NONE,
};


template<class t_real, class t_distr, class=void> struct distr_traits {};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::normal_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 2;
	static constexpr DistrType distr_type = DistrType::NORMAL;
	static constexpr const char* pcName = "Normal";
	static constexpr const char* pcParam1 = "mu";
	static constexpr const char* pcParam2 = "sigma";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::cauchy_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 2;
	static constexpr DistrType distr_type = DistrType::CAUCHY;
	static constexpr const char* pcName = "Cauchy";
	static constexpr const char* pcParam1 = "mu";
	static constexpr const char* pcParam2 = "sigma";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::poisson_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 1;
	static constexpr DistrType distr_type = DistrType::POISSON;
	static constexpr const char* pcName = "Poisson";
	static constexpr const char* pcParam1 = "lambda";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::binomial_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 2;
	static constexpr DistrType distr_type = DistrType::BINOMIAL;
	static constexpr const char* pcName = "Binomial";
	static constexpr const char* pcParam1 = "n";
	static constexpr const char* pcParam2 = "p";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::hypergeometric_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 3;
	static constexpr DistrType distr_type = DistrType::HYPERGEOMETRIC;
	static constexpr const char* pcName = "Hypergeometric";
	static constexpr const char* pcParam1 = "r";
	static constexpr const char* pcParam2 = "n";
	static constexpr const char* pcParam3 = "N";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::chi_squared_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 1;
	static constexpr DistrType distr_type = DistrType::CHI2;
	static constexpr const char* pcName = "Chi^2";
	static constexpr const char* pcParam1 = "dof";
};

template<class t_real, class t_distr>
struct distr_traits<t_real, t_distr, typename std::enable_if<std::is_same<t_distr, m::students_t_distribution<t_real>>::value>::type>
{
	using value_type = t_real;
	static constexpr std::size_t iNumArgs = 1;
	static constexpr DistrType distr_type = DistrType::STUDENT;
	static constexpr const char* pcName = "Student";
	static constexpr const char* pcParam1 = "dof";
};


template<class t_real=double>
class DistrBase
{
public:
	virtual t_real pdf(t_real x) const = 0;
	virtual t_real cdf(t_real x) const = 0;
};


template<class t_distr, class t_real=typename t_distr::value_type,
	std::size_t iParams=distr_traits<t_real, t_distr>::iNumArgs>
class Distr : public DistrBase<t_real>
{
public:
	using value_type = t_real;
	using distr_type = t_distr;
	using traits_type = distr_traits<t_real, t_distr>;

protected:
	t_distr distr;

public:
	template<std::size_t _iParams=iParams>
	Distr(t_real dParam1,
		typename std::enable_if<_iParams==1, void>::type* =nullptr)
		: distr(dParam1)
	{}

	template<std::size_t _iParams=iParams>
	Distr(t_real dParam1, t_real dParam2,
		typename std::enable_if<_iParams==2, void>::type* =nullptr)
		: distr(dParam1, dParam2)
	{}

	template<std::size_t _iParams=iParams>
	Distr(t_real dParam1, t_real dParam2, t_real dParam3,
		typename std::enable_if<_iParams==3, void>::type* =nullptr)
		: distr(dParam1, dParam2, dParam3)
	{}

	virtual t_real pdf(t_real x) const override
	{
		return m::pdf(distr, x);
	}

	virtual t_real cdf(t_real x) const override
	{
		return m::cdf(distr, x);
	}
};

template<class t_real> using t_normal_dist = Distr<m::normal_distribution<t_real>>;
template<class t_real> using t_cauchy_dist = Distr<m::cauchy_distribution<t_real>>;
template<class t_real> using t_poisson_dist = Distr<m::poisson_distribution<t_real>>;
template<class t_real> using t_binomial_dist = Distr<m::binomial_distribution<t_real>>;
template<class t_real> using t_hypergeo_dist = Distr<m::hypergeometric_distribution<t_real>>;
template<class t_real> using t_chi2_dist = Distr<m::chi_squared_distribution<t_real>>;
template<class t_real> using t_student_dist = Distr<m::students_t_distribution<t_real>>;

}
#endif
