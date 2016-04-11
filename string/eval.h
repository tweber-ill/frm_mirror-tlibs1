/**
 * minimalistic expression evaluator
 * @author tweber
 * @date apr-2016
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_EVAL_H__
#define __TLIBS_EVAL_H__

#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/phoenix.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <string>
#include <utility>
#include <unordered_map>
#include <type_traits>
#include "../log/log.h"
#include "string.h"

namespace tl
{
	// real functions with one parameter
	template<class t_str, class t_val,
		typename std::enable_if<std::is_floating_point<t_val>::value>::type* =nullptr>
	t_val call_func1(const t_str& strName, t_val t)
	{
		//std::cout << "calling " << strName << " with arg " << t << std::endl;
		static const std::unordered_map</*t_str*/std::string, t_val(*)(t_val)> s_funcs =
		{
			{ "sin", std::sin }, { "cos", std::cos }, { "tan", std::tan },
			{ "asin", std::asin }, { "acos", std::acos }, { "atan", std::atan },
			{ "sinh", std::sinh }, { "cosh", std::cosh }, { "tanh", std::tanh },
			{ "asinh", std::asinh }, { "acosh", std::acosh }, { "atanh", std::atanh },

			{ "sqrt", std::sqrt }, { "cbrt", std::cbrt },
			{ "exp", std::exp },
			{ "log", std::log }, { "log2", std::log2 }, { "log10", std::log10 },

			{ "erf", std::erf }, { "erfc", std::erfc }, { "erf_inv", boost::math::erf_inv },

			{ "round", std::round }, { "ceil", std::ceil }, { "floor", std::floor },
			{ "abs", std::abs },
		};

		return s_funcs.at(tl::wstr_to_str(strName))(t);
	}

	// real functions with two parameters
	template<class t_str, class t_val,
		typename std::enable_if<std::is_floating_point<t_val>::value>::type* =nullptr>
	t_val call_func2(const t_str& strName, t_val t1, t_val t2)
	{
		static const std::unordered_map</*t_str*/std::string, t_val(*)(t_val, t_val)> s_funcs =
		{
			{ "pow", std::pow }, { "atan2", std::atan2 },
			{ "mod", std::fmod },
		};

		return s_funcs.at(tl::wstr_to_str(strName))(t1, t2);
	}

	// real constants
	template<class t_str, class t_val,
		typename std::enable_if<std::is_floating_point<t_val>::value>::type* =nullptr>
	t_val get_const(const t_str& strName)
	{
		//std::cout << "requesting constant " << strName << std::endl;
		static const std::unordered_map</*t_str*/std::string, t_val> s_consts =
		{
			{ "pi", t_val(M_PI) }
		};

		return s_consts.at(tl::wstr_to_str(strName));
	}


	// alternative: int functions with one parameter
	template<class t_str, class t_val,
		typename std::enable_if<std::is_integral<t_val>::value>::type* =nullptr>
	t_val call_func1(const t_str& strName, t_val t)
	{
		static const std::unordered_map</*t_str*/std::string, t_val(*)(t_val)> s_funcs =
		{
			{ "abs", std::abs },
		};

		return s_funcs.at(tl::wstr_to_str(strName))(t);
	}

	// alternative: int functions with two parameters
	template<class t_str, class t_val,
		typename std::enable_if<std::is_integral<t_val>::value>::type* =nullptr>
	t_val call_func2(const t_str& strName, t_val t1, t_val t2)
	{
		static const std::unordered_map</*t_str*/std::string, std::function<t_val(t_val, t_val)>> s_funcs =
		{
			{ "pow", [t1, t2](t_val t1, t_val t2) -> t_val { return t_val(std::pow(t1, t2)); } },
			{ "mod", [t1, t2](t_val t1, t_val t2) -> t_val { return t1%t2; } },
		};

		return s_funcs.at(tl::wstr_to_str(strName))(t1, t2);
	}

	// alternative: int constants
	template<class t_str, class t_val,
		typename std::enable_if<std::is_integral<t_val>::value>::type* =nullptr>
	t_val get_const(const t_str& strName)
	{
		static const std::unordered_map</*t_str*/std::string, t_val> s_consts =
		{
		};

		return s_consts.at(tl::wstr_to_str(strName));
	}


	namespace qi = boost::spirit::qi;
	namespace asc = boost::spirit::ascii;
	namespace ph = boost::phoenix;

	template<class t_str, class t_val, class t_skip=asc::space_type>
	class ExprGrammar : public qi::grammar<
		typename t_str::const_iterator, t_val(), t_skip>
	{
		protected:
			using t_ch = typename t_str::value_type;
			using t_iter = typename t_str::const_iterator;
			using t_valparser = typename std::conditional<
				std::is_floating_point<t_val>::value,
				qi::real_parser<t_val>, qi::int_parser<t_val>>::type;

			qi::rule<t_iter, t_val(), t_skip> m_expr, m_term;
			qi::rule<t_iter, t_val(), t_skip> m_val, m_baseval, m_const, m_func;
			qi::rule<t_iter, t_val(), t_skip> m_pm, m_pm_opt;
			qi::rule<t_iter, t_str(), t_skip> m_ident;

		public:
			ExprGrammar() : ExprGrammar::base_type(m_expr)
			{
				// + or -
				m_expr = ((m_pm_opt >> m_term) [ qi::_val = qi::_1*qi::_2  ]
					>> *(m_pm >> m_term) [ qi::_val += qi::_1*qi::_2 ]
					);

				// * or /
				m_term = m_val [ qi::_val = qi::_1 ]
					>> *((t_ch('*') >> m_val) [ qi::_val *= qi::_1 ]
						| (t_ch('/') >> m_val) [ qi::_val /= qi::_1 ])
					| m_val [ qi::_val = qi::_1 ]
					;

				// + or -
				m_pm_opt = m_pm [ qi::_val = qi::_1 ]
					| qi::eps [ qi::_val = t_val(1) ]
					;
				m_pm = qi::char_(t_ch('+')) [ qi::_val = t_val(1) ]
					| qi::char_(t_ch('-')) [ qi::_val = t_val(-1) ]
					;

				// pow
				m_val = m_baseval [ qi::_val = qi::_1 ]
					>> *((t_ch('^') >> m_baseval)
					[ qi::_val = ph::bind([](t_val val1, t_val val2)->t_val
					{ return std::pow(val1, val2); }, qi::_val, qi::_1)]);

				m_baseval = t_valparser() | m_func | m_const
					| t_ch('(') >> m_expr >> t_ch(')')
					;

				// lazy evaluation of constants via phoenix bind
				m_const = m_ident [ qi::_val = ph::bind([](const t_str& strName) -> t_val
					{ return get_const<t_str, t_val>(strName); }, qi::_1) ];

				// lazy evaluation of functions via phoenix bind
				m_func = (m_ident >> t_ch('(') >> m_expr >> t_ch(',') >> m_expr >> t_ch(')'))
					[ qi::_val = ph::bind([](const t_str& strName, t_val val1, t_val val2) -> t_val
					{ return call_func2<t_str, t_val>(strName, val1, val2); },
					qi::_1, qi::_2, qi::_3) ]
					| (m_ident >> t_ch('(') >> m_expr >> t_ch(')'))
					[ qi::_val = ph::bind([](const t_str& strName, t_val val) -> t_val
					{ return call_func1<t_str, t_val>(strName, val); },
					qi::_1, qi::_2) ];

				m_ident = qi::lexeme[qi::char_("A-Za-z_") >> *qi::char_("0-9A-Za-z_")];
			}

			~ExprGrammar() {}
	};

	template<class t_str=std::string, class t_val=double>
	std::pair<bool, t_val> eval_expr(const t_str& str) noexcept
	{
		try
		{
			t_val valRes;
			ExprGrammar<t_str, t_val> gram;
			bool bOk = qi::phrase_parse(str.begin(), str.end(),
				gram, asc::space, valRes);

			return std::make_pair(bOk, valRes);
		}
		catch(const std::exception& ex)
		{
			log_err("Parsing failed with error: ", ex.what(), ".");
			return std::make_pair(false, t_val(0));
		}
	}
}

#endif
