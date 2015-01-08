/*
 * abstract function model base classes
 * @author tweber
 * @date 2013
 */

#ifndef __FUNC_MOD_H__
#define __FUNC_MOD_H__

// parametric function
template<class t_vec = boost::numeric::ublas::vector<double>, class T = typename t_vec::value_type>
class FunctionModel_param_gen
{
	public:
		virtual ~FunctionModel_param_gen() {}

		// t = 0..1
		virtual t_vec operator()(T t) const = 0;
		virtual const char* GetModelName() const = 0;
};

// explicit function
template<class T = double> class FunctionModel_gen
{
	public:
		virtual ~FunctionModel_gen() {}

		virtual T operator()(T x) const = 0;
		virtual const char* GetModelName() const = 0;
};

// interface for n dimensional function
template<class T = double> class FunctionModel_nd_gen
{
	protected:

	public:
		virtual ~FunctionModel_nd_gen() {}

		virtual T operator()(const T* px) const = 0;
		virtual const char* GetModelName() const = 0;
};



#include <boost/numeric/ublas/vector.hpp>

typedef class FunctionModel_param_gen<boost::numeric::ublas::vector<double>> FunctionModel_param;
typedef class FunctionModel_gen<double> FunctionModel;
typedef class FunctionModel_nd_gen<double> FunctionModel_nd;

#endif
