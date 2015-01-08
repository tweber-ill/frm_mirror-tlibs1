/*
 * abstract function model classes
 * @author tweber
 * @date 2013
 */

#ifndef __FUNC_MOD_H__
#define __FUNC_MOD_H__

// parametric function
class FunctionModel_param
{
	public:
		virtual ~FunctionModel_param() {}

		// t = 0..1
		virtual boost::numeric::ublas::vector<double> operator()(double t) const = 0;
		virtual const char* GetModelName() const = 0;
};

#endif
