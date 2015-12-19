/**
 * Debug helpers
 * @author tweber
 * @date 12-sep-2014
 * @license GPLv2 or GPLv3
 */

#ifndef __TL_DEBUG_H__
#define __TL_DEBUG_H__

#if BOOST_VERSION >= 105700
	#include <boost/type_index.hpp>
#endif

namespace tl{

#if BOOST_VERSION >= 105700
template<typename T>
std::string get_typename(bool bFull=1)
{
	boost::typeindex::type_index idx;

	if(bFull)
		idx = boost::typeindex::type_id_with_cvr<T>();
	else
		idx = boost::typeindex::type_id<T>();

	return idx.pretty_name();
}
#else
template<typename T>
std::string get_typename(bool bFull=1)
{
	return std::string("n/a");
}
#endif



extern void log_backtrace();

}
#endif
