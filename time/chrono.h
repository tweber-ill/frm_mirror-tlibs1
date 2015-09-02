/*
 * Chrono helpers
 * @author tweber
 * @date aug-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __CHRONO_HELPERS_H__
#define __CHRONO_HELPERS_H__

#include <chrono>


namespace tl {

template<typename T=double>
T epoch()
{
	namespace ch = std::chrono;

	using t_dur = ch::system_clock::duration;
	using t_per = t_dur::period;

	const T tSecsPerPeriod = T(t_per::num)/T(t_per::den);
	const t_dur durNow = ch::system_clock::now().time_since_epoch();

	return tSecsPerPeriod * T(durNow.count());
}

}
#endif
