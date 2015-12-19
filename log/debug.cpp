/**
 * Debug helpers
 * @author tweber
 * @date 12-sep-2014
 * @license GPLv2 or GPLv3
 */
#include "debug.h"


#ifdef DEBUG
#include <execinfo.h>

void tl::log_backtrace()
{
	void *pBtrBuf[512];
	int iBtr = backtrace(pBtrBuf, 512);
	char **ppcBtr = backtrace_symbols(pBtrBuf, iBtr);
	for(int _iBtr=0; _iBtr<iBtr; ++_iBtr)
		log_debug("Backtrace ", _iBtr, ": ", ppcBtr[_iBtr]);
}
#else
void tl::log_backtrace() {}
#endif
