/*
 * CL helpers
 * @author tweber
 * @date 20-jan-2015
 * @copyright GPLv2 or GPLv3
 */

#ifndef __CL_WRAP_H__
#define __CL_WRAP_H__

#include <vector>
#include <CL/cl.hpp>

namespace tl {

extern bool get_best_cl_dev(cl::Platform& plat, cl::Device& dev, bool bWantDouble=1);

}

#endif
