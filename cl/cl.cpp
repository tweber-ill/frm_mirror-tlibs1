/*
 * CL helpers
 * @author tweber
 * @date 20-jan-2015
 * @copyright GPLv2 or GPLv3
 */

#include "cl.h"
#include <algorithm>

namespace tl {

bool get_best_cl_dev(cl::Platform& platRet, cl::Device& devRet, bool bWantDouble)
{
	std::vector<cl::Platform> vecPlat;
	cl::Platform::get(&vecPlat);

	if(vecPlat.size() == 0)
		return false;


	struct _Dev
	{
		cl::Platform* pPlat = nullptr;
		cl::Device dev;

		cl_device_type devtype;
	};

	std::vector<_Dev> vecAllDevs;

	for(cl::Platform& plat : vecPlat)
	{
		std::vector<cl::Device> vecDevs;
		plat.getDevices(CL_DEVICE_TYPE_ALL, &vecDevs);

		for(cl::Device& dev : vecDevs)
		{
			_Dev _dev;
			_dev.pPlat = &plat;
			_dev.dev = dev;
			_dev.devtype = dev.getInfo<CL_DEVICE_TYPE>();

			std::string strExtensions = dev.getInfo<CL_DEVICE_EXTENSIONS>();

			if(bWantDouble)
			{
				bool bHasDouble = (strExtensions.find("cl_khr_fp64") != std::string::npos);
				if(!bHasDouble) continue;
			}

			vecAllDevs.push_back(_dev);
		}
	}

	if(vecAllDevs.size() == 0)
		return false;

	std::sort(vecAllDevs.begin(), vecAllDevs.end(),
		[](const _Dev& dev1, const _Dev& dev2) -> bool
		{
			int (*get_device_score)(cl_device_type ty) = [](cl_device_type ty) -> int
			{
				int iScore = 0;

				if(ty & CL_DEVICE_TYPE_GPU)
					iScore += 1000;
				if(ty & CL_DEVICE_TYPE_ACCELERATOR)
					iScore += 100;
				if(ty & CL_DEVICE_TYPE_CPU)
					iScore += 10;

				return iScore;
			};

			int iScore1 = get_device_score(dev1.devtype);
			int iScore2 = get_device_score(dev2.devtype);

			return iScore1 >= iScore2;
		});

	platRet = *vecAllDevs[0].pPlat;
	devRet = vecAllDevs[0].dev;

	return true;
}

}
