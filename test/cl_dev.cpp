// gcc -o cl_dev cl_dev.cpp ../cl/cl.cpp -std=c++11 -lstdc++ -lOpenCL

#include "../cl/cl.h"
#include <iostream>

int main()
{
	cl::Platform plat;
	cl::Device dev;
	if(!tl::get_best_cl_dev(plat, dev, 1))
	{
		std::cerr << "Cannot get devices." << std::endl;
		return -1;
	}

	std::cout << "Platform: " << plat.getInfo<CL_PLATFORM_NAME>() << std::endl;
	std::cout << "Device: " << dev.getInfo<CL_DEVICE_NAME>() << std::endl;

	return 0;
}
