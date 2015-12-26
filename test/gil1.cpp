// gcc -I . -o gil1 test/gil1.cpp log/log.cpp -std=c++11 -lstdc++ -ljpeg -lpng

#include "gfx/gil.h"
#include "log/debug.h"
#include <iostream>

namespace gil = tl::gil;

int main()
{
	auto pimg = tl::load_image<gil::rgb8_image_t>("/home/tweber/Pictures/rac.jpg");
	if(!pimg)
		return -1;

	std::size_t iChan = gil::view(*pimg).num_channels();
	std::cout << iChan << " channels." << std::endl;

	tl::for_each_in_view(gil::view(*pimg),
		[](gil::rgb8_pixel_t& pix)
		{
			std::swap(pix[0], pix[2]);
		});

	if(!tl::save_image("/tmp/tst.png", pimg/*.get()*/))
		return -1;

	return 0;
}
