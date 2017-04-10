/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

#include "../math/powder.h"
#include <iostream>

int main()
{
	tl::Powder<int> powder;

	powder.AddPeak(1,1,0);
	powder.AddPeak(1,-1,0);
	powder.AddPeak(-1,-1,0);

	powder.AddPeak(1,0,0);
	powder.AddPeak(0,1,0);

	powder.AddPeak(2,1,0);
	powder.AddPeak(-4,0,0);

	std::cout << powder << std::endl;

	return 0;
}
