/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

#include "../math/numint.h"
#include <iostream>

int main()
{
	double dRes = tl::newton<double>(
		[](double x)->double { return (x+2.)*(x+2.) + (x+2.); },
		[](double x)->double { return 2.*x + 5.; }, 
		0., 512, 0.00001);
	std::cout << dRes << std::endl;

	return 0;
}
