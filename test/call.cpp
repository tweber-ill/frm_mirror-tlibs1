/**
 * tlibs test file
 * @author Tobias Weber <tobias.weber@tum.de>
 * @license GPLv2 or GPLv3
 */

#include <iostream>
#include <vector>
#include <array>
#include "../helper/traits.h"


int tst(int a, int b)
{
	std::cout << "a = " << a << ", b = " << b << std::endl;

	return a+b;
}

int main()
{
	std::vector<int> vecArgs = {1, 2};
	std::array<int,2> arrArgs = {1, 2};

	tl::call<2>(tst, vecArgs);
	tl::call<2>(tst, arrArgs);

	return 0;
}
