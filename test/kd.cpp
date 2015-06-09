#include <iostream>
#include <vector>
#include <list>
#include "../math/kd.h"

int main()
{
	std::list<std::vector<int>> lst =
	{
		{1,9},
		{2,8},
		{3,7},
		{4,6},
		{5,5},
		{6,4},
		{7,3},
		{8,2},
		{9,1},
		{10,0},
	};

	Kd<int> kd(lst);
	kd.GetRootNode()->print(std::cout);
}
