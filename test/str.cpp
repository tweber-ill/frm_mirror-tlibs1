#include "../string/string.h"

int main()
{
	std::cout << tl::var_to_str<int>(1234567890, 10, 3) << std::endl;

	std::vector<int> vec({1,2,3,4,5,6,7,8});
	std::cout << tl::cont_to_str<decltype(vec)>(vec);
	std::cout << std::endl;

	return 0;
}
