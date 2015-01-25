#include "../string/string.h"

int main()
{
	std::cout << tl::var_to_str<int>(1234567890, 10, 3) << std::endl;

	return 0;
}
