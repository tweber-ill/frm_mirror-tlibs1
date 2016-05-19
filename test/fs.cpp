// gcc -o fs test/fs.cpp -lstdc++ -lboost_system -lboost_filesystem -std=c++11

#include <iostream>
#include "../file/file.h"

int main()
{
	std::cout << tl::get_file_size(std::string("test/fs.cpp")) << std::endl;
	return 0;
}
