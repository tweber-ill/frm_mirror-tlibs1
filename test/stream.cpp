// gcc -o stream stream.cpp ../helper/log.cpp -lstdc++ -std=c++11 -lboost_iostreams

#include "../file/comp.h"
#include <memory>

int main()
{
	std::ifstream ifstr("/home/tw/tmp/todo.txt");
	std::istream* pIstr = tl::create_autodecomp_istream(ifstr);
	std::unique_ptr<std::istream> ptrIstr(pIstr);

	while(!pIstr->eof())
	{
		std::string strLine;
		std::getline(*pIstr, strLine);
		std::cout << "Line: " << strLine << std::endl;
	}

	return 0;
}
