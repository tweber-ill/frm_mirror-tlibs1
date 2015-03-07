// gcc -o ostream ostream.cpp ../helper/log.cpp -lstdc++ -std=c++11 -lboost_iostreams

#include "../file/comp.h"
#include <memory>

int main()
{
	std::ofstream ofstr("tst.txt.gz");
	std::ostream *pOstr = tl::create_comp_ostream(ofstr, tl::Compressor::GZ);
	std::unique_ptr<std::ostream> ptrOstr(pOstr);

	for(int i=0; i<1000; ++i)
		(*pOstr) << i << "\n";

	return 0;
}
