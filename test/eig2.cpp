// gcc -I/usr/include/lapacke -o eig2 test/eig2.cpp math/linalg2.cpp log/log.cpp -lstdc++ -lm -llapacke -llapack -std=c++11

#include "../math/linalg.h"
#include "../math/linalg2.h"

namespace ublas = boost::numeric::ublas;

int main()
{
	ublas::matrix<double> M = tl::make_mat({{1.1,-2.2,3.3},
						{-2.2,5.5,8.8},
						{3.3,8.8,-9.9}});
	//ublas::matrix<double> M = tl::make_mat({{123.4,0},{0,567.8}});

	ublas::matrix<double> M_org = M;
	std::cout << M << std::endl;

	std::vector<ublas::vector<double>> evecs;
	std::vector<double> evals;
	tl::eigenvec_sym<double>(M, evecs, evals);
	for(int i=0; i<evals.size(); ++i)
		std::cout << "eval: " << evals[i] <<
		", evec: " << (evecs[i]/ublas::norm_2(evecs[i])) <<
		", len: " << ublas::norm_2(evecs[i]) << std::endl;
	std::cout << std::endl;

	std::vector<ublas::vector<double>> evecs2;
	std::vector<double> evals2;
	tl::eigenvec_sym_simple(M, evecs2, evals2);
	for(int i=0; i<evals2.size(); ++i)
		std::cout << "eval: " << evals2[i] <<
		", evec: " << (evecs2[i]/ublas::norm_2(evecs2[i])) <<
		", len: " << ublas::norm_2(evecs2[i]) << std::endl;
	std::cout << std::endl;

	return 0;
}
