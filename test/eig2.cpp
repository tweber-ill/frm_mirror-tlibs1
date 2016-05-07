// gcc -I/usr/include/lapacke -o eig2 test/eig2.cpp math/linalg2.cpp log/log.cpp -lstdc++ -lm -llapacke -llapack -std=c++11

#define TLIBS_INC_HDR_IMPLS
#include "../math/linalg.h"
#include "../math/linalg2.h"

namespace ublas = boost::numeric::ublas;
using T = double;

int main()
{
	ublas::matrix<T> M = tl::make_mat({
		{1.1,-2.2,3.3},
		{-2.2,5.5,8.8},
		{3.3,8.8,-9.9}	});
	//ublas::matrix<T> M = tl::make_mat({{123.4,0},{0,567.8}});

	ublas::matrix<T> M_org = M;
	std::cout << M << std::endl;

	std::vector<ublas::vector<T>> evecs;
	std::vector<T> evals;
	tl::eigenvec_sym<T>(M, evecs, evals);
	for(int i=0; i<evals.size(); ++i)
		std::cout << "eval: " << evals[i] <<
		", evec: " << (evecs[i]/ublas::norm_2(evecs[i])) <<
		", len: " << ublas::norm_2(evecs[i]) << std::endl;
	std::cout << std::endl;
	for(int i=0; i<evals.size(); ++i)
		std::cout <<
		"val*vec: " << evals[i]*evecs[i] <<
		"\nmat*vec:" << ublas::prod(M, evecs[i]) << std::endl;
	std::cout << std::endl;

	std::vector<ublas::vector<T>> evecs2;
	std::vector<T> evals2;
	tl::eigenvec_sym_simple(M, evecs2, evals2);
	for(int i=0; i<evals2.size(); ++i)
		std::cout << "eval: " << evals2[i] <<
		", evec: " << (evecs2[i]/ublas::norm_2(evecs2[i])) <<
		", len: " << ublas::norm_2(evecs2[i]) << std::endl;
	std::cout << std::endl;

	std::vector<ublas::vector<T>> evecs2_r, evecs2_i;
	std::vector<T> evals2_r, evals2_i;
	tl::eigenvec(M, evecs2_r, evecs2_i, evals2_r, evals2_i);
	for(int i=0; i<evals2_r.size(); ++i)
		std::cout << "eval r: " << evals2_r[i] <<
		", evec r: " << evecs2_r[i] << std::endl;
	for(int i=0; i<evals2_i.size(); ++i)
		std::cout << "eval i: " << evals2_i[i] <<
		", evec i: " << evecs2_i[i] << std::endl;
	std::cout << std::endl;


	// ----------------------------------------------------------------

	std::vector<ublas::vector<std::complex<T>>> evecs_c;
	std::vector<T> evals_c;
	ublas::matrix<std::complex<T>> Mc = tl::make_mat<ublas::matrix<std::complex<T>>>({
		{std::complex<T>(1., 0.), std::complex<T>(3., 1.)},
		{std::complex<T>(3., -1.), std::complex<T>(3., 0.)}
	});
	std::cout << Mc << std::endl;

	tl::eigenvec_herm<T>(Mc, evecs_c, evals_c);
	for(int i=0; i<evals_c.size(); ++i)
		std::cout << "eval: " << evals_c[i] <<
		", evec: " << evecs_c[i] << std::endl;
	std::cout << std::endl;
	for(int i=0; i<evals_c.size(); ++i)
		std::cout <<
		"val*vec: " << evals_c[i]*evecs_c[i] <<
		"\nmat*vec:" << ublas::prod(Mc, evecs_c[i]) << std::endl;
	std::cout << std::endl;


	std::vector<ublas::vector<std::complex<T>>> evecs_c2;
	std::vector<std::complex<T>> evals_c2;

	tl::eigenvec_cplx<T>(Mc, evecs_c2, evals_c2);
	for(int i=0; i<evals_c2.size(); ++i)
		std::cout << "eval: " << evals_c2[i] <<
		", evec: " << evecs_c2[i] << std::endl;
	std::cout << std::endl;

	return 0;
}
