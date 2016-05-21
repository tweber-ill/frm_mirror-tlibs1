#include <iostream>
#include "../math/term.h"

int main()
{
	double S,L,J;
	std::uint16_t l, es;

	std::cout << "l = "; std::cin >> l;
	std::cout << "Number of electrons: "; std::cin >> es;

	std::tie(S,L,J) = tl::hund(l, es);
	std::string strTerm = tl::get_termsymbol(S,L,J);

	std::cout << strTerm << "\n";
	std::cout << "S = " << S << "\n";
	std::cout << "L = " << L << "\n";
	std::cout << "J = " << J << "\n";
	std::cout.flush();

	return 0;
}
