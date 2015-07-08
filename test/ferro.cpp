// clang -o ferro test/ferro.cpp -lstdc++ -lm -std=c++11

#include "../math/disp.h"
#include "../gfx/gnuplot.h"
#include <iostream>

int main()
{
	std::vector<double> vecq = tl::linspace(0., 5., 128);
	std::vector<double> vecE;
	vecE.reserve(vecq.size());

	for(double q : vecq)
	{
		double E = tl::ferromag({
			{1., tl::make_vec({1., 0., 0.})},
			{1., tl::make_vec({-1., 0., 0.})},
			{1., tl::make_vec({0., 1., 0.})},
			{1., tl::make_vec({0., -1., 0.})},
			{1., tl::make_vec({0., 0., 1.})},
			{1., tl::make_vec({0., 0., -1.})},
				},
			tl::make_vec({q, 0., 0.}), 1.);

		vecE.push_back(E);
		//std::cout << q << " " << E << std::endl;
	}

	tl::GnuPlot plt;
	plt.Init();
	plt.SimplePlot(vecq, vecE, std::vector<double>(), std::vector<double>(), tl::STYLE_LINES_SOLID);
	return 0;
}
