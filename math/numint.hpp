/*
 * simple numerical integration
 * @author tweber
 * @date june-2015
 */

#ifndef __NUMINT_H__
#define __NUMINT_H__

#include <functional>

namespace tl {


template<class R=double, class A=double>
R numint_trap(const std::function<R(A)>& fkt,
	A x0, A x1)
{
	return 0.5*(x1-x0) * (fkt(x0) + fkt(x1));
}

template<class R=double, class A=double>
R numint_trapN(const std::function<R(A)>& fkt,
	A x0, A x1, unsigned int N)
{
	const A xstep = (x1-x0)/A(N);

	R xsum = fkt(x0) + fkt(x1);
	for(unsigned int i=1; i<N; ++i)
		xsum += 2.*fkt(x0 + i*xstep);

	xsum *= 0.5*xstep;
	return xsum;
}


template<class R=double, class A=double>
R numint_rect(const std::function<R(A)>& fkt,
	A x0, A x1, unsigned int N)
{
	const A xstep = (x1-x0)/A(N);

	R xsum = R(0);
	for(unsigned int i=0; i<N; ++i)
		xsum += fkt(x0 + i*xstep);

	xsum *= xstep;
	return xsum;
}


template<class R=double, class A=double>
R numint_simp(const std::function<R(A)>& fkt,
	A x0, A x1)
{
	return (fkt(x0) + 4.*fkt(0.5*(x0+x1)) + fkt(x1)) * (x1-x0)/6.;
}

template<class R=double, class A=double>
R numint_simpN(const std::function<R(A)>& fkt,
	A x0, A x1, unsigned int N)
{
	const A xstep = (x1-x0)/A(N);
	R xsum = fkt(x0) + fkt(x1);

	for(unsigned int i=1; i<=N/2; ++i)
	{
		xsum += 2.*fkt(x0 + 2.*i*xstep);
		xsum += 4.*fkt(x0 + (2.*i - 1.)*xstep);
	}
	xsum -= 2.*fkt(x0 + 2.*N/2*xstep);

	xsum *= xstep/3.;
	return xsum;
}


}
#endif
