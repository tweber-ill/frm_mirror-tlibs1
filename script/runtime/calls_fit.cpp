/*
 * external fit functions
 * @author tweber
 * @date jan 2014
 */

#include "calls_fit.h"
#include "../calls.h"


// --------------------------------------------------------------------------------
// fitting

static Symbol* fkt_fit(const std::vector<Symbol*>& vecSyms,
			ParseInfo& info, SymbolTable* pSymTab)
{
	// TODO
	return 0;
}


// --------------------------------------------------------------------------------



extern void init_ext_fit_calls()
{
	t_mapFkts mapFkts =
	{
		t_mapFkts::value_type("fit", fkt_fit),
//		t_mapFkts::value_type("fit_sin", fkt_fit_sin),
//		t_mapFkts::value_type("fit_gauss", fkt_fit_sin),

//		t_mapFkts::value_type("spline", fkt_spline),
//		t_mapFkts::value_type("bezier", fkt_bezier),
	};

	add_ext_calls(mapFkts);
}
