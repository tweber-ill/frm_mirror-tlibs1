/*
 * global symbols
 * @author tweber
 * @date 2013-2014
 */

#include "globals.h"
#include "helper/neutrons.hpp"

extern void init_global_syms(SymbolTable *pSymTab)
{
	pSymTab->InsertSymbol("pi", new SymbolDouble(M_PI));

	// hbar in eVs
	pSymTab->InsertSymbol("hbar_eVs", new SymbolDouble(co::hbar / one_eV / units::si::second));

	// hbar in Js
	pSymTab->InsertSymbol("hbar", new SymbolDouble(co::hbar / units::si::joule / units::si::second));

	// neutron mass
	pSymTab->InsertSymbol("m_n", new SymbolDouble(co::m_n / units::si::kilogram));
}
