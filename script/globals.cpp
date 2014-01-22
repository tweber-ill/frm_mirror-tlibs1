/*
 * global symbols
 * @author tweber
 * @date 2013-2014
 */

#include "globals.h"
#include "runtime/calls_plot.h"
#include "runtime/calls_math.h"
#include "runtime/calls_fit.h"
#include "runtime/calls_file.h"
#include "runtime/calls_thread.h"
#include "helper/neutrons.hpp"
#include <boost/units/systems/si/codata/electron_constants.hpp>

const char* g_pcVersion = "Hermelin script interpreter, version 0.4";

static inline void init_funcs()
{
	init_ext_thread_calls();
	init_ext_file_calls();
	init_ext_math_calls();
	init_ext_fit_calls();
	init_ext_plot_calls();
}


static inline void init_physics_const(SymbolTable *pSymTab)
{
	// hbar in eVs
	pSymTab->InsertSymbol("hbar_eVs", new SymbolDouble(co::hbar / one_eV / units::si::second));
	// hbar in Js
	pSymTab->InsertSymbol("hbar", new SymbolDouble(co::hbar / units::si::joule / units::si::second));
	// neutron mass
	pSymTab->InsertSymbol("m_n", new SymbolDouble(co::m_n / units::si::kilogram));
	// atomic mass unit
	pSymTab->InsertSymbol("m_u", new SymbolDouble(co::m_u / units::si::kilogram));
	// electron mass
	pSymTab->InsertSymbol("m_e", new SymbolDouble(co::m_e / units::si::kilogram));
	// Boltzmann const
	pSymTab->InsertSymbol("k_B", new SymbolDouble(co::k_B * units::si::kelvin/units::si::joules));
	// Boltzmann const in eV/K
	pSymTab->InsertSymbol("k_B_eVperK", new SymbolDouble(co::k_B * units::si::kelvin/one_eV));
	// Avogadro const
	pSymTab->InsertSymbol("N_A", new SymbolDouble(co::N_A * units::si::moles));
	// speed of light
	pSymTab->InsertSymbol("c_0", new SymbolDouble(co::c / units::si::meters*units::si::seconds));
}

extern void init_global_syms(SymbolTable *pSymTab)
{
	init_funcs();
	init_physics_const(pSymTab);

	pSymTab->InsertSymbol("pi", new SymbolDouble(M_PI));
}
