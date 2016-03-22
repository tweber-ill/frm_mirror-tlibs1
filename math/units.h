/*
 * wrapper for boost.units
 * @author Tobias Weber
 * @date dec-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_UNITS__
#define __TLIBS_UNITS__


#include <boost/units/unit.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/dimensionless_quantity.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/physical_dimensions.hpp>

#include <boost/units/systems/si.hpp>
#include <boost/units/systems/angle/degrees.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/neutron_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>


namespace tl {

namespace units = boost::units;
namespace co = boost::units::si::constants::codata;


// general quantities
template<class Sys, class T=double> using t_length =
	units::quantity<units::unit<units::length_dimension, Sys>, T>;
template<class Sys, class T=double> using t_momentum =
	units::quantity<units::unit<units::momentum_dimension, Sys>, T>;
template<class Sys, class T=double> using t_wavenumber =
	units::quantity<units::unit<units::wavenumber_dimension, Sys>, T>;
template<class Sys, class T=double> using t_velocity =
	units::quantity<units::unit<units::velocity_dimension, Sys>, T>;
template<class Sys, class T=double> using t_frequency =
	units::quantity<units::unit<units::frequency_dimension, Sys>, T>;
template<class Sys, class T=double> using t_energy =
	units::quantity<units::unit<units::energy_dimension, Sys>, T>;
template<class Sys, class T=double> using t_angle =
	units::quantity<units::unit<units::plane_angle_dimension, Sys>, T>;
template<class Sys, class T=double> using t_temperature =
	units::quantity<units::unit<units::temperature_dimension, Sys>, T>;
template<class Sys, class T=double> using t_mass =
	units::quantity<units::unit<units::mass_dimension, Sys>, T>;
template<class Sys, class T=double> using t_time =
	units::quantity<units::unit<units::time_dimension, Sys>, T>;
template<class Sys, class T=double> using t_flux =
	units::quantity<units::unit<units::magnetic_flux_density_dimension, Sys>, T>;
template<class Sys, class T=double> using t_area =
	units::quantity<units::unit<units::area_dimension, Sys>, T>;
template<class Sys, class T=double> using t_volume =
	units::quantity<units::unit<units::volume_dimension, Sys>, T>;

template<class Sys, class T=double> using t_length_inverse =
	units::quantity<units::unit<units::derived_dimension<units::length_base_dimension, -1>::type, Sys>, T>;
template<class Sys, class T=double> using t_length_square =
	units::quantity<units::unit<units::derived_dimension<units::length_base_dimension, 2>::type, Sys>, T>;
template<class Sys, class T=double> using t_momentum_square =
	units::quantity<units::unit<units::derived_dimension<units::momentum_dimension, 2>::type, Sys>, T>;
template<class Sys, class T=double> using t_action =
	units::quantity<units::unit<typename units::derived_dimension
	<units::mass_base_dimension,1, units::length_base_dimension,2, units::time_base_dimension,-1>::type, Sys>, T>;
template<class Sys, class T=double> using t_dimensionless =
	units::quantity<units::unit<units::dimensionless_type, Sys>, T>;

// synonyms
template<class Sys, class T=double> using t_freq = t_frequency<Sys, T>;
template<class Sys, class T=double> using t_temp = t_temperature<Sys,T>;


// constants
template<class Y=double> t_energy<units::si::system, Y> get_one_meV()
{ return Y(1e-3) * Y(co::e/units::si::coulombs)*units::si::coulombs*units::si::volts; }
template<class Y=double> t_length<units::si::system, Y> get_one_angstrom()
{ return Y(1e-10) * units::si::meters; }



// si quantities
typedef units::quantity<units::si::plane_angle> angle;
typedef units::quantity<units::si::wavenumber> wavenumber;
typedef units::quantity<units::si::energy> energy;
typedef units::quantity<units::si::momentum> momentum;
typedef units::quantity<units::si::velocity> velocity;
typedef units::quantity<units::si::length> length;
typedef decltype(1./length()) inv_length;
typedef units::quantity<units::si::area> area;
typedef units::quantity<units::si::time> time;
typedef units::quantity<units::si::magnetic_flux_density> flux;
typedef units::quantity<units::si::frequency> frequency;
typedef units::quantity<units::si::temperature> temperature;
typedef units::quantity<units::si::mass> mass;



// synonyms
typedef frequency freq;
typedef temperature temp;


// TODO: make C++14 template variables out of these
static const length meters = 1.*units::si::meters;
static const flux teslas = 1.*units::si::teslas;
static const time seconds = 1.*units::si::seconds;
static const angle radians = 1.*units::si::radians;
static const temp kelvins = 1.*units::si::kelvins;
static const mass amu = co::m_u;
static const area barns = 1e-28 * units::si::meters*units::si::meters;

static const energy one_meV = 1e-3 * co::e * units::si::volts;
static const energy one_eV = co::e * units::si::volts;
static const length angstrom = 1e-10 * meters;
static const length cm = meters/100.;
static const time ps = 1e-12 * seconds;

// synonyms
static const temp kelvin = kelvins;
static const length meter = meters;
static const time second = seconds;
static const energy meV = one_meV;
static const flux tesla = teslas;
static const area barn = barns;



// helper functions
template<class t_quant, class t_quant_sq>
t_quant my_units_sqrt(const t_quant_sq& val)
{
	t_quant one_quant;
	t_quant_sq one_quant_sq;

	using Y = typename t_quant::value_type;
	Y valsq = Y(val / one_quant_sq);
	return std::sqrt(valsq) * one_quant;
}

}
#endif
