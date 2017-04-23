/**
 * julia interface helpers
 * 
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 23-apr-2017
 * @license GPLv2 or GPLv3
 */

#ifndef __TL_JL_H__
#define __TL_JL_H__

#include <julia.h>


namespace tl
{

// ----------------------------------------------------------------------------
/**
 * Julia data type traits
 */
template<typename T> struct jl_traits {};

template<> struct jl_traits<int64_t>
{
	using value_type = int64_t;

	static jl_datatype_t* get_type() { return jl_int64_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_int64(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_int64(val); }
};

template<> struct jl_traits<uint64_t>
{
	using value_type = uint64_t;

	static jl_datatype_t* get_type() { return jl_uint64_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_uint64(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_uint64(val); }
};

template<> struct jl_traits<int32_t>
{
	using value_type = int32_t;

	static jl_datatype_t* get_type() { return jl_int32_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_int32(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_int32(val); }
};

template<> struct jl_traits<uint32_t>
{
	using value_type = uint32_t;

	static jl_datatype_t* get_type() { return jl_uint32_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_uint32(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_uint32(val); }
};

template<> struct jl_traits<int16_t>
{
	using value_type = int16_t;

	static jl_datatype_t* get_type() { return jl_int16_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_int16(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_int16(val); }
};

template<> struct jl_traits<uint16_t>
{
	using value_type = uint16_t;

	static jl_datatype_t* get_type() { return jl_uint16_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_uint16(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_uint16(val); }
};

template<> struct jl_traits<int8_t>
{
	using value_type = int8_t;

	static jl_datatype_t* get_type() { return jl_int8_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_int8(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_int8(val); }
};

template<> struct jl_traits<uint8_t>
{
	using value_type = uint8_t;

	static jl_datatype_t* get_type() { return jl_uint8_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_uint8(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_uint8(val); }
};

template<> struct jl_traits<float>
{
	using value_type = float;

	static jl_datatype_t* get_type() { return jl_float32_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_float32(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_float32(val); }
};

template<> struct jl_traits<double>
{
	using value_type = double;

	static jl_datatype_t* get_type() { return jl_float64_type; }
	static value_type unbox(jl_value_t *pVal) { return jl_unbox_float64(pVal); }
	static jl_value_t* box(value_type val) { return jl_box_float64(val); }
};
// ----------------------------------------------------------------------------


}
#endif
