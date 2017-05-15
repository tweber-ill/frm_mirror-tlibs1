/**
 * geometric primitives
 * @author: Tobias Weber <tobias.weber@tum.de>
 * @date: 14-may-2017
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_GEO_PRIM_H__
#define __TLIBS_GEO_PRIM_H__

#include "geo.h"


// TODO: Check polygon index order
namespace tl {

namespace ublas = boost::numeric::ublas;
namespace math = boost::math;


template<class t_vec = ublas::vector<double>>
class GeometricPrimitive
{
public:
	virtual std::size_t GetVertexCount() const = 0;
	virtual std::size_t GetPolyCount() const = 0;

	virtual const t_vec& GetVertex(std::size_t iVert) const = 0;
	virtual const std::vector<std::size_t>& GetPolyIndex(std::size_t iPoly) const = 0;

	virtual std::vector<t_vec> GetPoly(std::size_t iPoly) const
	{
		std::vector<t_vec> vecPoly;

		for(std::size_t iIdx : GetPolyIndex(iPoly))
			vecPoly.push_back(GetVertex(iIdx));

		return vecPoly;
	}
};


// ----------------------------------------------------------------------------


/**
 * Tetrahedron
 * see e.g.: https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec = ublas::vector<double>>
class Tetrahedron : public GeometricPrimitive<t_vec>
{
protected:
	// vertices
	std::vector<t_vec> m_vecVertices =
	{
		make_vec<t_vec>({  1,  1,  1 }),	// 0
		make_vec<t_vec>({  1, -1, -1 }),	// 1
		make_vec<t_vec>({ -1, -1,  1 }),	// 2
		make_vec<t_vec>({ -1,  1, -1 }),	// 3
	};

	// polygons
	std::vector<std::vector<std::size_t>> m_vecPolyIndices =
	{
		{ 0, 1, 2 }, { 0, 1, 3 },
		{ 0, 2, 3 }, { 1, 2, 3 },
	};

public:
	Tetrahedron() = default;
	~Tetrahedron() = default;

	virtual std::size_t GetVertexCount() const override
	{ return m_vecVertices.size(); }

	virtual std::size_t GetPolyCount() const override
	{ return m_vecPolyIndices.size(); }

	virtual const t_vec& GetVertex(std::size_t iVert) const override
	{ return m_vecVertices[iVert]; }

	virtual const std::vector<std::size_t>& GetPolyIndex(std::size_t iPoly) const override
	{ return m_vecPolyIndices[iPoly]; }
};


/**
 * Cube
 * see e.g.: https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec = ublas::vector<double>>
class Cube : public GeometricPrimitive<t_vec>
{
protected:
	// vertices
	std::vector<t_vec> m_vecVertices =
	{
		make_vec<t_vec>({  1,  1,  1 }),	// 0
		make_vec<t_vec>({  1,  1, -1 }),	// 1
		make_vec<t_vec>({  1, -1, -1 }),	// 2
		make_vec<t_vec>({  1, -1,  1 }),	// 3
		make_vec<t_vec>({ -1,  1,  1 }),	// 4
		make_vec<t_vec>({ -1,  1, -1 }),	// 5
		make_vec<t_vec>({ -1, -1, -1 }),	// 6
		make_vec<t_vec>({ -1, -1,  1 }),	// 7
	};

	// polygons
	std::vector<std::vector<std::size_t>> m_vecPolyIndices =
	{
		{ 0, 1, 2, 3 }, { 0, 3, 7, 4 }, { 3, 2, 6, 7 }, 
		{ 4, 5, 6, 7 }, { 1, 2, 6, 5 }, { 1, 0, 4, 5 },
	};

public:
	Cube() = default;
	~Cube() = default;

	virtual std::size_t GetVertexCount() const override
	{ return m_vecVertices.size(); }

	virtual std::size_t GetPolyCount() const override
	{ return m_vecPolyIndices.size(); }

	virtual const t_vec& GetVertex(std::size_t iVert) const override
	{ return m_vecVertices[iVert]; }

	virtual const std::vector<std::size_t>& GetPolyIndex(std::size_t iPoly) const override
	{ return m_vecPolyIndices[iPoly]; }
};


/**
 * Octahedron
 * see e.g.: https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec = ublas::vector<double>>
class Octahedron : public GeometricPrimitive<t_vec>
{
protected:
	// vertices
	std::vector<t_vec> m_vecVertices =
	{
		make_vec<t_vec>({  1,  0,  0 }),	// 0
		make_vec<t_vec>({ -1,  0,  0 }),	// 1
		make_vec<t_vec>({  0,  1,  0 }),	// 2
		make_vec<t_vec>({  0, -1,  0 }),	// 3
		make_vec<t_vec>({  0,  0,  1 }),	// 4
		make_vec<t_vec>({  0,  0, -1 }),	// 5
	};

	// polygons
	std::vector<std::vector<std::size_t>> m_vecPolyIndices =
	{
		{ 0, 2, 4 }, { 0, 2, 5 }, { 0, 3, 4 }, { 0, 3, 5 },
		{ 1, 2, 4 }, { 1, 2, 5 }, { 1, 3, 4 }, { 1, 3, 5 },
	};

public:
	Octahedron() = default;
	~Octahedron() = default;

	virtual std::size_t GetVertexCount() const override
	{ return m_vecVertices.size(); }

	virtual std::size_t GetPolyCount() const override
	{ return m_vecPolyIndices.size(); }

	virtual const t_vec& GetVertex(std::size_t iVert) const override
	{ return m_vecVertices[iVert]; }

	virtual const std::vector<std::size_t>& GetPolyIndex(std::size_t iPoly) const override
	{ return m_vecPolyIndices[iPoly]; }
};


/**
 * Icosahedron
 * see e.g.: https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec = ublas::vector<double>>
class Icosahedron : public GeometricPrimitive<t_vec>
{
protected:
	// Golden Ratio
	/*static constexpr*/ const typename t_vec::value_type s_g = 0.5 + 0.5*std::sqrt(5.);

	// vertices
	std::vector<t_vec> m_vecVertices =
	{
		make_vec<t_vec>({  0,  1,  s_g }),	// 0
		make_vec<t_vec>({  0,  1, -s_g }),	// 1
		make_vec<t_vec>({  0, -1,  s_g }),	// 2
		make_vec<t_vec>({  0, -1, -s_g }),	// 3

		make_vec<t_vec>({  1,  s_g,  0 }),	// 4
		make_vec<t_vec>({  1, -s_g,  0 }),	// 5
		make_vec<t_vec>({ -1,  s_g,  0 }),	// 6
		make_vec<t_vec>({ -1, -s_g,  0 }),	// 7

		make_vec<t_vec>({  s_g,  0,  1 }),	// 8
		make_vec<t_vec>({  s_g,  0, -1 }),	// 9
		make_vec<t_vec>({ -s_g,  0,  1 }),	// 10
		make_vec<t_vec>({ -s_g,  0, -1 }),	// 11
	};

	// polygons
	std::vector<std::vector<std::size_t>> m_vecPolyIndices =
	{
		{ 0, 10, 2 }, { 0, 2, 8 }, { 0, 8, 4 }, { 0,  4, 6 }, { 0,  6, 10 },	// upper cap
		{ 3,  5, 7 }, { 3, 9, 5 }, { 3, 1, 9 }, { 3, 11, 1 }, { 3,  7, 11 },	// lower cap
		{ 10, 7, 2 }, { 2, 5, 8 }, { 8, 9, 4 }, { 4,  1, 6 }, { 6, 11, 10 },	// sides
		{  7, 5, 2 }, { 5, 9, 8 }, { 9, 1, 4 }, { 1, 11, 6 }, { 11, 7, 10 },	// sides
	};

public:
	Icosahedron() = default;
	~Icosahedron() = default;

	virtual std::size_t GetVertexCount() const override
	{ return m_vecVertices.size(); }

	virtual std::size_t GetPolyCount() const override
	{ return m_vecPolyIndices.size(); }

	virtual const t_vec& GetVertex(std::size_t iVert) const override
	{ return m_vecVertices[iVert]; }

	virtual const std::vector<std::size_t>& GetPolyIndex(std::size_t iPoly) const override
	{ return m_vecPolyIndices[iPoly]; }
};


/**
 * Dodecahedron
 * see e.g.: https://en.wikipedia.org/wiki/Platonic_solid
 */
template<class t_vec = ublas::vector<double>>
class Dodecahedron : public GeometricPrimitive<t_vec>
{
protected:
	// Golden Ratio
	/*static constexpr*/ const typename t_vec::value_type s_g = 0.5 + 0.5*std::sqrt(5.);

	// vertices
	std::vector<t_vec> m_vecVertices =
	{
		make_vec<t_vec>({  0,  1./s_g,  s_g }),	// 0
		make_vec<t_vec>({  0,  1./s_g, -s_g }),	// 1
		make_vec<t_vec>({  0, -1./s_g,  s_g }),	// 2
		make_vec<t_vec>({  0, -1./s_g, -s_g }),	// 3

		make_vec<t_vec>({  1./s_g,  s_g,  0 }),	// 4
		make_vec<t_vec>({  1./s_g, -s_g,  0 }),	// 5
		make_vec<t_vec>({ -1./s_g,  s_g,  0 }),	// 6
		make_vec<t_vec>({ -1./s_g, -s_g,  0 }),	// 7

		make_vec<t_vec>({  s_g,  0,  1./s_g }),	// 8
		make_vec<t_vec>({  s_g,  0, -1./s_g }),	// 9
		make_vec<t_vec>({ -s_g,  0,  1./s_g }),	// 10
		make_vec<t_vec>({ -s_g,  0, -1./s_g }),	// 11

		make_vec<t_vec>({  1,  1,  1 }),	// 12
		make_vec<t_vec>({  1,  1, -1 }),	// 13
		make_vec<t_vec>({  1, -1, -1 }),	// 14
		make_vec<t_vec>({  1, -1,  1 }),	// 15
	
		make_vec<t_vec>({ -1,  1,  1 }),	// 16
		make_vec<t_vec>({ -1,  1, -1 }),	// 17
		make_vec<t_vec>({ -1, -1, -1 }),	// 18
		make_vec<t_vec>({ -1, -1,  1 }),	// 19
	};

	// polygons
	std::vector<std::vector<std::size_t>> m_vecPolyIndices =
	{
		{ 16, 10, 19,  2, 0 }, { 19, 7,  5, 15,  2 }, { 15, 8, 12, 0, 2 }, { 12, 4, 6, 16, 0 },	// top cap
		{ 18, 11, 17,  1, 3 }, { 17, 1, 13,  4,  6 }, { 13, 9, 14, 3, 1 }, { 14, 5, 7, 18, 3 },	// bottom cap
		{ 19, 10, 11, 18, 7 }, { 16, 6, 17, 11, 10 }, { 15, 5, 14, 9, 8 }, { 12, 8, 9, 13, 4 },	// sides
	};

public:
	Dodecahedron() = default;
	~Dodecahedron() = default;

	virtual std::size_t GetVertexCount() const override
	{ return m_vecVertices.size(); }

	virtual std::size_t GetPolyCount() const override
	{ return m_vecPolyIndices.size(); }

	virtual const t_vec& GetVertex(std::size_t iVert) const override
	{ return m_vecVertices[iVert]; }

	virtual const std::vector<std::size_t>& GetPolyIndex(std::size_t iPoly) const override
	{ return m_vecPolyIndices[iPoly]; }
};


// ----------------------------------------------------------------------------


/**
 * tessellated sphere
 */
template<class t_vec = ublas::vector<double>,
	template<class...> class t_underlying_solid = Icosahedron>
class TesselSphere : public t_underlying_solid<t_vec>
{
public:
	TesselSphere(typename t_vec::value_type dRad = 1.)
	{
		for(t_vec& vecVertex : t_underlying_solid<t_vec>::m_vecVertices)
		{
			vecVertex /= ublas::norm_2(vecVertex);
			vecVertex *= dRad;
		}
	}

	~TesselSphere() = default;
};

}
#endif
