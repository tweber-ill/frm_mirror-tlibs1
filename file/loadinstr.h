/*
 * Load instrument-specific data file
 * @author tweber
 * @date feb-2015
 * @copyright GPLv2 or GPLv3
 */

#ifndef __LOADINSTR_H__
#define __LOADINSTR_H__

#include <unordered_map>
#include <vector>
#include <array>
#include <iostream>
#include "../string/string.h"

namespace tl{

// psi & ill files
class FilePsi
{
	public:
		typedef std::unordered_map<std::string, std::string> t_mapParams;
		t_mapParams m_mapParams;

		// internal parameters in m_mapParams
		typedef std::unordered_map<std::string, double> t_mapIParams;
		t_mapIParams m_mapParameters, m_mapZeros, m_mapVariables, m_mapPosHkl, m_mapScanSteps;


		typedef std::vector<std::string> t_vecColNames;
		t_vecColNames m_vecColNames;

		typedef std::vector<double> t_vecVals;
		t_vecVals m_vecVals;

		std::vector<t_vecVals> m_vecData;

	protected:
		void ReadData(std::istream& istr);
		void GetInternalParams(const std::string& strAll, t_mapIParams& mapPara);

	public:
		FilePsi() = default;
		virtual ~FilePsi() = default;

		bool Load(const char* pcFile);

		void PrintParams(std::ostream& ostr) const;
		const t_mapParams& GetParams() const { return m_mapParams; }

		const std::string& GetColName(std::size_t iCol) const { return m_vecColNames[iCol]; }
		std::size_t GetColCount() const { return m_vecColNames.size(); }
		const std::vector<double>& GetCol(std::size_t iCol) const { return m_vecData[iCol]; }
		const std::vector<double>& GetCol(const std::string& strName) const;

	public:
		std::array<double, 3> GetSampleLattice() const;
		std::array<double, 3> GetSampleAngles() const;
		std::array<double, 2> GetMonoAnaD() const;

		std::array<bool, 3> GetScatterSenses() const;
		std::array<double, 3> GetScatterPlane0() const;
		std::array<double, 3> GetScatterPlane1() const;

		double GetKFix() const;

		std::array<double, 4> GetPosHKLE() const;	// zero pos.
		std::array<double, 4> GetStepsHKLE() const;	// scan steps
};

}

#endif
