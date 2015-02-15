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
#include <iostream>
#include "../string/string.h"

namespace tl{

// psi files
class FilePsi
{
	public:
		typedef std::unordered_map<std::string, std::string> t_mapParams;
		t_mapParams m_mapParams;

		typedef std::vector<std::string> t_vecColNames;
		t_vecColNames m_vecColNames;

		typedef std::vector<double> t_vecVals;
		t_vecVals m_vecVals;

		std::vector<t_vecVals> m_vecData;

	protected:
		void ReadData(std::istream& istr);

	public:
		FilePsi() = default;
		virtual ~FilePsi() = default;

		bool Load(const char* pcFile);
		void PrintParams(std::ostream& ostr) const;

		const std::string& GetColName(std::size_t iCol) const { return m_vecColNames[iCol]; }
		std::size_t GetColCount() const { return m_vecColNames.size(); }
		const std::vector<double>& GetCol(std::size_t iCol) const { return m_vecData[iCol]; }
		const std::vector<double>& GetCol(const std::string& strName) const;
};

}

#endif
