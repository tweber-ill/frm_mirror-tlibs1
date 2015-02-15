/*
 * Load instrument-specific data file
 * @author tweber
 * @date feb-2015
 * @copyright GPLv2 or GPLv3
 */

#include "loadinstr.h"
#include "../helper/log.h"
#include "../file/file.h"
#include <fstream>

namespace tl{

void FilePsi::ReadData(std::istream& istr)
{
	// header
	std::string strHdr;
	std::getline(istr, strHdr);
	tl::get_tokens<std::string, std::string, t_vecColNames>(strHdr, " \t", m_vecColNames);

	m_vecData.resize(m_vecColNames.size());

	// data
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		if(strLine.length() == 0)
			continue;

		std::vector<double> vecToks;
		tl::get_tokens<double, std::string>(strLine, " \t", vecToks);

		if(vecToks.size() != m_vecColNames.size())
		{
			tl::log_warn("Loader: Line size mismatch.");

			// add zeros
			while(m_vecColNames.size() > vecToks.size())
				vecToks.push_back(0.);
		}

		for(std::size_t iTok=0; iTok<vecToks.size(); ++iTok)
			m_vecData[iTok].push_back(vecToks[iTok]);
	}
}

bool FilePsi::Load(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return false;

	skip_after_line<char>(ifstr, "VVVV", true);

	while(!ifstr.eof())
	{
		std::string strLine;
		std::getline(ifstr, strLine);

		std::pair<std::string, std::string> pairLine =
				tl::split_first<std::string>(strLine, ":", 1);
		if(pairLine.first == "DATA_")
			ReadData(ifstr);
		else if(pairLine.first == "")
			continue;
		else
		{
			t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

			if(iter == m_mapParams.end())
				m_mapParams.insert(pairLine);
			else
				iter->second += "\n" + pairLine.second;
		}
	}

	return true;
}

const std::vector<double>& FilePsi::GetCol(const std::string& strName) const
{
	static const std::vector<double> vecNull;

	for(std::size_t i=0; i<m_vecColNames.size(); ++i)
	{
		if(m_vecColNames[i] == strName)
			return m_vecData[i];
	}

	return vecNull;
}

void FilePsi::PrintParams(std::ostream& ostr) const
{
	for(const t_mapParams::value_type& val : m_mapParams)
	{
		std::cout << "Param: " << val.first
					<< ", Val: " << val.second << "\n";
	}
}

}




/*
// test
// gcc -o 0 file/loadinstr.cpp helper/log.cpp string/string.cpp -std=c++11 -lstdc++
int main()
{
	tl::FilePsi psi;
	if(!psi.Load("/home/tw/Measurements/Eiger2014/eiger2014n001051.scn"))
	{
		tl::log_err("Cannot load data file.");
		return -1;
	}

	psi.PrintParams(std::cout);

	for(std::size_t iCol=0; iCol<psi.GetColCount(); ++iCol)
		std::cout << "Column: " << psi.GetColName(iCol) << std::endl;

	std::cout << "Cnts: ";
	for(std::size_t i=0; i<psi.GetCol("CNTS").size(); ++i)
		std::cout << psi.GetCol("CNTS")[i] << ", ";
	std::cout << std::endl;

	return 0;
}
*/
