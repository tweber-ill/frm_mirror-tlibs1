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

void FilePsi::GetInternalParams(const std::string& strAll, FilePsi::t_mapIParams& mapPara)
{
	std::vector<std::string> vecToks;
	tl::get_tokens<std::string, std::string>(strAll, ",\n", vecToks);
	//std::cout << strAll << std::endl;

	for(const std::string& strTok : vecToks)
	{
		std::pair<std::string, std::string> pair =
				tl::split_first<std::string>(strTok, "=", 1);

		if(pair.first == "")
			continue;

		double dVal = tl::str_to_var<double>(pair.second);
		mapPara.insert(t_mapIParams::value_type(pair.first, dVal));

		//std::cout << "Key: " << pair.first << ", Val: " << dVal << std::endl;
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

	t_mapParams::const_iterator iterParams = m_mapParams.find("PARAM"),
				iterZeros = m_mapParams.find("ZEROS"),
				iterVars = m_mapParams.find("VARIA"),
				iterPos = m_mapParams.find("POSQE"),
				iterSteps = m_mapParams.find("STEPS");

	if(iterParams!=m_mapParams.end()) GetInternalParams(iterParams->second, m_mapParameters);
	if(iterZeros!=m_mapParams.end()) GetInternalParams(iterZeros->second, m_mapParameters);
	if(iterVars!=m_mapParams.end()) GetInternalParams(iterVars->second, m_mapParameters);
	if(iterPos!=m_mapParams.end()) GetInternalParams(iterPos->second, m_mapParameters);
	if(iterSteps!=m_mapParams.end()) GetInternalParams(iterSteps->second, m_mapParameters);

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


std::array<double,3> FilePsi::GetSampleLattice() const
{
	t_mapIParams::const_iterator iterA = m_mapParameters.find("AS");
	t_mapIParams::const_iterator iterB = m_mapParameters.find("BS");
	t_mapIParams::const_iterator iterC = m_mapParameters.find("CS");

	double a = (iterA!=m_mapParameters.end() ? iterA->second : 0.);
	double b = (iterB!=m_mapParameters.end() ? iterB->second : 0.);
	double c = (iterC!=m_mapParameters.end() ? iterC->second : 0.);

	return std::array<double,3> {a,b,c};
}

std::array<double,3> FilePsi::GetSampleAngles() const
{
	t_mapIParams::const_iterator iterA = m_mapParameters.find("AA");
	t_mapIParams::const_iterator iterB = m_mapParameters.find("BB");
	t_mapIParams::const_iterator iterC = m_mapParameters.find("CC");

	double alpha = (iterA!=m_mapParameters.end() ? iterA->second/180.*M_PI : M_PI/2.);
	double beta = (iterB!=m_mapParameters.end() ? iterB->second/180.*M_PI : M_PI/2.);
	double gamma = (iterC!=m_mapParameters.end() ? iterC->second/180.*M_PI : M_PI/2.);

	return std::array<double,3> {alpha, beta, gamma};
}

std::array<double,2> FilePsi::GetMonoAnaD() const
{
	t_mapIParams::const_iterator iterM = m_mapParameters.find("DM");
	t_mapIParams::const_iterator iterA = m_mapParameters.find("DA");

	double m = (iterM!=m_mapParameters.end() ? iterM->second : 3.355);
	double a = (iterA!=m_mapParameters.end() ? iterA->second : 3.355);

	return std::array<double,2> {m, a};
}

std::array<bool, 3> FilePsi::GetScatterSenses() const
{
	t_mapIParams::const_iterator iterM = m_mapParameters.find("SM");
	t_mapIParams::const_iterator iterS = m_mapParameters.find("SS");
	t_mapIParams::const_iterator iterA = m_mapParameters.find("SA");

	bool m = (iterM!=m_mapParameters.end() ? iterM->second>0 : 0);
	bool s = (iterM!=m_mapParameters.end() ? iterS->second>0 : 1);
	bool a = (iterM!=m_mapParameters.end() ? iterA->second>0 : 0);

	return std::array<bool,3> {m, s, a};
}

std::array<double, 3> FilePsi::GetScatterPlane0() const
{
	t_mapIParams::const_iterator iterX = m_mapParameters.find("AX");
	t_mapIParams::const_iterator iterY = m_mapParameters.find("AY");
	t_mapIParams::const_iterator iterZ = m_mapParameters.find("AZ");

	double x = (iterX!=m_mapParameters.end() ? iterX->second : 1.);
	double y = (iterY!=m_mapParameters.end() ? iterY->second : 0.);
	double z = (iterZ!=m_mapParameters.end() ? iterZ->second : 0.);

	return std::array<double,3> {x,y,z};
}

std::array<double, 3> FilePsi::GetScatterPlane1() const
{
	t_mapIParams::const_iterator iterX = m_mapParameters.find("BX");
	t_mapIParams::const_iterator iterY = m_mapParameters.find("BY");
	t_mapIParams::const_iterator iterZ = m_mapParameters.find("BZ");

	double x = (iterX!=m_mapParameters.end() ? iterX->second : 0.);
	double y = (iterY!=m_mapParameters.end() ? iterY->second : 1.);
	double z = (iterZ!=m_mapParameters.end() ? iterZ->second : 0.);

	return std::array<double,3> {x,y,z};
}

double FilePsi::GetKFix() const
{
	t_mapIParams::const_iterator iterK = m_mapParameters.find("KFIX");
	double k = (iterK!=m_mapParameters.end() ? iterK->second : 0.);

	return k;
}

std::array<double, 4> FilePsi::GetPosHKLE() const
{
	t_mapIParams::const_iterator iterH = m_mapPosHkl.find("QH");
	t_mapIParams::const_iterator iterK = m_mapPosHkl.find("QK");
	t_mapIParams::const_iterator iterL = m_mapPosHkl.find("QL");
	t_mapIParams::const_iterator iterE = m_mapPosHkl.find("EN");

	double h = (iterH!=m_mapPosHkl.end() ? iterH->second : 0.);
	double k = (iterK!=m_mapPosHkl.end() ? iterK->second : 0.);
	double l = (iterL!=m_mapPosHkl.end() ? iterL->second : 0.);
	double E = (iterE!=m_mapPosHkl.end() ? iterE->second : 0.);

	return std::array<double,4> {h,k,l,E};
}

std::array<double, 4> FilePsi::GetStepsHKLE() const
{
        t_mapIParams::const_iterator iterH = m_mapScanSteps.find("DQH");
	if(iterH==m_mapScanSteps.end()) iterH = m_mapScanSteps.find("QH");

        t_mapIParams::const_iterator iterK = m_mapScanSteps.find("DQK");
	if(iterK==m_mapScanSteps.end()) iterK = m_mapScanSteps.find("QK");

        t_mapIParams::const_iterator iterL = m_mapScanSteps.find("DQL");
	if(iterL==m_mapScanSteps.end()) iterL = m_mapScanSteps.find("QL");

        t_mapIParams::const_iterator iterE = m_mapScanSteps.find("DEN");
	if(iterE==m_mapScanSteps.end()) iterE = m_mapScanSteps.find("EN");


        double h = (iterH!=m_mapScanSteps.end() ? iterH->second : 0.);
        double k = (iterK!=m_mapScanSteps.end() ? iterK->second : 0.);
        double l = (iterL!=m_mapScanSteps.end() ? iterL->second : 0.);
        double E = (iterE!=m_mapScanSteps.end() ? iterE->second : 0.);

        return std::array<double,4> {h,k,l,E};
}

}




/*
// test
// gcc -o 0 file/loadinstr.cpp helper/log.cpp -std=c++11 -lstdc++
int main()
{
	tl::FilePsi psi;
	if(!psi.Load("/home/tweber/tmp/tst.scn"))
	{
		tl::log_err("Cannot load data file.");
		return -1;
	}

	std::cout << "Lattice: ";
	for(int i=0; i<3; ++i)
		std::cout << psi.GetSampleLattice()[i] << ", ";
	std::cout << std::endl;

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
