/*
 * Load instrument-specific data file
 * @author tweber
 * @date feb-2015
 * @copyright GPLv2 or GPLv3
 */

#include "loadinstr.h"
#include "../helper/log.h"
#include "../helper/py.h"
#include "../file/file.h"
#include "../math/neutrons.hpp"
#include <fstream>


namespace tl{


// automatically choose correct instrument
FileInstr* FileInstr::LoadInstr(const char* pcFile)
{
	FileInstr* pDat = nullptr;

	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return nullptr;

	std::string strLine;
	std::getline(ifstr, strLine);
	ifstr.close();

	strLine = tl::str_to_lower(strLine);
	const std::string strNicos("nicos data file");

	if(strLine.find(strNicos) != std::string::npos)	// frm file
		pDat = new FileFrm();
	else											// psi or ill file
		pDat = new FilePsi();

	if(pDat && !pDat->Load(pcFile))
	{
		delete pDat;
		return nullptr;
	}

	return pDat;
}


std::array<double, 5> FileInstr::GetScanHKLKiKf(const char* pcH, const char* pcK,
											const char* pcL, const char* pcE,
											std::size_t i) const
{
	const t_vecVals& vecH = GetCol(pcH);
	const t_vecVals& vecK = GetCol(pcK);
	const t_vecVals& vecL = GetCol(pcL);
	const t_vecVals& vecE = GetCol(pcE);

	std::size_t minSize = tl::min4(vecH.size(), vecK.size(), vecL.size(), vecE.size());
	if(i>=minSize)
	{
		tl::log_err("Scan position ", i, " out of bounds. Size: ", minSize, ".");
		return std::array<double,5>{{0.,0.,0.,0.}};
	}

	double h = vecH[i];
	double k = vecK[i];
	double l = vecL[i];
	double E = vecE[i];

	bool bKiFix = IsKiFixed();
	double kfix = GetKFix();
	double kother = get_other_k(E*meV, kfix/angstrom, bKiFix) * angstrom;

	return std::array<double,5>{{h,k,l, bKiFix?kfix:kother, bKiFix?kother:kfix}};
}


bool FileInstr::MergeWith(const FileInstr* pDat)
{
	if(this->GetColNames().size() != pDat->GetColNames().size())
	{
		log_err("Cannot merge: Mismatching number of columns.");
		return false;
	}
	
	for(const std::string& strCol : GetColNames())
	{
		t_vecVals& col1 = this->GetCol(strCol);
		const t_vecVals& col2 = pDat->GetCol(strCol);
		
		if(col1.size() == 0 || col2.size() == 0)
		{
			log_err("Cannot merge: Column \"", strCol, "\" is empty.");
			return false;			
		}
		
		col1.insert(col1.end(), col2.begin(), col2.end());
	}
	
	return true;
}


// -----------------------------------------------------------------------------


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
	if(iterZeros!=m_mapParams.end()) GetInternalParams(iterZeros->second, m_mapZeros);
	if(iterVars!=m_mapParams.end()) GetInternalParams(iterVars->second, m_mapVariables);
	if(iterPos!=m_mapParams.end()) GetInternalParams(iterPos->second, m_mapPosHkl);
	if(iterSteps!=m_mapParams.end()) GetInternalParams(iterSteps->second, m_mapScanSteps);

	return true;
}

const FileInstr::t_vecVals& FilePsi::GetCol(const std::string& strName) const
{
	return const_cast<FilePsi*>(this)->GetCol(strName);
}

FileInstr::t_vecVals& FilePsi::GetCol(const std::string& strName)
{
	static std::vector<double> vecNull;

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

	return std::array<double,3>{{a,b,c}};
}

std::array<double,3> FilePsi::GetSampleAngles() const
{
	t_mapIParams::const_iterator iterA = m_mapParameters.find("AA");
	t_mapIParams::const_iterator iterB = m_mapParameters.find("BB");
	t_mapIParams::const_iterator iterC = m_mapParameters.find("CC");

	double alpha = (iterA!=m_mapParameters.end() ? iterA->second/180.*M_PI : M_PI/2.);
	double beta = (iterB!=m_mapParameters.end() ? iterB->second/180.*M_PI : M_PI/2.);
	double gamma = (iterC!=m_mapParameters.end() ? iterC->second/180.*M_PI : M_PI/2.);

	return std::array<double,3>{{alpha, beta, gamma}};
}

std::array<double,2> FilePsi::GetMonoAnaD() const
{
	t_mapIParams::const_iterator iterM = m_mapParameters.find("DM");
	t_mapIParams::const_iterator iterA = m_mapParameters.find("DA");

	double m = (iterM!=m_mapParameters.end() ? iterM->second : 3.355);
	double a = (iterA!=m_mapParameters.end() ? iterA->second : 3.355);

	return std::array<double,2>{{m, a}};
}

std::array<bool, 3> FilePsi::GetScatterSenses() const
{
	t_mapIParams::const_iterator iterM = m_mapParameters.find("SM");
	t_mapIParams::const_iterator iterS = m_mapParameters.find("SS");
	t_mapIParams::const_iterator iterA = m_mapParameters.find("SA");

	bool m = (iterM!=m_mapParameters.end() ? (iterM->second>0) : 0);
	bool s = (iterS!=m_mapParameters.end() ? (iterS->second>0) : 1);
	bool a = (iterA!=m_mapParameters.end() ? (iterA->second>0) : 0);

	return std::array<bool,3>{{m, s, a}};
}

std::array<double, 3> FilePsi::GetScatterPlane0() const
{
	t_mapIParams::const_iterator iterX = m_mapParameters.find("AX");
	t_mapIParams::const_iterator iterY = m_mapParameters.find("AY");
	t_mapIParams::const_iterator iterZ = m_mapParameters.find("AZ");

	double x = (iterX!=m_mapParameters.end() ? iterX->second : 1.);
	double y = (iterY!=m_mapParameters.end() ? iterY->second : 0.);
	double z = (iterZ!=m_mapParameters.end() ? iterZ->second : 0.);

	return std::array<double,3>{{x,y,z}};
}

std::array<double, 3> FilePsi::GetScatterPlane1() const
{
	t_mapIParams::const_iterator iterX = m_mapParameters.find("BX");
	t_mapIParams::const_iterator iterY = m_mapParameters.find("BY");
	t_mapIParams::const_iterator iterZ = m_mapParameters.find("BZ");

	double x = (iterX!=m_mapParameters.end() ? iterX->second : 0.);
	double y = (iterY!=m_mapParameters.end() ? iterY->second : 1.);
	double z = (iterZ!=m_mapParameters.end() ? iterZ->second : 0.);

	return std::array<double,3>{{x,y,z}};
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

	return std::array<double,4>{{h,k,l,E}};
}

std::array<double, 4> FilePsi::GetDeltaHKLE() const
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

        return std::array<double,4>{{h,k,l,E}};
}

bool FilePsi::MergeWith(const FileInstr* pDat)
{
	if(!FileInstr::MergeWith(pDat))
		return false;
		
	std::string strNr = pDat->GetScanNumber();
	if(strNr.length() != 0)
	{
		// include merged scan number
		t_mapParams::iterator iter = m_mapParams.find("FILE_");
		if(iter != m_mapParams.end())
			iter->second += std::string(" + ") + strNr;
	}

	return true;	
}

// TODO
bool FilePsi::IsKiFixed() const
{
	return 0;
}

std::size_t FilePsi::GetScanCount() const
{
	if(m_vecData.size() < 1)
		return 0;
	return m_vecData[0].size();
}

std::array<double, 5> FilePsi::GetScanHKLKiKf(std::size_t i) const
{
	return FileInstr::GetScanHKLKiKf("QH", "QK", "QL", "EN", i);
}


std::string FilePsi::GetTitle() const
{
	std::string strTitle;
	t_mapParams::const_iterator iter = m_mapParams.find("TITLE");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}

std::string FilePsi::GetScanNumber() const
{
	std::string strTitle;
	t_mapParams::const_iterator iter = m_mapParams.find("FILE_");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}

std::string FilePsi::GetSampleName() const
{
	return "";
}

std::string FilePsi::GetSpacegroup() const
{
	return "";
}


std::vector<std::string> FilePsi::GetScannedVars() const
{
	std::vector<std::string> vecVars;

	// steps parameter
	for(const t_mapIParams::value_type& pair : m_mapScanSteps)
	{
		if(!float_equal(pair.second, 0.) && pair.first.length())
		{
			if(std::tolower(pair.first[0]) == 'd')
				vecVars.push_back(pair.first.substr(1));
			else
				vecVars.push_back(pair.first);
		}
	}


	// nothing found yet -> try scan command instead
	if(!vecVars.size())
	{
		t_mapParams::const_iterator iter = m_mapParams.find("COMND");
		if(iter != m_mapParams.end())
		{
			std::vector<std::string> vecToks;
			tl::get_tokens<std::string, std::string>(iter->second, " \t", vecToks);
			std::transform(vecToks.begin(), vecToks.end(), vecToks.begin(), str_to_lower<std::string>);
			std::vector<std::string>::iterator iter = std::find(vecToks.begin(), vecToks.end(), "dqh");

			if(iter != vecToks.end())
			{
				double dh = str_to_var<double>(*(++iter));
				double dk = str_to_var<double>(*(++iter));
				double dl = str_to_var<double>(*(++iter));
				double dE = str_to_var<double>(*iter);

				if(!float_equal(dh, 0.)) vecVars.push_back("QH");
				if(!float_equal(dk, 0.)) vecVars.push_back("QK");
				if(!float_equal(dl, 0.)) vecVars.push_back("QL");
				if(!float_equal(dE, 0.)) vecVars.push_back("EN");
			}
		}
	}

	return vecVars;
}


// -----------------------------------------------------------------------------



void FileFrm::ReadHeader(std::istream& istr)
{
	istr.clear();
	istr.seekg(0, std::ios_base::beg);

	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		trim(strLine);
		if(strLine.length()==0 || strLine[0]!='#')
			continue;
		if(strLine.length()>=3 && strLine[0]=='#' && strLine[1]=='#' && strLine[2]=='#')
			continue;

		strLine = strLine.substr(1);
		//std::cout << strLine << std::endl;

		std::pair<std::string, std::string> pairLine =
				tl::split_first<std::string>(strLine, ":", 1);
		if(pairLine.first == "")
			continue;
		else
		{
			//std::cout << "Key: " << pairLine.first << ", Val: " << pairLine.second << std::endl;
			t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

			if(iter == m_mapParams.end())
				m_mapParams.insert(pairLine);
			else
				iter->second += "\n" + pairLine.second;
		}
	}
}

void FileFrm::ReadData(std::istream& istr)
{
	istr.clear();
	istr.seekg(0, std::ios_base::beg);
	skip_after_line<char>(istr, "### scan data", true, false);

	// column headers
	skip_after_char<char>(istr, '#');
	std::string strLineQuantities;
	std::getline(istr, strLineQuantities);
	tl::get_tokens<std::string, std::string, t_vecColNames>
		(strLineQuantities, " \t", m_vecQuantities);

	skip_after_char<char>(istr, '#');
	std::string strLineUnits;
	std::getline(istr, strLineUnits);
	tl::get_tokens<std::string, std::string, t_vecColNames>
		(strLineQuantities, " \t", m_vecUnits);


	m_vecData.resize(m_vecQuantities.size());

	// data
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		trim(strLine);
		if(strLine.length()==0 || strLine[0]=='#')
			continue;

		std::vector<double> vecToks;
		tl::get_tokens<double, std::string>(strLine, " \t", vecToks);

		if(vecToks.size() != m_vecQuantities.size())
		{
			tl::log_warn("Loader: Line size mismatch.");

			// add zeros
			while(m_vecQuantities.size() > vecToks.size())
				vecToks.push_back(0.);
		}

		for(std::size_t iTok=0; iTok<vecToks.size(); ++iTok)
			m_vecData[iTok].push_back(vecToks[iTok]);
	}
}


bool FileFrm::Load(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return false;

	ReadHeader(ifstr);
	ReadData(ifstr);

	return true;
}

const FileInstr::t_vecVals& FileFrm::GetCol(const std::string& strName) const
{
	return const_cast<FileFrm*>(this)->GetCol(strName);
}

FileInstr::t_vecVals& FileFrm::GetCol(const std::string& strName)
{
	static std::vector<double> vecNull;

	for(std::size_t i=0; i<m_vecQuantities.size(); ++i)
	{
		if(m_vecQuantities[i] == strName)
			return m_vecData[i];
	}

	return vecNull;
}


std::array<double, 3> FileFrm::GetSampleLattice() const
{
	t_mapParams::const_iterator iter = m_mapParams.find("Sample_lattice");
	std::vector<double> vec = tl::get_py_array(iter->second);
	if(vec.size() != 3)
	{
		tl::log_err("Invalid lattice array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}

	return std::array<double,3>{{vec[0],vec[1],vec[2]}};
}

std::array<double, 3> FileFrm::GetSampleAngles() const
{
	t_mapParams::const_iterator iter = m_mapParams.find("Sample_angles");
	std::vector<double> vec = tl::get_py_array(iter->second);
	if(vec.size() != 3)
	{
		tl::log_err("Invalid angle array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}

	return std::array<double,3>{{vec[0]/180.*M_PI, vec[1]/180.*M_PI, vec[2]/180.*M_PI}};
}

std::array<double, 2> FileFrm::GetMonoAnaD() const
{
	t_mapParams::const_iterator iterM = m_mapParams.find("mono_dvalue");
	t_mapParams::const_iterator iterA = m_mapParams.find("ana_dvalue");

	double m = (iterM!=m_mapParams.end() ? tl::str_to_var<double>(iterM->second) : 3.355);
	double a = (iterA!=m_mapParams.end() ? tl::str_to_var<double>(iterA->second) : 3.355);

	return std::array<double,2>{{m, a}};
}

std::array<bool, 3> FileFrm::GetScatterSenses() const
{
	std::vector<int> vec;

	t_mapParams::const_iterator iter;
	for(iter=m_mapParams.begin(); iter!=m_mapParams.end(); ++iter)
	{
		if(iter->first.find("scatteringsense") != std::string::npos)
		{
			vec = tl::get_py_array<std::string, std::vector<int>>(iter->second);
			break;
		}
	}

	if(vec.size() != 3)
	{
		vec.resize(3);
		vec[0] = 0; vec[1] = 1; vec[2] = 0;
	}

	return std::array<bool,3>{{vec[0]>0, vec[1]>0, vec[2]>0}};
}

std::array<double, 3> FileFrm::GetScatterPlane0() const
{
	t_mapParams::const_iterator iter = m_mapParams.find("Sample_orient1");
	std::vector<double> vec = tl::get_py_array(iter->second);
	if(vec.size() != 3)
	{
		tl::log_err("Invalid sample peak 1 array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}
	return std::array<double,3>{{vec[0],vec[1],vec[2]}};
}

std::array<double, 3> FileFrm::GetScatterPlane1() const
{
	t_mapParams::const_iterator iter = m_mapParams.find("Sample_orient2");
	std::vector<double> vec = tl::get_py_array(iter->second);
	if(vec.size() != 3)
	{
		tl::log_err("Invalid sample peak 2 array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}
	return std::array<double,3>{{-vec[0],-vec[1],-vec[2]}};	// LH -> RH
}

double FileFrm::GetKFix() const
{
	std::string strKey = (IsKiFixed() ? "ki_value" : "kf_value");

	t_mapParams::const_iterator iter = m_mapParams.find(strKey);
	return (iter!=m_mapParams.end() ? tl::str_to_var<double>(iter->second) : 0.);
}

bool FileFrm::IsKiFixed() const
{
	std::string strScanMode = "ckf";

	t_mapParams::const_iterator iter;
	for(iter=m_mapParams.begin(); iter!=m_mapParams.end(); ++iter)
	{
		if(iter->first.find("scanmode") != std::string::npos)
		{
			strScanMode = tl::str_to_lower(iter->second);
			tl::trim(strScanMode);
			break;
		}
	}

	if(strScanMode == "cki")
		return 1;
	return 0;
}

std::size_t FileFrm::GetScanCount() const
{
	if(m_vecData.size() < 1)
		return 0;
	return m_vecData[0].size();
}

std::array<double, 5> FileFrm::GetScanHKLKiKf(std::size_t i) const
{
	return FileInstr::GetScanHKLKiKf("h", "k", "l", "E", i);
}

bool FileFrm::MergeWith(const FileInstr* pDat)
{
	if(!FileInstr::MergeWith(pDat))
		return false;
		
	std::string strNr = pDat->GetScanNumber();
	if(strNr.length() != 0)
	{
		// include merged scan number
		t_mapParams::iterator iter = m_mapParams.find("number");
		if(iter != m_mapParams.end())
			iter->second += std::string(" + ") + strNr;
	}

	return true;	
}


std::string FileFrm::GetTitle() const
{
	std::string strTitle;
	t_mapParams::const_iterator iter = m_mapParams.find("Exp_title");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}

std::string FileFrm::GetScanNumber() const
{
	std::string strTitle;
	t_mapParams::const_iterator iter = m_mapParams.find("number");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;	
}

std::string FileFrm::GetSampleName() const
{
	std::string strName;
	t_mapParams::const_iterator iter = m_mapParams.find("Sample_samplename");
	if(iter != m_mapParams.end())
		strName = iter->second;
	return strName;
}

std::string FileFrm::GetSpacegroup() const
{
	std::string strSG;
	t_mapParams::const_iterator iter = m_mapParams.find("Sample_spacegroup");
	if(iter != m_mapParams.end())
		strSG = iter->second;
	return strSG;
}


std::vector<std::string> FileFrm::GetScannedVars() const
{
	std::vector<std::string> vecVars;

	// TODO

	return vecVars;
}

}




// -----------------------------------------------------------------------------


/*
// test
// gcc -o 0 file/loadinstr.cpp helper/log.cpp -std=c++11 -lstdc++
int main()
{
	tl::FileFrm dat;
	if(!dat.Load("/home/tweber/tmp/tst.dat"))
	{
		tl::log_err("Cannot load data file.");
		return -1;
	}

	std::array<double,3> latt = dat.GetSampleAngles();
	std::cout << latt[0] << ", " << latt[1] << ", " << latt[2] << std::endl;

	std::cout << "kfix = " << dat.GetKFix() << std::endl;

	return 0;
}
*/