/*
 * Loads instrument-specific data files
 * @author tweber
 * @date feb-2015
 * @license GPLv2 or GPLv3
 */

#include "loadinstr.h"
#include "../helper/log.h"
#include "../helper/py.h"
#include "../file/file.h"
#if !defined NO_IOSTR
	#include "../file/comp.h"
#endif
#include "../math/neutrons.hpp"
#include <fstream>


#ifndef USE_BOOST_REX
	#include <regex>
	namespace rex = ::std;
#else
	#include <boost/tr1/regex.hpp>
	namespace rex = ::boost;
#endif


namespace tl{


// automatically choose correct instrument
FileInstr* FileInstr::LoadInstr(const char* pcFile)
{
	FileInstr* pDat = nullptr;

	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return nullptr;

#if !defined NO_IOSTR
	std::istream* pIstr = create_autodecomp_istream(ifstr);
	std::unique_ptr<std::istream> ptrIstr(pIstr);
	if(!pIstr)
		return nullptr;
#else
	std::istream* pIstr = &ifstr;
#endif

	std::string strLine, strLine2;
	std::getline(*pIstr, strLine);
	std::getline(*pIstr, strLine2);
	//pIstr->close();


	trim(strLine);
	trim(strLine2);
	strLine = str_to_lower(strLine);
	strLine2 = str_to_lower(strLine2);

	if(strLine == "")
		return nullptr;

	const std::string strNicos("nicos data file");
	const std::string strMacs("ice");

	if(strLine.find(strNicos) != std::string::npos)		// frm file
		pDat = new FileFrm();
	else if(strLine.find(strMacs) != std::string::npos)	// macs file
		pDat = new FileMacs();
	else if(strLine2.find("scan start") != std::string::npos)
		pDat = new FileTrisp();
	else							// psi or ill file
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

	std::size_t minSize = min4(vecH.size(), vecK.size(), vecL.size(), vecE.size());
	if(i>=minSize)
	{
		log_err("Scan position ", i, " out of bounds. Size: ", minSize, ".");
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
	get_tokens<std::string, std::string, t_vecColNames>(strHdr, " \t", m_vecColNames);

	m_vecData.resize(m_vecColNames.size());

	// data
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		tl::trim(strLine);

		if(strLine.length() == 0)
			continue;
		if(strLine[0] == '#')
			continue;

		std::vector<double> vecToks;
		get_tokens<double, std::string>(strLine, " \t", vecToks);

		if(vecToks.size() != m_vecColNames.size())
		{
			log_warn("Loader: Line size mismatch.");

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
	get_tokens<std::string, std::string>(strAll, ",\n", vecToks);
	//std::cout << strAll << std::endl;

	for(const std::string& strTok : vecToks)
	{
		std::pair<std::string, std::string> pair =
				split_first<std::string>(strTok, "=", 1);

		if(pair.first == "")
			continue;

		double dVal = str_to_var<double>(pair.second);
		mapPara.insert(t_mapIParams::value_type(pair.first, dVal));

		//std::cout << "Key: " << pair.first << ", Val: " << dVal << std::endl;
	}
}

bool FilePsi::Load(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return false;

#if !defined NO_IOSTR
	std::istream* pIstr = create_autodecomp_istream(ifstr);
	if(!pIstr) return false;
	std::unique_ptr<std::istream> ptrIstr(pIstr);
#else
	std::istream* pIstr = &ifstr;
#endif

	while(!pIstr->eof())
	{
		std::string strLine;
		std::getline(*pIstr, strLine);

		if(strLine.substr(0,4) == "RRRR")
			skip_after_line<char>(*pIstr, "VVVV", true);

		std::pair<std::string, std::string> pairLine =
				split_first<std::string>(strLine, ":", 1);
		if(pairLine.first == "DATA_")
			ReadData(*pIstr);
		else if(pairLine.first == "")
			continue;
		else
		{
			t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

			if(iter == m_mapParams.end())
				m_mapParams.insert(pairLine);
			else
				iter->second += ", " + pairLine.second;
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
		ostr << "Param: " << val.first
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

	double alpha = (iterA!=m_mapParameters.end() ? tl::d2r(iterA->second) : M_PI/2.);
	double beta = (iterB!=m_mapParameters.end() ? tl::d2r(iterB->second) : M_PI/2.);
	double gamma = (iterC!=m_mapParameters.end() ? tl::d2r(iterC->second) : M_PI/2.);

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

std::string FilePsi::GetUser() const
{
	std::string strUser;
	t_mapParams::const_iterator iter = m_mapParams.find("USER_");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
}

std::string FilePsi::GetLocalContact() const
{
	std::string strUser;
	t_mapParams::const_iterator iter = m_mapParams.find("LOCAL");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
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
			//std::cout << iter->second << std::endl;
			get_tokens<std::string, std::string>(iter->second, " \t", vecToks);
			for(std::string& strTok : vecToks)
				tl::trim(strTok);

			std::transform(vecToks.begin(), vecToks.end(), vecToks.begin(), str_to_lower<std::string>);
			std::vector<std::string>::iterator iterTok = std::find(vecToks.begin(), vecToks.end(), "dqh");

			if(iterTok != vecToks.end())
			{
				double dh = str_to_var<double>(*(++iterTok));
				double dk = str_to_var<double>(*(++iterTok));
				double dl = str_to_var<double>(*(++iterTok));
				double dE = str_to_var<double>(*(++iterTok));

				if(!float_equal(dh, 0.)) vecVars.push_back("QH");
				if(!float_equal(dk, 0.)) vecVars.push_back("QK");
				if(!float_equal(dl, 0.)) vecVars.push_back("QL");
				if(!float_equal(dE, 0.)) vecVars.push_back("EN");
			}


			// still nothing found, try regex
			if(!vecVars.size())
			{
				const std::string strRegex = R"REX((SC|SCAN|sc|scan)[ \t]+([A-Za-z0-9]+)[ \t]+[0-9\.-]+[ \t]+[d|D]([A-Za-z0-9]+).*)REX";
				rex::regex rx(strRegex, rex::regex::ECMAScript);
				rex::smatch m;
				if(rex::regex_search(iter->second, m, rx) && m.size()>3)
				{
					const std::string& strSteps = m[3];
					vecVars.push_back(str_to_upper(strSteps));
				}
			}
		}
	}

	if(!vecVars.size())
	{
		tl::log_warn("Could not determine scan variable, using first column.");
		if(m_vecColNames.size() >= 1)
			vecVars.push_back(m_vecColNames[0]);
	}

	//for(std::string& strVar : vecVars)
	//	tl::log_info("Scan var: ", strVar);
	//tl::log_info("--------");
	return vecVars;
}

std::string FilePsi::GetCountVar() const
{
	// TODO
	return "CNTS";
}

std::string FilePsi::GetMonVar() const
{
	// TODO
	return "M1";
}

std::string FilePsi::GetScanCommand() const
{
	std::string strCmd;
	t_mapParams::const_iterator iter = m_mapParams.find("COMND");
	if(iter != m_mapParams.end())
		strCmd = iter->second;
	return strCmd;
}

std::string FilePsi::GetTimestamp() const
{
	std::string strDate;
	t_mapParams::const_iterator iter = m_mapParams.find("DATE_");
	if(iter != m_mapParams.end())
		strDate = iter->second;
	return strDate;
}

// -----------------------------------------------------------------------------



void FileFrm::ReadHeader(std::istream& istr)
{
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		//std::cout << strLine << std::endl;
		trim(strLine);
		if(strLine.length()==0 || strLine[0]!='#')
			continue;
		if(strLine.length()>=3 && strLine[0]=='#' && strLine[1]=='#' && strLine[2]=='#')
		{
			std::string strCreatedAt("created at");
			std::size_t iPosCreated = strLine.find(strCreatedAt);
			if(iPosCreated != std::string::npos)
			{
				iPosCreated += strCreatedAt.length();
				std::string strDate = strLine.substr(iPosCreated);
				tl::trim(strDate);

				m_mapParams["file_timestamp"] = strDate;
			}

			continue;
		}

		strLine = strLine.substr(1);
		//std::cout << strLine << std::endl;

		std::pair<std::string, std::string> pairLine =
				split_first<std::string>(strLine, ":", 1);
		if(pairLine.first == "")
			continue;
		else
		{
			//std::cout << "Key: " << pairLine.first << ", Val: " << pairLine.second << std::endl;
			t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

			if(iter == m_mapParams.end())
				m_mapParams.insert(pairLine);
			else
				iter->second += ", " + pairLine.second;
		}
	}
}

void FileFrm::ReadData(std::istream& istr)
{
	skip_after_line<char>(istr, "### scan data", true, false);

	// column headers
	skip_after_char<char>(istr, '#');
	std::string strLineQuantities;
	std::getline(istr, strLineQuantities);
	get_tokens<std::string, std::string, t_vecColNames>
		(strLineQuantities, " \t", m_vecQuantities);

	skip_after_char<char>(istr, '#');
	std::string strLineUnits;
	std::getline(istr, strLineUnits);
	get_tokens<std::string, std::string, t_vecColNames>
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
		get_tokens<double, std::string>(strLine, " \t", vecToks);

		if(vecToks.size() != m_vecQuantities.size())
		{
			log_warn("Loader: Line size mismatch.");

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
	for(int iStep : {0,1})
	{
		std::ifstream ifstr(pcFile);
		if(!ifstr.is_open())
			return false;

#if !defined NO_IOSTR
		std::istream* pIstr = create_autodecomp_istream(ifstr);
		if(!pIstr) return false;
		std::unique_ptr<std::istream> ptrIstr(pIstr);
#else
		std::ifstream *pIstr = &ifstr;
#endif

		//std::streampos posIstr = pIstr->tellg();
		//ReadHeader(*pIstr);
		//pIstr->seekg(posIstr, std::ios_base::beg);
		//pIstr->clear();
		//ReadData(*pIstr);

		if(iStep==0)
			ReadHeader(*pIstr);
		else if(iStep==1)
			ReadData(*pIstr);
	}

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
	if(iter == m_mapParams.end())
		return std::array<double,3>{{0.,0.,0.}};

	std::vector<double> vec = get_py_array(iter->second);
	if(vec.size() != 3)
	{
		log_err("Invalid lattice array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}

	return std::array<double,3>{{vec[0],vec[1],vec[2]}};
}

std::array<double, 3> FileFrm::GetSampleAngles() const
{
	t_mapParams::const_iterator iter = m_mapParams.find("Sample_angles");
	if(iter == m_mapParams.end())
		return std::array<double,3>{{0.,0.,0.}};

	std::vector<double> vec = get_py_array(iter->second);
	if(vec.size() != 3)
	{
		log_err("Invalid angle array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}

	return std::array<double,3>{{tl::d2r(vec[0]), tl::d2r(vec[1]), tl::d2r(vec[2])}};
}

std::array<double, 2> FileFrm::GetMonoAnaD() const
{
	t_mapParams::const_iterator iterM = m_mapParams.find("mono_dvalue");
	t_mapParams::const_iterator iterA = m_mapParams.find("ana_dvalue");

	double m = (iterM!=m_mapParams.end() ? str_to_var<double>(iterM->second) : 3.355);
	double a = (iterA!=m_mapParams.end() ? str_to_var<double>(iterA->second) : 3.355);

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
			vec = get_py_array<std::string, std::vector<int>>(iter->second);
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
	if(iter == m_mapParams.end())
		return std::array<double,3>{{0.,0.,0.}};

	std::vector<double> vec = get_py_array(iter->second);
	if(vec.size() != 3)
	{
		log_err("Invalid sample peak 1 array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}
	return std::array<double,3>{{vec[0],vec[1],vec[2]}};
}

std::array<double, 3> FileFrm::GetScatterPlane1() const
{
	t_mapParams::const_iterator iter = m_mapParams.find("Sample_orient2");
	if(iter == m_mapParams.end())
		return std::array<double,3>{{0.,0.,0.}};

	std::vector<double> vec = get_py_array(iter->second);
	if(vec.size() != 3)
	{
		log_err("Invalid sample peak 2 array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}
	return std::array<double,3>{{-vec[0],-vec[1],-vec[2]}};	// LH -> RH
}

double FileFrm::GetKFix() const
{
	std::string strKey = (IsKiFixed() ? "ki_value" : "kf_value");

	t_mapParams::const_iterator iter = m_mapParams.find(strKey);
	return (iter!=m_mapParams.end() ? str_to_var<double>(iter->second) : 0.);
}

bool FileFrm::IsKiFixed() const
{
	std::string strScanMode = "ckf";

	t_mapParams::const_iterator iter;
	for(iter=m_mapParams.begin(); iter!=m_mapParams.end(); ++iter)
	{
		if(iter->first.find("scanmode") != std::string::npos)
		{
			strScanMode = str_to_lower(iter->second);
			trim(strScanMode);
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

std::string FileFrm::GetUser() const
{
	std::string strUser;
	t_mapParams::const_iterator iter = m_mapParams.find("Exp_users");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
}

std::string FileFrm::GetLocalContact() const
{
	std::string strUser;
	t_mapParams::const_iterator iter = m_mapParams.find("Exp_localcontact");
	if(iter != m_mapParams.end())
		strUser = iter->second;
	return strUser;
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

	// scan command
	t_mapParams::const_iterator iter = m_mapParams.find("info");
	if(iter != m_mapParams.end())
	{
		const std::string& strInfo = iter->second;

		// try qscan/qcscan
		const std::string strRegex = R"REX((qscan|qcscan)\((\[.*\])[, ]+(\[.*\]).*\))REX";
		rex::regex rx(strRegex, rex::regex::ECMAScript);
		rex::smatch m;
		if(rex::regex_search(strInfo, m, rx) && m.size()>3)
		{
			const std::string& strSteps = m[3];
			std::vector<double> vecSteps = get_py_array(strSteps);

			if(vecSteps.size()>0 && !float_equal(vecSteps[0], 0.))
				vecVars.push_back("h");
			if(vecSteps.size()>1 && !float_equal(vecSteps[1], 0.))
				vecVars.push_back("k");
			if(vecSteps.size()>2 && !float_equal(vecSteps[2], 0.))
				vecVars.push_back("l");
			if(vecSteps.size()>3 && !float_equal(vecSteps[3], 0.))
				vecVars.push_back("E");
		}


		if(vecVars.size() == 0)
		{
			// try scan/cscan
			const std::string strRegexDevScan = R"REX((scan|cscan)\(([A-Za-z0-9_\.]+)[, ]+.*\))REX";
			rex::regex rxDev(strRegexDevScan, rex::regex::ECMAScript);
			rex::smatch mDev;
			if(rex::regex_search(strInfo, mDev, rxDev) && mDev.size()>2)
			{
				vecVars.push_back(mDev[2]);
			}
		}
	}

	if(!vecVars.size())
	{
		tl::log_warn("Could not determine scan variable, using first column.");
		if(m_vecQuantities.size() >= 1)
			vecVars.push_back(m_vecQuantities[0]);
	}

	return vecVars;
}

std::string FileFrm::GetCountVar() const
{
	// TODO
	return "ctr1";
}

std::string FileFrm::GetMonVar() const
{
	// TODO
	return "mon1";
}

std::string FileFrm::GetScanCommand() const
{
	std::string strCmd;
	t_mapParams::const_iterator iter = m_mapParams.find("info");
	if(iter != m_mapParams.end())
		strCmd = iter->second;
	return strCmd;
}

std::string FileFrm::GetTimestamp() const
{
	std::string strDate;
	t_mapParams::const_iterator iter = m_mapParams.find("file_timestamp");
	if(iter != m_mapParams.end())
		strDate = iter->second;
	return strDate;
}


// -----------------------------------------------------------------------------



void FileMacs::ReadHeader(std::istream& istr)
{
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);

		trim(strLine);
		if(strLine.length()==0 || strLine[0]!='#')
			continue;

		strLine = strLine.substr(1);

		std::pair<std::string, std::string> pairLine =
				split_first<std::string>(strLine, " \t", 1);
		//std::cout << "key: " << pairLine.first << ", val: " << pairLine.second << std::endl;

		if(pairLine.first == "")
			continue;
		else if(pairLine.first == "Columns")
		{
			tl::get_tokens<std::string, std::string>(pairLine.second, " \t", m_vecQuantities);
			//for(const std::string& strCol : m_vecQuantities)
			//	std::cout << "column: \"" << strCol << "\"" << std::endl;
			continue;
		}
		else
		{
			t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

			if(iter == m_mapParams.end())
				m_mapParams.insert(pairLine);
			else
				iter->second += ", " + pairLine.second;
		}
	}
}

void FileMacs::ReadData(std::istream& istr)
{
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
		get_tokens<double, std::string>(strLine, " \t", vecToks);

		if(vecToks.size() != m_vecQuantities.size())
		{
			log_warn("Loader: Line size mismatch.");

			// add zeros
			while(m_vecQuantities.size() > vecToks.size())
				vecToks.push_back(0.);
		}

		for(std::size_t iTok=0; iTok<vecToks.size(); ++iTok)
			m_vecData[iTok].push_back(vecToks[iTok]);
	}
}


bool FileMacs::Load(const char* pcFile)
{
	for(int iStep : {0,1})
	{
		std::ifstream ifstr(pcFile);
		if(!ifstr.is_open())
			return false;

#if !defined NO_IOSTR
		std::istream* pIstr = create_autodecomp_istream(ifstr);
		if(!pIstr) return false;
		std::unique_ptr<std::istream> ptrIstr(pIstr);
#else
		std::ifstream *pIstr = &ifstr;
#endif

		if(iStep==0)
			ReadHeader(*pIstr);
		else if(iStep==1)
			ReadData(*pIstr);
	}

	return true;
}

const FileInstr::t_vecVals& FileMacs::GetCol(const std::string& strName) const
{
	return const_cast<FileMacs*>(this)->GetCol(strName);
}

FileInstr::t_vecVals& FileMacs::GetCol(const std::string& strName)
{
	static std::vector<double> vecNull;

	for(std::size_t i=0; i<m_vecQuantities.size(); ++i)
	{
		if(m_vecQuantities[i] == strName)
			return m_vecData[i];
	}

	return vecNull;
}


std::array<double, 3> FileMacs::GetSampleLattice() const
{
	t_mapParams::const_iterator iter = m_mapParams.find("Lattice");
	if(iter == m_mapParams.end())
		return std::array<double,3>{{0.,0.,0.}};

	std::vector<double> vecToks;
	tl::get_tokens<double, std::string>(iter->second, " \t", vecToks);
	if(vecToks.size() != 6)
	{
		log_err("Invalid sample lattice array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}

	return std::array<double,3>{{vecToks[0], vecToks[1], vecToks[2]}};
}

std::array<double, 3> FileMacs::GetSampleAngles() const
{
	t_mapParams::const_iterator iter = m_mapParams.find("Lattice");
	if(iter == m_mapParams.end())
		return std::array<double,3>{{0.,0.,0.}};

	std::vector<double> vecToks;
	tl::get_tokens<double, std::string>(iter->second, " \t", vecToks);
	if(vecToks.size() != 6)
	{
		log_err("Invalid sample lattice array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}

	return std::array<double,3>{{tl::d2r(vecToks[3]), tl::d2r(vecToks[4]), tl::d2r(vecToks[5])}};
}

std::array<double, 2> FileMacs::GetMonoAnaD() const
{
	t_mapParams::const_iterator iterM = m_mapParams.find("MonoSpacing");
	t_mapParams::const_iterator iterA = m_mapParams.find("AnaSpacing");

	double m = (iterM!=m_mapParams.end() ? str_to_var<double>(iterM->second) : 3.355);
	double a = (iterA!=m_mapParams.end() ? str_to_var<double>(iterA->second) : 3.355);

	return std::array<double,2>{{m, a}};
}

std::array<bool, 3> FileMacs::GetScatterSenses() const
{
	return std::array<bool,3>{{0, 1, 0}};
}

std::array<double, 3> FileMacs::GetScatterPlane0() const
{
	t_mapParams::const_iterator iter = m_mapParams.find("Orient");
	if(iter == m_mapParams.end())
		return std::array<double,3>{{0.,0.,0.}};

	std::vector<double> vecToks;
	tl::get_tokens<double, std::string>(iter->second, " \t", vecToks);
	if(vecToks.size() != 6)
	{
		log_err("Invalid sample orientation array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}

	return std::array<double,3>{{vecToks[0],vecToks[1],vecToks[2]}};
}

std::array<double, 3> FileMacs::GetScatterPlane1() const
{
	t_mapParams::const_iterator iter = m_mapParams.find("Orient");
	if(iter == m_mapParams.end())
		return std::array<double,3>{{0.,0.,0.}};

	std::vector<double> vecToks;
	tl::get_tokens<double, std::string>(iter->second, " \t", vecToks);
	if(vecToks.size() != 6)
	{
		log_err("Invalid sample orientation array size.");
		return std::array<double,3>{{0.,0.,0.}};
	}

	return std::array<double,3>{{vecToks[3],vecToks[4],vecToks[5]}};
}

double FileMacs::GetKFix() const
{
	// 1) look in data columns
	const std::string strKey = (IsKiFixed() ? "Ei" : "Ef");
	const t_vecVals& vecVals = GetCol(strKey);
	if(vecVals.size() != 0)
	{
		bool bImag;
		double k = tl::E2k(vecVals[0] * tl::meV, bImag) * tl::angstrom;
		return k;
	}


	// 2) look in header
	t_mapParams::const_iterator iter = m_mapParams.find("FixedE");
	if(iter==m_mapParams.end())
	{
		tl::log_err("Cannot determine kfix.");
		return 0.;
	}

	std::vector<std::string> vecToks;
	tl::get_tokens<std::string, std::string>(iter->second, " \t", vecToks);

	if(vecToks.size()<2)
	{
		tl::log_err("Cannot determine kfix.");
		return 0.;
	}

	double dEfix = tl::str_to_var<double>(vecToks[1]);
	bool bImag;
	double k = tl::E2k(dEfix * tl::meV, bImag) * tl::angstrom;
	return k;
}

bool FileMacs::IsKiFixed() const
{
	t_mapParams::const_iterator iter = m_mapParams.find("FixedE");
	if(iter==m_mapParams.end())
		return 0;	// assume ckf

	std::vector<std::string> vecToks;
	tl::get_tokens<std::string, std::string>(iter->second, " \t", vecToks);

	if(vecToks.size()==0)
		return 0;	// assume ckf

	std::string strFixedE = vecToks[0];
	tl::trim(strFixedE);

	if(strFixedE == "Ef")
		return 0;
	else if(strFixedE == "Ei")
		return 1;

	return 0;		// assume ckf
}

std::size_t FileMacs::GetScanCount() const
{
	if(m_vecData.size() < 1)
		return 0;
	return m_vecData[0].size();
}

std::array<double, 5> FileMacs::GetScanHKLKiKf(std::size_t i) const
{
	return FileInstr::GetScanHKLKiKf("QX", "QY", "QZ", "E", i);
}

bool FileMacs::MergeWith(const FileInstr* pDat)
{
	if(!FileInstr::MergeWith(pDat))
		return false;

	std::string strNr = pDat->GetScanNumber();
	if(strNr.length() != 0)
	{
		// include merged scan number
		t_mapParams::iterator iter = m_mapParams.find("Filename");
		if(iter != m_mapParams.end())
			iter->second += std::string(" + ") + strNr;
	}

	return true;
}


std::string FileMacs::GetTitle() const
{
	std::string strTitle;
	t_mapParams::const_iterator iter = m_mapParams.find("ExptID");
	if(iter != m_mapParams.end())
		strTitle = iter->second;

	iter = m_mapParams.find("ExptName");
	if(iter != m_mapParams.end() && iter->second != "")
		strTitle += " - " + iter->second;

	return strTitle;
}

std::string FileMacs::GetUser() const
{
	std::string str;
	t_mapParams::const_iterator iter = m_mapParams.find("User");
	if(iter != m_mapParams.end())
		str = iter->second;
	return str;
}

std::string FileMacs::GetLocalContact() const
{
	// TODO
	return "";
}


std::string FileMacs::GetScanNumber() const
{
	std::string strTitle;
	t_mapParams::const_iterator iter = m_mapParams.find("Filename");
	if(iter != m_mapParams.end())
		strTitle = iter->second;
	return strTitle;
}

std::string FileMacs::GetSampleName() const
{
	return "";
}

std::string FileMacs::GetSpacegroup() const
{
	return "";
}

std::vector<std::string> FileMacs::GetScannedVars() const
{
	std::vector<std::string> vecScan;

	t_mapParams::const_iterator iter = m_mapParams.find("Scan");
	if(iter != m_mapParams.end())
	{
		std::vector<std::string> vecToks;
		tl::get_tokens<std::string, std::string>(iter->second, " \t", vecToks);

		if(vecToks.size() >= 2)
			vecScan.push_back(vecToks[1]);
	}

	if(!vecScan.size())
	{
		tl::log_warn("Could not determine scan variable, using first column.");
		if(m_vecQuantities.size() >= 1)
			vecScan.push_back(m_vecQuantities[0]);
	}

	return vecScan;
}

std::string FileMacs::GetCountVar() const
{
	// TODO
	return "SPEC";
}

std::string FileMacs::GetMonVar() const
{
	// TODO
	return "Monitor";
}

std::string FileMacs::GetScanCommand() const
{
	// TODO
	return "";
}

std::string FileMacs::GetTimestamp() const
{
	std::string str;
	t_mapParams::const_iterator iter = m_mapParams.find("Date");
	if(iter != m_mapParams.end())
		str = iter->second;
	return str;
}



// -----------------------------------------------------------------------------



void FileTrisp::ReadHeader(std::istream& istr)
{
	bool bInVarSection = 0;
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);

		trim(strLine);
		if(strLine.length()==0)
			continue;

		if(str_contains<std::string>(strLine, "----", 0))	// new variable section beginning
		{
			bInVarSection = 1;
			//std::cout << "Section: " << strLine << std::endl;

			if(str_contains<std::string>(strLine, "steps", 0))
				break;
			continue;
		}

		if(bInVarSection)
		{
			std::pair<std::string, std::string> pairLine =
					split_first<std::string>(strLine, " \t", 1);

			if(pairLine.first == "")
				continue;
			//std::cout << "key: " << pairLine.first << ", val: " << pairLine.second << std::endl;

			t_mapParams::iterator iter = m_mapParams.find(pairLine.first);

			if(iter == m_mapParams.end())
				m_mapParams.insert(pairLine);
			else
				iter->second += ", " + pairLine.second;
		}
		else
		{
			if(begins_with<std::string>(str_to_lower(strLine), "scan start:"))
				m_mapParams["scan_start_timestamp"] = trimmed(strLine.substr(11));
			else if(begins_with<std::string>(str_to_lower(strLine), "sc"))
				m_mapParams["scan_command"] = strLine;
		}
	}
}

void FileTrisp::ReadData(std::istream& istr)
{
	bool bAtStepsBeginning = 0;
	bool bInFooter = 0;

	// data
	while(!istr.eof())
	{
		std::string strLine;
		std::getline(istr, strLine);
		trim(strLine);

		if(!bAtStepsBeginning)
		{
			if(begins_with<std::string>(str_to_lower(strLine), "pnt"))
			{
				get_tokens<std::string, std::string>(strLine, " \t", m_vecQuantities);
				//for(const std::string& strCol : m_vecQuantities)
				//	std::cout << "col: " << strCol << std::endl;

				bAtStepsBeginning = 1;
				m_vecData.resize(m_vecQuantities.size());
			}
			continue;
		}

		if(strLine.length()==0 || strLine[0]=='#')
			continue;


		// character in scan data -> beginning of footer
		for(std::string::value_type c : split_first<std::string>(strLine, " \t", 1).first)
		{
			if(std::isalpha(c))
			{
				if(begins_with<std::string>(str_to_lower(strLine), "scan end:"))
					m_mapParams["scan_finish_timestamp"] = trimmed(strLine.substr(9));
				else if(begins_with<std::string>(str_to_lower(strLine), "scan"))
				{
					std::pair<std::string, std::string> pairLine = 
						split_first<std::string>(strLine, " \t", 1);

					m_mapParams["scan_vars"] = trimmed(pairLine.second);
				}

				bInFooter = 1;
			}
		}


		if(!bInFooter)
		{
			std::vector<double> vecToks;
			get_tokens<double, std::string>(strLine, " \t", vecToks);

			if(vecToks.size() != m_vecQuantities.size())
			{
				log_warn("Loader: Line size mismatch.");

				// add zeros
				while(m_vecQuantities.size() > vecToks.size())
					vecToks.push_back(0.);
			}

			for(std::size_t iTok=0; iTok<vecToks.size(); ++iTok)
				m_vecData[iTok].push_back(vecToks[iTok]);
		}
	}
}


bool FileTrisp::Load(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return false;

#if !defined NO_IOSTR
	std::istream* pIstr = create_autodecomp_istream(ifstr);
	if(!pIstr) return false;
	std::unique_ptr<std::istream> ptrIstr(pIstr);
#else
	std::ifstream *pIstr = &ifstr;
#endif

	ReadHeader(*pIstr);
	ReadData(*pIstr);

	return true;
}

const FileInstr::t_vecVals& FileTrisp::GetCol(const std::string& strName) const
{
	return const_cast<FileTrisp*>(this)->GetCol(strName);
}

FileInstr::t_vecVals& FileTrisp::GetCol(const std::string& strName)
{
	static std::vector<double> vecNull;

	for(std::size_t i=0; i<m_vecQuantities.size(); ++i)
	{
		if(m_vecQuantities[i] == strName)
			return m_vecData[i];
	}

	return vecNull;
}

std::array<double,3> FileTrisp::GetSampleLattice() const
{
	t_mapParams::const_iterator iterA = m_mapParams.find("AS");
	t_mapParams::const_iterator iterB = m_mapParams.find("BS");
	t_mapParams::const_iterator iterC = m_mapParams.find("CS");

	double a = (iterA!=m_mapParams.end() ? str_to_var<double>(iterA->second) : 0.);
	double b = (iterB!=m_mapParams.end() ? str_to_var<double>(iterB->second) : 0.);
	double c = (iterC!=m_mapParams.end() ? str_to_var<double>(iterC->second) : 0.);

	return std::array<double,3>{{a,b,c}};
}

std::array<double,3> FileTrisp::GetSampleAngles() const
{
	t_mapParams::const_iterator iterA = m_mapParams.find("AA");
	t_mapParams::const_iterator iterB = m_mapParams.find("BB");
	t_mapParams::const_iterator iterC = m_mapParams.find("CC");

	double alpha = (iterA!=m_mapParams.end() ? tl::d2r(str_to_var<double>(iterA->second)) : M_PI/2.);
	double beta = (iterB!=m_mapParams.end() ? tl::d2r(str_to_var<double>(iterB->second)) : M_PI/2.);
	double gamma = (iterC!=m_mapParams.end() ? tl::d2r(str_to_var<double>(iterC->second)) : M_PI/2.);

	return std::array<double,3>{{alpha, beta, gamma}};
}

std::array<double,2> FileTrisp::GetMonoAnaD() const
{
	t_mapParams::const_iterator iterM = m_mapParams.find("DM");
	t_mapParams::const_iterator iterA = m_mapParams.find("DA");

	double m = (iterM!=m_mapParams.end() ? str_to_var<double>(iterM->second) : 3.355);
	double a = (iterA!=m_mapParams.end() ? str_to_var<double>(iterA->second) : 3.355);

	return std::array<double,2>{{m, a}};
}

std::array<bool, 3> FileTrisp::GetScatterSenses() const
{
	t_mapParams::const_iterator iterM = m_mapParams.find("SM");
	t_mapParams::const_iterator iterS = m_mapParams.find("SS");
	t_mapParams::const_iterator iterA = m_mapParams.find("SA");

	bool m = (iterM!=m_mapParams.end() ? (str_to_var<int>(iterM->second)>0) : 0);
	bool s = (iterS!=m_mapParams.end() ? (str_to_var<int>(iterS->second)>0) : 1);
	bool a = (iterA!=m_mapParams.end() ? (str_to_var<int>(iterA->second)>0) : 0);

	return std::array<bool,3>{{m, s, a}};
}

std::array<double, 3> FileTrisp::GetScatterPlane0() const
{
	t_mapParams::const_iterator iterX = m_mapParams.find("AX");
	t_mapParams::const_iterator iterY = m_mapParams.find("AY");
	t_mapParams::const_iterator iterZ = m_mapParams.find("AZ");

	double x = (iterX!=m_mapParams.end() ? str_to_var<double>(iterX->second) : 1.);
	double y = (iterY!=m_mapParams.end() ? str_to_var<double>(iterY->second) : 0.);
	double z = (iterZ!=m_mapParams.end() ? str_to_var<double>(iterZ->second) : 0.);

	return std::array<double,3>{{x,y,z}};
}

std::array<double, 3> FileTrisp::GetScatterPlane1() const
{
	t_mapParams::const_iterator iterX = m_mapParams.find("BX");
	t_mapParams::const_iterator iterY = m_mapParams.find("BY");
	t_mapParams::const_iterator iterZ = m_mapParams.find("BZ");

	double x = (iterX!=m_mapParams.end() ? str_to_var<double>(iterX->second) : 0.);
	double y = (iterY!=m_mapParams.end() ? str_to_var<double>(iterY->second) : 1.);
	double z = (iterZ!=m_mapParams.end() ? str_to_var<double>(iterZ->second) : 0.);

	return std::array<double,3>{{x,y,z}};
}

double FileTrisp::GetKFix() const
{
	const std::string strKey = IsKiFixed() ? "KI" : "KF";
	
	t_mapParams::const_iterator iter = m_mapParams.find(strKey);
	if(iter==m_mapParams.end())
	{
		tl::log_err("Cannot determine kfix.");
		return 0.;
	}

	return tl::str_to_var<double>(iter->second);
}

bool FileTrisp::IsKiFixed() const
{
	return 0;		// assume ckf
}

std::size_t FileTrisp::GetScanCount() const
{
	if(m_vecData.size() < 1)
		return 0;
	return m_vecData[0].size();
}

std::array<double, 5> FileTrisp::GetScanHKLKiKf(std::size_t i) const
{
	return FileInstr::GetScanHKLKiKf("QH", "QK", "QL", "E", i);
}

bool FileTrisp::MergeWith(const FileInstr* pDat) 
{ 
	return FileInstr::MergeWith(pDat); 
}

std::string FileTrisp::GetTitle() const { return ""; }
std::string FileTrisp::GetUser() const { return ""; }
std::string FileTrisp::GetLocalContact() const { return ""; }
std::string FileTrisp::GetScanNumber() const { return ""; }
std::string FileTrisp::GetSampleName() const { return ""; }
std::string FileTrisp::GetSpacegroup() const { return ""; }

std::vector<std::string> FileTrisp::GetScannedVars() const
{
	std::vector<std::string> vecScan;

	t_mapParams::const_iterator iter = m_mapParams.find("scan_vars");
	if(iter != m_mapParams.end())
		tl::get_tokens<std::string, std::string>(iter->second, " \t", vecScan);

	if(!vecScan.size())
	{
		tl::log_warn("Could not determine scan variable, using first column.");
		if(m_vecQuantities.size() >= 1)
			vecScan.push_back(m_vecQuantities[0]);
	}

	return vecScan;
}

std::string FileTrisp::GetCountVar() const { return "c1"; }
std::string FileTrisp::GetMonVar() const { return "mon"; }

std::string FileTrisp::GetScanCommand() const
{
	std::string str;
	t_mapParams::const_iterator iter = m_mapParams.find("scan_command");
	if(iter != m_mapParams.end())
		str = iter->second;
	return str;
}

std::string FileTrisp::GetTimestamp() const
{
	std::string str;
	t_mapParams::const_iterator iter = m_mapParams.find("scan_start_timestamp");
	if(iter != m_mapParams.end())
		str = iter->second;
	return str;
}


}


// -----------------------------------------------------------------------------


/*
// test
// gcc -DNO_IOSTR -o 0 file/loadinstr.cpp helper/log.cpp -std=c++11 -lstdc++ -lm
int main()
{
	//tl::FileFrm dat;
	//tl::FileMacs dat;
	tl::FileTrisp dat;
	//if(!dat.Load("/home/tweber/tmp/tst.dat"))
	//if(!dat.Load("/home/tweber/Messdaten/MACS_2014/data/Escan_31896.ng0"))
	if(!dat.Load("/home/tweber/Messdaten/trisp-15/data/sc77087.log"))
	{
		tl::log_err("Cannot load data file.");
		return -1;
	}

	std::array<double,3> latt = dat.GetSampleLattice();
	std::cout << latt[0] << ", " << latt[1] << ", " << latt[2] << std::endl;

	std::cout << "kfix = " << dat.GetKFix() << std::endl;

	return 0;
}
*/
