/*
 * Load .dat files
 *
 * @author: Tobias Weber
 * @date: April 2012
 */

#include "loadtxt.h"
#include "../helper/file.h"
#include "../helper/string.h"
#include "../helper/misc.h"
#include "../helper/math.h"
#include "../helper/log.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <boost/algorithm/minmax_element.hpp>
#include <limits>

static bool CleanString(std::string& strLine)
{
	bool bModified = 0;

	while(strLine.find("nan")!=std::string::npos)
	{
		find_and_replace(strLine, std::string("nan"), std::string("0"));
		bModified = 1;
	}
	while(strLine.find("inf")!=std::string::npos)
	{
		find_and_replace(strLine, std::string("inf"), std::string("0"));
		bModified = 1;
	}

	return bModified;
}

static void get_limits_from_str(const std::string& str, double &dMin, double &dMax, bool &bLog)
{
        std::istringstream istr(str);

        std::vector<double> vec;
        while(!istr.eof())
        {
                double dval;
                istr >> dval;

                vec.push_back(dval);
        }

        std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator>
                        minmax = boost::minmax_element(vec.begin(), vec.end());
        dMin = *minmax.first;
        dMax = *minmax.second;

        bLog = 0;

        //std::cout << "min: " << dMin << ", max: " << dMax << std::endl;
}

LoadTxt::LoadTxt(const char* pcFile, bool bOnlyHeader, bool bVerbose)
					: m_uiColLen(0), m_bLoadOnlyHeader(bOnlyHeader),
					  m_bVerbose(bVerbose)
{ Load(pcFile); }

LoadTxt::~LoadTxt() { Unload(); }

const TxtType LoadTxt::GetFileType() const
{
	for(const std::string& strComm : m_vecComments)
	{
		if(strComm.find("NICOS data file") != std::string::npos)
			return NICOS_DATA;
	}

	return MCSTAS_DATA;
}

void LoadTxt::StrTrim(std::string& str)
{
	if(str.length()==0) return;

	bool bIsTripleComment = false;
	const char* pcTripleComment = "###";
	std::size_t iFirstTripleComment = str.find(pcTripleComment);
	if(iFirstTripleComment!=std::string::npos)
	{
		std::string strTripleComment = str.substr(iFirstTripleComment+3, std::string::npos);
		trim(strTripleComment);

		bIsTripleComment = true;
		m_vecComments.push_back(strTripleComment);
	}

	const char* pcComment = "#";
	std::size_t iFirstComment = str.find(pcComment);
	if(iFirstComment!=std::string::npos)
	{
		// first: get any information in the comment
		// Mcstas header parameters look like this: "Param: use_div=0"
		std::string strComm = str.substr(iFirstComment+1, std::string::npos);
		std::size_t iSeparator = strComm.find(":");
		if(iSeparator != std::string::npos && iSeparator < strComm.length())
		{
			std::string strKey = strComm.substr(0, iSeparator);
			std::string strVal = strComm.substr(iSeparator+1, std::string::npos);

			trim(strKey);
			trim(strVal);

			// if key is "Param" then the value itself is a "Key = Value" pair!
			if(is_equal(strKey, std::string("param"), false))
			{
				std::size_t iParamSep = strVal.find("=");
				if(iParamSep != std::string::npos && iParamSep < strVal.length())
				{
					std::string strParamKey = strVal.substr(0, iParamSep);
					std::string strParamVal = strVal.substr(iParamSep+1, std::string::npos);

					trim(strParamKey);
					trim(strParamVal);

					//std::cout << "Param Key=\"" << strParamKey << "\", Value=\"" << strParamVal << "\"." << std::endl;

					strKey = std::string("param_") + strParamKey;
					strVal = strParamVal;
				}
			}

			//std::cout << "Key=\"" << strKey << "\", Value=\"" << strVal << "\"." << std::endl;

			if(strKey!="" && strVal!="")
				m_mapComm[strKey].push_back(strVal);
		}
		else
		{
			if(!bIsTripleComment)
			{
				// Strings not fitting anywhere
				trim(strComm);
				m_vecAuxStrings.push_back(strComm);
			}
		}

		// second: remove the comment from the string
		str = str.substr(0, iFirstComment);
	}
	::trim(str);
}


bool LoadTxt::Load(const char* pcFile)
{
	Unload();
	if(!pcFile) return 0;

	bool bHasAxesLines = false;
	bool bTranspose = false;
	bool bIsResData = false;


	std::string strExt = get_fileext(std::string(pcFile));
	if(is_equal(strExt, std::string("sqw")))
	{
		m_mapComm["type"].push_back("array_2d");
		m_mapComm["subtype"].push_back("tobisown");		// 1 data block
		m_mapComm["zlabel"].push_back("S(q,w)");

		bHasAxesLines = true;
		bTranspose = true;
	}
	else if(is_equal(strExt, std::string("res")))
	{
		bIsResData = true;
	}


	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
		return 0;

	unsigned int uiFileSize = get_file_size(ifstr);

	std::vector<std::vector<double> > vecLines;
	unsigned int uiLine = 0;
	unsigned int uiLineWithoutComment = 0;

	std::vector<double> vecLine;
	int iLineSize=-1;
	int iLastPercentLoaded=-1;
	while(!ifstr.eof())
	{
		std::string strLine;
		std::getline(ifstr, strLine);
		++uiLine;

		if(m_bVerbose)
		{
			unsigned int uiFilePos = get_file_pos(ifstr);
			int iPercentLoaded = int(100.*double(uiFilePos)/double(uiFileSize));
			if(ifstr.eof())
				iPercentLoaded = 100;

			if(iPercentLoaded != iLastPercentLoaded)
			{
				std::cout << "\rLoading... " << iPercentLoaded << "%";
				iLastPercentLoaded = iPercentLoaded;

				if(iPercentLoaded==100)
					std::cout << "\n";

				std::cout.flush();
			}
		}

		StrTrim(strLine);
		if(strLine.size()==0) continue;
		++uiLineWithoutComment;

		// end of head comment section?
		if(strLine.size()!=0 && m_bLoadOnlyHeader)
			break;

		// first two lines describe axes
		if(bHasAxesLines && (uiLineWithoutComment==1 || uiLineWithoutComment==2))
		{
			const char* pcLabel[] = {"ylabel", "xlabel"};
			const char* pcLabelQty[] = {"q [1/A]", "E [meV]"};
			const char* pcLimits[] = {"ylimits", "xlimits"};
			const char* pcLog[] = {"ylog", "xlog"};

			unsigned int iIdx = uiLineWithoutComment-1;
			if(bTranspose)
				iIdx = (iIdx==0)?1:0;

			m_mapComm[pcLabel[iIdx]].push_back(pcLabelQty[uiLineWithoutComment-1]);

			double dMin=0., dMax=0.;
			bool bLog=0;

			get_limits_from_str(strLine, dMin, dMax, bLog);

			std::ostringstream ostrLimits;
			ostrLimits << dMin << " " << dMax;
			m_mapComm[pcLimits[iIdx]].push_back(ostrLimits.str());

			std::ostringstream ostrLog;
			ostrLog << bLog;
			m_mapComm[pcLog[iIdx]].push_back(ostrLog.str());

			continue;
		}

		if(CleanString(strLine))
			log_warn("Replaced \"inf\"/\"nan\" with \"0\" in line ",uiLine);

		vecLine.clear();
		get_tokens<double>(strLine, std::string(" \t"), vecLine);

		if(iLineSize==-1)
		{
			iLineSize = vecLine.size();
			vecLine.reserve(iLineSize);
		}
		if(vecLine.size() != (unsigned int)(iLineSize))
		{
			log_warn("Line ", uiLine, " has wrong size!",
					  " Expected ", iLineSize,
					  " tokens, got ", vecLine.size(), ".");

			//return 0;

			if(int(vecLine.size()) > iLineSize)
			{
				log_warn("Removing last elements.");
				vecLine.resize(iLineSize);
			}
			else if(int(vecLine.size()) < iLineSize)
			{
				log_warn("Inserting zeros.");

				unsigned int iOldSize = vecLine.size();
				vecLine.resize(iLineSize);

				for(unsigned int iElem=iOldSize; int(iElem)<iLineSize; ++iElem)
					vecLine[iElem] = 0.;
			}
		}

		vecLines.push_back(vecLine);
	}


	if(!m_bLoadOnlyHeader)
	{
		if(iLineSize<=0)
		{
			log_warn("No data in file \"", pcFile, "\".");
			return 0;
		}


		if(bTranspose)
		{
			m_vecColumns.reserve(vecLines.size());

			for(unsigned int iLine=0; iLine<vecLines.size(); ++iLine)
			{
				double *pd = new double[iLineSize];
				for(unsigned int iCol=0; iCol<(unsigned int)iLineSize; ++iCol)
					pd[iCol] = vecLines[iLine][iCol];

				m_vecColumns.push_back(pd);
			}

			m_uiColLen = iLineSize;
		}
		else
		{
			m_uiColLen = vecLines.size();
			m_vecColumns.reserve(iLineSize);
			for(unsigned int iCol=0; iCol<(unsigned int)iLineSize; ++iCol)
			{
				double *pd = new double[m_uiColLen];
				for(unsigned int iLine=0; iLine<m_uiColLen; ++iLine)
					pd[iLine] = vecLines[iLine][iCol];

				m_vecColumns.push_back(pd);
			}
		}
	}


	// data is from resolution sample/monitor
	std::string strComp;
	// TODO: find a better way to identify resolution data
	if(GetMapString("ylabel", strComp) &&
		(strComp=="ki_x ki_y ki_z kf_x kf_y kf_z x y z p_i p_f" ||
		 strComp=="ki_x ki_y ki_z kf_x kf_y kf_z x y z p_i p_f sx sy sz"))
		bIsResData = true;

	if(bIsResData)
		SetMapString("type", "resdata");

	m_strFileName = pcFile;
	return 1;
}

bool LoadTxt::Save(const char* pcFile, bool bSaveComments) const
{
	std::ofstream ofstr(pcFile);
	if(!ofstr.is_open())
		return false;

	ofstr << std::scientific;
	//ofstr.precision(16);
	
	if(bSaveComments)
	{
		for(const auto& pairComm : m_mapComm)
		{
			for(unsigned int uiVal=0; uiVal<pairComm.second.size(); ++uiVal)
				ofstr << "# " << pairComm.first
						<< ": " << pairComm.second[uiVal] << "\n";
		}
	}

	for(unsigned int iY=0; iY<GetColLen(); ++iY)
	{
		for(unsigned int iX=0; iX<GetColCnt(); ++iX)
			ofstr << GetColumn(iX)[iY] << " ";
		ofstr << "\n";
	}
	return true;
}

void LoadTxt::Unload()
{
	m_strFileName = "";
	m_mapComm.clear();

	for(double*& pdCol : m_vecColumns)
	{
		if(pdCol)
		{
			delete[] pdCol;
			pdCol = 0;
		}
	}

	m_vecColumns.clear();
	m_uiColLen=0;
}

std::map<std::string, std::string> LoadTxt::GetCommMapSingle() const
{
	typedef std::map<std::string, std::string> tmap;
	tmap mapRet;

	for(const t_mapComm::value_type& val : m_mapComm)
	{
		std::string strKey;
		if(val.second.size() > 0)
			strKey = val.second[0];

		mapRet.insert(tmap::value_type(val.first, strKey));
	}

	return mapRet;
}

void LoadTxt::GetMinMax(double& dMin, double& dMax) const
{
	dMin = std::numeric_limits<double>::max();
	dMax = -dMin;
	
	for(unsigned int uiCol=0; uiCol<GetColCnt(); ++uiCol)
	{
		const double* pCol = GetColumn(uiCol);
		std::pair<const double*, const double*> minmax = boost::minmax_element(pCol, pCol+GetColLen());
		
		double dColMin = *minmax.first;
		double dColMax = *minmax.second;
		
		dMin = std::min(dMin, dColMin);
		dMax = std::max(dMax, dColMax);
	}
}

void LoadTxt::for_each(void (*pfkt)(double*))
{
	for(unsigned int uiCol=0; uiCol<GetColCnt(); ++uiCol)
		for(unsigned int uiVal=0; uiVal<GetColLen(); ++uiVal)
			pfkt(&m_vecColumns[uiCol][uiVal]);
}

bool LoadTxt::GetMapString(const std::string& strKey, std::string& strVal) const
{
	bool bOk = true;

	LoadTxt::t_mapComm::const_iterator iter = m_mapComm.find(strKey);
	if(iter!=m_mapComm.end() && (*iter).second.size()>=1)
		strVal = (*iter).second[0];
	else
		bOk = false;

	convert_mcstas_string(strVal);
	return bOk;
}

void LoadTxt::SetMapString(const std::string& strKey, const std::string& strVal)
{
	LoadTxt::t_mapComm::iterator iter = m_mapComm.find(strKey);
	if(iter!=m_mapComm.end() && (*iter).second.size()>=1)
		(*iter).second[0] = strVal;
	else
		m_mapComm[strKey].push_back(strVal);

	/*
	std::cout << strKey << " -> " << strVal << std::endl;
	std::string strVal_new;
	GetMapString(strKey, strVal_new);
	std::cout << strKey << " = " << strVal_new << std::endl;
	*/
}

std::ostream& operator<<(std::ostream& ostr, LoadTxt& txt)
{
	for(unsigned int uiCol=0; uiCol<txt.GetColCnt(); ++uiCol)
	{
		ostr << "Column " << uiCol+1 << ": ";
		for(unsigned int uiVal=0; uiVal<txt.GetColLen(); ++uiVal)
			ostr << txt.GetColumn(uiCol)[uiVal] << " ";
		ostr << "\n";
	}
	return ostr;
}


//------------------------------------------------------------------------------

void convert_mcstas_string(std::string& str)
{
	find_and_replace(str, std::string("\\gf"), std::string("phi"));
	find_and_replace(str, std::string("\\gh"), std::string("theta"));
	find_and_replace(str, std::string("\\gm"), std::string("u"));
}

//------------------------------------------------------------------------------


bool McData::GetTitle(std::string& strTitle) const
{
	return GetString("title", strTitle);
}

// not a mcstas file, but a personal one?
bool McData::IsOwnDataFile() const
{
	std::string strVal;
	bool bExists = GetString("subtype", strVal);
	if(bExists && strVal == "tobisown")
		return true;
		
	return false;
}

bool McData::GetLogScale(bool &bXLog, bool &bYLog) const
{
	// default to no log if field not found
	bXLog = bYLog = 0;
	
	std::string strVal;
	if(GetString("xylog", strVal))
	{
		std::istringstream istr(strVal);
		istr >> bXLog >> bYLog;

		return true;
	}

	bool bOk = false;

	// if not found look for these fields
	if(GetString("xlog", strVal))
	{
		std::istringstream istr(strVal);
		istr >> bXLog;
		bOk = true;
	}

	if(GetString("ylog", strVal))
	{
		std::istringstream istr(strVal);
		istr >> bYLog;
		bOk = true;
	}
	return bOk;
}

bool McData::GetFit(std::string& strFit) const
{
	return GetString("fit", strFit);
}


//------------------------------------------------------------------------------


DataRes::DataRes(const LoadTxt& data) : McData(data)
{
	const unsigned int iResColCnt = 11;
	m_bOk = true;

	m_iDim = m_data.GetColLen();
	m_iCols = m_data.GetColCnt();

	if(m_data.GetColCnt()!=iResColCnt && m_data.GetColCnt()!=(iResColCnt+3))
	{
		log_err("Wrong number of columns (need ", iResColCnt, " or ", iResColCnt+3,  ").");
		m_bOk = false;
	}
}

const double* DataRes::GetColumn(unsigned int iIdx) const
{
	if(iIdx < m_iCols)
		return m_data.GetColumn(iIdx);
	else
		return 0;
}


//------------------------------------------------------------------------------


Data1D::Data1D(const LoadTxt& data) : McData(data)
{
	m_iDim = m_data.GetColLen();

	std::string strVal;
	if(GetString("type", strVal))		// e.g.: "type: array_1d(32)"
	{
		bool bHasArrayDim = 0;
		std::size_t idx = strVal.find('(');
		if(idx != std::string::npos)
			bHasArrayDim = 1;

		m_iDim = m_data.GetColLen();

		if(bHasArrayDim)
		{
			std::string strDim = strVal.substr(idx+1);
			trim(strDim);

			if(strDim.size() != 0)
			{
				std::istringstream istr(strDim);
				istr >> m_iDim;

				//std::cout << "dim: " << m_iDim << std::endl;
			}
		}

		if(m_data.GetColLen() != m_iDim)
		{
			m_iDim = m_data.GetColLen();
			log_warn("Mismatch in dimension, assuming ", m_iDim, ".");
		}

		m_bOk = true;
	}
	else
	{
		// we can still continue even without the header
		//std::cerr << "Warning: No dimension specified in header." << std::endl;
		m_bOk = false;
	}
}

const double* Data1D::GetColumn(unsigned int iIdx) const
{
	return m_data.GetColumn(iIdx);
}

double Data1D::GetMin(unsigned int iColIdx, unsigned int *pMinIdx) const
{
	const double* pCol =m_data.GetColumn(iColIdx);
	const double *pMin = std::min_element(pCol, pCol + m_data.GetColLen());

	if(pMinIdx)
		*pMinIdx = pMin - pCol;

	return *pMin;
}

double Data1D::GetMax(unsigned int iColIdx, unsigned int *pMaxIdx) const
{
	const double* pCol =m_data.GetColumn(iColIdx);
	const double *pMax = std::max_element(pCol, pCol + m_data.GetColLen());

	if(pMaxIdx)
		*pMaxIdx = pMax - pCol;

	return *pMax;
}

bool Data1D::GetXLimits(double& dXMin, double& dXMax) const
{
	std::string strVal;
	if(GetString("xlimits", strVal))
	{
		std::istringstream istr(strVal);
		istr >> dXMin;
		istr >> dXMax;

		return true;
	}
	return false;
}

bool Data1D::GetYLimits(double& dYMin, double& dYMax) const
{
	std::string strVal;
	if(GetString("ylimits", strVal))
	{
		std::istringstream istr(strVal);
		istr >> dYMin;
		istr >> dYMax;

		return true;
	}
	return false;
}

void Data1D::CalcLimits(double& dXMin, double& dXMax, int iCol) const
{
	// no column given, use x
	if(iCol<0)
	{
		int iX=0, iY=1, iYErr=2;
		GetColumnIndices(iX, iY, iYErr);

		iCol = iX;
	}

	const double *pX = GetColumn(iCol);
	unsigned int uiDim = GetDim();
	
	std::pair<const double*, const double*> pair_minmax =
								boost::minmax_element(pX, pX+uiDim);
	
	dXMin = *pair_minmax.first;
	dXMax = *pair_minmax.second;
}

bool Data1D::GetColLabel(unsigned int iCol, std::string& strLabel) const
{
	std::ostringstream strColLab;
	strColLab << "col-" << iCol << "-label";
	return GetString(strColLab.str(), strLabel);
}

bool Data1D::GetLabels(std::string& strLabelX, std::string& strLabelY) const
{	
	// first see if explicit column labels exist
	int iX=0, iY=1, iYErr=2;
	GetColumnIndices(iX,iY,iYErr);
	
	bool bHasX = GetColLabel(iX, strLabelX);
	bool bHasY = GetColLabel(iY, strLabelY);


	bool bOk = true;

	// otherwise use the standard labels
	if(!bHasX)
		if(!GetString("xlabel", strLabelX))
			bOk = false;

	if(!bHasY)
		if(!GetString("ylabel", strLabelY))
			bOk = false;
	
	return bOk;
}

bool Data1D::GetColumnIndices(int& iX, int& iY, int& iYErr) const
{
	std::string strVal;
	bool bHasKeyX = GetString("which-col-is-x", strVal);
	if(bHasKeyX)
		std::istringstream(strVal) >> iX;

	bool bHasKeyY = GetString("which-col-is-y", strVal);
	if(bHasKeyY)
		std::istringstream(strVal) >> iY;

	bool bHasKeyYErr = GetString("which-col-is-yerr", strVal);
	if(bHasKeyYErr)
		std::istringstream(strVal) >> iYErr;

	return bHasKeyX || bHasKeyY || bHasKeyYErr;
}

void Data1D::NormY(double dNorm, double dNormErr)
{
	int iX=0, iY=1, iYErr=2;
	GetColumnIndices(iX, iY, iYErr);

	double* pY = const_cast<double*>(GetColumn(iY));
	double* pYErr = const_cast<double*>(GetColumn(iYErr));

	if(!pY || GetDim()==0)
		return;

	double dY0 = dNorm;
	double dY0_err = dNormErr;
	for(unsigned int i=0; i<GetDim(); ++i)
	{
		pY[i] /= dY0;
		pYErr[i] = sqrt( (pYErr[i]/dY0)*(pYErr[i]/dY0) +
						 (-pY[i]/(dY0*dY0)*dY0_err)*(-pY[i]/(dY0*dY0)*dY0_err) );
	}
}

/*
 * for intensity and intensity error see chapter 4.2 in the mcstas manual
 */
double Data1D::GetIntensity(unsigned int iICol, unsigned int iIErrCol, unsigned int iCountCol) const
{
	const double* pI = GetColumn(iICol);
	const double* pI_err = GetColumn(iIErrCol);
	const double* pCount = GetColumn(iCountCol);

	double dInt = 0.;
	for(unsigned int i=0; i<GetDim(); ++i)
		dInt += pI[i];
	return dInt;
}

double Data1D::GetIntensityError(unsigned int iICol, unsigned int iIErrCol, unsigned int iCountCol) const
{
	const double* pI = GetColumn(iICol);
	const double* pI_err = GetColumn(iIErrCol);
	const double* pCount = GetColumn(iCountCol);
	
	double dErr = 0.;
	for(unsigned int i=0; i<GetDim(); ++i)
		dErr += pI_err[i]*pI_err[i];
	dErr = sqrt(dErr);
	return dErr;
}

void Data1D::ForEach(unsigned int iColIdx, double (*pFkt)(double))
{
	double* pCol = const_cast<double*>(GetColumn(iColIdx));
	if(!pCol)
		return;

	for(unsigned int i=0; i<GetDim(); ++i)
		pCol[i] = pFkt(pCol[i]);
}

//------------------------------------------------------------------------------


Data2D::Data2D(const LoadTxt& data)
		: McData(data), m_iNumBlocks(3), m_bRecalcZLimits(0)
{
	if(IsOwnDataFile())
		m_iNumBlocks = 1;

	m_iXDim = m_data.GetColCnt();
	m_iYDim = m_data.GetColLen()/m_iNumBlocks;

	std::string strVal;
	if(GetString("type", strVal) && strVal.length()>8)
	{
		const std::string strDim = strVal.substr(9);
		std::istringstream istr(strDim);

		char c;
		// data now looks like this: "123, 234, 345"
		istr >> m_iXDim >> c >> m_iYDim;
		//std::cout << m_iXDim << " " << m_iYDim << " " << m_iTDim << std::endl;

		if(m_data.GetColCnt()!=m_iXDim)
		{
			m_iXDim = m_data.GetColCnt();
			log_warn("Mismatch in x dimension, assuming ", m_iXDim, ".");
		}
		if(m_data.GetColLen()/m_iNumBlocks!=m_iYDim)
		{
			m_iYDim = m_data.GetColLen()/m_iNumBlocks;
			log_warn("Mismatch in y dimension, assuming ", m_iYDim, ".");
		}

		m_bOk = true;
	}
	else
	{
		// we can still continue even without the header
		//std::cerr << "Warning: Unknown dimensions of 2d data." << std::endl;
		m_bOk = false;
	}
}

void Data2D::ForEach(unsigned int iBlockIdx, double (*pfkt)(double))
{
	if(iBlockIdx>=m_iNumBlocks)
		return;
	
	for(unsigned int iY=0; iY<GetYDim(); ++iY)
		for(unsigned int iX=0; iX<GetXDim(); ++iX)
		{
			double* pdval = const_cast<double*>(m_data.GetColumn(iX) + m_iYDim*iBlockIdx + iY);
			*pdval = pfkt(*pdval);
		}

	m_bRecalcZLimits = 1;
}

void Data2D::Add(unsigned int iBlockIdx, double d)
{
	if(iBlockIdx>=m_iNumBlocks)
		return;

	for(unsigned int iY=0; iY<GetYDim(); ++iY)
		for(unsigned int iX=0; iX<GetXDim(); ++iX)
		{
			double* pdval = const_cast<double*>(m_data.GetColumn(iX) + m_iYDim*iBlockIdx + iY);
			*pdval += d;
		}

	m_bRecalcZLimits = 1;
}

void Data2D::Mod0to2Pi(unsigned int uiBlock)
{
	if(uiBlock>=m_iNumBlocks)
		return;
	
	for(unsigned int iY=0; iY<GetYDim(); ++iY)
		for(unsigned int iX=0; iX<GetXDim(); ++iX)
		{
			double* pdval = const_cast<double*>(m_data.GetColumn(iX) + m_iYDim*uiBlock + iY);
			*pdval = fmod(*pdval, 2.*M_PI);

			if(*pdval < 0.)
				*pdval += 2.*M_PI;
		}
}

double Data2D::GetCenterVal(unsigned int iBlock) const
{
	unsigned int iX = GetXDim()/2;
	unsigned int iY = GetYDim()/2;

	return GetBlockVal(iBlock, iX, iY);
}

double Data2D::GetBlockVal(unsigned int iBlock, unsigned int iX, unsigned int iY) const
{
	if(iBlock>=m_iNumBlocks)
		return 0.;
	if(iX>=m_iXDim || iY>=m_iYDim)
		return 0.;
	
	const double* pcol = m_data.GetColumn(iX) + m_iYDim*iBlock;
	return pcol[iY];
}

void Data2D::GetBlockMinMax(unsigned int iBlock, double& dMin, double& dMax) const
{
	dMin = std::numeric_limits<double>::max();
	dMax = -dMin;
	
	for(unsigned int uiCol=0; uiCol<m_iXDim; ++uiCol)
	{
		const double* pCol = m_data.GetColumn(uiCol) + m_iYDim*iBlock;
		std::pair<const double*, const double*> minmax = boost::minmax_element(pCol, pCol+m_iYDim);
		
		double dColMin = *minmax.first;
		double dColMax = *minmax.second;
		
		dMin = std::min(dMin, dColMin);
		dMax = std::max(dMax, dColMax);
	}
}

double Data2D::GetVal(unsigned int iX, unsigned int iY) const
{ return GetBlockVal(0, iX, iY); }
double Data2D::GetErr(unsigned int iX, unsigned int iY) const
{ return GetBlockVal(1, iX, iY); }
double Data2D::GetCount(unsigned int iX, unsigned int iY) const
{ return GetBlockVal(2, iX, iY); }

void Data2D::GetValMinMax(double& dMin, double& dMax) const
{ GetBlockMinMax(0, dMin, dMax); }
void Data2D::GetErrMinMax(double& dMin, double& dMax) const
{ GetBlockMinMax(1, dMin, dMax); }


bool Data2D::GetLimits(double& dXMin, double& dXMax,
			   double& dYMin, double& dYMax,
			   double& dZMin, double& dZMax) const
{
	std::string strVal;
	if(GetString("xylimits", strVal))
	{
		std::istringstream istr(strVal);
		istr >> dXMin >> dXMax >> dYMin >> dYMax >> dZMin >> dZMax;
		
		return true;
	}

	bool bOk = false;
	
	// if not found look for these fields
	if(GetString("xlimits", strVal))
	{
		std::istringstream istr(strVal);
		istr >> dXMin;
		istr >> dXMax;

		bOk = true;
	}
	
	if(GetString("ylimits", strVal))
	{
		std::istringstream istr(strVal);
		istr >> dYMin;
		istr >> dYMax;

		bOk = true;
	}
	
	if(GetString("zlimits", strVal))
	{
		std::istringstream istr(strVal);
		istr >> dZMin;
		istr >> dZMax;

		bOk = true;

		if(m_bRecalcZLimits)
			GetValMinMax(dZMin, dZMax);
	}
	
	return bOk;
}

bool Data2D::GetLabels(std::string& strLabelX, std::string& strLabelY, std::string& strLabelZ) const
{
	bool bOk = true;

	if(!GetString("xlabel", strLabelX))
		bOk = false;

	if(!GetString("ylabel", strLabelY))
		bOk = false;

	if(!GetString("zlabel", strLabelZ))
		bOk = false;

	return bOk;
}

bool Data2D::ExtractColOrRow(int iCol, int iRow, const char* pcOutFile, const Data2D* pErrors)
{
	std::ofstream ofstr(pcOutFile);
	if(!ofstr.is_open())
	{
		log_err("Could not open \"", pcOutFile, "\".");
		return 0;
	}

	bool bWantRow=0, bWantCol=0;
	if(iCol>=0) bWantCol=1;
	if(iRow>=0) bWantRow=1;

	if(bWantRow==0 && bWantCol==0)
	{
		log_err("Neither column nor row selected.");
		return 0;
	}
	if(bWantRow==1 && bWantCol==1)
	{
		log_err("Both column and row selected.");
		return 0;
	}


	std::string strLabX="x", strLabY="y", strLabZ="z";
	GetLabels(strLabX, strLabY, strLabZ);
	
	double dXMin, dXMax, dYMin, dYMax, dZMin, dZMax;
	GetLimits(dXMin, dXMax, dYMin, dYMax, dZMin, dZMax);

	bool bLogX=0, bLogY=0;
	GetLogScale(bLogX, bLogY);

	unsigned int uiXDim = GetXDim(),
				 uiYDim = GetYDim();

	unsigned int iDim = (bWantRow ? uiXDim : uiYDim);
	double dMin = (bWantRow ? dXMin : dYMin);
	double dMax = (bWantRow ? dXMax : dYMax);
	bool bLog = (bWantRow ? bLogX : bLogY );

	unsigned int iDim_other = (bWantRow ? uiYDim : uiXDim);
	double dMin_other = (bWantRow ? dYMin : dXMin);
	double dMax_other = (bWantRow ? dYMax : dXMax);
	bool bLog_other = (bWantRow ? bLogY : bLogX);

	std::ostringstream ostrTitle;
	ostrTitle << (bWantRow ? strLabY : strLabX);
	ostrTitle << " = "
			  << tic_trafo(iDim_other, dMin_other, dMax_other, bLog_other, double(bWantRow?iRow:iCol))
			  << "\n";

	ofstr << "# title: " << ostrTitle.str();
	ofstr << "# xlabel: " << (bWantRow?strLabX:strLabY) << "\n";
	ofstr << "# ylabel: " << strLabZ << "\n";
	ofstr << "# which-col-is-x: 0\n";
	ofstr << "# which-col-is-y: 1\n";
	ofstr << "# which-col-is-yerr: 2\n";

	double dValMin=0., dValMax=1.;
	GetValMinMax(dValMin, dValMax);

	for(unsigned int uiPix=0; uiPix<iDim; ++uiPix)
	{
		double dXval = tic_trafo(iDim, dMin, dMax, bLog, double(uiPix));

		double dYval, dYerr;
		if(bWantRow)
		{
			dYval = GetVal(uiPix, iRow);
			dYerr = GetErr(uiPix, iRow);

			if(pErrors)
				dYerr = pErrors->GetVal(uiPix, iRow);
		}
		else
		{
			dYval = GetVal(iCol, uiPix);
			dYerr = GetErr(iCol, uiPix);

			if(pErrors)
				dYerr = pErrors->GetVal(iCol, uiPix);
		}

		// assume an error
		if(pErrors==0 && IsOwnDataFile())
		{
			dYerr = dValMax / 100.;
			//std::cerr << "Warning: Assuming an error of " << dYerr << "."
			//			<< std::endl;
		}

		ofstr << dXval << " " << dYval << " " << dYerr << "\n";
	}

	ofstr.flush();
	ofstr.close();

	return 1;
}


//------------------------------------------------------------------------------

Data3D::Data3D(const LoadTxt& data)
		: McData(data), m_iNumBlocks(3),
		  m_iXDim(0), m_iYDim(0), m_iTDim(0)
{
	if(IsOwnDataFile())
		m_iNumBlocks = 1;

	std::string strVal;
	if(GetString("type", strVal) && strVal.length()>8)
	{
		const std::string strDim = strVal.substr(9);
		std::istringstream istr(strDim);

		char c;
		// data now looks like this: "123, 234, 345"
		istr >> m_iXDim >> c >> m_iYDim >> c >> m_iTDim;
		//std::cout << m_iXDim << " " << m_iYDim << " " << m_iTDim << std::endl;

		if(m_data.GetColCnt()!=m_iXDim)
		{
			m_iXDim = m_data.GetColCnt();
			log_warn("Mismatch in x dimension, assuming ", m_iXDim, ".");
		}
		if(m_data.GetColLen()/m_iNumBlocks/m_iTDim!=m_iYDim)
		{
			m_iYDim = m_data.GetColLen()/m_iNumBlocks/m_iTDim;
			log_warn("Mismatch in y dimension, assuming ", m_iYDim, ".");
		}

		m_bOk = true;
	}
	else
	{
		log_err("Unknown dimensions of 3d data.");
		m_bOk = false;
	}
}

void Data3D::ForEach(unsigned int iBlockIdx, double (*pfkt)(double))
{
	if(iBlockIdx>=m_iNumBlocks)
		return;

	for(unsigned int iY=0; iY<GetYDim(); ++iY)
		for(unsigned int iX=0; iX<GetXDim(); ++iX)
			for(unsigned int iT=0; iT<GetTDim(); ++iT)
			{
				double* pdval = const_cast<double*>(m_data.GetColumn(iX)
														+ m_iYDim*m_iTDim*iBlockIdx
														+ iY*m_iTDim + iT);
				*pdval = pfkt(*pdval);
			}
}

double Data3D::GetBlockVal(unsigned int iBlock, unsigned int iX, unsigned int iY, unsigned int iT) const
{
	if(iBlock>=m_iNumBlocks)
		return 0.;
	if(iX>=m_iXDim || iY>=m_iYDim || iT>=m_iTDim)
		return 0.;

	const double* pcol = m_data.GetColumn(iX) + m_iYDim*m_iTDim*iBlock;
	return pcol[iY*m_iTDim + iT];
}

double Data3D::GetVal(unsigned int iX, unsigned int iY, unsigned int iT) const
{ return GetBlockVal(0, iX, iY, iT); }
double Data3D::GetErr(unsigned int iX, unsigned int iY, unsigned int iT) const
{ return GetBlockVal(1, iX, iY, iT); }
double Data3D::GetCount(unsigned int iX, unsigned int iY, unsigned int iT) const
{ return GetBlockVal(2, iX, iY, iT); }


const double* Data3D::GetBlockT(unsigned int iBlock, unsigned int iX, unsigned int iY) const
{
	if(iBlock>=m_iNumBlocks)
		return 0;
	if(iX>=m_iXDim || iY>=m_iYDim)
		return 0;
	
	const double* pcol = m_data.GetColumn(iX) + m_iYDim*m_iTDim*iBlock;
	return pcol + iY*m_iTDim;
}

const double* Data3D::GetT(unsigned int iX, unsigned int iY) const
{ return GetBlockT(0, iX, iY); }
const double* Data3D::GetTErr(unsigned int iX, unsigned int iY) const
{ return GetBlockT(1, iX, iY); }
const double* Data3D::GetTCount(unsigned int iX, unsigned int iY) const
{ return GetBlockT(2, iX, iY); }


bool Data3D::GetLimits(double& dXMin, double& dXMax,
						double& dYMin, double& dYMax,
						double& dTMin, double& dTMax) const
{
	std::string strVal;
	if(GetString("xylimits", strVal))
	{
		std::istringstream istr(strVal);
		istr >> dXMin >> dXMax >> dYMin >> dYMax >> dTMin >> dTMax;

		return true;
	}

	bool bOk = false;

	// if not found look for these fields
	if(GetString("xlimits", strVal))
	{
		std::istringstream istr(strVal);
		istr >> dXMin;
		istr >> dXMax;

		bOk = true;
	}

	if(GetString("ylimits", strVal))
	{
		std::istringstream istr(strVal);
		istr >> dYMin;
		istr >> dYMax;

		bOk = true;
	}

	if(GetString("zlimits", strVal))
	{
		std::istringstream istr(strVal);
		istr >> dTMin;
		istr >> dTMax;

		bOk = true;
	}
	
	return bOk;
}

bool Data3D::GetLabels(std::string& strLabelX, std::string& strLabelY, std::string& strLabelT) const
{
	bool bOk = true;
	
	if(!GetString("xlabel", strLabelX))
		bOk = false;
	
	if(!GetString("ylabel", strLabelY))
		bOk = false;

	if(!GetString("zlabel", strLabelT))
		bOk = false;
	
	return bOk;
}

std::string Data3D::GetBlockName(unsigned int uiBlock) const
{
	switch(uiBlock)
	{
		case 0: return "intensities";
		case 1: return "intensity errors";
		case 2: return "counts";
		default: return "<unknown block>";
	}
}

bool Data3D::SaveTChan(const char* pcFile, unsigned int iX, unsigned int iY) const
{
	double dXMin, dXMax, dYMin, dYMax, dTMin, dTMax;
	bool bHasLimits = GetLimits(dXMin, dXMax, dYMin, dYMax, dTMin, dTMax);
	
	double *px = new double[GetTDim()];
	autodeleter<double> _adel(px, true);

	if(bHasLimits)
	{
		for(unsigned int i=0; i<GetTDim(); ++i)
			px[i] = dTMin + double(i+0.5)/double(GetTDim()) * (dTMax-dTMin);
	}
	else
	{
		for(unsigned int i=0; i<GetTDim(); ++i)
			px[i] = double(i);
	}

	std::string strLabX, strLabY, strLabT;
	GetLabels(strLabX, strLabY, strLabT);

	std::string strTitle;
	GetTitle(strTitle);

	const double* py = GetT(iX, iY);
	const double* pyerr = GetTErr(iX, iY);

	if(!py || !pyerr)
	{
		log_err("Invalid pixel.");
		return false;
	}

	std::ofstream ofstr(pcFile);
	if(!ofstr.is_open())
	{
		log_err("Could not open file \"", pcFile, "\".");
		return false;
	}

	ofstr << std::scientific;
	ofstr << "# title: " << strTitle << ", pixel (" << iX << ", " << iY << ")\n" ;
	ofstr << "# xlabel: " << strLabT << "\n";
	ofstr << "# ylabel: Intensity" << "\n";
	ofstr << "# xlimits: " << dTMin << " " << dTMax << "\n";
	ofstr << "# type: array_1d(" << GetTDim() << ")\n";

	for(unsigned int i=0; i<GetTDim(); ++i)
		ofstr << px[i] << " " << py[i] << " " << pyerr[i] << "\n";

	ofstr.flush();
	return true;
}

bool Data3D::SaveXY(const char* pcFile, unsigned int uiBlock, int iT) const
{
	std::ofstream ofstr(pcFile);
	if(!ofstr.is_open())
	{
		log_err("Could not open file \"", pcFile, "\".");
		return false;
	}

	ofstr << std::scientific;

	std::string strTitle;
	if(GetTitle(strTitle))
	{
		ofstr << "# title: " << strTitle;
		if(iT<0)
			ofstr << ", all timechannels";
		else
			ofstr << ", timechannel " << iT << "";

		ofstr << ", showing " << GetBlockName(uiBlock) << "\n";
	}

	ofstr << "# type: array_2d\n";
	ofstr << "# subtype: tobisown\n";

	std::string strLabX, strLabY, strLabT;
	if(GetLabels(strLabX, strLabY, strLabT))
	{
		ofstr << "# xlabel: " << strLabX << "\n";
		ofstr << "# ylabel: " << strLabY << "\n";
	}

	double dXMin, dXMax, dYMin, dYMax, dTMin, dTMax;
	if(GetLimits(dXMin, dXMax, dYMin, dYMax, dTMin, dTMax))
	{
		ofstr << "# xlimits: " << dXMin << " " << dXMax << "\n";
		ofstr << "# ylimits: " << dYMin << " " << dYMax << "\n";
		ofstr << "# zlimits: 0 0\n";
	}

	for(unsigned int iY=0; iY<GetYDim(); ++iY)
	{
		for(unsigned int iX=0; iX<GetXDim(); ++iX)
		{
			double dVal = 0.;

			if(iT>=0)		// single time channel
			{
				dVal = GetBlockVal(uiBlock, iX, iY, iT);
			}
			else			// sum all time channels
			{
				for(unsigned int iTime=0; iTime<GetTDim(); ++iTime)
					dVal += GetBlockVal(uiBlock, iX, iY, iTime);
			}

			ofstr << dVal << " ";
		}
		ofstr << "\n";
	}

	ofstr.flush();
	return true;
}

double Data3D::GetIntensity() const
{
	double dInt = 0.;
	for(unsigned int iY=0; iY<GetYDim(); ++iY)
		for(unsigned int iX=0; iX<GetXDim(); ++iX)
			for(unsigned int iT=0; iT<GetTDim(); ++iT)
				dInt += GetVal(iX, iY, iT);
	return dInt;
}

double Data3D::GetIntensityError() const
{
	double dErr = 0.;
	for(unsigned int iY=0; iY<GetYDim(); ++iY)
		for(unsigned int iX=0; iX<GetXDim(); ++iX)
			for(unsigned int iT=0; iT<GetTDim(); ++iT)
				dErr += GetErr(iX, iY, iT)*GetErr(iX, iY, iT);
	dErr = sqrt(dErr);
	return dErr;
}
