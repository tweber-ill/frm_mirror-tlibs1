/*
 * Load .dat files
 *
 * @author: Tobias Weber
 * @date: April 2012
 */

/*
	available options
	-----------------

	Mcstas compatible = [x]


	[x] title

	[x] type
	subtype

	[x] xlabel
	[x] ylabel
	[?] zlabel

	[x] xlimits
	[x] ylimits
	[?] zlimits
	[x] xylimits

	xlog
	ylog
	xylog

	which-col-is-x
	which-col-is-y
	which-col-is-yerr

	col-0-label, col-1-label, etc.
*/

#ifndef __LOADTXT__
#define __LOADTXT__

#include <vector>
#include <map>
#include <ostream>
#include <sstream>
#include <string>

enum TxtType
{
	MCSTAS_DATA = 0,
	NICOS_DATA,

	UNKNOWN_DATA
};

class LoadTxt
{
	public:
		typedef std::map<std::string, std::vector<std::string> > t_mapComm;

	protected:
		typedef std::vector<double*> t_vecColumns;

		t_vecColumns m_vecColumns;
		unsigned int m_uiColLen;

		t_mapComm m_mapComm;
		std::vector<std::string> m_vecAuxStrings, m_vecComments;

		// trim and collect McStas comments
		void StrTrim(std::string& str);

		bool m_bLoadOnlyHeader;
		bool m_bVerbose;

		std::string m_strFileName;

	public:
		bool Load(const char* pcFile);
		void Unload();

		bool Save(const char* pcFile, bool bSaveComments=1) const;

		LoadTxt(const char* pcFile=0, bool bOnlyHeader=false, bool bVerbose=false);
		virtual ~LoadTxt();

		unsigned int GetColLen() const { return m_uiColLen; }
		unsigned int GetColCnt() const { return m_vecColumns.size(); }

		const double* GetColumn(unsigned int uiCol) const
		{
			if(uiCol < m_vecColumns.size())
				return m_vecColumns[uiCol];
			return 0;
		}
		double* GetColumn(unsigned int uiCol)
		{
			return const_cast<double*>(const_cast<const LoadTxt*>(this)->GetColumn(uiCol));
		}

		const t_mapComm& GetCommMap() const { return m_mapComm; }
		t_mapComm& GetCommMap() { return m_mapComm; }

		// only first value in map's "second" vector
		std::map<std::string, std::string> GetCommMapSingle() const;

		void GetMinMax(double& dMin, double& dMax) const;

		void for_each(void (*)(double*));

		bool GetMapString(const std::string& strKey, std::string& strVal) const;
		template<typename T> bool GetMapVal(const std::string& strKey, T& tVal) const
		{
			std::string strVal;
			bool bOk = GetMapString(strKey, strVal);

			if(!bOk)
				return false;

			std::istringstream istr(strVal);
			istr >> tVal;

			return true;
		}

		void SetMapString(const std::string& strKey, const std::string& strVal);
		template<typename T> void SetMapVal(const std::string& strKey, const T& t)
		{
			std::ostringstream ostr;
			ostr << t;
			SetMapString(strKey, ostr.str());
		}

		const std::string& GetFileName() const { return m_strFileName; }
		const TxtType GetFileType() const;

		const std::vector<std::string>& GetAuxStrings() const { return m_vecAuxStrings; }

		friend std::ostream& operator<<(std::ostream& ostr, LoadTxt& txt);
};


void convert_mcstas_string(std::string& str);

// base class for Mcstas data
class McData
{
	protected:
		const LoadTxt& m_data;
		bool m_bOk;

	public:
		McData(const LoadTxt& data) : m_data(data), m_bOk(0) {}
		virtual ~McData() {}

		bool IsOk() const { return m_bOk; }

		bool GetTitle(std::string& strTitle) const;
		bool IsOwnDataFile() const;
		bool GetFit(std::string& strFit) const;

		bool GetString(const std::string& strKey, std::string& strVal) const
		{ return m_data.GetMapString(strKey, strVal); }
		bool GetDouble(const std::string& strKey, double& dVal) const
		{ return m_data.GetMapVal<double>(strKey, dVal); }

		const LoadTxt& GetRawData() const { return m_data; }

		bool GetLogScale(bool &bXLog, bool &bYLog) const;

		const std::string& GetFileName() const { return m_data.GetFileName(); }
};

// resolution data
// layout:
//			1st column: ki_x
//			2nd column: ki_y
//			3rd column: ki_z
//			4th column: kf_x
//			5th column: kf_y
//			6th column: kf_z
//			7th column: x
//			8th column: y
//			9th column: z
//			10th column: p_i
//			11th column: p_f
class DataRes : public McData
{
	protected:
		unsigned int m_iDim;
		unsigned int m_iCols;

	public:
		DataRes(const LoadTxt& data);
		virtual ~DataRes() {}

		unsigned int GetDim(void) const { return m_iDim; }

		const double* GetColumn(unsigned int iIdx) const;

		const double* Getki_x() const { return GetColumn(0); }
		const double* Getki_y() const { return GetColumn(1); }
		const double* Getki_z() const { return GetColumn(2); }

		const double* Getkf_x() const { return GetColumn(3); }
		const double* Getkf_y() const { return GetColumn(4); }
		const double* Getkf_z() const { return GetColumn(5); }

		const double* Getx() const { return GetColumn(6); }
		const double* Gety() const { return GetColumn(7); }
		const double* Getz() const { return GetColumn(8); }

		const double* Getp_i() const { return GetColumn(9); }
		const double* Getp_f() const { return GetColumn(10); }
};

// class for interpreting Mcstas 1d arrays
// layout:
//			first column: x
//			second column: y
//			third column: y err
//			forth column: counts
class Data1D : public McData
{
	protected:
		unsigned int m_iDim;

	public:
		Data1D(const LoadTxt& data);
		virtual ~Data1D() {}

		void ForEach(unsigned int iColIdx, double (*)(double));

		unsigned int GetDim(void) const { return m_iDim; }
		unsigned int GetColCnt() const { return m_data.GetColCnt(); }

		bool GetColumnIndices(int& iX, int& iY, int&iYErr) const;

		const double* GetColumn(unsigned int iIdx) const;

		double GetMin(unsigned int iColIdx, unsigned int *pMinIdx=0) const;
		double GetMax(unsigned int iColIdx, unsigned int *pMaxIdx=0) const;

		void CalcLimits(double& dMin, double& dMax, int iCol=-1) const;

		bool GetXLimits(double& dXMin, double& dXMax) const;
		bool GetYLimits(double& dYMin, double& dYMax) const;

		bool GetLabels(std::string& strLabelX, std::string& strLabelY) const;
		bool GetColLabel(unsigned int iCol, std::string& strLabel) const;

		void NormY(double dNorm, double dNormErr);

		double GetIntensity(unsigned int iICol, unsigned int iIErrCol, unsigned int iCountCol) const;
		double GetIntensityError(unsigned int iICol, unsigned int iIErrCol, unsigned int iCountCol) const;
};

// class for interpreting Mcstas 2d arrays
// layout:
//			first block: data
//			second block: errors
//			third block: counts
class Data2D;
class Data2D : public McData
{
	protected:
		unsigned int m_iNumBlocks;
		unsigned int m_iXDim, m_iYDim;

		bool m_bRecalcZLimits;

	public:
		Data2D(const LoadTxt& data);
		virtual ~Data2D() {}

		void ForEach(unsigned int iBlockIdx, double (*)(double));
		void Add(unsigned int iBlockIdx, double d);

		double GetCenterVal(unsigned int iBlock) const;
		double GetBlockVal(unsigned int iBlock, unsigned int iX, unsigned int iY) const;
		double GetVal(unsigned int iX, unsigned int iY) const;
		double GetErr(unsigned int iX, unsigned int iY) const;
		double GetCount(unsigned int iX, unsigned int iY) const;

		void GetBlockMinMax(unsigned int iBlock, double& dMin, double& dMax) const;
		void GetValMinMax(double& dMin, double& dMax) const;
		void GetErrMinMax(double& dMin, double& dMax) const;

		unsigned int GetXDim(void) const { return m_iXDim; }
		unsigned int GetYDim(void) const { return m_iYDim; }

		bool GetLimits(double& dXMin, double& dXMax,
						double& dYMin, double& dYMax,
						double& dZMin, double& dZMax) const;
		void CalcLimits(double& dZMin, double& dZMax, unsigned int iBlock=0) const;

		bool GetLabels(std::string& strLabelX, std::string& strLabelY, std::string& strLabelZ) const;

		bool ExtractColOrRow(int iCol, int iRow, const char* pcOutFile, const Data2D* pErrors=0);

		void Mod0to2Pi(unsigned int uiBlock);
};

// class for interpreting Mcstas 3d arrays
// layout:
//			first block: data
//			second block: errors
//			third block: counts
class Data3D : public McData
{
	protected:
		unsigned int m_iNumBlocks;
		unsigned int m_iXDim, m_iYDim, m_iTDim;

	public:
		Data3D(const LoadTxt& data);
		virtual ~Data3D() {}

		void ForEach(unsigned int iBlockIdx, double (*)(double));

		std::string GetBlockName(unsigned int uiBlock) const;

		double GetBlockVal(unsigned int iBlock, unsigned int iX, unsigned int iY, unsigned int iT) const;
		double GetVal(unsigned int iX, unsigned int iY, unsigned int iT) const;
		double GetErr(unsigned int iX, unsigned int iY, unsigned int iT) const;
		double GetCount(unsigned int iX, unsigned int iY, unsigned int iT) const;

		const double* GetBlockT(unsigned int iBlock, unsigned int iX, unsigned int iY) const;
		const double* GetT(unsigned int iX, unsigned int iY) const;
		const double* GetTErr(unsigned int iX, unsigned int iY) const;
		const double* GetTCount(unsigned int iX, unsigned int iY) const;

		unsigned int GetXDim(void) const { return m_iXDim; }
		unsigned int GetYDim(void) const { return m_iYDim; }
		unsigned int GetTDim(void) const { return m_iTDim; }

		bool GetLimits(double& dXMin, double& dXMax,
					   double& dYMin, double& dYMax,
					   double& dTMin, double& dTMax) const;
		bool GetLabels(std::string& strLabelX, std::string& strLabelY, std::string& strLabelT) const;

		bool SaveTChan(const char* pcFile, unsigned int iX, unsigned int iY) const;
		bool SaveXY(const char* pcFile, unsigned int uiBlock, int iT=-1) const;

		double GetIntensity() const;
		double GetIntensityError() const;
};

#endif
