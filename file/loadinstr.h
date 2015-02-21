/*
 * Load instrument-specific data file
 * @author tweber
 * @date feb-2015
 * @copyright GPLv2 or GPLv3
 */

#ifndef __LOADINSTR_H__
#define __LOADINSTR_H__

#include <unordered_map>
#include <map>
#include <vector>
#include <array>
#include <iostream>
#include "../string/string.h"

namespace tl{

// interface for instrument-specific data files
class FileInstr
{
	public:
		typedef std::unordered_map<std::string, std::string> t_mapParams;
		typedef std::vector<std::string> t_vecColNames;
		typedef std::vector<double> t_vecVals;
		typedef std::vector<t_vecVals> t_vecDat;

	protected:
		std::array<double, 5> GetScanHKLKiKf(const char* pcH, const char* pcK,
											const char* pcL, const char* pcE,
											std::size_t i) const;

	public:
		FileInstr() = default;
		virtual ~FileInstr() = default;

		virtual bool Load(const char* pcFile) = 0;

		virtual std::array<double, 3> GetSampleLattice() const = 0;
		virtual std::array<double, 3> GetSampleAngles() const = 0;
		virtual std::array<double, 2> GetMonoAnaD() const = 0;

		virtual std::array<bool, 3> GetScatterSenses() const = 0;
		virtual std::array<double, 3> GetScatterPlane0() const = 0;
		virtual std::array<double, 3> GetScatterPlane1() const = 0;

		virtual double GetKFix() const = 0;
		virtual bool IsKiFixed() const = 0;

		virtual const std::vector<double>& GetCol(const std::string& strName) const = 0;

		virtual std::size_t GetScanCount() const = 0;
		virtual std::array<double, 5> GetScanHKLKiKf(std::size_t i) const = 0;

		virtual std::string GetTitle() const = 0;
		virtual std::string GetSampleName() const = 0;
		virtual std::string GetSpacegroup() const = 0;

		virtual const t_vecDat& GetData() const = 0;
		virtual const t_vecColNames& GetColNames() const = 0;
		virtual const t_mapParams& GetAllParams() const = 0;

		virtual std::vector<std::string> GetScannedVars() const = 0;

	public:
		static FileInstr* LoadInstr(const char* pcFile);
};


// psi & ill files
class FilePsi : public FileInstr
{
	public:
		// internal parameters in m_mapParams
		typedef std::map<std::string, double> t_mapIParams;

	protected:
		t_mapParams m_mapParams;
		t_mapIParams m_mapParameters, m_mapZeros, m_mapVariables, m_mapPosHkl, m_mapScanSteps;
		t_vecColNames m_vecColNames;
		t_vecDat m_vecData;

	protected:
		void ReadData(std::istream& istr);
		void GetInternalParams(const std::string& strAll, t_mapIParams& mapPara);

	public:
		FilePsi() = default;
		virtual ~FilePsi() = default;

		virtual bool Load(const char* pcFile) override;

		void PrintParams(std::ostream& ostr) const;
		const t_mapParams& GetParams() const { return m_mapParams; }

		const std::string& GetColName(std::size_t iCol) const { return m_vecColNames[iCol]; }
		std::size_t GetColCount() const { return m_vecColNames.size(); }

		const std::vector<double>& GetCol(std::size_t iCol) const { return m_vecData[iCol]; }
		virtual const std::vector<double>& GetCol(const std::string& strName) const override;

	public:
		virtual std::array<double, 3> GetSampleLattice() const override;
		virtual std::array<double, 3> GetSampleAngles() const override;
		virtual std::array<double, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<double, 3> GetScatterPlane0() const override;
		virtual std::array<double, 3> GetScatterPlane1() const override;

		virtual double GetKFix() const override;
		virtual bool IsKiFixed() const override;

		std::array<double, 4> GetPosHKLE() const;	// zero pos.
		std::array<double, 4> GetDeltaHKLE() const;	// scan steps

		virtual std::size_t GetScanCount() const override;
		virtual std::array<double, 5> GetScanHKLKiKf(std::size_t i) const override;

		virtual std::string GetTitle() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecColNames; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
};


// frm/nicos files
class FileFrm : public FileInstr
{
	protected:
		t_mapParams m_mapParams;
		t_vecColNames m_vecQuantities, m_vecUnits;
		t_vecDat m_vecData;

	public:
		FileFrm() = default;
		virtual ~FileFrm() = default;

	protected:
		void ReadHeader(std::istream& istr);
		void ReadData(std::istream& istr);

	public:
		virtual bool Load(const char* pcFile) override;

		virtual std::array<double, 3> GetSampleLattice() const override;
		virtual std::array<double, 3> GetSampleAngles() const override;
		virtual std::array<double, 2> GetMonoAnaD() const override;

		virtual std::array<bool, 3> GetScatterSenses() const override;
		virtual std::array<double, 3> GetScatterPlane0() const override;
		virtual std::array<double, 3> GetScatterPlane1() const override;

		virtual double GetKFix() const override;
		virtual bool IsKiFixed() const override;

		virtual std::size_t GetScanCount() const override;
		virtual std::array<double, 5> GetScanHKLKiKf(std::size_t i) const override;

		virtual const std::vector<double>& GetCol(const std::string& strName) const override;

		virtual std::string GetTitle() const override;
		virtual std::string GetSampleName() const override;
		virtual std::string GetSpacegroup() const override;

		virtual const t_vecDat& GetData() const override { return m_vecData; }
		virtual const t_vecColNames& GetColNames() const override { return m_vecQuantities; }
		virtual const t_mapParams& GetAllParams() const override { return m_mapParams; }

		virtual std::vector<std::string> GetScannedVars() const override;
};

}

#endif
