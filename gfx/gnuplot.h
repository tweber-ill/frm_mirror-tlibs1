/*
 * invoke gnuplot
 * @autor tweber
 * @date 24-dec-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __M_PLOTTER_GPL__
#define __M_PLOTTER_GPL__

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <iostream>
#include <string>
#include <vector>

namespace tl {
enum LineStyle
{
	STYLE_POINTS,
	STYLE_LINES_SOLID,
	STYLE_LINES_DASHED
};

struct PlotObj
{
	std::vector<double> vecX, vecY;
	std::vector<double> vecErrX, vecErrY;

	std::string strLegend;

	LineStyle linestyle;

	bool bHasSize;
	double dSize;

	bool bHasColor;
	unsigned int iColor;

	PlotObj() : linestyle(STYLE_POINTS), bHasSize(0), bHasColor(0)
	{}
};

class GnuPlot
{
protected:
	FILE *m_pipe = 0;
	boost::iostreams::file_descriptor_sink *m_pfds = 0;
	boost::iostreams::stream_buffer<boost::iostreams::file_descriptor_sink> *m_psbuf = 0;
	std::ostream *m_postr = 0;

	std::vector<PlotObj> m_vecObjs;
	// has to be 0 to show plot
	int m_iStartCounter = 0;

	std::string BuildCmd();
	std::string BuildTable(const std::vector<double>& vecX, const std::vector<double>& vecY,
			const std::vector<double>& vecYErr, const std::vector<double>& vecXErr);

	bool m_bTermLocked = false;
	bool m_bHasLegend = false;

	std::string m_strLegendOpts;
	std::string m_strLegendPlacement = "default";

	void RefreshVars();

public:
	GnuPlot() = default;
	virtual ~GnuPlot();

	void Init();
	void DeInit();

	bool IsReady() const;
	std::ostream& GetStream();

	void SetTerminal(int iWnd=0, const char* pcBackend="x11");
	void SetFileTerminal(const char* pcFile);

	void StartPlot();
	void AddLine(const PlotObj& obj);
	void FinishPlot();

	void SimplePlot(const std::vector<double>& vecX, const std::vector<double>& vecY,
			const std::vector<double>& vecYErr, const std::vector<double>& vecXErr,
			LineStyle style=STYLE_POINTS);
	void SimplePlot2d(const std::vector<std::vector<double> >& vec,
			double dMinX=0., double dMaxX=-1., double dMinY=0., double dMaxY=-1.);

	void SetXLabel(const char* pcLab);
	void SetYLabel(const char* pcLab);
	void SetTitle(const char* pcTitle);
	void SetGrid(bool bOn);
	void AddArrow(double dX0, double dY0, double dX1, double dY1, bool bHead=1);

	void SetXRange(double dMin, double dMax);
	void SetYRange(double dMin, double dMax);

	void SetColorBarRange(double dMin, double dMax, bool bCyclic=0);

	void LockTerminal() { m_bTermLocked = 1; }
	void UnlockTerminal() { m_bTermLocked = 0; }

	void SetLegendOpts(const std::string& strOpts) { m_strLegendOpts = strOpts; }
	void SetLegendPlace(const std::string& strPlace) { m_strLegendPlacement = strPlace; }
};
}

#endif
