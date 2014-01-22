/*
 * invoke gnuplot
 * @autor tweber
 * @date 24-dec-2013
 */

#ifndef __M_PLOTTER_GPL__
#define __M_PLOTTER_GPL__

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <iostream>
#include <string>
#include <vector>

struct PlotObj
{
	std::vector<double> vecX, vecY;
	std::vector<double> vecErrX, vecErrY;

	std::string strLegend;

	bool bConnectLines;

	PlotObj() : bConnectLines(0)
	{}
};

class GnuPlot
{
protected:
	FILE *m_pipe;
	boost::iostreams::file_descriptor_sink *m_pfds;
	boost::iostreams::stream_buffer<boost::iostreams::file_descriptor_sink> *m_psbuf;
	std::ostream *m_postr;

	std::vector<PlotObj> m_vecObjs;
	// has to be 0 to show plot
	int m_iStartCounter;

	std::string BuildCmd();
	std::string BuildTable(const std::vector<double>& vecX, const std::vector<double>& vecY,
			const std::vector<double>& vecYErr, const std::vector<double>& vecXErr);

	bool m_bTermLocked;
	bool m_bHasLegend;

	void RefreshVars();

public:
	GnuPlot();
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
			const std::vector<double>& vecYErr, const std::vector<double>& vecXErr);
	void SimplePlot2d(const std::vector<std::vector<double> >& vec,
			double dMinX=0., double dMaxX=-1., double dMinY=0., double dMaxY=-1.);

	void SetXLabel(const char* pcLab);
	void SetYLabel(const char* pcLab);
	void SetTitle(const char* pcTitle);

	void SetXRange(double dMin, double dMax);
	void SetYRange(double dMin, double dMax);

	void SetColorBarRange(double dMin, double dMax, bool bCyclic=0);

	void LockTerminal() { m_bTermLocked = 1; }
	void UnlockTerminal() { m_bTermLocked = 0; }
};


#endif
