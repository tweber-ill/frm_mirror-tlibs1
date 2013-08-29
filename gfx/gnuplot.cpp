/*
 * invoke gnuplot
 * @autor tweber
 * @date 24-dec-2013
 */

#include "flags.h"
#include "gnuplot.h"
#include "misc.h"
#include "string.h"
#include "log.h"

#include <cstdio>
#include <sstream>
#include <iomanip>

namespace ios = boost::iostreams;

GnuPlot::~GnuPlot() { DeInit(); }

void GnuPlot::DeInit()
{
	if(m_postr) delete m_postr;
	if(m_psbuf) delete m_psbuf;
	if(m_pfds) delete m_pfds;
	if(m_pipe) ::pclose(m_pipe);

	m_postr = 0;
	m_psbuf = 0;
	m_pfds = 0;
	m_pipe = 0;

	m_iStartCounter = 0;
}

void GnuPlot::Init()
{
	if(IsReady()) return;
	DeInit();

	m_pipe = (FILE*)::/*my_*/popen("gnuplot -p 2>/dev/null 1>/dev/null", "w");
	//m_pipe = popen("gnuplot -persist 2>gnuplot.err 1>gnuplot.out", "w");
	if(!m_pipe)
	{
		log_err("Could not load gnuplot.");
		return;
	}

	m_pfds = new ios::file_descriptor_sink(fileno(m_pipe), ios::close_handle /*ios::never_close_handle*/);
	m_psbuf = new ios::stream_buffer<ios::file_descriptor_sink>(*m_pfds);
	m_postr = new std::ostream(m_psbuf);

	(*m_postr) << "set grid\n";
	(*m_postr) << "set nokey\n";
	//(*m_postr) << "set noborder\n";
	(*m_postr) << "set size 1,1\n";
	(*m_postr) << "set palette rgbformulae 33,13,10\n";
	//(*m_postr) << "set termopt dashed\n";
	//(*m_postr) << "set palette defined (0 \"blue\", 0.3333 \"cyan\", 0.6666 \"yellow\", 1 \"red\")\n";
}

void GnuPlot::SetTerminal(int iWnd, const char* pcBackend)
{
	if(m_bTermLocked) return;

	(*m_postr) << "set output\n";
	(*m_postr) << "set obj 1 rectangle behind fillcolor rgbcolor \"white\" from screen 0,0 to screen 1,1\n";

	(*m_postr) << "set term " << pcBackend <<  " " << iWnd << " "
			<< "size 640,480 "
			<< "enhanced "
			<< "font 'Helvetica,11' "
//			<< "title \"" << "Plot " << (iWnd+1) << "\" " 
			<< "persist "
			<< "dashed "
			<<  "\n";
}

void GnuPlot::SetFileTerminal(const char* pcFile)
{
	if(m_bTermLocked) return;

	std::string strFile = pcFile;
	std::string strExt = get_fileext(strFile);

	//std::cout << "File: " << strFile << "\nExtension: " << strExt << std::endl;

	if(str_is_equal(strExt, std::string("pdf"), 0))
	{
		//(*m_postr) << "set term pdfcairo enhanced color font \"Helvetica\"\n";
		(*m_postr) << "set term pdf enhanced color\n";
	}
	else if(str_is_equal(strExt, std::string("ps"), 0))
	{
		//(*m_postr) << "set term postscript eps enhanced color \"Helvetica\" 24\n";
		(*m_postr) << "set term postscript eps enhanced color\n";
	}
	else
	{
		log_err("Unknown file extension \"", 
			strExt, "\" for output terminal.");
		return;
	}

	(*m_postr) << "set output \"" << strFile << "\"\n";
}

void GnuPlot::RefreshVars()
{
	if(m_bHasLegend)
		(*m_postr) << "set key on default box 0 " << m_strLegendOpts << "\n";
	else
		(*m_postr) << "set nokey\n";
}

void GnuPlot::SimplePlot(const std::vector<double>& vecX, const std::vector<double>& vecY,
		const std::vector<double>& vecYErr, const std::vector<double>& vecXErr)
{
	if(!IsReady()) return;
	(*m_postr) << "plot '-'\n";

	std::ostream* postr = m_postr;
	std::string strTable = BuildTable(vecX, vecY, vecYErr, vecXErr);
	(*postr) << strTable;
	postr->flush();
}


void GnuPlot::SimplePlot2d(const std::vector<std::vector<double> >& vec,
		double dMinX, double dMaxX, double dMinY, double dMaxY)
{
	if(!IsReady()) return;

	std::vector<unsigned int> vecSizes;
	vecSizes.reserve(vec.size());

	for(unsigned int iY=0; iY<vec.size(); ++iY)
		vecSizes.push_back(vec[iY].size());

	std::vector<unsigned int>::iterator iterMin =
		std::min_element(vecSizes.begin(), vecSizes.end());
	unsigned int iXCntMin = *iterMin;


	unsigned int iYDim = vec.size();
	unsigned int iXDim = iXCntMin;

	// invalid values select image dimensions
	if(dMinX > dMaxX)
	{
		dMinX = 0.;
		dMaxX = iXDim-1;
	}
	if(dMinY > dMaxY)
	{
		dMinY = 0.;
		dMaxY = iYDim-1;
	}

	// ----------------------------------------
	// ranges
	(*m_postr) << "set tics out scale 0.8\n";

	double dRangeMinX = tic_trafo(iXDim, dMinX, dMaxX, 0, -0.5);
	double dRangeMaxX = tic_trafo(iXDim, dMinX, dMaxX, 0, double(iXDim)-0.5);
	double dRangeMinY = tic_trafo(iYDim, dMinY, dMaxY, 0, -0.5);
	double dRangeMaxY = tic_trafo(iYDim, dMinY, dMaxY, 0, double(iYDim)-0.5);

	(*m_postr) << "set xrange [" << dRangeMinX << ":" << dRangeMaxX << "]\n";
	(*m_postr) << "set yrange [" << dRangeMinY << ":" << dRangeMaxY << "]\n";
	// ----------------------------------------

	// ----------------------------------------
	// tics
	std::ostringstream ostrTicsX, ostrTicsY;
	ostrTicsX << "(" << dMinX << " + " << "($1)/" << iXDim
			<< " * (" << dMaxX << "-" << dMinX << "))";
	ostrTicsY << "(" << dMinY << " + " << "($2)/" << iYDim
			<< " * (" << dMaxY << "-" << dMinY << "))";

	std::string strTics = "using " + ostrTicsX.str() + ":" + ostrTicsY.str() + ":3";
	// ----------------------------------------

	(*m_postr) << "plot '-' " << strTics << " matrix with image\n";


	for(unsigned int iY=0; iY<vec.size(); ++iY)
	{
		for(unsigned int iX=0; iX<iXCntMin; ++iX)
			(*m_postr) << vec[iY][iX] << " ";
		(*m_postr) << "\n";
	}

	(*m_postr) << "e\ne\n";


	m_postr->flush();
}

void GnuPlot::AddLine(const PlotObj& obj)
{
	m_vecObjs.push_back(obj);
}

void GnuPlot::StartPlot()
{
	if(!IsReady()) return;

	if(m_iStartCounter == 0)
		m_vecObjs.clear();

	++m_iStartCounter;
}

void GnuPlot::FinishPlot()
{
	if(!IsReady()) return;

	if(--m_iStartCounter == 0)
	{
		std::string strCmd = BuildCmd();
		RefreshVars();

		//std::cout << "Plot cmd: " << strCmd << std::endl;
		(*m_postr) << strCmd;
		//(*m_postr) << "replot\n";
		m_postr->flush();
		m_vecObjs.clear();
	}
}

std::string GnuPlot::BuildCmd()
{
	m_bHasLegend = 0;

	std::ostringstream ostr;
	ostr << "plot ";

	for(const PlotObj& obj : m_vecObjs)
	{
		const bool bConnectLines = (obj.linestyle != STYLE_POINTS);
		const bool bHasXErr = (obj.vecErrX.size() != 0);
		const bool bHasYErr = (obj.vecErrY.size() != 0);
		const bool bHasSize = obj.bHasSize;
		double dSize = obj.dSize;

		if(!bHasSize)
		{
			if(bConnectLines)
				dSize = 1.;
			else
				dSize = 1.25;
		}

		std::ostringstream ostrTmp;
		std::string strPointStyle;

		if(bHasXErr && bHasYErr)
			ostrTmp << "with xyerrorbars";
		else if(bHasXErr && !bHasYErr)
			ostrTmp << "with xerrorbars";
		else if(!bHasXErr && bHasYErr)
			ostrTmp << "with yerrorbars";
		else if(!bHasXErr && !bHasYErr)
			ostrTmp << "with points";

		ostrTmp << " pt 7 ps " << dSize;
		//std::cout << "*** point size: " << dSize << std::endl;
		strPointStyle = ostrTmp.str();


		ostr << "'-' ";
		switch(obj.linestyle)
		{
			case STYLE_LINES_SOLID:
				ostr << "with lines lt 1 lw " << dSize;
				break;
			case STYLE_LINES_DASHED:
				ostr << "with lines lt 2 lw " << dSize;
				break;
			default:
				log_warn("Unknown line style.");
			case STYLE_POINTS:
				ostr << strPointStyle;
				break;
		}

		if(obj.bHasColor)
		{
			char chFill = ostr.fill();
			ostr << std::setfill('0');
			ostr << "lc rgb \"#" << std::hex
				<< std::setw(2) << ((obj.iColor & 0xff0000) >> 16)
				<< std::setw(2) << ((obj.iColor & 0x00ff00) >> 8)
				<< std::setw(2) << ((obj.iColor & 0x0000ff))
				<< "\" " << std::dec;
			ostr << std::setfill(chFill);
		}

		if(obj.strLegend != "")
		{
			m_bHasLegend = 1;
			ostr << " title \"" << obj.strLegend << "\" ";
		}

		if(&obj != &(*m_vecObjs.rbegin()))
			ostr << ", ";
	}
	ostr << "\n";

	for(const PlotObj& obj : m_vecObjs)
	{
		std::string strTab = BuildTable(obj.vecX, obj.vecY, obj.vecErrY, obj.vecErrX);
		ostr << strTab;
	}

	return ostr.str();
}

std::string GnuPlot::BuildTable(const std::vector<double>& vecX, const std::vector<double>& vecY,
								const std::vector<double>& vecYErr, const std::vector<double>& vecXErr)
{
	std::ostringstream ostr;

	const unsigned int iSize = std::min(vecX.size(), vecY.size());
	const bool bHasXErr = (vecXErr.size() != 0);
	const bool bHasYErr = (vecYErr.size() != 0);

	for(unsigned int iDat=0; iDat<iSize; ++iDat)
	{
		ostr << vecX[iDat] << " " << vecY[iDat];

		if(bHasXErr) ostr << " " << vecXErr[iDat];
		if(bHasYErr) ostr << " " << vecYErr[iDat];

		ostr << "\n";
	}
	ostr << "e\n";

	return ostr.str();
}

void GnuPlot::SetXLabel(const char* pcLab)
{
	if(!IsReady()) return;
	(*m_postr) << "set xlabel \"" << pcLab << "\"\n";
	//m_postr->flush();
}
void GnuPlot::SetYLabel(const char* pcLab)
{
	if(!IsReady()) return;
	(*m_postr) << "set ylabel \"" << pcLab << "\"\n";
	//m_postr->flush();
}
void GnuPlot::SetTitle(const char* pcTitle)
{
	if(!IsReady()) return;
	(*m_postr) << "set title \"" << pcTitle << "\"\n";
	//m_postr->flush();
}

void GnuPlot::SetXRange(double dMin, double dMax)
{
	if(!IsReady()) return;
	//std::cout << "xmin: "  << dMin << ", xmax: " << dMax << std::endl;
	(*m_postr) << "set xrange [" << dMin << ":" << dMax << "]\n";
	m_postr->flush();
}

void GnuPlot::SetYRange(double dMin, double dMax)
{
	if(!IsReady()) return;
	(*m_postr) << "set yrange [" << dMin << ":" << dMax << "]\n";
	m_postr->flush();
}

void GnuPlot::SetColorBarRange(double dMin, double dMax, bool bCyclic)
{
	if(!IsReady()) return;

	(*m_postr) << "set cbrange [" << dMin << ":" << dMax << "]\n";

	if(bCyclic)
		(*m_postr) << "set palette defined (0 \"blue\", 0.25 \"cyan\", 0.5 \"yellow\", 0.75 \"red\", 1 \"blue\")\n";
	else
		(*m_postr) << "set palette defined (0 \"blue\", 0.3333 \"cyan\", 0.6666 \"yellow\", 1 \"red\")\n";

	m_postr->flush();
}

bool GnuPlot::IsReady() const { return m_postr!=0; }
std::ostream& GnuPlot::GetStream() { return *m_postr; }