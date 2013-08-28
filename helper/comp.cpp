/*
 * Compression
 * @author tweber
 * @date 15-aug-2013
 */

#include "comp.h"
#include "string.h"

#include <fstream>
#include <list>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/copy.hpp>

namespace ios = boost::iostreams;


static Compressor comp_from_magic(const unsigned char* pcMagic, unsigned int iLen)
{
	if(iLen<2)
		return COMP_INVALID;

	if(pcMagic[0]==0x1f && pcMagic[1]==0x8b)
		return COMP_GZ;

	if(pcMagic[0]==0x78 && (pcMagic[1]==0xda || pcMagic[1]==0x9c || pcMagic[1]==0x01))
		return COMP_Z;


	if(iLen<3)
		return COMP_INVALID;

	if(pcMagic[0]=='B' && pcMagic[1]=='Z' && pcMagic[2]=='h')
		return COMP_BZ2;


	return COMP_INVALID;
}


//--------------------------------------------------------------------------------


bool decomp_stream_to_stream(std::istream& istr, std::ostream& ostr, Compressor comp)
{
	if(comp == COMP_AUTO)
	{
		unsigned char pcMagic[3] = {0,0,0};
		istr.read((char*)pcMagic, 3);
		istr.seekg(0,std::ios::beg);

		comp = comp_from_magic(pcMagic, 3);
	}


	ios::filtering_streambuf<ios::input> input;

	if(comp == COMP_GZ)
		input.push(ios::gzip_decompressor());
	else if(comp == COMP_BZ2)
		input.push(ios::bzip2_decompressor());
	else if(comp == COMP_Z)
		input.push(ios::zlib_decompressor());
	else
	{
		std::cerr << "Error: Unknown decompression selected." << std::endl;
		return false;
	}

	input.push(istr);
	ios::copy(input, ostr);

	return true;
}

bool comp_stream_to_stream(std::istream& istr, std::ostream& ostr, Compressor comp)
{
	if(comp == COMP_AUTO)
		comp = COMP_GZ;


	ios::filtering_streambuf<ios::input> input;

	if(comp == COMP_GZ)
		input.push(ios::gzip_compressor());
	else if(comp == COMP_BZ2)
		input.push(ios::bzip2_compressor());
	else if(comp == COMP_Z)
		input.push(ios::zlib_compressor());
	else
	{
		std::cerr << "Error: Unknown compression selected." << std::endl;
		return false;
	}

	input.push(istr);
	ios::copy(input, ostr);

	return true;
}


//--------------------------------------------------------------------------------


static inline bool __comp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp, bool bDecomp=0)
{
	std::ifstream ifstr(pcFileIn, std::ios_base::binary);
	if(!ifstr.is_open())
	{
		std::cerr << "Error: Cannot open \"" << pcFileIn << "\"." << std::endl;
		return false;
	}


	std::ofstream ofstr(pcFileOut, std::ios_base::binary);
	if(!ofstr.is_open())
	{
		std::cerr << "Error: Cannot open \"" << pcFileOut << "\"." << std::endl;
		return false;
	}


	if(comp == COMP_AUTO)
	{
		std::string strExt;
		if(bDecomp)
			strExt = get_fileext(pcFileIn);
		else
			strExt = get_fileext(pcFileOut);

		if(strExt == "gz")
			comp = COMP_GZ;
		else if(strExt == "bz2")
			comp = COMP_BZ2;
	}

	if(bDecomp)
		return decomp_stream_to_stream(ifstr, ofstr, comp);
	else
		return comp_stream_to_stream(ifstr, ofstr, comp);
}

bool comp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp)
{
	return __comp_file_to_file(pcFileIn, pcFileOut, comp, 0);
}

bool decomp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp)
{
	return __comp_file_to_file(pcFileIn, pcFileOut, comp, 1);
}


//--------------------------------------------------------------------------------


static inline bool __comp_mem_to_mem(void* pvIn, unsigned int iLenIn, void*& pvOut, unsigned int& iLenOut, Compressor comp, bool bDecomp=0)
{
	char *pcIn = (char*)pvIn;
	ios::stream<ios::basic_array_source<char> > istr(pcIn, iLenIn);

	std::list<char> lstOut;
	ios::filtering_ostream arrOut(ios::back_inserter(lstOut));


	bool bOk=0;
	if(bDecomp)
		bOk = decomp_stream_to_stream(istr, arrOut, comp);
	else
		bOk =comp_stream_to_stream(istr, arrOut, comp);


	iLenOut = lstOut.size();
	pvOut = new char[iLenOut];
	char *pcOut = (char*)pvOut;

	unsigned int iIdx = 0;
	for(char c : lstOut)
		pcOut[iIdx++] = c;

	return bOk;
}

/*
 * example:
 *	char pc[] = "123456\nABCDEF\n\n";
 *	void *pvTst;
 *	unsigned int iLenOut=0;
 *	::comp_mem_to_mem(pc, strlen(pc), pvTst,iLenOut,COMP_BZ2);
 *	std::ofstream ofstr("tst.txt.bz2");
 *	ofstr.write((char*)pvTst, iLenOut);
 *	ofstr.close();
 *	delete[] pvTst;
 */
bool comp_mem_to_mem(void* pvIn, unsigned int iLenIn, void*& pvOut, unsigned int& iLenOut, Compressor comp)
{
	return __comp_mem_to_mem(pvIn, iLenIn, pvOut, iLenOut, comp, 0);
}

bool decomp_mem_to_mem(void* pvIn, unsigned int iLenIn, void*& pvOut, unsigned int& iLenOut, Compressor comp)
{
	return __comp_mem_to_mem(pvIn, iLenIn, pvOut, iLenOut, comp, 1);
}
