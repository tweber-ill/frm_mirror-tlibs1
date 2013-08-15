/*
 * Compression
 * @author tweber
 * @date 15-aug-2013
 */

#include "comp.h"
#include "string.h"

#include <fstream>
#include <list>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/copy.hpp>

namespace ios = boost::iostreams;

bool decomp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp)
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


	ios::filtering_streambuf<ios::input> input;

	if(comp == COMP_AUTO)
	{
		std::string strExt = get_fileext(pcFileIn);
		if(strExt == "gz")
			comp = COMP_GZ;
		else if(strExt == "bz2")
			comp = COMP_BZ2;
	}

	if(comp == COMP_GZ)
		input.push(ios::gzip_decompressor());
	else if(comp == COMP_BZ2)
		input.push(ios::bzip2_decompressor());
	else
	{
		std::cerr << "Error: Unknown decompression for \"" << pcFileIn << "\"." << std::endl;
		return false;
	}

	input.push(ifstr);

	ios::copy(input, ofstr);

	return true;
}


bool decomp_mem_to_mem(void* pvIn, unsigned int iLenIn, void*& pvOut, unsigned int& iLenOut, Compressor comp)
{
	char *pcIn = (char*)pvIn;
	ios::basic_array_source<char> arrIn(pcIn, pcIn+iLenIn);

	std::list<char> lstOut;
	ios::filtering_ostream arrOut(ios::back_inserter(lstOut));


	ios::filtering_streambuf<ios::input> input;

	if(comp == COMP_GZ)
		input.push(ios::gzip_decompressor());
	else if(comp == COMP_BZ2)
		input.push(ios::bzip2_decompressor());
	else
	{
		std::cerr << "Error: Unknown decompression selected." << std::endl;
		return false;
	}

	input.push(arrIn);

	ios::copy(input, arrOut);


	iLenOut = lstOut.size();
	pvOut = new char[iLenOut];
	char *pcOut = (char*)pvOut;

	unsigned int iIdx = 0;
	for(char c : lstOut)
		pcOut[iIdx++] = c;

	return true;
}
