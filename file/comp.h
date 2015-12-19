/*
 * Compression
 * @author tweber
 * @date 15-aug-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_COMP_H__
#define __TLIBS_COMP_H__

#include <iostream>
#include <type_traits>
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

#include "../string/string.h"
#include "../log/log.h"


namespace tl {

namespace ios = boost::iostreams;


enum class Compressor
{
	GZ,
	BZ2,
	Z,
	XZ,

	AUTO,
	INVALID
};



// -----------------------------------------------------------------------------



template<class t_char=unsigned char>
Compressor comp_from_magic(const t_char* pcMagic, std::size_t iLen)
{
	if(iLen<2)
		return Compressor::INVALID;
	if(pcMagic[0]==0x1f && pcMagic[1]==0x8b)
		return Compressor::GZ;
	if(pcMagic[0]==0x78 && (pcMagic[1]==0xda ||
		pcMagic[1]==0x9c || pcMagic[1]==0x01))
		return Compressor::Z;


	if(iLen<3)
		return Compressor::INVALID;
	if(pcMagic[0]=='B' && pcMagic[1]=='Z' && pcMagic[2]=='h')
		return Compressor::BZ2;


	if(iLen<6)
		return Compressor::INVALID;
	if(pcMagic[0]==0xfd && pcMagic[1]=='7' && pcMagic[2]=='z' &&
		pcMagic[3]=='X' && pcMagic[4]=='Z' && pcMagic[5]==0x00)
		return Compressor::XZ;


	return Compressor::INVALID;
}



// -----------------------------------------------------------------------------



template<class t_char=char>
bool decomp_stream_to_stream(std::basic_istream<t_char>& istr,
							std::basic_ostream<t_char>& ostr,
							Compressor comp=Compressor::AUTO)
{
	typedef typename std::make_unsigned<t_char>::type t_uchar;

	if(comp == Compressor::AUTO)
	{
		t_uchar pcMagic[] = {0,0,0,0,0,0};
		istr.read((t_char*)pcMagic, 6);
		istr.seekg(0, std::ios::beg);

		comp = comp_from_magic(pcMagic, 6);
	}


	ios::filtering_streambuf<ios::input> input;

	if(comp == Compressor::GZ)
		input.push(ios::gzip_decompressor());
	else if(comp == Compressor::BZ2)
		input.push(ios::bzip2_decompressor());
	else if(comp == Compressor::Z)
		input.push(ios::zlib_decompressor());
	else if(comp == Compressor::XZ)
	{
		log_err("XZ decompression not yet supported.");
		return false;
	}
	else
	{
		log_err("Unknown decompression selected.");
		return false;
	}

	input.push(istr);
	ios::copy(input, ostr);

	return true;
}

template<class t_char=char>
bool comp_stream_to_stream(std::basic_istream<t_char>& istr,
						std::basic_ostream<t_char>& ostr,
						Compressor comp=Compressor::GZ)
{
	if(comp == Compressor::AUTO)
		comp = Compressor::GZ;


	ios::filtering_streambuf<ios::input> input;

	if(comp == Compressor::GZ)
		input.push(ios::gzip_compressor(/*ios::gzip_params(9)*/));
	else if(comp == Compressor::BZ2)
		input.push(ios::bzip2_compressor());
	else if(comp == Compressor::Z)
		input.push(ios::zlib_compressor(/*ios::zlib_params(9)*/));
	else if(comp == Compressor::XZ)
	{
		log_err("XZ compression not yet supported.");
		return false;
	}
	else
	{
		log_err("Unknown compression selected.");
		return false;
	}

	input.push(istr);
	ios::copy(input, ostr);

	return true;
}



// -----------------------------------------------------------------------------



template<class t_char=char>
inline bool __comp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp, bool bDecomp=0)
{
	std::basic_ifstream<t_char> ifstr(pcFileIn, std::ios_base::binary);
	if(!ifstr.is_open())
	{
		log_err("Cannot open \"", pcFileIn, "\".");
		return false;
	}


	std::basic_ofstream<t_char> ofstr(pcFileOut, std::ios_base::binary);
	if(!ofstr.is_open())
	{
		log_err("Cannot open \"", pcFileOut, "\".");
		return false;
	}


	if(comp == Compressor::AUTO)
	{
		std::string strExt;
		if(bDecomp)
			strExt = get_fileext(std::string(pcFileIn));
		else
			strExt = get_fileext(std::string(pcFileOut));

		if(strExt == "gz")
			comp = Compressor::GZ;
		else if(strExt == "bz2")
			comp = Compressor::BZ2;
	}

	if(bDecomp)
		return decomp_stream_to_stream<t_char>(ifstr, ofstr, comp);
	else
		return comp_stream_to_stream<t_char>(ifstr, ofstr, comp);
}

template<class t_char=char>
bool comp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp=Compressor::AUTO)
{
	return __comp_file_to_file<t_char>(pcFileIn, pcFileOut, comp, 0);
}

template<class t_char=char>
bool decomp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp=Compressor::AUTO)
{
	return __comp_file_to_file<t_char>(pcFileIn, pcFileOut, comp, 1);
}


//--------------------------------------------------------------------------------



// TODO: find better way: this is too slow
template<class t_char=char>
inline bool __comp_mem_to_mem(const void* pvIn, std::size_t iLenIn,
					void*& pvOut, std::size_t& iLenOut,
					Compressor comp, bool bDecomp=0)
{
	t_char *pcIn = (char*)pvIn;
	ios::stream<ios::basic_array_source<t_char>> istr(pcIn, iLenIn);

	std::list<t_char> lstOut;
	ios::filtering_ostream arrOut(ios::back_inserter(lstOut));


	bool bOk=0;
	if(bDecomp)
		bOk = decomp_stream_to_stream<t_char>(istr, arrOut, comp);
	else
		bOk = comp_stream_to_stream<t_char>(istr, arrOut, comp);


	iLenOut = lstOut.size();
	pvOut = new t_char[iLenOut];
	t_char *pcOut = (t_char*)pvOut;

	unsigned int iIdx = 0;
	for(t_char c : lstOut)
		pcOut[iIdx++] = c;

	return bOk;
}


/*
 * example:
 *	char pc[] = "123456\nABCDEF\n\n";
 *	void *pvTst;
 *	unsigned int iLenOut=0;
 *	comp_mem_to_mem(pc, strlen(pc), pvTst,iLenOut,COMP_BZ2);
 *	std::ofstream ofstr("tst.txt.bz2");
 *	ofstr.write((char*)pvTst, iLenOut);
 *	ofstr.close();
 *	delete[] pvTst;
 */
template<class t_char=char>
bool comp_mem_to_mem(const void* pvIn, std::size_t iLenIn, void*& pvOut, std::size_t& iLenOut, Compressor comp=Compressor::GZ)
{
	return __comp_mem_to_mem<t_char>(pvIn, iLenIn, pvOut, iLenOut, comp, 0);
}

template<class t_char=char>
bool decomp_mem_to_mem(const void* pvIn, std::size_t iLenIn, void*& pvOut, std::size_t& iLenOut, Compressor comp=Compressor::AUTO)
{
	return __comp_mem_to_mem<t_char>(pvIn, iLenIn, pvOut, iLenOut, comp, 1);
}



//--------------------------------------------------------------------------------


template<class t_char=char>
inline bool __comp_mem_to_stream(const void* pvIn, std::size_t iLenIn,
					std::basic_ostream<t_char>& ostr, Compressor comp,
					bool bDecomp=0)
{
	ios::stream<ios::basic_array_source<t_char>> istr((t_char*)pvIn, iLenIn);

	bool bOk=0;
	if(bDecomp)
		bOk = decomp_stream_to_stream<t_char>(istr, ostr, comp);
	else
		bOk = comp_stream_to_stream<t_char>(istr, ostr, comp);

	return bOk;
}

template<class t_char=char>
bool comp_mem_to_stream(const void* pvIn, std::size_t iLenIn, std::basic_ostream<t_char>& ostr, Compressor comp=Compressor::GZ)
{
	return __comp_mem_to_stream<t_char>(pvIn, iLenIn, ostr, comp, 0);
}

template<class t_char=char>
bool decomp_mem_to_stream(const void* pvIn, std::size_t iLenIn, std::ostream& ostr, Compressor comp=Compressor::AUTO)
{
	return __comp_mem_to_stream<t_char>(pvIn, iLenIn, ostr, comp, 1);
}



//--------------------------------------------------------------------------------


template<class t_char=char>
inline bool __comp_mem_to_mem_fix(const void* pvIn, std::size_t iLenIn,
								void* pvOut, std::size_t iLenOut,
								Compressor comp, bool bDecomp=0)
{
	t_char *pcIn = (t_char*)pvIn;
	t_char *pcOut = (t_char*)pvOut;
	ios::stream<ios::basic_array_source<t_char>> istr(pcIn, iLenIn);
	ios::stream<ios::basic_array_sink<t_char>> ostr(pcOut, iLenOut);

	bool bOk = 0;
	if(bDecomp)
		bOk = decomp_stream_to_stream<t_char>(istr, ostr, comp);
	else
		bOk = comp_stream_to_stream<t_char>(istr, ostr, comp);

	return bOk;
}

template<class t_char=char>
bool comp_mem_to_mem_fix(const void* pvIn, std::size_t iLenIn, void* pvOut, std::size_t iLenOut, Compressor comp=Compressor::GZ)
{
	return __comp_mem_to_mem_fix<t_char>(pvIn, iLenIn, pvOut, iLenOut, comp, 0);
}

template<class t_char=char>
bool decomp_mem_to_mem_fix(const void* pvIn, std::size_t iLenIn, void* pvOut, std::size_t iLenOut, Compressor comp=Compressor::AUTO)
{
	return __comp_mem_to_mem_fix<t_char>(pvIn, iLenIn, pvOut, iLenOut, comp, 1);
}


//--------------------------------------------------------------------------------


template<class t_char=char>
std::basic_istream<t_char>* create_autodecomp_istream(std::basic_istream<t_char>& istr)
{
	typedef typename std::make_unsigned<t_char>::type t_uchar;

	t_uchar pcMagic[] = {0,0,0,0,0,0};
	std::streampos pos = istr.tellg();
	istr.read((t_char*)pcMagic, 6);
	istr.seekg(pos, std::ios::beg);
	Compressor comp = comp_from_magic(pcMagic, 6);

	ios::filtering_istream *pIstr = new ios::filtering_istream();

	if(comp == Compressor::GZ)
		pIstr->push(ios::gzip_decompressor());
	else if(comp == Compressor::BZ2)
		pIstr->push(ios::bzip2_decompressor());
	else if(comp == Compressor::Z)
		pIstr->push(ios::zlib_decompressor());
	else if(comp == Compressor::XZ)
	{
		log_err("XZ decompression not yet supported.");
		delete pIstr;
		return nullptr;
	}

	pIstr->push(istr);
	return pIstr;
}

template<class t_char=char>
std::basic_ostream<t_char>* create_comp_ostream(std::basic_ostream<t_char>& ostr,
		Compressor comp = Compressor::INVALID)
{
	ios::filtering_ostream *pOstr = new ios::filtering_ostream();

	if(comp == Compressor::GZ)
		pOstr->push(ios::gzip_compressor());
	else if(comp == Compressor::BZ2)
		pOstr->push(ios::bzip2_compressor());
	else if(comp == Compressor::Z)
		pOstr->push(ios::zlib_compressor());
	else if(comp == Compressor::XZ)
		log_err("XZ compression not yet supported.");

	pOstr->push(ostr);
	return pOstr;
}


}
#endif
