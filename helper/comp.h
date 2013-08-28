/*
 * Compression
 * @author tweber
 * @date 15-aug-2013
 */

#ifndef __M_COMP_H__
#define __M_COMP_H__

#include <iostream>

enum Compressor
{
	COMP_GZ,
	COMP_BZ2,
	COMP_Z,

	COMP_AUTO,
	COMP_INVALID
};

extern bool decomp_stream_to_stream(std::istream& istr, std::ostream& ostr, Compressor comp=COMP_AUTO);
extern bool decomp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp=COMP_AUTO);
extern bool decomp_mem_to_mem(void* pvIn, unsigned int iLenIn, void*& pvOut, unsigned int& iLenOut, Compressor comp=COMP_AUTO);
extern bool decomp_mem_to_stream(void* pvIn, unsigned int iLenIn, std::ostream& ostr, Compressor=COMP_AUTO);

extern bool comp_stream_to_stream(std::istream& istr, std::ostream& ostr, Compressor comp=COMP_GZ);
extern bool comp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp=COMP_AUTO);
extern bool comp_mem_to_mem(void* pvIn, unsigned int iLenIn, void*& pvOut, unsigned int& iLenOut, Compressor comp=COMP_GZ);
extern bool comp_mem_to_stream(void* pvIn, unsigned int iLenIn, std::ostream& ostr, Compressor=COMP_GZ);
#endif
