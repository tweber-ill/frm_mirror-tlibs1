/*
 * Compression
 * @author tweber
 * @date 15-aug-2013
 */

#ifndef __M_COMP_H__
#define __M_COMP_H__

enum Compressor
{
	COMP_GZ,
	COMP_BZ2,

	COMP_AUTO
};

extern bool decomp_file_to_file(const char* pcFileIn, const char* pcFileOut, Compressor comp=COMP_AUTO);
extern bool decomp_mem_to_mem(void* pvIn, unsigned int iLenIn, void*& pvOut, unsigned int& iLenOut, Compressor comp=COMP_GZ);

#endif
