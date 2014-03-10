/*
 * Simple Script
 * @author tweber
 * @date 2013
 */

#include "script_helper.h"
#include <fstream>
#include <iostream>

t_char* load_file(const char* pcFile)
{
	t_ifstream ifstr(pcFile/*, std::ios::binary*/);
	//ifstr.imbue(std::locale(""));

	if(!ifstr.is_open())
	{
		G_CERR << "Error: Cannot open \"" << pcFile << "\"." << std::endl;
		return 0;
	}

	ifstr.seekg(0, std::ios::end);
	std::size_t iFileLen = ifstr.tellg();
	ifstr.seekg(0, std::ios::beg);

	t_char* pcInput = new t_char[iFileLen+1];
	ifstr.read(pcInput, iFileLen);
	pcInput[iFileLen] = 0;

	if(ifstr.fail())
		G_CERR << "Error: Failed loading file \"" << pcFile << "\"." << std::endl;

	return pcInput;
}
