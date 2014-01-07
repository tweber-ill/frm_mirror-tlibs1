/*
 * Simple Script
 * @author tweber
 * @date 2013
 */

#include "script_helper.h"
#include <fstream>
#include <iostream>

char* load_file(const char* pcFile)
{
	std::ifstream ifstr(pcFile);
	if(!ifstr.is_open())
	{
		std::cerr << "Error: Cannot open \"" << pcFile << "\"." << std::endl;
		return 0;
	}

	ifstr.seekg(0, std::ios::end);
	std::size_t iFileLen = ifstr.tellg();
	ifstr.seekg(0, std::ios::beg);

	char* pcInput = new char[iFileLen+1];
	ifstr.read(pcInput, iFileLen);
	pcInput[iFileLen] = 0;

	return pcInput;
}
