/*
 * file helper
 * @author tweber
 * @date 07-mar-2013
 */

#ifndef __MIEZE_FILE_HELPER__
#define __MIEZE_FILE_HELPER__

#include<iostream>

inline unsigned int get_file_size(std::istream& istr)
{
        std::streampos iPos = istr.tellg();

        istr.seekg(0, std::ios_base::end);
        std::streampos iSize = istr.tellg();

        istr.seekg(iPos, std::ios_base::beg);

        return (unsigned int)iSize;
}

inline unsigned int get_file_pos(std::istream& istr)
{
        std::streampos iPos = istr.tellg();

        if(iPos < 0) return 0;
        return (unsigned int) iPos;
}

#endif
