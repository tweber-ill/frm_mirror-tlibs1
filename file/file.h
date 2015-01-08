/*
 * file helper
 * @author tweber
 * @date 07-mar-2013
 */

#ifndef __MIEZE_FILE_HELPER__
#define __MIEZE_FILE_HELPER__

#include <iostream>
#include <string>

template<typename t_char=char>
std::streampos get_file_size(std::basic_istream<t_char>& istr)
{
        std::streampos iPos = istr.tellg();

        istr.seekg(0, std::ios_base::end);
        std::streampos iSize = istr.tellg();
        istr.seekg(iPos, std::ios_base::beg);

        return iSize;
}

template<typename t_char=char>
std::streampos get_file_pos(std::basic_istream<t_char>& istr)
{
        std::streampos iPos = istr.tellg();

        if(iPos < 0) return 0;
        return iPos;
}


extern bool dir_exists(const char* pcDir);


class TmpFile
{
protected:
	std::string m_strFile;
	std::string m_strPrefix;
	int m_iHandle;

public:
	TmpFile();
	virtual ~TmpFile();

	bool open();
	void close();
	const std::string& GetFileName() const;

	void SetPrefix(const char* pcStr);

	static int mkstemp(std::string& strFile);
};

#endif
