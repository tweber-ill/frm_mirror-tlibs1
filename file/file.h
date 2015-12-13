/*
 * file helper
 * @author tweber
 * @date 07-mar-2013
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIB_FILE_HELPER__
#define __TLIB_FILE_HELPER__

#include <iostream>
#include <list>
#include <string>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>


namespace tl {

// -----------------------------------------------------------------------------

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


template<typename t_char=char>
bool dir_exists(const t_char* pcDir)
{
	boost::filesystem::path path(pcDir);
	bool bExists = boost::filesystem::exists(path);
	bool bIsDir = boost::filesystem::is_directory(path);

	return bExists && bIsDir;
}

template<typename t_char=char>
bool file_exists(const t_char* pcDir)
{
	boost::filesystem::path path(pcDir);
	bool bExists = boost::filesystem::exists(path);
	bool bIsDir = boost::filesystem::is_directory(path);
	bool bIsFile = boost::filesystem::is_regular_file(path);
	bool bIsLink = boost::filesystem::is_symlink(path);

	return bExists && (bIsFile || bIsLink) && !bIsDir;
}


// -----------------------------------------------------------------------------


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


// -----------------------------------------------------------------------------

}

#endif
