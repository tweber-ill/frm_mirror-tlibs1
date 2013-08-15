/*
 * file helper
 * @author tweber
 * @date 07-mar-2013
 */

#ifndef __MIEZE_FILE_HELPER__
#define __MIEZE_FILE_HELPER__

#include <iostream>
#include <string>

extern unsigned int get_file_size(std::istream& istr);
extern unsigned int get_file_pos(std::istream& istr);

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
};

#endif
