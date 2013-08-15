/*
 * file helper
 * @author tweber
 * @date 07-mar-2013
 */

#include "file.h"
#include <unistd.h>
#include <stdio.h>
#include <QtCore/QDir>

unsigned int get_file_size(std::istream& istr)
{
        std::streampos iPos = istr.tellg();

        istr.seekg(0, std::ios_base::end);
        std::streampos iSize = istr.tellg();

        istr.seekg(iPos, std::ios_base::beg);

        return (unsigned int)iSize;
}

unsigned int get_file_pos(std::istream& istr)
{
        std::streampos iPos = istr.tellg();

        if(iPos < 0) return 0;
        return (unsigned int) iPos;
}

bool dir_exists(const char* pcDir)
{
	return QDir(pcDir).exists();
}




TmpFile::TmpFile() : m_strPrefix("cattus"), m_iHandle(-1)
{}

TmpFile::~TmpFile()
{
	close();
}

bool TmpFile::open()
{
	char pcTmpFile[256];
	const char pcMemDir[] = "/dev/shm";

	if(dir_exists(pcMemDir))
		strcpy(pcTmpFile, pcMemDir);
	else
		strcpy(pcTmpFile, QDir::tempPath().toStdString().c_str());


	if(pcTmpFile[strlen(pcTmpFile)-1] != '/')
		strcat(pcTmpFile, "/");
	strcat(pcTmpFile, m_strPrefix.c_str());
	strcat(pcTmpFile, "_tmp.XXXXXX");

	m_iHandle = mkstemp(pcTmpFile);
	if(m_iHandle == -1)
		return false;

	m_strFile = pcTmpFile;
	//std::cout << "temp file: " << m_strFile << std::endl;
	return true;
}

void TmpFile::close()
{
	if(m_iHandle != -1)
		::close(m_iHandle);

	if(m_strFile != "")
		remove(m_strFile.c_str());
}

const std::string& TmpFile::GetFileName() const
{
	return m_strFile;
}

void TmpFile::SetPrefix(const char* pcStr)
{
	m_strPrefix = pcStr;
}
