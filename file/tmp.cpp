/*
 * file helper
 * @author tweber
 * @date 07-mar-2013
 * @license GPLv2 or GPLv3
 */

#include "tmp.h"
#include "file.h"
#include "../math/rand.h"
#include "../string/string.h"

#include <unistd.h>
#include <fcntl.h>
#include <cstdio>

#ifndef NO_QT
	#include <QDir>
#endif


namespace tl {

TmpFile::TmpFile() : m_strPrefix("tlibs_tmp"), m_iHandle(-1)
{}

TmpFile::~TmpFile()
{
	close();
}

// cygwin does not seem to have a ::mkstemp...
int TmpFile::mkstemp(std::string& strFile)
{
	static const std::string strChars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_";
	static const unsigned int iLenChars = strChars.length();

	std::string strRnd;
	strRnd.reserve(6);
	for(unsigned int iRnd=0; iRnd<6; ++iRnd)
		strRnd.push_back(strChars[simple_rand(iLenChars)]);

	if(!find_and_replace(strFile, std::string("XXXXXX"), strRnd))
		return -1;

	int iFile = ::open(strFile.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);

	/*std::cout << "random temp file: " << strFile
			<< ", file: " << iFile
			<< std::endl;*/
	return iFile;
}

bool TmpFile::open()
{
	const std::string strMemDir = "/dev/shm";
	std::string strTmpDir = "/tmp";
#ifndef NO_QT
		strTmpDir = QDir::tempPath().toStdString();
#endif

	if(dir_exists(strMemDir.c_str()))
		m_strFile =  strMemDir;
	else
		m_strFile = strTmpDir;

	if(m_strFile[m_strFile.length()-1] != '/')
		m_strFile += "/";
	m_strFile += m_strPrefix;
	m_strFile += "_tmp.XXXXXX";

	m_iHandle = mkstemp(m_strFile);
	if(m_iHandle == -1)
		return false;

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


}
