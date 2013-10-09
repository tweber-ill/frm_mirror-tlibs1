/*
 * Handles
 * @author tweber
 * @date 09-oct-2013
 */

#include "handles.h"

HandleFile::HandleFile(FILE *pFile) : m_pFile(pFile)
{}

HandleFile::~HandleFile()
{
	if(m_pFile)
	{
		fclose(m_pFile);
		m_pFile = 0;
	}
}





HandleThread::HandleThread(std::thread* pThread) : m_pThread(pThread)
{}

HandleThread::~HandleThread()
{
	if(m_pThread)
	{
		delete m_pThread;
		m_pThread = 0;
	}
}





HandleManager::HandleManager()
{
	m_vecHandles.reserve(64);
}

HandleManager::~HandleManager()
{
	for(Handle*& pHandle : m_vecHandles)
	{
		if(pHandle)
		{
			delete pHandle;
			pHandle = 0;
		}
	}

	m_vecHandles.clear();
}

Handle* HandleManager::GetHandle(unsigned int iIdx)
{
	if(iIdx >= m_vecHandles.size())
		return 0;

	return m_vecHandles[iIdx];
}

unsigned int HandleManager::AddHandle(Handle* pHandle)
{
	m_vecHandles.push_back(pHandle);
	return m_vecHandles.size()-1;
}

void HandleManager::CloseHandle(unsigned int iIdx)
{
	if(iIdx >= m_vecHandles.size())
		return;

	if(m_vecHandles[iIdx])
	{
		delete m_vecHandles[iIdx];
		m_vecHandles[iIdx] = 0;
	}
}
