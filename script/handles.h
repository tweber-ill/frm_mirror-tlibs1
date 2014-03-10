/*
 * Handles
 * @author tweber
 * @date 09-oct-2013
 */

#ifndef __SCRIPT_HANDLES__
#define __SCRIPT_HANDLES__

#include "types.h"

#include <vector>
#include <thread>
#include <stdio.h>


enum HandleType
{
	HANDLE_FILE,

	HANDLE_THREAD,
	HANDLE_MUTEX
};




class Handle
{
public:
	virtual ~Handle() {}

	virtual HandleType GetType() = 0;
};


class HandleFile : public Handle
{
protected:
	FILE *m_pFile;

public:
	HandleFile(FILE *pFile);
	virtual ~HandleFile();

	virtual HandleType GetType() { return HANDLE_FILE; }
};


class HandleThread : public Handle
{
protected:
	std::thread* m_pThread;

public:
	HandleThread(std::thread* pThread);
	virtual ~HandleThread();

	virtual HandleType GetType() { return HANDLE_THREAD; }
	std::thread* GetInternalHandle() { return m_pThread; }
};



class HandleManager
{
protected:
	std::vector<Handle*> m_vecHandles;

public:
	HandleManager();
	virtual ~HandleManager();

	Handle* GetHandle(unsigned int iIdx);
	unsigned int AddHandle(Handle* pHandle);
	void CloseHandle(unsigned int iIdx);
};


#endif
