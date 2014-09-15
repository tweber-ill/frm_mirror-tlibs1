/*
 * Handles
 * @author tweber
 * @date 09-oct-2013
 */

#ifndef __SCRIPT_HANDLES__
#define __SCRIPT_HANDLES__

#include "types.h"
#include "symbol.h"

#include <vector>
#include <thread>
#include <future>
#include <mutex>
#include <stdio.h>


enum HandleType : unsigned int
{
	HANDLE_FILE,

	HANDLE_THREAD,
	HANDLE_TASK,
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


template<class HANDLE, HandleType HANDLE_TYPE>
class GenericHandle : public Handle
{
protected:
	HANDLE* m_pHandle;

public:
	GenericHandle(HANDLE* pHandle) : m_pHandle(pHandle)
	{}
	virtual ~GenericHandle()
	{
		if(m_pHandle)
		{
			delete m_pHandle;
			m_pHandle = 0;
		}
	}


	virtual HandleType GetType() { return HANDLE_TYPE; }
	HANDLE* GetInternalHandle() { return m_pHandle; }
};


using HandleThread = GenericHandle<std::thread, HANDLE_THREAD>;
using HandleTask = GenericHandle<std::future<Symbol*>, HANDLE_TASK>;
using HandleMutex = GenericHandle<std::mutex, HANDLE_MUTEX>;



class HandleManager
{
protected:
	std::vector<Handle*> m_vecHandles;

public:
	HandleManager();
	virtual ~HandleManager();

	Handle* GetHandle(t_int iIdx);
	t_int AddHandle(Handle* pHandle);
	void CloseHandle(t_int iIdx);
};


#endif
