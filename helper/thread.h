/*
 * Thread helpers
 * @author tweber
 * @date aug-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_THREAD_H__
#define __TLIBS_THREAD_H__

#include <future>
#include <thread>
#include <mutex>
#include <list>
#include <functional>
#include <type_traits>
#include <condition_variable>
//#include <iostream>

namespace tl {

template<class t_func>
class ThreadPool
{
	public:
		using t_ret = typename std::result_of<t_func&()>::type;
		using t_fut = std::list<std::future<t_ret>>;
		using t_task = std::list<std::packaged_task<t_ret()>>;

	protected:
		std::list<std::thread*> m_lstThreads;
		t_task m_lstTasks;
		t_fut m_lstFutures;

		std::mutex m_mtx, m_mtxStart;
		std::condition_variable m_signalStart;
		bool m_bStart = 0;

	public:
		ThreadPool(unsigned int iNumThreads = std::thread::hardware_concurrency())
		{
			for(unsigned int i=0; i<iNumThreads; ++i)
			{
				m_lstThreads.push_back(new std::thread([this]()
				{
					std::unique_lock<std::mutex> lock(m_mtxStart);
					m_signalStart.wait(lock, [this]()->bool { return m_bStart; });
					lock.unlock();

					while(1)
					{
						std::unique_lock<std::mutex> lock0(m_mtx);
						if(m_lstTasks.size() > 0)
						{
							std::packaged_task<t_ret()> task = std::move(m_lstTasks.front());
							m_lstTasks.pop_front();
							lock0.unlock();

							task();
						}
						else
							break;
					}
				}));
			}
		}

		virtual ~ThreadPool()
		{
			JoinAll();

			for(std::thread *pThread : m_lstThreads)
				delete pThread;
			m_lstThreads.clear();
		}

		void AddTask(const std::function<t_func>& fkt)
		{
			std::packaged_task<t_ret()> task(fkt);
			std::future<t_ret> fut = task.get_future();

			std::lock_guard<std::mutex> lock(m_mtx);
			m_lstTasks.push_back(std::move(task));
			m_lstFutures.push_back(std::move(fut));
		}

		void StartTasks()
		{
			std::lock_guard<std::mutex> lock(m_mtxStart);
			m_bStart = 1;
			m_signalStart.notify_all();
		}

		t_fut& GetFutures() { return m_lstFutures; }
		t_task& GetTasks() { return m_lstTasks; }

		void JoinAll()
		{
			for(std::thread *pThread : m_lstThreads)
				pThread->join();
		}
};

}
#endif
