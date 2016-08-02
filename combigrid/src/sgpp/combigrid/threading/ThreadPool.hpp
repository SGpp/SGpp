/*
 * ThreadPool.hpp
 *
 *  Created on: 20.06.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_THREADING_THREADPOOL_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_THREADING_THREADPOOL_HPP_

#include <sgpp/globaldef.hpp>

#include <thread>
#include <mutex>
#include <memory>
#include <deque>
#include <functional>
#include <vector>
#include <condition_variable>

namespace SGPP {
namespace combigrid {

class ThreadPool {
public:

	typedef std::function<void()> Task;
	typedef std::function<void(ThreadPool&)> IdleCallback;

private:

	size_t numThreads;
	std::vector<std::shared_ptr<std::thread>> threads;
	std::deque<Task> tasks;
	std::mutex poolMutex;
	std::mutex idleMutex;
	bool terminateFlag;
	bool useIdleCallback;
	IdleCallback idleCallback;

public:
	ThreadPool(size_t numThreads);
	ThreadPool(size_t numThreads, IdleCallback idleCallback);
	~ThreadPool();

	void addTask(Task const &task);
	void addTasks(std::vector<Task> const &newTasks);

	void start();

	void triggerTermination();
	void join();

	/**
	 * This function can be used as a parameter to the constructor if the threads shall terminate as soon as there are no tasks left.
	 */
	static void terminateWhenIdle(ThreadPool &tp);
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_THREADING_THREADPOOL_HPP_ */
