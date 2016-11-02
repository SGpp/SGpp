// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_THREADING_THREADPOOL_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_THREADING_THREADPOOL_HPP_

#include <sgpp/globaldef.hpp>

#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This implements a thread-pool with a pre-specified number of threads that process a list of
 * tasks.
 */
class ThreadPool {
 public:
  typedef std::function<void()> Task;
  typedef std::function<void(ThreadPool &)> IdleCallback;

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
  /**
   * Creates a ThreadPool that processes available tasks. When no more tasks are available, the
   * threads terminate. Another way to terminate earlier is using triggerTermination().
   * The ThreadPool starts its computation only when start() is called.
   */
  explicit ThreadPool(size_t numThreads);

  /**
   * Creates a ThreadPool that processes available tasks. When no more tasks are available, the
   * callback function idleCallback is called. This callback should either add more tasks or call
   * triggerTermination(), which will cause the threads to terminate. The ThreadPool is passed as a
   * parameter to the callback.
   * The ThreadPool starts its computation only when start() is called.
   */
  ThreadPool(size_t numThreads, IdleCallback idleCallback);
  ~ThreadPool();

  /**
   * Adds a single task to the task list (thread-safe).
   */
  void addTask(Task const &task);

  /**
   * Adds a list of tasks to the task list (thread-safe).
   */
  void addTasks(std::vector<Task> const &newTasks);

  /**
   * Starts the threads.
   */
  void start();

  /**
   * Sets a termination flag (thread-safe). When a thread completes a task, it will check the
   * termination flag and terminate.
   */
  void triggerTermination();

  /**
   * Waits until all threads are finished.
   */
  void join();

  /**
   * This function can be used as a parameter to the constructor if the threads shall terminate as
   * soon as there are no tasks left.
   */
  static void terminateWhenIdle(ThreadPool &tp);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_THREADING_THREADPOOL_HPP_ */
