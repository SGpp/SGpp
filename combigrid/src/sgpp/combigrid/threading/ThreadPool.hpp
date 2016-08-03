// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_THREADING_THREADPOOL_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_THREADING_THREADPOOL_HPP_

#include <sgpp/globaldef.hpp>

#include <thread>
#include <mutex>
#include <memory>
#include <deque>
#include <functional>
#include <vector>

namespace sgpp {
namespace combigrid {

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
  explicit ThreadPool(size_t numThreads);
  ThreadPool(size_t numThreads, IdleCallback idleCallback);
  ~ThreadPool();

  void addTask(Task const &task);
  void addTasks(std::vector<Task> const &newTasks);

  void start();

  void triggerTermination();
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
