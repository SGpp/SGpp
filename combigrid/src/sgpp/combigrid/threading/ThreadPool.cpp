// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/threading/ThreadPool.hpp>
#include <sgpp/combigrid/definitions.hpp>

#include <chrono>

#include <vector>

namespace sgpp {
namespace combigrid {

ThreadPool::ThreadPool(size_t numThreads)
    : numThreads(numThreads),
      threads(),
      tasks(),
      poolMutex(),
      terminateFlag(false),
      useIdleCallback(false),
      idleCallback() {}

ThreadPool::ThreadPool(size_t numThreads, IdleCallback idleCallback)
    : numThreads(numThreads),
      threads(),
      tasks(),
      poolMutex(),
      terminateFlag(false),
      useIdleCallback(true),
      idleCallback(idleCallback) {}

ThreadPool::~ThreadPool() {
  triggerTermination();
  join();
}

void ThreadPool::addTask(const Task& task) {
  CGLOG_SURROUND(std::lock_guard<std::mutex> guard(poolMutex));
  tasks.push_back(task);
}

void ThreadPool::addTasks(const std::vector<Task>& newTasks) {
  CGLOG_SURROUND(std::lock_guard<std::mutex> guard(poolMutex));
  tasks.insert(tasks.cend(), newTasks.begin(), newTasks.end());
}

void ThreadPool::start() {
  for (size_t i = 0; i < numThreads; ++i) {
    threads.push_back(std::make_shared<std::thread>([=]() {
      while (true) {
        Task nextTask;

        // wait for terminate or next task
        while (true) {
          {
            CGLOG_SURROUND(std::lock_guard<std::mutex> guard(this->poolMutex));
            if (this->terminateFlag) {
              return;
            }
            if (!this->tasks.empty()) {
              nextTask = tasks.front();
              tasks.pop_front();
              break;
            }
          }

          // no tasks, so acquire tasks

          if (useIdleCallback) {
            CGLOG_SURROUND(std::lock_guard<std::mutex> idleLock(idleMutex));

            {
              CGLOG_SURROUND(std::lock_guard<std::mutex> guard(this->poolMutex));
              if (this->terminateFlag || !this->tasks.empty()) {
                continue;
              }
            }

            idleCallback(*this);
          } else {
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
          }
        }

        // execute next task
        nextTask();
      }
    }));
  }
}

void ThreadPool::triggerTermination() {
  CGLOG_SURROUND(std::lock_guard<std::mutex> guard(poolMutex));
  terminateFlag = true;
}

void ThreadPool::join() {
  for (auto thread_ptr : threads) {
    thread_ptr->join();
  }

  threads.clear();
}

// static
void ThreadPool::terminateWhenIdle(ThreadPool& tp) { tp.triggerTermination(); }

} /* namespace combigrid */
} /* namespace sgpp*/
