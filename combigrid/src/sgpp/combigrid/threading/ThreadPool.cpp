// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>

#include <chrono>
#include <vector>

namespace sgpp {
namespace combigrid {

ThreadPool::IdleCallback ThreadPool::terminateWhenIdle((ThreadPool::doTerminateWhenIdle));

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
  CGLOG_SURROUND(std::lock_guard<std::recursive_mutex> guard(poolMutex));
  tasks.push_back(task);
}

void ThreadPool::addTasks(const std::vector<Task>& newTasks) {
  CGLOG_SURROUND(std::lock_guard<std::recursive_mutex> guard(poolMutex));
  tasks.insert(tasks.end(), newTasks.begin(), newTasks.end());
}

void ThreadPool::start() {
  for (size_t i = 0; i < numThreads; ++i) {
    threads.push_back(std::make_shared<std::thread>([this]() {
      while (true) {
        Task nextTask;

        // wait for terminate or next task
        while (true) {
          {
            CGLOG_SURROUND(std::lock_guard<std::recursive_mutex> guard(this->poolMutex));
            if (this->terminateFlag) {
              return;
            }
            if (!this->tasks.empty()) {
              nextTask = tasks.front();
              tasks.pop_front();
              break;
            } else if (!useIdleCallback) {
              return;
            }
            CGLOG("leave outer guard(this->poolMutex)");
          }

          // no tasks, so acquire tasks

          if (useIdleCallback) {
            CGLOG_SURROUND(std::lock_guard<std::recursive_mutex> idleLock(idleMutex));

            {
              CGLOG_SURROUND(std::lock_guard<std::recursive_mutex> guard(this->poolMutex));
              if (this->terminateFlag || !this->tasks.empty()) {
                CGLOG("leave inner guard(this->poolMutex)");
                continue;
              }
              CGLOG("leave inner guard(this->poolMutex)");
            }

            idleCallback(*this);
            CGLOG("leave idleLock(idleMutex)");
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
  CGLOG_SURROUND(std::lock_guard<std::recursive_mutex> guard(poolMutex));
  terminateFlag = true;
}

void ThreadPool::join() {
  for (auto thread_ptr : threads) {
    thread_ptr->join();
  }

  threads.clear();
}

// static
void ThreadPool::doTerminateWhenIdle(ThreadPool& tp) { tp.triggerTermination(); }

} /* namespace combigrid */
} /* namespace sgpp*/
