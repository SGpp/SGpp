// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <boost/test/unit_test.hpp>
#include <sgpp/combigrid/integration/MCIntegrator.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

using sgpp::base::DataVector;
using sgpp::combigrid::ThreadPool;

static int counter = 0;
static std::vector<int> data;
static std::recursive_mutex dataMutex;

void idleCallback(ThreadPool &tp) {
  if (counter < 100) {
    int local_counter = counter;
    tp.addTask(ThreadPool::Task([local_counter]() {
      std::lock_guard<std::recursive_mutex> guard(dataMutex);
      data.push_back(local_counter);
    }));
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    ++counter;
  } else {
    tp.triggerTermination();
  }
}

void checkCorrectness() {
  std::sort(data.begin(), data.end());
  bool correct = true;

  for (size_t i = 0; i < data.size(); ++i) {
    if (data[i] != static_cast<int>(i)) {
      correct = false;
      std::cout << i << ", " << data[i] << "\n";
    }
  }

  BOOST_CHECK(correct);
}

BOOST_AUTO_TEST_CASE(testThreading) {
  auto tp = std::make_shared<ThreadPool>(4);
  data.clear();

  for (int i = 0; i < 100; ++i) {
    tp->addTask(ThreadPool::Task([i]() {
      std::lock_guard<std::recursive_mutex> guard(dataMutex);
      data.push_back(i);
    }));
  }

  tp->start();

  std::this_thread::sleep_for(std::chrono::milliseconds(100));

  tp->triggerTermination();
  tp->join();

  checkCorrectness();
}

BOOST_AUTO_TEST_CASE(testThreadingWithCallback) {
  auto tp = std::make_shared<ThreadPool>(4, ThreadPool::IdleCallback(idleCallback));
  data.clear();
  tp->start();
  tp->join();

  checkCorrectness();
}
