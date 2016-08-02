// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STOPWATCH_HPP_
#define STOPWATCH_HPP_

#include <chrono>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace combigrid {

class Stopwatch {
  std::chrono::time_point<std::chrono::high_resolution_clock> startTime;

 public:
  Stopwatch();

  void start();

  double elapsedSeconds();
  void log();
};
}  // namespace combigrid
} /* namespace sgpp*/

#endif /* STOPWATCH_HPP_ */
