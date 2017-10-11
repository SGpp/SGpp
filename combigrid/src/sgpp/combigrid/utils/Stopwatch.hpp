// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STOPWATCH_HPP_
#define STOPWATCH_HPP_

#include <sgpp/globaldef.hpp>

#include <chrono>

namespace sgpp {
namespace combigrid {

/**
 * Simple stopwatch implementation.
 */
class Stopwatch {
  std::chrono::time_point<std::chrono::high_resolution_clock> startTime;

 public:
  /**
   * Starts the stopwatch.
   */
  Stopwatch();

  /**
   * Re-starts the stopwatch.
   */
  void start();

  /**
   * @return the number of seconds that have passed since the last start.
   */
  double elapsedSeconds();

  /**
   * Prints the result of elapsedSeconds() to the console.
   */
  void log();
};
}  // namespace combigrid
} /* namespace sgpp*/

#endif /* STOPWATCH_HPP_ */
