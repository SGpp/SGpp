// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/utils/Stopwatch.hpp>

#include <iostream>

namespace sgpp {
namespace combigrid {

Stopwatch::Stopwatch() : startTime(std::chrono::high_resolution_clock::now()) {}

void Stopwatch::start() { startTime = std::chrono::high_resolution_clock::now(); }

double Stopwatch::elapsedSeconds() {
  std::chrono::duration<double> diff = std::chrono::high_resolution_clock::now() - startTime;
  return diff.count();
}

void Stopwatch::log() { std::cout << "Time: " << elapsedSeconds() << "s." << std::endl; }
}  // namespace combigrid
} /* namespace sgpp*/
