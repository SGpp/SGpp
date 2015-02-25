// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    SGppStopwatch::SGppStopwatch() :
      start_time(std::chrono::high_resolution_clock::now()) {
    }

    SGppStopwatch::~SGppStopwatch() {
    }

    void SGppStopwatch::start() {
      start_time = std::chrono::high_resolution_clock::now();
    }

    float_t SGppStopwatch::stop() {
      std::chrono::time_point<std::chrono::high_resolution_clock> end_time =
        std::chrono::high_resolution_clock::now();
      std::chrono::duration<float_t> elapsed_seconds = end_time - start_time;
      return elapsed_seconds.count();
    }

  }
}