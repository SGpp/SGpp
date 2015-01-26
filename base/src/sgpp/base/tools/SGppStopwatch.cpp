/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Julian Valentin (Julian.Valentin@ipvs.uni-stuttgart.de)

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

    double SGppStopwatch::stop() {
      std::chrono::time_point<std::chrono::high_resolution_clock> end_time =
        std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed_seconds = end_time - start_time;
      return elapsed_seconds.count();
    }

  }
}
