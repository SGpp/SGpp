/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    SGppStopwatch::SGppStopwatch() {
#ifdef _WIN32
      QueryPerformanceFrequency(&ticksPerSecond);
#endif
#ifndef _WIN32

#endif
    }

    SGppStopwatch::~SGppStopwatch() {
    }

    void SGppStopwatch::start() {
#ifdef _WIN32
      QueryPerformanceCounter(&begin);
#endif
#ifndef _WIN32
      gettimeofday(&begin, (struct timezone*)0);
#endif
    }

    double SGppStopwatch::stop() {
#ifdef _WIN32
      LARGE_INTEGER end;
      QueryPerformanceCounter(&end);

      double ret, ticksps;

      end.QuadPart -= begin.QuadPart;
      ret = (double)(end.QuadPart);
      ticksps = (double)(ticksPerSecond.QuadPart);
      ret /= ticksps;

      return ret;
#endif
#ifndef _WIN32
      timeval end;
      gettimeofday(&end, (struct timezone*)0);
      double seconds, useconds;
      double ret, tmp;

      if (end.tv_usec >= begin.tv_usec) {
        seconds = (double)end.tv_sec - (double)begin.tv_sec;
        useconds = (double)end.tv_usec - (double)begin.tv_usec;
      } else {
        seconds = (double)end.tv_sec - (double)begin.tv_sec;
        seconds -= 1;         // Correction
        useconds = (double)end.tv_usec - (double)begin.tv_usec;
        useconds += 1000000;      // Correction
      }

      // get time in seconds
      tmp = (double)useconds;
      ret = (double)seconds;
      tmp /= 1000000;
      ret += tmp;

      return ret;
#endif
    }

  }
}
