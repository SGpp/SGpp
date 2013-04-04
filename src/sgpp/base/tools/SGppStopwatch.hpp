/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SGPPSTOPWATCH_H
#define SGPPSTOPWATCH_H

#ifdef _WIN32
#include <windows.h>
#endif
#ifndef _WIN32
#include <sys/time.h>
#endif

namespace sg {
  namespace base {

    /**
     *  OS-independent (per Preprocessor) version of a stopwatch
     *
     *  Part of SGpp, so you can easily calculate the needed time of SGpp computations with a high precision
     */
    class SGppStopwatch {
      private:
#ifdef _WIN32
        LARGE_INTEGER ticksPerSecond;
        LARGE_INTEGER begin;
#endif
#ifndef _WIN32
        timeval begin;
#endif

      public:
        /**
         *  Constructor
         *
         *  resets the Stopwatch
         */
        SGppStopwatch();

        /**
         *  Destructor
         */
        ~SGppStopwatch();

        /**
         *  starts the time measuring
         */
        void start();

        /**
         *  stops time measuring
         *
         *  \return measured time in seconds
         */
        double stop();
    };

  }
}

#endif  /* SGPPSTOPWATCH_H */
