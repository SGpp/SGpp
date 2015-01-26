// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef SGPPSTOPWATCH_H
#define SGPPSTOPWATCH_H

#include <chrono>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     *  OS-independent version of a stop watch (using std::chrono).
     *
     *  Part of SGpp, so you can easily calculate the needed time of SGpp computations with a high precision
     */
    class SGppStopwatch {
      public:
        /**
         * Constructor. Resets the stop watch.
         */
        SGppStopwatch();

        /**
         * Destructor.
         */
        ~SGppStopwatch();

        /**
         * Starts the stop watch.
         */
        void start();

        /**
         * Stops the stop watch.
         *
         * @return elapsed time since the last call to start in seconds
         */
        double stop();

      protected:
        std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    };

  }
}

#endif  /* SGPPSTOPWATCH_H */