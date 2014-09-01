/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include <algorithm>
#include <cmath>
#include <ctime>

#include "opt/tools/RNG.hpp"

namespace sg {
  namespace opt {
    namespace tools {

      RNG rng;

      RNG::RNG() :
        seed(0) {
        initialize_mt();
      }

// source: https://en.wikipedia.org/w/index.php?title=Mersenne_twister&oldid=617192530
      void RNG::initialize_mt() {
        mt_index = 0;
        mt_state[0] = seed;

        for (size_t i = 1; i < 624; i++) {
          mt_state[i] = 0xffffffff & (i + 0x6c078965 * (mt_state[i-1] ^ (mt_state[i-1] >> 30)));
        }
      }

// source: https://en.wikipedia.org/w/index.php?title=Mersenne_twister&oldid=617192530
      size_t RNG::run_mt() {
        if (mt_index == 0) {
          for (size_t i = 0; i < 624; i++) {
            size_t y = (0x80000000 & mt_state[i]) | (0x7fffffff & mt_state[(i+1) % 624]);
            mt_state[i] = mt_state[(i+397) % 624] ^ (y >> 1);

            if (y % 2 == 1) {
              mt_state[i] ^= 0x9908b0df;
            }
          }
        }

        size_t y = mt_state[mt_index];
        y ^= (y >> 11);
        y ^= 0x9d2c5680 & (y << 7);
        y ^= 0xefc60000 & (y << 15);
        y ^= (y >> 18);

        mt_index = (mt_index + 1) % 624;

        return y;
      }

      double RNG::getUniformRN(double a, double b) {
        run_mt();
        return a + (b-a) * static_cast<double>(run_mt()) / static_cast<double>(0xffffffff);
      }

      size_t RNG::getUniformIndexRN(size_t size) {
        return std::min(size - 1,
                        static_cast<size_t>(RNG::getUniformRN(0.0, static_cast<double>(size))));
      }

// source: /usr/include/c++/4.8/bits/random.h, /usr/include/c++/4.8/bits/random.tcc
      double RNG::getGaussianRN(double std_dev, double mean) {
        static bool has_saved = false;
        static double saved;

        double standard_gaussian;

        if (has_saved) {
          has_saved = false;
          standard_gaussian = saved;
        } else {
          double x, y, r2;

          do {
            x = getUniformRN(-1.0, 1.0);
            y = getUniformRN(-1.0, 1.0);
            r2 = x*x + y*y;
          } while ((r2 > 1.0) || (r2 == 0.0));

          double mult = std::sqrt(-2.0 * std::log(r2) / r2);
          saved = x * mult;
          has_saved = true;
          standard_gaussian = y * mult;
        }

        return mean + std_dev * standard_gaussian;
      }

      RNG::SeedType RNG::getSeed() {
        return seed;
      }

      void RNG::setSeed() {
        setSeed(static_cast<SeedType>(std::time(0)));
      }

      void RNG::setSeed(RNG::SeedType seed) {
        this->seed = seed;
        initialize_mt();
      }

    }
  }
}
