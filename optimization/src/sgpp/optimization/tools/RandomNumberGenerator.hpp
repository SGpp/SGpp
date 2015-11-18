// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TOOLS_RANDOMNUMBERGENERATOR_HPP
#define SGPP_OPTIMIZATION_TOOLS_RANDOMNUMBERGENERATOR_HPP

#include <sgpp/globaldef.hpp>
#include <random>

namespace SGPP {
  namespace optimization {

    /**
     * Singleton class for generating pseudo-random numbers
     * (wrapper around std::mt19937 from &lt;random&gt;).
     */
    class RandomNumberGenerator {
      public:
        /// type of the seed
        typedef std::mt19937::result_type SeedType;

        /**
         * @return singleton instance
         */
        inline static RandomNumberGenerator& getInstance() {
          static RandomNumberGenerator rng;
          return rng;
        }

        /**
         * Generate a uniform pseudo-random number.
         *
         * @param a lower bound
         * @param b upper bound
         * @return  uniform pseudo-random number in \f$[a, b]\f$
         */
        float_t getUniformRN(float_t a = 0.0, float_t b = 1.0);

        /**
         * Generate a uniform pseudo-random array index.
         *
         * @param size  size of the array
         * @return      discrete uniform pseudo-random number in
         *              \f$\{0, \dotsc, \text{\texttt{size}} - 1\}\f$
         */
        size_t getUniformIndexRN(size_t size);

        /**
         * Generate a Gaussian pseudo-random number.
         *
         * @param stdDev    standard deviation of the Gaussian distribution
         * @param mean      mean of the Gaussian distribution
         * @return          Gaussian pseudo-random number
         */
        float_t getGaussianRN(float_t stdDev = 1.0, float_t mean = 0.0);

        /**
         * @return      seed
         */
        SeedType getSeed() const;

        /**
         * Reseeds with time-dependent seed.
         */
        void setSeed();

        /**
         * Reseeds.
         *
         * @param seed  seed to be used
         */
        void setSeed(SeedType seed);

      protected:
        /// random number generator
        std::mt19937 generator;
        /// seed
        SeedType seed;

      private:
        /**
         * Constructor, initializes with time-dependent seed.
         */
        RandomNumberGenerator();

        /**
         * Deleted copy constructor.
         */
        RandomNumberGenerator(const RandomNumberGenerator&) = delete;

        /**
         * Deleted assignment operator.
         */
        void operator=(const RandomNumberGenerator&) = delete;
    };
  }
}

#endif /* SGPP_OPTIMIZATION_TOOLS_RANDOMNUMBERGENERATOR_HPP */
