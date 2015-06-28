// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORSOO_HPP
#define SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORSOO_HPP

#include <cstddef>

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * TODO
     */
    class IterativeGridGeneratorSOO : public IterativeGridGenerator {
      public:
        /// default adaptivity
        //static constexpr float_t DEFAULT_GAMMA = 0.85;
        /// default maximal level of grid points
        static const size_t DEFAULT_MAX_LEVEL = 20;

        /**
         * Constructor.
         * Do not destruct the grid before this object!
         *
         * @param f             objective function
         * @param grid          grid (should be empty)
         * @param N             maximal number of grid points
         * @param maxLevel      maximal level of grid points
         */
        IterativeGridGeneratorSOO(ObjectiveFunction& f,
                                  base::Grid& grid,
                                  size_t N,
                                  //float_t gamma = DEFAULT_GAMMA,
                                  size_t maxLevel = DEFAULT_MAX_LEVEL);

        /**
         * Generate the grid.
         *
         * @return true on success, otherwise false
         */
        bool generate();

        /*
         * @return      adaptivity
         */
        //float_t getGamma() const;

        /*
         * @param gamma adaptivity
         */
        //void setGamma(float_t gamma);

        /**
         * @return          maximal level of grid points
         */
        size_t getMaxLevel() const;

        /**
         * @param maxLevel  maximal level of grid points
         */
        void setMaxLevel(size_t maxLevel);

      protected:
        /// adaptivity
        //float_t gamma;
        /// maximal level of grid points
        size_t maxLevel;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORSOO_HPP */
