// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORRITTERNOVAK_HPP
#define SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORRITTERNOVAK_HPP

#include <vector>
#include <cstddef>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Iterative grid generation based on Ritter/Novak's refinement criterion.
     * Caution: This class uses HashRefinementMultiple, so it generates grids
     * that don't meet the "hierarchical ancestors" requirement!
     *
     * Literature: Erich Novak, Klaus Ritter: Global Optimization Using
     * Hyperbolic Cross Points.
     * In: Christodoulos A. Floudas, Panos M. Pardalos (eds.): State of the
     * Art in Global Optimization,
     * Computational Methods and Applications, Vol. 7. Springer 1996.
     * DOI: 10.1007/978-1-4613-3437-8_2
     *
     * @see HashRefinementMultiple
     */
    class IterativeGridGeneratorRitterNovak : public IterativeGridGenerator {
      public:
        /// default adaptivity
        static constexpr float_t DEFAULT_GAMMA = 0.85;
        /// default maximal level of grid points
        static const size_t DEFAULT_MAX_LEVEL = 15;

        /// exponentiation methods
        enum PowMethod {
          STD_POW, FAST_POW
        };

        /**
         * Constructor.
         * Do not destruct the grid before this object!
         *
         * @param f             objective function
         * @param grid          grid (should be empty)
         * @param N             maximal number of grid points
         * @param gamma         adaptivity
         * @param maxLevel      maximal level of grid points
         * @param powMethod     exponentiation method
         *                      (fastPow is faster than std::pow,
         *                      but only approximative)
         */
        IterativeGridGeneratorRitterNovak(ObjectiveFunction& f,
                                          base::Grid& grid,
                                          size_t N,
                                          float_t gamma = DEFAULT_GAMMA,
                                          size_t maxLevel = DEFAULT_MAX_LEVEL,
                                          PowMethod powMethod = STD_POW);

        /**
         * Generate the grid.
         *
         * @return true on success, otherwise false
         */
        bool generate();

        /**
         * @return      adaptivity
         */
        float_t getGamma() const;

        /**
         * @param gamma adaptivity
         */
        void setGamma(float_t gamma);

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
        float_t gamma;
        /// maximal level of grid points
        size_t maxLevel;
        /// exponentiation method
        PowMethod powMethod;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORRITTERNOVAK_HPP */
