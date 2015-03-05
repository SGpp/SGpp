// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP
#define SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP

#include <vector>
#include <cstddef>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/ObjectiveFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Abstract base class for iterative grid generation methods.
     */
    class IterativeGridGenerator {
      public:
        /**
         * Constructor.
         * Do not destruct the grid before this object!
         *
         * @param f     objective function
         * @param grid  grid (should be empty)
         * @param N     maximal number of grid points
         */
        IterativeGridGenerator(ObjectiveFunction& f,
                               base::Grid& grid, size_t N) :
          f(f), grid(grid), N(N) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~IterativeGridGenerator() {
        }

        /**
         * Pure virtual method for iterative grid generation.
         *
         * @return true on success, otherwise false
         */
        virtual bool generate() = 0;

        /**
         * @return underlying grid
         */
        base::Grid& getGrid() const {
          return grid;
        }

        /**
         * @return vector of function values at the grid points
         */
        const std::vector<float_t>& getFunctionValues() const {
          return functionValues;
        }

      protected:
        /// objective function
        ObjectiveFunction& f;
        /// underlying grid
        base::Grid& grid;
        /// maximal number of grid points
        size_t N;
        /// vector of function values at the grid points
        std::vector<float_t> functionValues;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP */
