// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP
#define SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP

#include <cstddef>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>

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
        IterativeGridGenerator(
          ScalarFunction& f, base::Grid& grid, size_t N);

        /**
         * Destructor.
         */
        virtual ~IterativeGridGenerator();

        /**
         * Pure virtual method for iterative grid generation.
         *
         * @return true on success, otherwise false
         */
        virtual bool generate() = 0;

        /**
         * @return underlying grid
         */
        base::Grid& getGrid() const;

        /**
         * @return vector of function values at the grid points
         */
        const base::DataVector& getFunctionValues() const;

      protected:
        /// objective function
        ScalarFunction& f;
        /// underlying grid
        base::Grid& grid;
        /// maximal number of grid points
        size_t N;
        /// vector of function values at the grid points
        base::DataVector functionValues;

        /**
         * Removes grid points with indices
         * [oldGridSize, oldGridSize + 1, ..., grid.getStorage()->size() - 1]
         * from the grid.
         *
         * @param oldGridSize   number of grid points after removal
         */
        void undoRefinement(size_t oldGridSize);

        /**
         * Evaluates the objective function at grid points with indices
         * [oldGridSize, oldGridSize + 1, ..., grid.getStorage()->size() - 1]
         * and saves values in functionValues.
         *
         * @param oldGridSize   number of grid points already evaluated
         */
        void evalFunction(size_t oldGridSize = 0);
    };

  }
}

#endif /* SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP */
