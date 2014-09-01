/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP
#define SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATOR_HPP

#include <vector>
#include <cstddef>

#include "opt/function/Objective.hpp"
#include "base/grid/Grid.hpp"

namespace sg {
  namespace opt {
    namespace gridgen {

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
          IterativeGridGenerator(function::Objective& f, base::Grid& grid, size_t N) :
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
          const std::vector<double>& getFunctionValues() const {
            return function_values;
          }

        protected:
          /// objective function
          function::Objective& f;
          /// underlying grid
          base::Grid& grid;
          /// maximal number of grid points
          size_t N;
          /// vector of function values at the grid points
          std::vector<double> function_values;
      };

    }
  }
}

#endif
