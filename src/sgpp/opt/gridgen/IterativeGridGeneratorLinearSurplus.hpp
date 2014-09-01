/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP
#define SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP

#include <cstddef>

#include "opt/gridgen/IterativeGridGenerator.hpp"

#include "base/basis/linear/noboundary/LinearBasis.hpp"
#include "base/basis/linear/boundary/LinearBoundaryBasis.hpp"
#include "base/basis/linear/clenshawcurtis/LinearClenshawCurtisBasis.hpp"
#include "base/basis/linear/modified/ModLinearBasis.hpp"

namespace sg {
  namespace opt {
    namespace gridgen {

      /**
       * Iterative grid generation based on linear surplusses.
       * In each iteration, the fraction of \f$alpha\f$ (e.g. \f$\alpha = 0.2\f$ means 20%)
       * of the grid points with the largest hierarchical linear surplusses are refined.
       */
      class IterativeGridGeneratorLinearSurplus : public IterativeGridGenerator {
        public:
          /// default adaptivity
          static const double DEFAULT_ALPHA;

          /**
           * Constructor.
           * Do not destruct the grid before this object!
           *
           * @param f         objective function
           * @param grid      grid (should be empty)
           * @param N         maximal number of grid points
           * @param alpha     adaptivity
           */
          IterativeGridGeneratorLinearSurplus(function::Objective& f, base::Grid& grid,
                                              size_t N, double alpha = DEFAULT_ALPHA);

          /**
           * Generate the grid.
           *
           * @return true on success, otherwise false
           */
          bool generate();

          /**
           * @return      adaptivity
           */
          double getAlpha() const;

          /**
           * @param alpha adaptivity
           */
          void setAlpha(double alpha);

        protected:
          /// pointer to a linear grid (of the same type as the "real" grid)
          tools::SmartPointer<base::Grid> linear_grid;
          /// adaptivity
          double alpha;
      };

    }
  }
}

#endif
