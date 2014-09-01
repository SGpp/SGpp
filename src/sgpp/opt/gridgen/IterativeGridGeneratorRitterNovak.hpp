/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATORRITTERNOVAK_HPP
#define SGPP_OPT_GRIDGEN_ITERATIVEGRIDGENERATORRITTERNOVAK_HPP

#include <vector>
#include <cstddef>

#include "opt/gridgen/IterativeGridGenerator.hpp"

namespace sg {
  namespace opt {
    namespace gridgen {

      /**
       * Iterative grid generation based on Ritter/Novak's refinement criterion.
       * Caution: This class uses HashRefinementMultiple, so it generates grids that don't meet
       * the "hierarchical ancestors" requirement!
       *
       * Literature: Erich Novak, Klaus Ritter: Global Optimization Using Hyperbolic Cross Points.
       * In: Christodoulos A. Floudas, Panos M. Pardalos (eds.): State of the Art in Global Optimization,
       * Computational Methods and Applications, Vol. 7. Springer 1996. DOI: 10.1007/978-1-4613-3437-8_2
       *
       * @see HashRefinementMultiple
       */
      class IterativeGridGeneratorRitterNovak : public IterativeGridGenerator {
        public:
          /// default adaptivity
          static const double DEFAULT_ALPHA;
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
           * @param alpha         adaptivity
           * @param max_level     maximal level of grid points
           * @param pow_method    exponentiation method
           *                      (fastPow is faster than std::pow, but only approximative)
           */
          IterativeGridGeneratorRitterNovak(function::Objective& f, base::Grid& grid, size_t N,
                                            double alpha = DEFAULT_ALPHA,
                                            size_t max_level = DEFAULT_MAX_LEVEL,
                                            PowMethod pow_method = STD_POW);

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

          /**
           * @return          maximal level of grid points
           */
          size_t getMaxLevel() const;

          /**
           * @param max_level maximal level of grid points
           */
          void setMaxLevel(size_t max_level);

        protected:
          /// adaptivity
          double alpha;
          /// maximal level of grid points
          size_t max_level;
          /// exponentiation method
          PowMethod pow_method;
      };

    }
  }
}

#endif
