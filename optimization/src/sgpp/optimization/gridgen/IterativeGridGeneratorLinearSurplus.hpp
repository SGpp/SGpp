// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP
#define SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP

#include <cstddef>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Iterative grid generation based on linear surplusses.
     * In each iteration, the fraction of \f$alpha\f$
     * (e.g. \f$\alpha = 0.2\f$ means 20%)
     * of the grid points with the largest hierarchical linear
     * surplusses are refined.
     */
    class IterativeGridGeneratorLinearSurplus : public IterativeGridGenerator {
      public:
        /// default adaptivity
        static const float_t DEFAULT_ALPHA;

        /**
         * Constructor.
         * Do not destruct the grid before this object!
         *
         * @param f         objective function
         * @param grid      grid (should be empty)
         * @param N         maximal number of grid points
         * @param alpha     adaptivity
         */
        IterativeGridGeneratorLinearSurplus(ObjectiveFunction& f,
                                            base::Grid& grid,
                                            size_t N,
                                            float_t alpha = DEFAULT_ALPHA);

        /**
         * Generate the grid.
         *
         * @return true on success, otherwise false
         */
        bool generate();

        /**
         * @return      adaptivity
         */
        float_t getAlpha() const;

        /**
         * @param alpha adaptivity
         */
        void setAlpha(float_t alpha);

      protected:
        /// pointer to a linear grid (of the same type as the "real" grid)
        std::unique_ptr<base::Grid> linearGrid;
        /// adaptivity
        float_t alpha;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORLINEARSURPLUS_HPP */
