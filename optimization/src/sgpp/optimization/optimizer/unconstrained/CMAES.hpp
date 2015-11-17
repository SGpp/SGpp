// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_CMAES_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_CMAES_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Gradient-free CMA-ES method.
       */
      class CMAES : public UnconstrainedOptimizer {
        public:
          /// default maximal number of function evaluations
          static const size_t DEFAULT_MAX_FCN_EVAL_COUNT = 1000;

          /**
           * Constructor.
           * The starting point is set to
           * \f$(0.5, \dotsc, 0.5)^{\mathrm{T}}\f$.
           *
           * @param f                     objective function
           * @param maxFcnEvalCount       maximal number of
           *                              function evaluations
           */
          CMAES(ScalarFunction& f,
                size_t maxFcnEvalCount = DEFAULT_MAX_FCN_EVAL_COUNT);

          /**
           * Destructor.
           */
          virtual ~CMAES() override;

          virtual void optimize() override;

          /**
           * @param[out] clone pointer to cloned object
           */
          virtual void clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const override;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_CMAES_HPP */
