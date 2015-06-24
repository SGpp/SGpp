// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_CONSTRAINEDOPTIMIZER_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_CONSTRAINEDOPTIMIZER_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>
#include <sgpp/optimization/function/ObjectiveFunction.hpp>
#include <sgpp/optimization/function/ConstraintFunction.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <cstddef>
#include <memory>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Abstract class for solving constrained optimization problems.
       */
      class ConstrainedOptimizer : public UnconstrainedOptimizer {
        public:
          /**
           * Constructor.
           * The starting point is set to
           * \f$(0.5, \dotsc, 0.5)^{\mathrm{T}}\f$.
           * Depending on the implementation $g$ and/or $h$ may be ignored
           * (if only equality or inequality constraints can be handled
           * by the underlying algorithm).
           *
           * @param f     function to optimize
           * @param g     inequality constraint function
           *              (\f$g(\vec{x}) \le 0\f$)
           * @param h     equality constraint function
           *              (\f$h(\vec{x}) = 0\f$)
           * @param N     maximal number of iterations or
           *              objective function evaluations
           *              (depending on the implementation)
           */
          ConstrainedOptimizer(ObjectiveFunction& f,
                               ConstraintFunction& g,
                               ConstraintFunction& h,
                               size_t N = DEFAULT_N) :
            UnconstrainedOptimizer(f, N),
            g(g),
            h(h) {
          }

          /**
           * Virtual destructor.
           */
          virtual ~ConstrainedOptimizer() {
          }

          /**
           * @return inequality constraint function
           */
          ConstraintFunction& getInequalityConstraintFunction() const {
            return g;
          }

          /**
           * @return equality constraint function
           */
          ConstraintFunction& getEqualityConstraintFunction() const {
            return h;
          }

        protected:
          /// inequality constraint function
          ConstraintFunction& g;
          /// equality constraint function
          ConstraintFunction& h;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_CONSTRAINED_CONSTRAINEDOPTIMIZER_HPP */
