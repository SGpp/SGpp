// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NELDERMEAD_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NELDERMEAD_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Gradient-free Nelder-Mead method.
       */
      class NelderMead : public UnconstrainedOptimizer {
        public:
          /// default reflection coefficient
          static constexpr float_t DEFAULT_ALPHA = 1.0;
          /// default expansion coefficient
          static constexpr float_t DEFAULT_BETA = 2.0;
          /// default contraction coefficient
          static constexpr float_t DEFAULT_GAMMA = 0.5;
          /// default shrinking coefficient
          static constexpr float_t DEFAULT_DELTA = 0.5;
          /// default maximal number of function evaluations
          static const size_t DEFAULT_MAX_FCN_EVAL_COUNT = 1000;
          /// edge length of starting simplex
          static constexpr float_t STARTING_SIMPLEX_EDGE_LENGTH = 0.4;

          /**
           * Constructor.
           * The starting point is set to
           * \f$(0.5, \dotsc, 0.5)^{\mathrm{T}}\f$.
           *
           * @param f                     objective function
           * @param maxFcnEvalCount       maximal number of
           *                              function evaluations
           * @param alpha                 reflection coefficient
           * @param beta                  expansion coefficient
           * @param gamma                 contraction coefficient
           * @param delta                 shrinking coefficient
           */
          NelderMead(ScalarFunction& f,
                     size_t maxFcnEvalCount = DEFAULT_MAX_FCN_EVAL_COUNT,
                     float_t alpha = DEFAULT_ALPHA,
                     float_t beta = DEFAULT_BETA,
                     float_t gamma = DEFAULT_GAMMA,
                     float_t delta = DEFAULT_DELTA);

          void optimize();

          /**
           * @return          reflection coefficient
           */
          float_t getAlpha() const;

          /**
           * @param alpha     reflection coefficient
           */
          void setAlpha(float_t alpha);

          /**
           * @return          expansion coefficient
           */
          float_t getBeta() const;

          /**
           * @param beta      expansion coefficient
           */
          void setBeta(float_t beta);

          /**
           * @return          contraction coefficient
           */
          float_t getGamma() const;

          /**
           * @param gamma     contraction coefficient
           */
          void setGamma(float_t gamma);

          /**
           * @return          shrinking coefficient
           */
          float_t getDelta() const;

          /**
           * @param delta     shrinking coefficient
           */
          void setDelta(float_t delta);

          /**
           * @param[out] clone pointer to cloned object
           */
          void clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const;

        protected:
          /// reflection coefficient
          float_t alpha;
          /// expansion coefficient
          float_t beta;
          /// contraction coefficient
          float_t gamma;
          /// shrinking coefficient
          float_t delta;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_NELDERMEAD_HPP */
