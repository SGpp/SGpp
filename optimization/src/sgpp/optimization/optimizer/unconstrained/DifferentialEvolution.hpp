// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_DIFFERENTIALEVOLUTION_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_DIFFERENTIALEVOLUTION_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Gradient-free Differential Evolution method.
       */
      class DifferentialEvolution : public UnconstrainedOptimizer {
        public:
          /// default crossover probability
          static constexpr float_t DEFAULT_CROSSOVER_PROBABILITY = 0.5;
          /// default crossover scaling factor
          static constexpr float_t DEFAULT_SCALING_FACTOR = 0.6;
          /// default stopping criterion parameter 1
          static const size_t DEFAULT_IDLE_GENERATIONS_COUNT = 20;
          /// default stopping criterion parameter 2
          static constexpr float_t DEFAULT_AVG_IMPROVEMENT_THRESHOLD = 1e-6;
          /// default stopping criterion parameter 3
          static constexpr float_t DEFAULT_MAX_DISTANCE_THRESHOLD = 1e-4;

          /**
           * Constructor.
           *
           * @param f                         objective function
           * @param maxFcnEvalCount           maximal number of
           *                                  function evaluations
           * @param populationSize            number of individuals
           *                                  (default: \f$10d\f$)
           * @param crossoverProbability      crossover probability
           * @param scalingFactor             crossover scaling factor
           * @param idleGenerationsCount      stopping criterion parameter 1
           * @param avgImprovementThreshold   stopping criterion parameter 2
           * @param maxDistanceThreshold      stopping criterion parameter 3
           */
          DifferentialEvolution(ScalarFunction& f,
                                size_t maxFcnEvalCount =
                                  DEFAULT_N,
                                size_t populationSize =
                                  0,
                                float_t crossoverProbability =
                                  DEFAULT_CROSSOVER_PROBABILITY,
                                float_t scalingFactor =
                                  DEFAULT_SCALING_FACTOR,
                                size_t idleGenerationsCount =
                                  DEFAULT_IDLE_GENERATIONS_COUNT,
                                float_t avgImprovementThreshold =
                                  DEFAULT_AVG_IMPROVEMENT_THRESHOLD,
                                float_t maxDistanceThreshold =
                                  DEFAULT_MAX_DISTANCE_THRESHOLD);

          /**
           * @param[out] xOpt optimal point
           * @return          optimal objective function value
           */
          float_t optimize(base::DataVector& xOpt);

          /**
           * @return                  number of individuals
           */
          size_t getPopulationSize() const;

          /**
           * @param populationSize    number of individuals
           */
          void setPopulationSize(size_t populationSize);

        protected:
          /// number of individuals
          size_t populationSize;
          /// crossover probability
          float_t crossoverProbability;
          /// crossover scaling factor
          float_t scalingFactor;
          /// stopping criterion parameter 1
          size_t idleGenerationsCount;
          /// stopping criterion parameter 2
          float_t avgImprovementThreshold;
          /// stopping criterion parameter 3
          float_t maxDistanceThreshold;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_DIFFERENTIALEVOLUTION_HPP */
