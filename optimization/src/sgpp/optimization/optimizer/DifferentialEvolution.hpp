// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_DIFFERENTIALEVOLUTION_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_DIFFERENTIALEVOLUTION_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/Optimizer.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Gradient-free Differential Evolution method.
       */
      class DifferentialEvolution : public Optimizer {
        public:
          /// default crossover probability
          static const float_t DEFAULT_CROSSOVER_PROBABILITY;
          /// default crossover scaling factor
          static const float_t DEFAULT_SCALING_FACTOR;
          /// default stopping criterion parameter 1
          static const size_t DEFAULT_IDLE_GENERATIONS_COUNT = 20;
          /// default stopping criterion parameter 2
          static const float_t DEFAULT_AVG_IMPROVEMENT_THRESHOLD;
          /// default stopping criterion parameter 3
          static const float_t DEFAULT_MAX_DISTANCE_THRESHOLD;

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
          DifferentialEvolution(const ObjectiveFunction& f,
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
          float_t optimize(std::vector<float_t>& xOpt);

          /**
           * @param[out] clone pointer to cloned object
           */
          void clone(std::unique_ptr<Optimizer>& clone) const;

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

          void initialize(
            size_t populationSize,
            float_t crossoverProbability,
            float_t scalingFactor,
            size_t idleGenerationsCount,
            float_t avgImprovementThreshold,
            float_t maxDistanceThreshold);
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_DIFFERENTIALEVOLUTION_HPP */
