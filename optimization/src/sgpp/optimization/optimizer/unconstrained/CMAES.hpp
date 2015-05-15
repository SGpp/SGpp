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

      class CMAES : public UnconstrainedOptimizer {
        public:
          // TODO
          /// default initial step size
          //static constexpr float_t DEFAULT_INITIAL_STEP_SIZE = 0.3;
          static constexpr float_t DEFAULT_INITIAL_STEP_SIZE = 0.01;

          /**
           * Constructor.
           *
           * @param f                 objective function
           * @param maxFcnEvalCount   maximal number of
           *                          function evaluations
           * @param populationSize    number of individuals
           * @param recombinationSize number of individuals used to recombine
           * @param initialStepSize   initial step size
           */
          CMAES(ObjectiveFunction& f,
                size_t maxFcnEvalCount = DEFAULT_N,
                size_t populationSize = 0,
                size_t recombinationSize = 0,
                float_t initialStepSize = DEFAULT_INITIAL_STEP_SIZE);

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

          /**
           * @return                  number of individuals used to recombine
           */
          size_t getRecombinationSize() const;

          /**
           * @param recombinationSize number of individuals used to recombine
           */
          void setRecombinationSize(size_t recombinationSize);

          /**
           * @return                  initial step size
           */
          float_t getInitialStepSize() const;

          /**
           * @param initialStepSize   initial step size
           */
          void setInitialStepSize(float_t initialStepSize);

        protected:
          /// number of individuals
          size_t lambda;
          /// number of individuals used to recombine
          size_t mu;
          /// initial step size
          float_t sigma0;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_CMAES_HPP */
