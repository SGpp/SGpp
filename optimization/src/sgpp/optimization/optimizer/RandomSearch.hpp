// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_RANDOMSEARCH_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_RANDOMSEARCH_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/Optimizer.hpp>
#include <sgpp/optimization/optimizer/NelderMead.hpp>

namespace SGPP {
  namespace optimization {
    namespace optimizer {

      /**
       * Gradient-free random search.
       */
      class RandomSearch : public Optimizer {
        public:
          /// default maximal number of function evaluations
          static const size_t DEFAULT_MAX_FCN_EVAL_COUNT = 1000;

          /**
           * Constructor.
           * By default, Nelder-Mead is used as optimization algorithm.
           *
           * @param f               objective function
           * @param maxFcnEvalCount maximal number of function evaluations
           * @param populationSize  number of individual points
           *                        (default: \f$\min(10d, 100)\f$)
           */
          RandomSearch(ObjectiveFunction& f,
                       size_t maxFcnEvalCount = DEFAULT_MAX_FCN_EVAL_COUNT,
                       size_t populationSize = 0);

          /**
           * Constructor with custom optimization algorithm.
           *
           * @param optimizer        optimization algorithm and
           *                         objective function
           * @param maxFcnEvalCount  maximal number of function evaluations
           * @param populationSize   number of individual points
           *                         (default: \f$\min(10d, 100)\f$)
           */
          RandomSearch(Optimizer& optimizer,
                       size_t maxFcnEvalCount = DEFAULT_MAX_FCN_EVAL_COUNT,
                       size_t populationSize = 0);

          /**
           * @param[out] xOpt optimal point
           * @return          optimal objective function value
           */
          float_t optimize(std::vector<float_t>& xOpt);

          /**
           * @return                  number of individual points
           */
          size_t getPopulationSize() const;

          /**
           * @param populationSize    number of individual points
           */
          void setPopulationSize(size_t populationSize);

        protected:
          /// default optimization algorithm
          NelderMead defaultOptimizer;
          /// optimization algorithm
          Optimizer& optimizer;
          /// number of individual points
          size_t populationSize;
      };

    }
  }
}

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_RANDOMSEARCH_HPP */
