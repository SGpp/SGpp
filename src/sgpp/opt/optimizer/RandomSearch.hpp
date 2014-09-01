/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPTIMIZER_RANDOMSEARCH_HPP
#define SGPP_OPT_OPTIMIZER_RANDOMSEARCH_HPP

#include "opt/optimizer/Optimizer.hpp"
#include "opt/optimizer/NelderMead.hpp"

namespace sg {
  namespace opt {
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
           * @param f                     objective function
           * @param max_fcn_eval_count    maximal number of function evaluations
           * @param population_size       number of individual points (default: \f$\min(10d, 100)\f$)
           */
          RandomSearch(function::Objective& f,
                       size_t max_fcn_eval_count = DEFAULT_MAX_FCN_EVAL_COUNT,
                       size_t population_size = 0);

          /**
           * Constructor with custom optimization algorithm.
           *
           * @param optimizer             optimization algorithm and objective function
           * @param max_fcn_eval_count    maximal number of function evaluations
           * @param population_size       number of individual points (default: \f$\min(10d, 100)\f$)
           */
          RandomSearch(Optimizer& optimizer,
                       size_t max_fcn_eval_count = DEFAULT_MAX_FCN_EVAL_COUNT,
                       size_t population_size = 0);

          /**
           * @return smart pointer to cloned object
           */
          tools::SmartPointer<Optimizer> clone();

          /**
           * @param[out] xopt optimal point
           * @return          optimal objective function value
           */
          double optimize(std::vector<double>& xopt);

          /**
           * @return                  number of individual points
           */
          size_t getPopulationSize() const;

          /**
           * @param population_size   number of individual points
           */
          void setPopulationSize(size_t population_size);

        protected:
          /// default optimization algorithm
          NelderMead default_optimizer;
          /// optimization algorithm
          Optimizer& optimizer;
          /// number of individual points
          size_t population_size;

          /**
           * Internal function for initializing the member variables.
           *
           * @param population_size   number of individual points
           */
          void initialize(size_t population_size);
      };

    }
  }
}

#endif
