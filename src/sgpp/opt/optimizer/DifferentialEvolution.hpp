/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPTIMIZER_DIFFERENTIALEVOLUTION_HPP
#define SGPP_OPT_OPTIMIZER_DIFFERENTIALEVOLUTION_HPP

#include "opt/optimizer/Optimizer.hpp"

namespace sg {
  namespace opt {
    namespace optimizer {

      /**
       * Gradient-free Differential Evolution method.
       */
      class DifferentialEvolution : public Optimizer {
        public:
          /// default crossover probability
          static const double DEFAULT_CROSSOVER_PROBABILITY;
          /// default crossover scaling factor
          static const double DEFAULT_SCALING_FACTOR;
          /// default stopping criterion parameter 1
          static const size_t DEFAULT_IDLE_GENERATIONS_COUNT = 20;
          /// default stopping criterion parameter 2
          static const double DEFAULT_AVG_IMPROVEMENT_THRESHOLD;
          /// default stopping criterion parameter 3
          static const double DEFAULT_MAX_DISTANCE_THRESHOLD;

          /**
           * Constructor.
           *
           * @param f                         objective function
           * @param max_fcn_eval_count        maximal number of function evaluations
           * @param population_size           number of individuals (default: \f$10d\f$)
           * @param crossover_probability     crossover probability
           * @param scaling_factor            crossover scaling factor
           * @param idle_generations_count    stopping criterion parameter 1
           * @param avg_improvement_threshold stopping criterion parameter 2
           * @param max_distance_threshold    stopping criterion parameter 3
           */
          DifferentialEvolution(function::Objective& f,
                                size_t max_fcn_eval_count = DEFAULT_N,
                                size_t population_size = 0,
                                double crossover_probability = DEFAULT_CROSSOVER_PROBABILITY,
                                double scaling_factor = DEFAULT_SCALING_FACTOR,
                                size_t idle_generations_count = DEFAULT_IDLE_GENERATIONS_COUNT,
                                double avg_improvement_threshold = DEFAULT_AVG_IMPROVEMENT_THRESHOLD,
                                double max_distance_threshold = DEFAULT_MAX_DISTANCE_THRESHOLD);

          /**
           * @param[out] xopt optimal point
           * @return          optimal objective function value
           */
          double optimize(std::vector<double>& xopt);

          /**
           * @return smart pointer to cloned object
           */
          tools::SmartPointer<Optimizer> clone();

          /**
           * @return                  number of individuals
           */
          size_t getPopulationSize() const;

          /**
           * @param population_count  number of individuals
           */
          void setPopulationSize(size_t population_size);

        protected:
          /// number of individuals
          size_t population_size;
          /// crossover probability
          double crossover_probability;
          /// crossover scaling factor
          double scaling_factor;
          /// stopping criterion parameter 1
          size_t idle_generations_count;
          /// stopping criterion parameter 2
          double avg_improvement_threshold;
          /// stopping criterion parameter 3
          double max_distance_threshold;

          void initialize(
            size_t population_size,
            double crossover_probability,
            double scaling_factor,
            size_t idle_generations_count,
            double avg_improvement_threshold,
            double max_distance_threshold);
      };

    }
  }
}

#endif
