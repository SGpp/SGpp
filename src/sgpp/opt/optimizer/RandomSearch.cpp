/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/optimizer/RandomSearch.hpp"
#include "opt/tools/Printer.hpp"
#include "opt/tools/RNG.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>

namespace sg {
  namespace opt {
    namespace optimizer {

      RandomSearch::RandomSearch(function::Objective& f,
                                 size_t max_fcn_eval_count,
                                 size_t population_size) :
        Optimizer(f, max_fcn_eval_count),
        default_optimizer(NelderMead(f)),
        optimizer(default_optimizer) {
        initialize(population_size);
      }

      RandomSearch::RandomSearch(Optimizer& optimizer,
                                 size_t max_fcn_eval_count,
                                 size_t population_size) :
        Optimizer(*optimizer.getObjectiveFunction(), max_fcn_eval_count),
        default_optimizer(NelderMead(*f)),
        optimizer(optimizer) {
        initialize(population_size);
      }

      void RandomSearch::initialize(size_t population_size) {
        if (population_size == 0) {
          this->population_size = std::min(10*f->getDimension(), static_cast<size_t>(100));
        } else {
          this->population_size = population_size;
        }
      }

      double RandomSearch::optimize(std::vector<double>& xopt) {
        tools::printer.printStatusBegin("Optimizing (random search)...");

        size_t d = f->getDimension();
        std::vector<std::vector<double> > x0(population_size, std::vector<double>(d, 0.0));
        std::vector<size_t> round_N(population_size, 0);
        size_t remaining_N = N;

        // split the number of function evaluations evenly up for all points,
        // generate pseudorandom starting points
        for (size_t k = 0; k < population_size; k++) {
          round_N[k] = static_cast<size_t>(std::ceil(static_cast<double>(remaining_N) /
                                           static_cast<double>(population_size - k)));
          remaining_N -= round_N[k];

          for (size_t t = 0; t < d; t++) {
            x0[k][t] = sg::opt::tools::rng.getUniformRN();
          }
        }

        double fopt = INFINITY;

        tools::printer.disableStatusPrinting();

        #pragma omp parallel shared(d, x0, round_N, xopt, fopt, tools::printer) default(none)
        {
          Optimizer* cur_optimizer_ptr = &optimizer;
#ifdef _OPENMP
          tools::SmartPointer<Optimizer> cur_optimizer;
          if (omp_get_max_threads() > 1) {
            cur_optimizer = tools::SmartPointer<Optimizer>(optimizer.clone());
            cur_optimizer_ptr = cur_optimizer.get();
          }
#endif

          std::vector<double> cur_xopt(d, 0.0);
          double cur_fopt;

          #pragma omp for ordered schedule(dynamic)
          for (size_t k = 0; k < population_size; k++) {
            // optimize with k-th starting point
            cur_optimizer_ptr->setStartingPoint(x0[k]);
            cur_optimizer_ptr->setN(round_N[k]);
            cur_optimizer_ptr->optimize(cur_xopt);

            cur_fopt = cur_optimizer_ptr->getObjectiveFunction()->eval(cur_xopt);

            #pragma omp critical
            {
              if (cur_fopt < fopt) {
                // this point is the best so far
                fopt = cur_fopt;
                xopt = cur_xopt;
              }
            }

            // status printing
            #pragma omp ordered
            {
              char str[10];
              snprintf(str, 10, "%.1f%%",
                       static_cast<double>(k) / static_cast<double>(population_size) * 100.0);
              tools::printer.getMutex().lock();
              tools::printer.enableStatusPrinting();
              tools::printer.printStatusUpdate(std::string(str) +
                                               ", f(x) = " + toString(fopt));
              tools::printer.disableStatusPrinting();
              tools::printer.getMutex().unlock();
            }
          }
        }

        tools::printer.enableStatusPrinting();
        tools::printer.printStatusUpdate("100.0%, f(x) = " + toString(fopt));
        tools::printer.printStatusEnd();

        return fopt;
      }

      tools::SmartPointer<Optimizer> RandomSearch::clone() {
        tools::SmartPointer<Optimizer> result(new RandomSearch(optimizer, N, population_size));
        result->setStartingPoint(x0);
        return result;
      }

      size_t RandomSearch::getPopulationSize() const {
        return population_size;
      }

      void RandomSearch::setPopulationSize(size_t population_size) {
        this->population_size = population_size;
      }

    }
  }
}
