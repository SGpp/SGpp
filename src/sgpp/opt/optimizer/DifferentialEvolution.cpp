/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/optimizer/DifferentialEvolution.hpp"
#include "opt/tools/Printer.hpp"
#include "opt/tools/RNG.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>

namespace sg {
  namespace opt {
    namespace optimizer {

      const double DifferentialEvolution::DEFAULT_CROSSOVER_PROBABILITY = 0.5;
      const double DifferentialEvolution::DEFAULT_SCALING_FACTOR = 0.6;
      const double DifferentialEvolution::DEFAULT_AVG_IMPROVEMENT_THRESHOLD = 1e-6;
      const double DifferentialEvolution::DEFAULT_MAX_DISTANCE_THRESHOLD = 1e-4;

      DifferentialEvolution::DifferentialEvolution(function::Objective& f,
          size_t max_fcn_eval_count, size_t population_size,
          double crossover_probability, double scaling_factor,
          size_t idle_generations_count, double avg_improvement_threshold,
          double max_distance_threshold) :
        Optimizer(f, max_fcn_eval_count) {
        initialize(population_size, crossover_probability, scaling_factor,
                   idle_generations_count, avg_improvement_threshold, max_distance_threshold);
      }

      void DifferentialEvolution::initialize(
        size_t population_size,
        double crossover_probability,
        double scaling_factor,
        size_t idle_generations_count,
        double avg_improvement_threshold,
        double max_distance_threshold) {
        if (population_size == 0) {
          this->population_size = 10 * f->getDimension();
        } else {
          this->population_size = population_size;
        }

        this->crossover_probability = crossover_probability;
        this->scaling_factor = scaling_factor;
        this->idle_generations_count = idle_generations_count;
        this->avg_improvement_threshold = avg_improvement_threshold;
        this->max_distance_threshold = max_distance_threshold;
      }

      double DifferentialEvolution::optimize(std::vector<double>& xopt) {
        tools::printer.printStatusBegin("Optimizing (differential evolution)...");

        size_t d = f->getDimension();

        // vector of individuals
        std::vector<std::vector<double> > x1(population_size, std::vector<double>(d, 0.0));
        // another vector for the new population
        std::vector<std::vector<double> > x2 = x1;

        // pointers for swapping both vectors at the end of each iterations
        std::vector<std::vector<double> >* x_old = &x1;
        std::vector<std::vector<double> >* x_new = &x2;

        // function values at the points of the populations (no need to swape those)
        std::vector<double> fx(population_size, 0.0);

        // initial pseudorandom points
        for (size_t i = 0; i < population_size; i++) {
          for (size_t t = 0; t < d; t++) {
            (*x_old)[i][t] = sg::opt::tools::rng.getUniformRN();
          }

          fx[i] = f->eval((*x_old)[i]);
        }

        // smallest function value in the population
        double fopt = INFINITY;
        // index of the point with value fopt
        size_t xopt_index = 0;
        // iteration number of the last iteration with significant improvement
        size_t last_nonidle_k = 0;
        // average of all function values
        double avg = 0.0;
        // average in the previous round
        double last_avg = 0.0;
        // number of iterations
        size_t max_k = std::max(static_cast<size_t>(2), N / population_size) - 1;

        std::vector<std::vector<size_t> > a(max_k, std::vector<size_t>(population_size, 0)),
            b = a, c = a, j = a;
        std::vector<std::vector<std::vector<double> > > prob(max_k,
            std::vector<std::vector<double> >(population_size, std::vector<double>(d, 0)));

        // pregenerate all pseudorandom numbers because the real algorithm is parallelized
        // (for comparability of results, and maybe the RNG isn't thread-safe)
        for (size_t k = 0; k < max_k; k++) {
          for (size_t i = 0; i < population_size; i++) {
            do {
              a[k][i] = sg::opt::tools::rng.getUniformIndexRN(population_size);
            } while (a[k][i] == i);

            do {
              b[k][i] = sg::opt::tools::rng.getUniformIndexRN(population_size);
            } while ((b[k][i] == i) || (b[k][i] == a[k][i]));

            do {
              c[k][i] = sg::opt::tools::rng.getUniformIndexRN(population_size);
            } while ((c[k][i] == i) || (c[k][i] == a[k][i]) || (c[k][i] == b[k][i]));

            j[k][i] = sg::opt::tools::rng.getUniformIndexRN(d);

            for (size_t t = 0; t < d; t++) {
              if (t != j[k][i]) {
                prob[k][i][t] = sg::opt::tools::rng.getUniformRN();
              }
            }
          }
        }

        // "real" algorithm loop
        for (size_t k = 0; k < max_k; k++) {
          // abbreviations
          const std::vector<size_t>& a_k = a[k], &b_k = b[k], &c_k = c[k];
          const std::vector<size_t>& j_k = j[k];
          const std::vector<std::vector<double> >& prob_k = prob[k];

          #pragma omp parallel shared(k, a_k, b_k, c_k, j_k, prob_k, \
          x_old, d, fx, fopt, xopt_index, x_new) default(none)
          {
            std::vector<double> y(d, 0.0);
            tools::SmartPointer<function::Objective> cur_f(f->clone());

            // for each point in the population
            #pragma omp for schedule(dynamic)
            for (size_t i = 0; i < population_size; i++) {
              const size_t& cur_a = a_k[i], &cur_b = b_k[i], &cur_c = c_k[i];
              const size_t& cur_j = j_k[i];
              const std::vector<double>& prob_ki = prob_k[i];
              bool in_domain = true;

              // for each dimension
              for (size_t t = 0; t < d; t++) {
                const double& cur_prob = prob_ki[t];

                if ((t == cur_j) || (cur_prob < crossover_probability)) {
                  // mutate point in this dimension
                  y[t] = (*x_old)[cur_a][t] +
                         scaling_factor * ((*x_old)[cur_b][t] - (*x_old)[cur_c][t]);
                } else {
                  // don't mutate point in this dimension
                  y[t] = (*x_old)[i][t];
                }

                // mutated point is out of bounds ==> discard
                if ((y[t] < 0.0) || (y[t] > 1.0)) {
                  in_domain = false;
                  break;
                }
              }

              // evaluate mutated point (if not out of bounds)
              double fy = (in_domain ? cur_f->eval(y) : INFINITY);

              if (fy < fx[i]) {
                // function_value is better ==> replace point with mutated one
                #pragma omp critical
                {
                  fx[i] = fy;

                  if (fy < fopt) {
                    xopt_index = i;
                    fopt = fy;
                  }
                }

                for (size_t t = 0; t < d; t++) {
                  (*x_new)[i][t] = y[t];
                }
              } else {
                // function value not better ==> keep old point
                for (size_t t = 0; t < d; t++) {
                  (*x_new)[i][t] = (*x_old)[i][t];
                }
              }
            }
          }

          // swap populations
          std::swap(x_old, x_new);
          avg = 0.0;

          // calculate average function value
          for (size_t i = 0; i < population_size; i++) {
          avg += fx[i];
          }

          avg /= static_cast<double>(population_size);

          if (last_avg - avg >= avg_improvement_threshold) {
          // significant improvement
          last_nonidle_k = k;
        } else if (k - last_nonidle_k >= idle_generations_count) {
          // last significant improvement too long ago
          // ==> calculate maximum squared distance of all points to the best one
          double max_distance2 = 0.0;

          for (size_t i = 0; i < population_size; i++) {
              if (i == xopt_index) {
                continue;
              }

              double distance2 = 0.0;

              for (size_t t = 0; t < d; t++) {
                distance2 += ((*x_old)[i][t] - (*x_old)[xopt_index][t]) *
                             ((*x_old)[i][t] - (*x_old)[xopt_index][t]);
              }

              if (distance2 > max_distance2) {
                max_distance2 = distance2;
              }
            }

            // stopping criterion
            if (std::sqrt(max_distance2) < max_distance_threshold) {
              break;
            }
          }

          // save average in last_avg
          last_avg = avg;

          // status message
          if (k % 10 == 0) {
          tools::printer.printStatusUpdate(toString(k) + " steps, f(x) = " + toString(fopt));
          }
        }

        // optimal point
        xopt = (*x_old)[xopt_index];

        tools::printer.printStatusUpdate(toString(max_k) + " steps, f(x) = " + toString(fopt));
        tools::printer.printStatusEnd();

        return fopt;
      }

      tools::SmartPointer<Optimizer> DifferentialEvolution::clone() {
        tools::SmartPointer<Optimizer> result(new DifferentialEvolution(
                                                *f, N, population_size, crossover_probability, scaling_factor,
                                                idle_generations_count, avg_improvement_threshold, max_distance_threshold));
        result->setStartingPoint(x0);
        return result;
      }

      size_t DifferentialEvolution::getPopulationSize() const {
        return population_size;
      }

      void DifferentialEvolution::setPopulationSize(size_t population_size) {
        this->population_size = population_size;
      }

    }
  }
}
