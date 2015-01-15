/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include <cstring>

#include "opt/gridgen/IterativeGridGeneratorRitterNovak.hpp"
#include "opt/gridgen/HashRefinementMultiple.hpp"
#include "opt/tools/Permuter.hpp"
#include "opt/tools/Printer.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"

namespace sg {
  namespace opt {
    namespace gridgen {

      const double IterativeGridGeneratorRitterNovak::DEFAULT_ALPHA = 0.85;

      /**
       * Fast and approximative version of std::pow.
       * Source: http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
       *
       * @param a     base
       * @param b     exponent
       * @return      approximation to \f$a^b\f$
       */
      inline double fastPow(double a, double b) {
        union {
          double d;
          int x[2];
        } u = {a};

        u.x[1] = static_cast<int>(b * (u.x[1] - 1072632447) + 1072632447);
        u.x[0] = 0;

        return u.d;
      }

      IterativeGridGeneratorRitterNovak::IterativeGridGeneratorRitterNovak(
        function::Objective& f, base::Grid& grid, size_t N,
        double alpha, size_t max_level, PowMethod pow_method) :
        IterativeGridGenerator(f, grid, N),
        alpha(alpha),
        max_level(max_level),
        pow_method(pow_method) {
      }

      double IterativeGridGeneratorRitterNovak::getAlpha() const {
        return alpha;
      }

      void IterativeGridGeneratorRitterNovak::setAlpha(double alpha) {
        this->alpha = alpha;
      }

      size_t IterativeGridGeneratorRitterNovak::getMaxLevel() const {
        return max_level;
      }

      void IterativeGridGeneratorRitterNovak::setMaxLevel(size_t max_level) {
        this->max_level = max_level;
      }

      bool IterativeGridGeneratorRitterNovak::generate() {
        tools::printer.printStatusBegin("Adaptive grid generation (Ritter-Novak)...");

        bool result = true;
        base::GridIndex::PointDistribution distr = base::GridIndex::Normal;
        base::GridStorage* grid_storage = grid.getStorage();
        const size_t d = f.getDimension();

        HashRefinementMultiple refinement;

        if ((std::strcmp(grid.getType(), "BsplineClenshawCurtis") == 0) ||
            (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0)) {
          // Clenshaw-Curtis grid
          distr = base::GridIndex::ClenshawCurtis;
        }

        // generate initial grid
        {
          tools::SmartPointer<base::GridGenerator> grid_gen(grid.createGridGenerator());
	  // TODO: was changed for SGaA 2014 from 3 to 1
          grid_gen->regular(1);
        }

        size_t current_N = grid_storage->size();

        // abbreviation (function_values is a member variable of IterativeGridGenerator)
        std::vector<double>& fX = function_values;
        // fX_sorted is fX sorted ascendingly
        std::vector<double> fX_sorted(current_N, 0);
        // fX_order fulfills fX_sorted[i] = fX[fX_order[i]]
        std::vector<size_t> fX_order(current_N, 0);

        fX.assign(std::max(N, current_N), 0.0);
        // degree[i] is the number of times the i-th grid point was chosen for refinement
        std::vector<size_t> degree(fX.size(), 0);
        // level_sum[i] is the 1-norm of the level vector of the i-th grid point
        std::vector<size_t> level_sum(fX.size(), 0);
        // rank fulfills rank[i] = #{j | fX[j] <= fX[i]}
        std::vector<size_t> rank(fX.size(), 0);
        // for those grid points with ignore[i] == true the refinement criterion won't be evaluated
        std::vector<bool> ignore(fX.size(), false);

        // refinement_alpha will be a standard basis vector
        // (e.g. refinement_alpha[i] == 1.0 for exactly one i) to refine exactly one grid point
        base::DataVector refinement_alpha(current_N);
        refinement_alpha.setAll(0.0);

        for (size_t i = 0; i < current_N; i++) {
          base::GridIndex* gp = grid_storage->get(i);
          gp->setPointDistribution(distr);
          // prepare fX_order and rank
          fX_order[i] = i;
          rank[i] = i;

          // calculate sum of levels
          for (size_t t = 0; t < d; t++) {
            level_sum[i] += gp->getLevel(t);
          }
        }

        // parallel evaluation of f in the initial grid points
        #pragma omp parallel shared(fX, current_N, grid_storage) default(none)
        {
          std::vector<double> x(d, 0.0);
          function::Objective* cur_f_ptr = &f;
#ifdef _OPENMP
          tools::SmartPointer<function::Objective> cur_f;
          if (omp_get_max_threads() > 1) {
            cur_f = tools::SmartPointer<function::Objective>(f.clone());
            cur_f_ptr = cur_f.get();
          }
#endif

          #pragma omp for
          for (size_t i = 0; i < current_N; i++) {
            // convert grid point to coordinate vector
            #pragma omp critical
            {
              base::GridIndex* gp = grid_storage->get(i);

              for (size_t t = 0; t < d; t++) {
                x[t] = gp->abs(t);
              }
            }

            double fx = cur_f_ptr->eval(x);

            #pragma omp critical
            fX[i] = fx;
          }
        }

        // determine fX_order and rank (prepared above)
        // (C++11 would eliminate the need of tools::Permuter.)
        {
          tools::Permuter<double> permuter1(fX);
          std::sort(fX_order.begin(), fX_order.end(), permuter1);

          tools::Permuter<size_t> permuter2(fX_order);
          std::sort(rank.begin(), rank.end(), permuter2);
        }

        // determine fX_sorted
        for (size_t i = 0; i < current_N; i++) {
          fX_sorted[i] = fX[fX_order[i]];
        }

        // iteration counter
        size_t k = 0;

        while (current_N < N) {
          // status printing
          {
            char str[10];
            snprintf(str, 10, "%.1f%%",
                     static_cast<double>(current_N) / static_cast<double>(N) * 100.0);
            tools::printer.printStatusUpdate(std::string(str) +
                                             " (N = " + toString(current_N) +
                                             ", k = " + toString(k) + ")");
          }

          // determine the best i (i.e. i_best = argmin_i g_i)
          size_t i_best = 0;
          double g_best = INFINITY;

          for (size_t i = 0; i < current_N; i++) {
            if (ignore[i]) {
              continue;
            }

            // refinement criterion
            double g;

            if (pow_method == STD_POW) {
              g = std::pow(static_cast<double>(level_sum[i] + degree[i]) + 1.0, alpha) *
                  std::pow(static_cast<double>(rank[i]) + 1.0, 1.0 - alpha);
            } else {
              g = fastPow(static_cast<double>(level_sum[i] + degree[i]) + 1.0, alpha) *
                  fastPow(static_cast<double>(rank[i]) + 1.0, 1.0 - alpha);
            }

            if (g < g_best) {
              // so far the best value
              // ==> check if a refinement of this point would generate children with
              // a level greater than max_level (in one coordinate), if yes ignore the point
              base::GridIndex* gp = grid_storage->get(i);

              {
                base::GridIndex::index_type source_index, child_index;
                base::GridIndex::level_type source_level, child_level;

                // for each dimension
                for (size_t t = 0; t < d; t++) {
                  gp->get(t, source_level, source_index);

                  // inspect the left child to be generated
                  if ((source_level > 0) || (source_index == 1)) {
                    child_index = source_index;
                    child_level = source_level;

                    while (grid_storage->has_key(gp)) {
                      child_index *= 2;
                      child_level++;
                      gp->set(t, child_level, child_index - 1);
                    }

                    gp->set(t, source_level, source_index);

                    if (child_level > max_level) {
                      ignore[i] = true;
                      break;
                    }
                  }

                  // inspect the right child to be generated
                  if ((source_level > 0) || (source_index == 0)) {
                    child_index = source_index;
                    child_level = source_level;

                    while (grid_storage->has_key(gp)) {
                      child_index *= 2;
                      child_level++;
                      gp->set(t, child_level, child_index + 1);
                    }

                    gp->set(t, source_level, source_index);

                    if (child_level > max_level) {
                      ignore[i] = true;
                      break;
                    }
                  }
                }
              }

              // children were too "deep" ==> ignore the point
              if (ignore[i]) {
                continue;
              }

              // no ignore ==> new candidate for the point to be refined
              i_best = i;
              g_best = g;
            }
          }

          // refine point no. i_best
          degree[i_best]++;
          refinement_alpha[i_best] = 1.0;
          base::SurplusRefinementFunctor refine_func(&refinement_alpha, 1);
          refinement.free_refine(grid_storage, &refine_func);

          // new grid size
          size_t new_N = grid_storage->size();

          if (new_N == current_N) {
            // size unchanged ==> point not refined (should not happen)
            tools::printer.printStatusEnd(
              "error: size unchanged in IterativeGridGeneratorRitterNovak");
            result = false;
            break;
          }

          if (new_N > N) {
            // too many new points ==> undo refinement and exit
            std::list<size_t> indices_to_remove;

            for (size_t i = current_N; i < new_N; i++) {
              indices_to_remove.push_back(i);
            }

            grid_storage->deletePoints(indices_to_remove);
            break;
          }

          // resize refinement vector and set to all zeros (in the following loop)
          refinement_alpha.resize(new_N);
          refinement_alpha[i_best] = 0.0;

          for (size_t i = current_N; i < new_N; i++) {
            base::GridIndex* gp = grid_storage->get(i);
            // set point distribution accordingly to normal/Clenshaw-Curtis grids
            gp->setPointDistribution(distr);
            refinement_alpha[i] = 0.0;

            // calculate sum of levels
            for (size_t t = 0; t < d; t++) {
              level_sum[i] += gp->getLevel(t);
            }
          }

          // parallel evaluation of f in the new grid points
          #pragma omp parallel shared(fX, current_N, new_N, grid_storage) default(none)
          {
            std::vector<double> x(d, 0.0);
            function::Objective* cur_f_ptr = &f;
#ifdef _OPENMP
            tools::SmartPointer<function::Objective> cur_f;
            if (omp_get_max_threads() > 1) {
              cur_f = tools::SmartPointer<function::Objective>(f.clone());
              cur_f_ptr = cur_f.get();
            }
#endif

            #pragma omp for
            for (size_t i = current_N; i < new_N; i++) {
              // convert grid point to coordinate vector
              #pragma omp critical
              {
                base::GridIndex* gp = grid_storage->get(i);

                for (size_t t = 0; t < d; t++) {
                  x[t] = gp->abs(t);
                }
              }

              double fx = cur_f_ptr->eval(x);

              #pragma omp critical
              fX[i] = fx;
            }
          }

          for (size_t i = current_N; i < new_N; i++) {
            // update rank and fX_order by insertion sort
            // ==> go through fX from biggest entry to lowest
            for (size_t j = i - 1; j-- > 0; ) {
              if (fX_sorted[j] < fX[i]) {
                // new function value is bigger as current one ==> insert here
                fX_order.insert(fX_order.begin() + (j+1), i);
                fX_sorted.insert(fX_sorted.begin() + (j+1), fX[i]);
                rank[i] = j+1;
                break;
              } else {
                // new function value value not bigger yet ==> increase rank
                rank[fX_order[j]]++;
              }
            }

            if (fX_order.size() == i) {
              // this happens if the new function value is smaller than all of the previous ones
              // ==> insert at the beginning of fX_order, rank = 0
              fX_order.insert(fX_order.begin(), i);
              fX_sorted.insert(fX_sorted.begin(), fX[i]);
              rank[i] = 0;
            }
          }

          // next round
          current_N = new_N;
          k++;
        }

        // delete superfluous entries in fX
        fX.erase(fX.begin() + current_N, fX.end());

        if (result) {
          tools::printer.printStatusUpdate("100.0% (N = " + toString(current_N) +
                                           ", k = " + toString(k) + ")");
          tools::printer.printStatusEnd();
          return true;
        } else {
          return false;
        }
      }

    }
  }
}
