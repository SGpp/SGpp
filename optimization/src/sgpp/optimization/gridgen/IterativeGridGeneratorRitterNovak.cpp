// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <cstring>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/gridgen/HashRefinementMultiple.hpp>
#include <sgpp/optimization/tools/Permuter.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Fast and approximative version of std::pow.
     * Source: http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
     *
     * @param a     base
     * @param b     exponent
     * @return      approximation to \f$a^b\f$
     */
    inline float_t fastPow(float_t a, float_t b) {
      union {
        float_t d;
        int x[2];
      } u = {a};

      u.x[1] = static_cast<int>(b * (u.x[1] - 1072632447) + 1072632447);
      u.x[0] = 0;

      return u.d;
    }

    IterativeGridGeneratorRitterNovak::IterativeGridGeneratorRitterNovak(
      ObjectiveFunction& f, base::Grid& grid, size_t N,
      float_t gamma, size_t maxLevel, PowMethod powMethod) :
      IterativeGridGenerator(f, grid, N),
      gamma(gamma),
      maxLevel(maxLevel),
      powMethod(powMethod) {
    }

    float_t IterativeGridGeneratorRitterNovak::getGamma() const {
      return gamma;
    }

    void IterativeGridGeneratorRitterNovak::setGamma(float_t gamma) {
      this->gamma = gamma;
    }

    size_t IterativeGridGeneratorRitterNovak::getMaxLevel() const {
      return maxLevel;
    }

    void IterativeGridGeneratorRitterNovak::setMaxLevel(size_t maxLevel) {
      this->maxLevel = maxLevel;
    }

    bool IterativeGridGeneratorRitterNovak::generate() {
      printer.printStatusBegin("Adaptive grid generation (Ritter-Novak)...");

      bool result = true;
      base::GridIndex::PointDistribution distr = base::GridIndex::Normal;
      base::GridStorage& gridStorage = *grid.getStorage();
      const size_t d = f.getDimension();

      HashRefinementMultiple refinement;

      if ((std::strcmp(grid.getType(), "bsplineClenshawCurtis") == 0) ||
          (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0)) {
        // Clenshaw-Curtis grid
        distr = base::GridIndex::ClenshawCurtis;
      }

      // generate initial grid
      {
        std::unique_ptr<base::GridGenerator> gridGen(
          grid.createGridGenerator());
        gridGen->regular(3);
      }

      size_t currentN = gridStorage.size();

      // abbreviation (functionValues is a member variable of
      // IterativeGridGenerator)
      std::vector<float_t>& fX = functionValues;
      // fX_sorted is fX sorted ascendingly
      std::vector<float_t> fXSorted(currentN, 0);
      // fX_order fulfills fX_sorted[i] = fX[fX_order[i]]
      std::vector<size_t> fXOrder(currentN, 0);

      fX.assign(std::max(N, currentN), 0.0);
      // degree[i] is the number of times the i-th grid point was
      // chosen for refinement
      std::vector<size_t> degree(fX.size(), 0);
      // level_sum[i] is the 1-norm of the level vector of the i-th grid point
      std::vector<size_t> levelSum(fX.size(), 0);
      // rank fulfills rank[i] = #{j | fX[j] <= fX[i]}
      std::vector<size_t> rank(fX.size(), 0);
      // for those grid points with ignore[i] == true the refinement
      // criterion won't be evaluated
      std::vector<bool> ignore(fX.size(), false);

      // refinementAlpha will be a standard basis vector
      // (e.g. refinementAlpha[i] == 1.0 for exactly one i) to refine
      // exactly one grid point
      base::DataVector refinementAlpha(currentN);
      refinementAlpha.setAll(0.0);

      for (size_t i = 0; i < currentN; i++) {
        base::GridIndex& gp = *gridStorage.get(i);
        gp.setPointDistribution(distr);
        // prepare fX_order and rank
        fXOrder[i] = i;
        rank[i] = i;

        // calculate sum of levels
        for (size_t t = 0; t < d; t++) {
          levelSum[i] += gp.getLevel(t);
        }
      }

      // parallel evaluation of f in the initial grid points
      #pragma omp parallel shared(fX, currentN, gridStorage) default(none)
      {
        std::vector<float_t> x(d, 0.0);
        ObjectiveFunction* curFPtr = &f;
#ifdef _OPENMP
        std::unique_ptr<ObjectiveFunction> curF;

        if (omp_get_max_threads() > 1) {
          f.clone(curF);
          curFPtr = curF.get();
        }

#endif /* _OPENMP */

        #pragma omp for

        for (size_t i = 0; i < currentN; i++) {
          // convert grid point to coordinate vector
          #pragma omp critical
          {
            base::GridIndex& gp = *gridStorage.get(i);

            for (size_t t = 0; t < d; t++) {
              x[t] = gp.getCoord(t);
            }
          }

          float_t fx = curFPtr->eval(x);

          #pragma omp critical
          fX[i] = fx;
        }
      }

      // determine fX_order and rank (prepared above)
      // (C++11 would eliminate the need of Permuter.)
      {
        Permuter<float_t> permuter1(fX);
        std::sort(fXOrder.begin(), fXOrder.end(), permuter1);

        Permuter<size_t> permuter2(fXOrder);
        std::sort(rank.begin(), rank.end(), permuter2);
      }

      // determine fX_sorted
      for (size_t i = 0; i < currentN; i++) {
        fXSorted[i] = fX[fXOrder[i]];
      }

      // iteration counter
      size_t k = 0;

      while (currentN < N) {
        // status printing
        {
          char str[10];
          snprintf(str, 10, "%.1f%%",
                   static_cast<float_t>(currentN) /
                   static_cast<float_t>(N) * 100.0);
          printer.printStatusUpdate(std::string(str) +
                                    " (N = " + std::to_string(currentN) +
                                    ", k = " + std::to_string(k) + ")");
        }

        // determine the best i (i.e. i_best = argmin_i g_i)
        size_t iBest = 0;
        float_t gBest = INFINITY;

        for (size_t i = 0; i < currentN; i++) {
          if (ignore[i]) {
            continue;
          }

          // refinement criterion
          float_t g;

          if (powMethod == STD_POW) {
            g = std::pow(static_cast<float_t>(levelSum[i] + degree[i]) + 1.0,
                         gamma) *
                std::pow(static_cast<float_t>(rank[i]) + 1.0, 1.0 - gamma);
          } else {
            g = fastPow(static_cast<float_t>(levelSum[i] + degree[i]) + 1.0,
                        gamma) *
                fastPow(static_cast<float_t>(rank[i]) + 1.0, 1.0 - gamma);
          }

          if (g < gBest) {
            // so far the best value
            // ==> check if a refinement of this point would generate
            // children with a level greater than max_level
            // (in one coordinate), if yes ignore the point
            base::GridIndex& gp = *gridStorage.get(i);

            {
              base::GridIndex::index_type sourceIndex, childIndex;
              base::GridIndex::level_type sourceLevel, childLevel;

              // for each dimension
              for (size_t t = 0; t < d; t++) {
                gp.get(t, sourceLevel, sourceIndex);

                // inspect the left child to be generated
                if ((sourceLevel > 0) || (sourceIndex == 1)) {
                  childIndex = sourceIndex;
                  childLevel = sourceLevel;

                  while (gridStorage.has_key(&gp)) {
                    childIndex *= 2;
                    childLevel++;
                    gp.set(t, childLevel, childIndex - 1);
                  }

                  gp.set(t, sourceLevel, sourceIndex);

                  if (childLevel > maxLevel) {
                    ignore[i] = true;
                    break;
                  }
                }

                // inspect the right child to be generated
                if ((sourceLevel > 0) || (sourceIndex == 0)) {
                  childIndex = sourceIndex;
                  childLevel = sourceLevel;

                  while (gridStorage.has_key(&gp)) {
                    childIndex *= 2;
                    childLevel++;
                    gp.set(t, childLevel, childIndex + 1);
                  }

                  gp.set(t, sourceLevel, sourceIndex);

                  if (childLevel > maxLevel) {
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
            iBest = i;
            gBest = g;
          }
        }

        // refine point no. i_best
        degree[iBest]++;
        refinementAlpha[iBest] = 1.0;
        base::SurplusRefinementFunctor refineFunc(&refinementAlpha, 1);
        refinement.free_refine(&gridStorage, &refineFunc);

        // new grid size
        size_t newN = gridStorage.size();

        if (newN == currentN) {
          // size unchanged ==> point not refined (should not happen)
          printer.printStatusEnd(
            "error: size unchanged in IterativeGridGeneratorRitterNovak");
          result = false;
          break;
        }

        if (newN > N) {
          // too many new points ==> undo refinement and exit
          std::list<size_t> indicesToRemove;

          for (size_t i = currentN; i < newN; i++) {
            indicesToRemove.push_back(i);
          }

          gridStorage.deletePoints(indicesToRemove);
          break;
        }

        // resize refinement vector and set to all zeros
        // (in the following loop)
        refinementAlpha.resize(newN);
        refinementAlpha[iBest] = 0.0;

        for (size_t i = currentN; i < newN; i++) {
          base::GridIndex& gp = *gridStorage.get(i);
          // set point distribution accordingly to normal/Clenshaw-Curtis grids
          gp.setPointDistribution(distr);
          refinementAlpha[i] = 0.0;

          // calculate sum of levels
          for (size_t t = 0; t < d; t++) {
            levelSum[i] += gp.getLevel(t);
          }
        }

        // parallel evaluation of f in the new grid points
        #pragma omp parallel shared(fX, currentN, newN, gridStorage) \
        default(none)
        {
          std::vector<float_t> x(d, 0.0);
          ObjectiveFunction* curFPtr = &f;
#ifdef _OPENMP
          std::unique_ptr<ObjectiveFunction> curF;

          if (omp_get_max_threads() > 1) {
            f.clone(curF);
            curFPtr = curF.get();
          }

#endif /* _OPENMP */

          #pragma omp for

          for (size_t i = currentN; i < newN; i++) {
            // convert grid point to coordinate vector
            #pragma omp critical
            {
              base::GridIndex& gp = *gridStorage.get(i);

              for (size_t t = 0; t < d; t++) {
                x[t] = gp.getCoord(t);
              }
            }

            float_t fx = curFPtr->eval(x);

            #pragma omp critical
            fX[i] = fx;
          }
        }

        for (size_t i = currentN; i < newN; i++) {
          // update rank and fX_order by insertion sort
          // ==> go through fX from biggest entry to lowest
          for (size_t j = i - 1; j-- > 0; ) {
            if (fXSorted[j] < fX[i]) {
              // new function value is bigger as current one ==> insert here
              fXOrder.insert(fXOrder.begin() + (j + 1), i);
              fXSorted.insert(fXSorted.begin() + (j + 1), fX[i]);
              rank[i] = j + 1;
              break;
            } else {
              // new function value value not bigger yet ==> increase rank
              rank[fXOrder[j]]++;
            }
          }

          if (fXOrder.size() == i) {
            // this happens if the new function value is smaller than all
            // of the previous ones
            // ==> insert at the beginning of fX_order, rank = 0
            fXOrder.insert(fXOrder.begin(), i);
            fXSorted.insert(fXSorted.begin(), fX[i]);
            rank[i] = 0;
          }
        }

        // next round
        currentN = newN;
        k++;
      }

      // delete superfluous entries in fX
      fX.erase(fX.begin() + currentN, fX.end());

      if (result) {
        printer.printStatusUpdate("100.0% (N = " + std::to_string(currentN) +
                                  ", k = " + std::to_string(k) + ")");
        printer.printStatusEnd();
        return true;
      } else {
        return false;
      }
    }

  }
}
