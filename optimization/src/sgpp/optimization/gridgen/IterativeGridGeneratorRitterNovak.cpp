// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/gridgen/HashRefinementMultiple.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>

#include <algorithm>
#include <cstring>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

namespace sgpp {
namespace optimization {

namespace {
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

  u.x[1] = static_cast<int>(b * static_cast<double>(u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;

  return u.d;
}
}  // namespace

IterativeGridGeneratorRitterNovak::IterativeGridGeneratorRitterNovak(
    base::ScalarFunction& f, base::Grid& grid, size_t N, double adaptivity,
    base::level_t initialLevel, base::level_t maxLevel, PowMethod powMethod)
    : IterativeGridGenerator(f, grid, N),
      gamma(adaptivity),
      initialLevel(initialLevel),
      maxLevel(maxLevel),
      powMethod(powMethod) {}

IterativeGridGeneratorRitterNovak::~IterativeGridGeneratorRitterNovak() {}

double IterativeGridGeneratorRitterNovak::getAdaptivity() const { return gamma; }

void IterativeGridGeneratorRitterNovak::setAdaptivity(double adaptivity) {
  this->gamma = adaptivity;
}

base::level_t IterativeGridGeneratorRitterNovak::getInitialLevel() const { return initialLevel; }

void IterativeGridGeneratorRitterNovak::setInitialLevel(base::level_t initialLevel) {
  this->initialLevel = initialLevel;
}

base::level_t IterativeGridGeneratorRitterNovak::getMaxLevel() const { return maxLevel; }

void IterativeGridGeneratorRitterNovak::setMaxLevel(base::level_t maxLevel) {
  this->maxLevel = maxLevel;
}

IterativeGridGeneratorRitterNovak::PowMethod IterativeGridGeneratorRitterNovak::getPowMethod()
    const {
  return powMethod;
}

void IterativeGridGeneratorRitterNovak::setPowMethod(
    IterativeGridGeneratorRitterNovak::PowMethod powMethod) {
  this->powMethod = powMethod;
}

bool IterativeGridGeneratorRitterNovak::generate() {
  sgpp::base::Printer::getInstance().printStatusBegin("Adaptive grid generation (Ritter-Novak)...");

  bool result = true;
  base::GridStorage& gridStorage = grid.getStorage();
  const size_t d = f.getNumberOfParameters();

  HashRefinementMultiple refinement;

  // generate initial grid
  grid.getGenerator().regular(initialLevel);

  size_t currentN = gridStorage.getSize();

  // abbreviation (functionValues is a member variable of
  // IterativeGridGenerator)
  base::DataVector& fX = functionValues;
  // fXSorted is fX sorted ascendingly
  base::DataVector fXSorted(currentN);
  // fXOrder fulfills fXSorted[i] = fX[fXOrder[i]]
  std::vector<size_t> fXOrder(currentN);

  fX.resize(std::max(N, currentN));
  fX.setAll(0.0);
  // degree[i] is the number of times the i-th grid point was
  // chosen for refinement
  std::vector<size_t> degree(fX.getSize(), 0);
  // level_sum[i] is the 1-norm of the level vector of the i-th grid point
  std::vector<size_t> levelSum(fX.getSize(), 0);
  // rank fulfills rank[i] = #{j | fX[j] <= fX[i]}
  std::vector<size_t> rank(fX.getSize(), 0);
  // for those grid points with ignore[i] == true the refinement
  // criterion won't be evaluated
  std::vector<bool> ignore(fX.getSize(), false);

  // refinementAlpha will be a standard basis vector
  // (e.g. refinementAlpha[i] == 1.0 for exactly one i) to refine
  // exactly one grid point
  base::DataVector refinementAlpha(currentN);

  for (size_t i = 0; i < currentN; i++) {
    base::GridPoint& gp = gridStorage[i];
    // prepare fXOrder and rank
    fXOrder[i] = i;
    rank[i] = i + 1;

    // calculate sum of levels
    for (size_t t = 0; t < d; t++) {
      levelSum[i] += gp.getLevel(t);
    }
  }

  // evaluation of f in the initial grid points
  evalFunction();

  // determine fXOrder and rank (prepared above)
  std::sort(fXOrder.begin(), fXOrder.begin() + currentN,
            [&fX](size_t a, size_t b) { return (fX[a] < fX[b]); });
  std::sort(rank.begin(), rank.begin() + currentN,
            [&fXOrder](size_t a, size_t b) { return (fXOrder[a - 1] < fXOrder[b - 1]); });

  // determine fXSorted
  for (size_t i = 0; i < currentN; i++) {
    fXSorted[i] = fX[fXOrder[i]];
  }

  // iteration counter
  size_t k = 0;

  while (currentN < N) {
    // status printing
    {
      char str[10];
      snprintf(str, sizeof(str), "%.1f%%",
               static_cast<double>(currentN) / static_cast<double>(N) * 100.0);
      sgpp::base::Printer::getInstance().printStatusUpdate(std::string(str) +
                                                           " (N = " + std::to_string(currentN) +
                                                           ", k = " + std::to_string(k) + ")");
    }

    // determine the best i (i.e., iBest = argmin_i g_i)
    size_t iBest = 0;
    double gBest = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < currentN; i++) {
      if (ignore[i]) {
        continue;
      }

      // refinement criterion
      double g;

      if (powMethod == PowMethod::STD_POW) {
        g = std::pow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
            std::pow(static_cast<double>(rank[i]) + 1.0, 1.0 - gamma);
      } else {
        g = fastPow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
            fastPow(static_cast<double>(rank[i]) + 1.0, 1.0 - gamma);
      }

      if (g < gBest) {
        // so far the best value
        // ==> check if a refinement of this point would generate
        // children with a level greater than max_level
        // (in one coordinate), if yes ignore the point
        base::GridPoint& gp = gridStorage[i];

        {
          base::index_t sourceIndex, childIndex;
          base::level_t sourceLevel, childLevel;

          // for each dimension
          for (size_t t = 0; t < d; t++) {
            gp.get(t, sourceLevel, sourceIndex);

            // inspect the left child to be generated
            if ((sourceLevel > 0) || (sourceIndex == 1)) {
              childIndex = sourceIndex;
              childLevel = sourceLevel;

              while (gridStorage.isContaining(gp)) {
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

              while (gridStorage.isContaining(gp)) {
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

    // refine point no. iBest
    degree[iBest]++;
    refinementAlpha[iBest] = 1.0;
    base::SurplusRefinementFunctor refineFunc(refinementAlpha, 1);
    refinement.free_refine(gridStorage, refineFunc);
    // new grid size
    const size_t newN = gridStorage.getSize();

    if (newN == currentN) {
      // size unchanged ==> point not refined (should not happen)
      sgpp::base::Printer::getInstance().printStatusEnd(
          "error: size unchanged in IterativeGridGeneratorRitterNovak");
      result = false;
      break;
    }

    if (newN > N) {
      // too many new points ==> undo refinement and exit
      undoRefinement(currentN);
      break;
    }

    // resize refinement vector and set to all zeros
    // (in the following loop)
    refinementAlpha.resize(newN);
    refinementAlpha[iBest] = 0.0;

    for (size_t i = currentN; i < newN; i++) {
      base::GridPoint& gp = gridStorage[i];
      refinementAlpha[i] = 0.0;

      // calculate sum of levels
      for (size_t t = 0; t < d; t++) {
        levelSum[i] += gp.getLevel(t);
      }
    }

    // evaluation of f in the new grid points
    evalFunction(currentN);

    for (size_t i = currentN; i < newN; i++) {
      const double fXi = fX[i];

      // update rank and fXOrder by insertion sort
      // ==> go through fX from greatest entry to lowest
      for (size_t j = i; j-- > 0;) {
        if (fXSorted[j] < fXi) {
          // new function value is greater than current one ==> insert here
          fXOrder.insert(fXOrder.begin() + (j + 1), i);
          fXSorted.insert(j + 1, fXi);
          rank[i] = j + 1;
          break;
        } else {
          // new function value value not greater yet ==> increase rank
          rank[fXOrder[j]]++;
        }
      }

      if (fXOrder.size() == i) {
        // this happens if the new function value is smaller than all
        // of the previous ones
        // ==> insert at the beginning of fXOrder, rank = 0
        fXOrder.insert(fXOrder.begin(), i);
        fXSorted.insert(0, fXi);
        rank[i] = 0;
      }
    }

    // next round
    currentN = newN;
    k++;
  }

  // delete superfluous entries in fX
  fX.resize(currentN);

  if (result) {
    sgpp::base::Printer::getInstance().printStatusUpdate("100.0% (N = " + std::to_string(currentN) +
                                                         ", k = " + std::to_string(k) + ")");
    sgpp::base::Printer::getInstance().printStatusEnd();
    return true;
  } else {
    return false;
  }
}
}  // namespace optimization
}  // namespace sgpp
