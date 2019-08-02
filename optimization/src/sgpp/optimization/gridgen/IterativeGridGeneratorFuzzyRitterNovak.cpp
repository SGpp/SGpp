// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/gridgen/IterativeGridGeneratorFuzzyRitterNovak.hpp>
#include <sgpp/optimization/gridgen/HashRefinementMultiple.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

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

IterativeGridGeneratorFuzzyRitterNovak::IterativeGridGeneratorFuzzyRitterNovak(
    ScalarFunction& f, base::Grid& grid, size_t N,
    const std::vector<const FuzzyInterval*>& xFuzzy, size_t numberOfAlphaSegments,
    double adaptivity, base::level_t initialLevel,
    base::level_t maxLevel, PowMethod powMethod)
    : IterativeGridGeneratorRitterNovak(f, grid, N, adaptivity, initialLevel, maxLevel, powMethod),
      xFuzzy(xFuzzy),
      numberOfAlphaSegments(numberOfAlphaSegments) {}

IterativeGridGeneratorFuzzyRitterNovak::~IterativeGridGeneratorFuzzyRitterNovak() {}

const std::vector<const FuzzyInterval*>&
IterativeGridGeneratorFuzzyRitterNovak::getXFuzzy() const {
  return xFuzzy;
}

void IterativeGridGeneratorFuzzyRitterNovak::setXFuzzy(
    const std::vector<const FuzzyInterval*>& xFuzzy) {
  this->xFuzzy = xFuzzy;
}

size_t IterativeGridGeneratorFuzzyRitterNovak::getNumberOfAlphaSegments() const {
  return numberOfAlphaSegments;
}

void IterativeGridGeneratorFuzzyRitterNovak::setNumberOfAlphaSegments(
    size_t numberOfAlphaSegments) {
  this->numberOfAlphaSegments = numberOfAlphaSegments;
}

bool IterativeGridGeneratorFuzzyRitterNovak::generate() {
  Printer::getInstance().printStatusBegin("Adaptive grid generation (fuzzy Ritter-Novak)...");

  bool result = true;
  base::GridStorage& gridStorage = grid.getStorage();
  const size_t d = f.getNumberOfParameters();

  HashRefinementMultiple refinement;

  grid.getGenerator().regular(initialLevel);

  size_t currentN = gridStorage.getSize();
  base::DataVector& fX = functionValues;
  fX.resize(std::max(N, currentN));
  fX.setAll(0.0);
  std::vector<size_t> degree(fX.getSize(), 0);
  std::vector<size_t> levelSum(fX.getSize(), 0);
  std::vector<bool> ignore(fX.getSize(), false);

  base::DataVector refinementAlpha(currentN);
  refinementAlpha.setAll(0.0);

  for (size_t i = 0; i < currentN; i++) {
    base::GridPoint& gp = gridStorage[i];

    for (size_t t = 0; t < d; t++) {
      levelSum[i] += gp.getLevel(t);
    }
  }

  evalFunction();

  const size_t numberOfBounds = numberOfAlphaSegments + 2;
  base::DataMatrix lowerBounds(numberOfBounds, d);
  base::DataMatrix upperBounds(numberOfBounds, d);

  for (size_t t = 0; t < d; t++) {
    lowerBounds(0, t) = 0.0;
    upperBounds(0, t) = 1.0;
  }

  for (size_t j = 0; j <= numberOfAlphaSegments; j++) {
    const double alphaLevel = static_cast<double>(j) / static_cast<double>(numberOfAlphaSegments);

    for (size_t t = 0; t < d; t++) {
      double lowerBound = xFuzzy[t]->evaluateConfidenceIntervalLowerBound(alphaLevel);
      double upperBound = xFuzzy[t]->evaluateConfidenceIntervalUpperBound(alphaLevel);
      double boundsSize = upperBound - lowerBound;
      const double minimumBoundsSize = 0.05;

      if (boundsSize < minimumBoundsSize) {
        const double boundsCenter = (lowerBound + upperBound) / 2.0;
        lowerBound = boundsCenter - minimumBoundsSize / 2.0;
        upperBound = boundsCenter + minimumBoundsSize / 2.0;
        boundsSize = minimumBoundsSize;
      }

      lowerBound = std::max(lowerBound - 0.05 * boundsSize, 0.0);
      upperBound = std::min(upperBound + 0.05 * boundsSize, 1.0);

      lowerBounds(j + 1, t) = lowerBound;
      upperBounds(j + 1, t) = upperBound;
    }
  }

  size_t k = 0;
  base::DataVector x(d);

  while (currentN < N) {
    {
      char str[10];
      snprintf(str, sizeof(str),
               "%.1f%%", static_cast<double>(currentN) / static_cast<double>(N) * 100.0);
      Printer::getInstance().printStatusUpdate(std::string(str) + " (N = " +
                                               std::to_string(currentN) + ", k = " +
                                               std::to_string(k) + ")");
    }

    std::vector<size_t> pointsToBeRefined;
    std::vector<bool> pointAlreadyChecked(currentN, false);

    for (size_t j = 0; j < numberOfBounds; j++) {
      std::vector<size_t> pointsRestricted;

      for (size_t i = 0; i < currentN; i++) {
        if (ignore[i]) {
          continue;
        }

        gridStorage[i].getStandardCoordinates(x);
        bool inBounds = true;

        for (size_t t = 0; t < d; t++) {
          if ((x[t] < lowerBounds(j, t)) || (x[t] > upperBounds(j, t))) {
            inBounds = false;
            break;
          }
        }

        if (inBounds) {
          pointsRestricted.push_back(i);
        }
      }

      if (pointsRestricted.empty()) {
        continue;
      }

      const size_t NRestricted = pointsRestricted.size();
      base::DataVector fXRestricted(NRestricted);
      std::vector<size_t> fXOrderRestricted(NRestricted);
      std::vector<size_t> rankRestricted(NRestricted);

      for (size_t k = 0; k < NRestricted; k++) {
        fXRestricted[k] = fX[pointsRestricted[k]];
        fXOrderRestricted[k] = k;
        rankRestricted[k] = k + 1;
      }

      std::sort(fXOrderRestricted.begin(), fXOrderRestricted.end(),
                [&fXRestricted](size_t a, size_t b) {
                  return (fXRestricted[a] < fXRestricted[b]);
                });
      std::sort(rankRestricted.begin(), rankRestricted.end(),
                [&fXOrderRestricted](size_t a, size_t b) {
                  return (fXOrderRestricted[a - 1] < fXOrderRestricted[b - 1]);
                });

      double gMinBest = std::numeric_limits<double>::infinity();
      size_t iMinBest;
      double gMaxBest = std::numeric_limits<double>::infinity();
      size_t iMaxBest;

      for (size_t k = 0; k < NRestricted; k++) {
        const size_t i = pointsRestricted[k];
        double gMin;
        double gMax;

        if (powMethod == STD_POW) {
          gMin = std::pow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                 std::pow(static_cast<double>(rankRestricted[k]) + 1.0, 1.0 - gamma);
          gMax = std::pow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                 std::pow(static_cast<double>(NRestricted - rankRestricted[k]) + 2.0, 1.0 - gamma);
        } else {
          gMin = fastPow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                 fastPow(static_cast<double>(rankRestricted[k]) + 1.0, 1.0 - gamma);
          gMax = fastPow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                 fastPow(static_cast<double>(NRestricted - rankRestricted[k]) + 2.0, 1.0 - gamma);
        }

        if ((gMin >= gMinBest) && (gMax >= gMaxBest)) {
          continue;
        }

        if (ignore[i]) {
          continue;
        }

        if (!pointAlreadyChecked[i]) {
          pointAlreadyChecked[i] = true;

          base::GridPoint& gp = gridStorage[i];
          base::index_t sourceIndex, childIndex;
          base::level_t sourceLevel, childLevel;

          for (size_t t = 0; t < d; t++) {
            gp.get(t, sourceLevel, sourceIndex);

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

        if (ignore[i]) {
          continue;
        }

        if (gMin < gMinBest) {
          gMinBest = gMin;
          iMinBest = i;
        }

        if (gMax < gMaxBest) {
          gMaxBest = gMax;
          iMaxBest = i;
        }
      }

      if (gMinBest < std::numeric_limits<double>::infinity()) {
        pointsToBeRefined.push_back(iMinBest);
      }

      if (gMaxBest < std::numeric_limits<double>::infinity()) {
        pointsToBeRefined.push_back(iMaxBest);
      }
    }

    std::sort(pointsToBeRefined.begin(), pointsToBeRefined.end());
    pointsToBeRefined.erase(std::unique(pointsToBeRefined.begin(), pointsToBeRefined.end()),
                            pointsToBeRefined.end());

    bool exit = false;

    for (size_t iBest : pointsToBeRefined) {
      degree[iBest]++;
      refinementAlpha[iBest] = 1.0;
      base::SurplusRefinementFunctor refineFunc(refinementAlpha, 1);
      refinement.free_refine(gridStorage, refineFunc);

      const size_t newN = gridStorage.getSize();

      if (newN == currentN) {
        Printer::getInstance().printStatusEnd(
            "error: size unchanged in IterativeGridGeneratorFuzzyRitterNovak");
        result = false;
        exit = true;
        break;
      }

      if (newN > N) {
        undoRefinement(currentN);
        exit = true;
        break;
      }

      refinementAlpha.resize(newN);
      refinementAlpha[iBest] = 0.0;

      for (size_t i = currentN; i < newN; i++) {
        base::GridPoint& gp = gridStorage[i];
        refinementAlpha[i] = 0.0;

        for (size_t t = 0; t < d; t++) {
          levelSum[i] += gp.getLevel(t);
        }
      }

      evalFunction(currentN);

      currentN = newN;
    }

    if (exit) {
      break;
    }

    k++;
  }

  fX.resize(currentN);

  if (result) {
    Printer::getInstance().printStatusUpdate(
        "100.0% (N = " + std::to_string(currentN) + ", k = " + std::to_string(k) + ")");
    Printer::getInstance().printStatusEnd();
    return true;
  } else {
    return false;
  }
}

/*bool IterativeGridGeneratorFuzzyRitterNovak::generate() {
  Printer::getInstance().printStatusBegin("Adaptive grid generation (fuzzy Ritter-Novak)...");

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
  refinementAlpha.setAll(0.0);

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

  const size_t numberOfBounds = numberOfAlphaSegments + 2;
  base::DataMatrix lowerBounds(numberOfBounds, d);
  base::DataMatrix upperBounds(numberOfBounds, d);

  for (size_t t = 0; t < d; t++) {
    lowerBounds(0, t) = 0.0;
    upperBounds(0, t) = 1.0;
  }

  for (size_t j = 0; j <= numberOfAlphaSegments; j++) {
    const double alphaLevel = static_cast<double>(j) / static_cast<double>(numberOfAlphaSegments);

    for (size_t t = 0; t < d; t++) {
      lowerBounds(j + 1, t) = xFuzzy[t]->evaluateConfidenceIntervalLowerBound(alphaLevel);
      upperBounds(j + 1, t) = xFuzzy[t]->evaluateConfidenceIntervalUpperBound(alphaLevel);
    }
  }

  base::DataVector x(d);
  std::vector<size_t> iBests(numberOfBounds);
  base::DataVector gBests(numberOfBounds);
  std::vector<size_t> refineIndices;

  // iteration counter
  size_t k = 0;

  bool exit;

  while (currentN < N) {
    // status printing
    {
      char str[10];
      snprintf(str, sizeof(str),
               "%.1f%%", static_cast<double>(currentN) / static_cast<double>(N) * 100.0);
      Printer::getInstance().printStatusUpdate(std::string(str) + " (N = " +
                                               std::to_string(currentN) + ", k = " +
                                               std::to_string(k) + ")");
    }

    for (size_t q = 0; q < 2; q++) {
      // determine the best i (i.e. i_best = argmin_i g_i)
      std::fill(iBests.begin(), iBests.end(), 0);
      gBests.setAll(INFINITY);

      for (size_t i = 0; i < currentN; i++) {
        if (ignore[i]) {
          continue;
        }

        // refinement criterion
        double g;

        if (q == 0) {
          if (powMethod == STD_POW) {
            g = std::pow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                std::pow(static_cast<double>(rank[i]) + 1.0, 1.0 - gamma);
          } else {
            g = fastPow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                fastPow(static_cast<double>(rank[i]) + 1.0, 1.0 - gamma);
          }
        } else {
          if (powMethod == STD_POW) {
            g = std::pow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                std::pow(static_cast<double>(currentN - rank[i]) + 2.0, 1.0 - gamma);
          } else {
            g = fastPow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                fastPow(static_cast<double>(currentN - rank[i]) + 2.0, 1.0 - gamma);
          }
        }

        base::GridPoint& gp = gridStorage[i];
        gp.getStandardCoordinates(x);

        bool alreadyCheckedForIgnore = false;

        for (size_t j = 0; j < numberOfBounds; j++) {
          bool inBounds = true;

          for (size_t t = 0; t < d; t++) {
            if ((x[t] < lowerBounds(j, t)) || (x[t] > upperBounds(j, t))) {
              inBounds = false;
              break;
            }
          }

          if (inBounds && (g < gBests[j])) {
            // so far the best value
            if (!alreadyCheckedForIgnore) {
              // ==> check if a refinement of this point would generate
              // children with a level greater than max_level
              // (in one coordinate), if yes ignore the point
              alreadyCheckedForIgnore = true;

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
              break;
            }

            // no ignore ==> new candidate for the point to be refined
            iBests[j] = i;
            gBests[j] = g;
          }
        }
      }

      refineIndices.clear();

      for (size_t j = 0; j < numberOfBounds; j++) {
        if (gBests[j] != INFINITY) {
          refineIndices.push_back(iBests[j]);
        }
      }

      std::sort(refineIndices.begin(), refineIndices.end());
      refineIndices.erase(std::unique(refineIndices.begin(), refineIndices.end()),
                          refineIndices.end());

      exit = false;

      for (size_t iBest : refineIndices) {
        // refine point no. iBest
        degree[iBest]++;
        refinementAlpha[iBest] = 1.0;
        base::SurplusRefinementFunctor refineFunc(refinementAlpha, 1);
        refinement.free_refine(gridStorage, refineFunc);

        // new grid size
        const size_t newN = gridStorage.getSize();

        if (newN == currentN) {
          // size unchanged ==> point not refined (should not happen)
          Printer::getInstance().printStatusEnd(
              "error: size unchanged in IterativeGridGeneratorRitterNovak");
          result = false;
          exit = true;
          break;
        }

        if (newN > N) {
          // too many new points ==> undo refinement and exit
          undoRefinement(currentN);
          exit = true;
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
      }

      if (exit) {
        break;
      }
    }

    if (exit) {
      break;
    }

    k++;
  }

  // delete superfluous entries in fX
  fX.resize(currentN);

  if (result) {
    Printer::getInstance().printStatusUpdate("100.0% (N = " + std::to_string(currentN) + ", k = " +
                                             std::to_string(k) + ")");
    Printer::getInstance().printStatusEnd();
    return true;
  } else {
    return false;
  }
}*/
}  // namespace optimization
}  // namespace sgpp
