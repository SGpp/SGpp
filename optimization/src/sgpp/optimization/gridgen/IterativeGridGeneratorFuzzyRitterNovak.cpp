// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorFuzzyRitterNovak.hpp>
#include <sgpp/optimization/gridgen/HashRefinementMultiple.hpp>

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
    base::ScalarFunction& f, base::Grid& grid, size_t N,
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
  base::Printer::getInstance().printStatusBegin("Adaptive grid generation (fuzzy Ritter-Novak)...");

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
  fX.resize(std::max(N, currentN));
  fX.setAll(0.0);
  // degree[i] is the number of times the i-th grid point was
  // chosen for refinement
  std::vector<size_t> degree(fX.getSize(), 0);
  // level_sum[i] is the 1-norm of the level vector of the i-th grid point
  std::vector<size_t> levelSum(fX.getSize(), 0);
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

    // calculate sum of levels
    for (size_t t = 0; t < d; t++) {
      levelSum[i] += gp.getLevel(t);
    }
  }

  // evaluation of f in the initial grid points
  evalFunction();

  // tall matrices containing the bounding boxes for the
  // confidence intervals for the different alpha values; the first row
  // corresponds to the whole domain [0, 1]^d, while the j-th row
  // corresponds to alpha_{j-2} = (j-2)/n
  // (j = 2, ..., n+2, where n := numberOfAlphaSegments)
  const size_t numberOfBounds = numberOfAlphaSegments + 2;
  base::DataMatrix lowerBounds(numberOfBounds, d);
  base::DataMatrix upperBounds(numberOfBounds, d);

  // initialize first row to unit hyper-cube
  for (size_t t = 0; t < d; t++) {
    lowerBounds(0, t) = 0.0;
    upperBounds(0, t) = 1.0;
  }

  // initialize remaining rows
  for (size_t j = 0; j <= numberOfAlphaSegments; j++) {
    const double alphaLevel = static_cast<double>(j) / static_cast<double>(numberOfAlphaSegments);

    for (size_t t = 0; t < d; t++) {
      // determine confidence interval for alpha value
      double lowerBound = xFuzzy[t]->evaluateConfidenceIntervalLowerBound(alphaLevel);
      double upperBound = xFuzzy[t]->evaluateConfidenceIntervalUpperBound(alphaLevel);
      double boundsSize = upperBound - lowerBound;
      const double minimumBoundsSize = 0.05;

      // maintain some minimum size for the confidence interval
      // if too small
      if (boundsSize < minimumBoundsSize) {
        const double boundsCenter = (lowerBound + upperBound) / 2.0;
        lowerBound = boundsCenter - minimumBoundsSize / 2.0;
        upperBound = boundsCenter + minimumBoundsSize / 2.0;
        boundsSize = minimumBoundsSize;
      }

      // enlarge confidence interval by 10%
      // (optimum might be near the boundary)
      lowerBound = std::max(lowerBound - 0.05 * boundsSize, 0.0);
      upperBound = std::min(upperBound + 0.05 * boundsSize, 1.0);

      lowerBounds(j + 1, t) = lowerBound;
      upperBounds(j + 1, t) = upperBound;
    }
  }

  // iteration counter
  size_t k = 0;
  base::DataVector x(d);

  while (currentN < N) {
    // status printing
    {
      char str[10];
      snprintf(str, sizeof(str),
               "%.1f%%", static_cast<double>(currentN) / static_cast<double>(N) * 100.0);
      base::Printer::getInstance().printStatusUpdate(std::string(str) + " (N = " +
                                               std::to_string(currentN) + ", k = " +
                                               std::to_string(k) + ")");
    }

    // list of points to be refined (usually two per confidence interval)
    std::vector<size_t> pointsToBeRefined;
    // for each point to be refined, this says if we already checked
    // whether the children (if inserted) of this point would be too deep
    std::vector<bool> pointAlreadyChecked(currentN, false);

    // for each confidence interval
    for (size_t j = 0; j < numberOfBounds; j++) {
      // vector of indices of the grid point which are contained
      // in the confidence interval
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
        // confidence interval doesn't contain any grid points
        continue;
      }

      // number of grid points in the confidence interval
      const size_t NRestricted = pointsRestricted.size();
      // fXRestricted fulfills fXRestricted[j] = fX[pointsRestricted[j]]
      base::DataVector fXRestricted(NRestricted);
      // fXOrderRestricted fulfills
      // fXOrderRestricted[j] = fXRestricted[fXOrderRestricted[i]]
      std::vector<size_t> fXOrderRestricted(NRestricted);
      // rank fulfills
      // rankRestricted[j] = #{a | fXRestricted[a] <= fXRestricted[j]}
      std::vector<size_t> rankRestricted(NRestricted);

      for (size_t j = 0; j < NRestricted; j++) {
        fXRestricted[j] = fX[pointsRestricted[j]];
        fXOrderRestricted[j] = j;
        rankRestricted[j] = j + 1;
      }

      // calculate order and rank
      std::sort(fXOrderRestricted.begin(), fXOrderRestricted.end(),
                [&fXRestricted](size_t a, size_t b) {
                  return (fXRestricted[a] < fXRestricted[b]);
                });
      std::sort(rankRestricted.begin(), rankRestricted.end(),
                [&fXOrderRestricted](size_t a, size_t b) {
                  return (fXOrderRestricted[a - 1] < fXOrderRestricted[b - 1]);
                });

      // determine the best "i"s (i.e.,
      // iMinBest = argmin_i g_i and iMaxBest = argmax_i g_i)
      double gMinBest = std::numeric_limits<double>::infinity();
      size_t iMinBest;
      double gMaxBest = std::numeric_limits<double>::infinity();
      size_t iMaxBest;

      for (size_t j = 0; j < NRestricted; j++) {
        const size_t i = pointsRestricted[j];

        // refinement criterion
        double gMin;
        double gMax;

        if (powMethod == STD_POW) {
          gMin = std::pow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                 std::pow(static_cast<double>(rankRestricted[j]) + 1.0, 1.0 - gamma);
          gMax = std::pow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                 std::pow(static_cast<double>(NRestricted - rankRestricted[j]) + 2.0, 1.0 - gamma);
        } else {
          gMin = fastPow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                 fastPow(static_cast<double>(rankRestricted[j]) + 1.0, 1.0 - gamma);
          gMax = fastPow(static_cast<double>(levelSum[i] + degree[i]) + 1.0, gamma) *
                 fastPow(static_cast<double>(NRestricted - rankRestricted[j]) + 2.0, 1.0 - gamma);
        }

        // skip if not new minimum or new maximum
        if ((gMin >= gMinBest) && (gMax >= gMaxBest)) {
          continue;
        }

        // children were too "deep" ==> ignore the point
        if (ignore[i]) {
          continue;
        }

        // check if a refinement of this point would generate
        // children with a level greater than max_level
        // (in one coordinate), if yes ignore the point
        if (!pointAlreadyChecked[i]) {
          pointAlreadyChecked[i] = true;

          base::GridPoint& gp = gridStorage[i];
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

        if (gMin < gMinBest) {
          // no ignore ==> new candidate for the point to be refined
          gMinBest = gMin;
          iMinBest = i;
        }

        if (gMax < gMaxBest) {
          // no ignore ==> new candidate for the point to be refined
          gMaxBest = gMax;
          iMaxBest = i;
        }
      }

      // add minimal point to pointsToBeRefined
      if (gMinBest < std::numeric_limits<double>::infinity()) {
        pointsToBeRefined.push_back(iMinBest);
      }

      // add maximal point to pointsToBeRefined
      if (gMaxBest < std::numeric_limits<double>::infinity()) {
        pointsToBeRefined.push_back(iMaxBest);
      }
    }

    // sort pointsToBeRefined and eliminate duplicates
    std::sort(pointsToBeRefined.begin(), pointsToBeRefined.end());
    pointsToBeRefined.erase(std::unique(pointsToBeRefined.begin(), pointsToBeRefined.end()),
                            pointsToBeRefined.end());

    bool exit = false;

    for (size_t iBest : pointsToBeRefined) {
      // refine point no. iBest
      degree[iBest]++;
      refinementAlpha[iBest] = 1.0;
      base::SurplusRefinementFunctor refineFunc(refinementAlpha, 1);
      refinement.free_refine(gridStorage, refineFunc);

      // new grid size
      const size_t newN = gridStorage.getSize();

      if (newN == currentN) {
        // size unchanged ==> point not refined (should not happen)
        base::Printer::getInstance().printStatusEnd(
            "error: size unchanged in IterativeGridGeneratorFuzzyRitterNovak");
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

      currentN = newN;
    }

    if (exit) {
      break;
    }

    // next round
    k++;
  }

  // delete superfluous entries in fX
  fX.resize(currentN);

  if (result) {
    base::Printer::getInstance().printStatusUpdate(
        "100.0% (N = " + std::to_string(currentN) + ", k = " + std::to_string(k) + ")");
    base::Printer::getInstance().printStatusEnd();
    return true;
  } else {
    return false;
  }
}

}  // namespace optimization
}  // namespace sgpp
