// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/gridgen/HashRefinementMultiple.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorSOO.hpp>

#include <algorithm>
#include <cstring>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

namespace sgpp {
namespace optimization {

IterativeGridGeneratorSOO::IterativeGridGeneratorSOO(base::ScalarFunction& f, base::Grid& grid,
                                                     size_t N, double adaptivity)
    : IterativeGridGenerator(f, grid, N) {
  setAdaptivity(adaptivity);
}

IterativeGridGeneratorSOO::~IterativeGridGeneratorSOO() {}

IterativeGridGeneratorSOO::AdaptivityFunction IterativeGridGeneratorSOO::getAdaptivity() const {
  return hMax;
}

void IterativeGridGeneratorSOO::setAdaptivity(double adaptivity) {
  hMax = [adaptivity](size_t n) { return static_cast<size_t>(std::pow(n, adaptivity)); };
}

void IterativeGridGeneratorSOO::setAdaptivity(
    IterativeGridGeneratorSOO::AdaptivityFunction adaptivity) {
  hMax = adaptivity;
}

bool IterativeGridGeneratorSOO::generate() {
  sgpp::base::Printer::getInstance().printStatusBegin("Adaptive grid generation (SOO)...");

  bool result = true;
  base::GridStorage& gridStorage = grid.getStorage();
  const size_t d = f.getNumberOfParameters();

  HashRefinementMultiple refinement;

  // generate initial grid
  {
    base::GridPoint gp(d);

    for (size_t t = 0; t < d; t++) {
      gp.set(t, 1, 1);
    }

    gridStorage.insert(gp);
  }

  size_t currentN = 1;

  // depthMap[depth] is a list of indices j with delta(x[j]) = depth
  std::vector<std::vector<size_t>> depthMap{{0}};
  // if the i-th grid point is refinable
  std::vector<bool> refinable(N, true);

  // abbreviation (functionValues is a member variable of
  // IterativeGridGenerator)
  base::DataVector& fX = functionValues;
  fX.resize(N);

  {
    base::DataVector x(gridStorage.getCoordinates(gridStorage[0]));
    fX[0] = f.eval(x);
  }

  base::DataVector refinementAlpha(1, 0.0);

  size_t depthBoundOffset = 0;
  size_t n = 0;
  bool breakLoop = false;

  // iteration counter
  size_t k = 0;

  while (!breakLoop) {
    // status printing
    {
      char str[10];
      snprintf(str, sizeof(str), "%.1f%%",
               static_cast<double>(currentN) / static_cast<double>(N) * 100.0);
      sgpp::base::Printer::getInstance().printStatusUpdate(std::string(str) +
                                                           " (N = " + std::to_string(currentN) +
                                                           ", k = " + std::to_string(k) + ")");
    }

    const size_t curDepthBound =
        std::min(depthMap.size() - 1, static_cast<size_t>(hMax(n)) + depthBoundOffset);
    double nuMin = std::numeric_limits<double>::infinity();

    for (size_t depth = 0; depth <= curDepthBound; depth++) {
      double fBest = std::numeric_limits<double>::infinity();
      size_t iBest = 0;

      for (size_t i : depthMap[depth]) {
        if (refinable[i] && (fX[i] < fBest)) {
          fBest = fX[i];
          iBest = i;
        }
      }

      if (fBest < nuMin) {
        refinementAlpha[iBest] = 1.0;
        base::SurplusRefinementFunctor refineFunc(refinementAlpha, 1);
        refinement.free_refine(gridStorage, refineFunc);

        // new grid size
        const size_t newN = gridStorage.getSize();

        if (newN == currentN) {
          // size unchanged ==> point not refined (should not happen)
          sgpp::base::Printer::getInstance().printStatusEnd(
              "error: size unchanged in IterativeGridGeneratorSOO");
          result = false;
          breakLoop = true;
          break;
        }

        if (newN > N) {
          // too many new points ==> undo refinement and exit
          undoRefinement(currentN);
          breakLoop = true;
          break;
        }

        // resize refinement vector and set to all zeros
        // (in the following loop)
        refinementAlpha.resize(newN);
        refinementAlpha[iBest] = 0.0;

        for (size_t i = currentN; i < newN; i++) {
          base::GridPoint& gp = gridStorage[i];
          refinementAlpha[i] = 0.0;
          size_t depth = 0;

          // calculate sum of levels
          for (size_t t = 0; t < d; t++) {
            depth += gp.getLevel(t);
          }

          depth -= d;

          for (size_t depth2 = depthMap.size(); depth2 <= depth; depth2++) {
            depthMap.push_back(std::vector<size_t>());
          }

          depthMap[depth].push_back(i);
        }

        // evaluation of f in the new grid points
        evalFunction(currentN);

        refinable[iBest] = false;
        n++;
        nuMin = fBest;

        currentN = newN;
      }
    }

    if (std::isinf(nuMin)) {
      depthBoundOffset++;
    }

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
