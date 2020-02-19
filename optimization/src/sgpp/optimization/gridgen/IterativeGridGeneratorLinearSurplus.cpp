// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/gridgen/IterativeGridGeneratorLinearSurplus.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/solver/BiCGStab.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisBoundaryGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <string>

namespace sgpp {
namespace optimization {

IterativeGridGeneratorLinearSurplus::IterativeGridGeneratorLinearSurplus(base::ScalarFunction& f,
                                                                         base::Grid& grid, size_t N,
                                                                         double adaptivity,
                                                                         base::level_t initialLevel)
    : IterativeGridGenerator(f, grid, N), gamma(adaptivity), initialLevel(initialLevel) {
  if ((grid.getType() == base::GridType::Bspline) || (grid.getType() == base::GridType::Wavelet) ||
      (grid.getType() == base::GridType::Linear) ||
      (grid.getType() == base::GridType::NakBspline) ||
      (grid.getType() == base::GridType::FundamentalSpline)) {
    linearGrid = std::unique_ptr<base::Grid>(new base::LinearGrid(f.getNumberOfParameters()));
  } else if ((grid.getType() == base::GridType::BsplineBoundary) ||
             (grid.getType() == base::GridType::WaveletBoundary) ||
             (grid.getType() == base::GridType::LinearBoundary) ||
             (grid.getType() == base::GridType::FundamentalNakSplineBoundary) ||
             (grid.getType() == base::GridType::FundamentalSplineBoundary) ||
             (grid.getType() == base::GridType::WeaklyFundamentalNakSplineBoundary) ||
             (grid.getType() == base::GridType::WeaklyFundamentalSplineBoundary) ||
             (grid.getType() == base::GridType::NaturalBsplineBoundary) ||
             (grid.getType() == base::GridType::NakBsplineBoundary)) {
    linearGrid =
        std::unique_ptr<base::Grid>(new base::LinearBoundaryGrid(f.getNumberOfParameters()));
  } else if ((grid.getType() == base::GridType::BsplineClenshawCurtis) ||
             (grid.getType() == base::GridType::LinearClenshawCurtisBoundary)) {
    linearGrid = std::unique_ptr<base::Grid>(
        new base::LinearClenshawCurtisBoundaryGrid(f.getNumberOfParameters()));
  } else if (grid.getType() == base::GridType::LinearClenshawCurtis) {
    linearGrid = std::unique_ptr<base::Grid>(
        new base::LinearClenshawCurtisGrid(f.getNumberOfParameters()));
  } else if ((grid.getType() == base::GridType::ModBspline) ||
             (grid.getType() == base::GridType::ModWavelet) ||
             (grid.getType() == base::GridType::ModLinear) ||
             (grid.getType() == base::GridType::ModFundamentalSpline) ||
             (grid.getType() == base::GridType::ModWeaklyFundamentalNakSpline) ||
             (grid.getType() == base::GridType::ModNakBspline)) {
    linearGrid = std::unique_ptr<base::Grid>(new base::ModLinearGrid(f.getNumberOfParameters()));
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    linearGrid = std::unique_ptr<base::Grid>(
        new base::ModBsplineClenshawCurtisGrid(f.getNumberOfParameters(), 1));
  } else if (grid.getType() == base::GridType::NakBsplineExtended) {
    linearGrid =
        std::unique_ptr<base::Grid>(new base::NakBsplineExtendedGrid(f.getNumberOfParameters(), 1));
  } else {
    throw std::invalid_argument("Grid type not supported.");
  }
}

IterativeGridGeneratorLinearSurplus::~IterativeGridGeneratorLinearSurplus() {}

double IterativeGridGeneratorLinearSurplus::getAdaptivity() const { return gamma; }

void IterativeGridGeneratorLinearSurplus::setAdaptivity(double adaptivity) {
  this->gamma = adaptivity;
}

base::level_t IterativeGridGeneratorLinearSurplus::getInitialLevel() const { return initialLevel; }

void IterativeGridGeneratorLinearSurplus::setInitialLevel(base::level_t initialLevel) {
  this->initialLevel = initialLevel;
}

bool IterativeGridGeneratorLinearSurplus::generate() {
  base::Printer::getInstance().printStatusBegin("Adaptive grid generation (linear surplus)...");

  std::unique_ptr<base::AbstractRefinement> abstractRefinement;

  if ((grid.getType() == base::GridType::BsplineBoundary) ||
      (grid.getType() == base::GridType::WaveletBoundary) ||
      (grid.getType() == base::GridType::LinearBoundary) ||
      (grid.getType() == base::GridType::BsplineClenshawCurtis) ||
      (grid.getType() == base::GridType::LinearClenshawCurtisBoundary) ||
      (grid.getType() == base::GridType::WeaklyFundamentalNakSplineBoundary) ||
      (grid.getType() == base::GridType::WeaklyFundamentalSplineBoundary) ||
      (grid.getType() == base::GridType::NaturalBsplineBoundary) ||
      (grid.getType() == base::GridType::NakBsplineBoundary)) {
    // grid with boundaries
    abstractRefinement =
        std::unique_ptr<base::AbstractRefinement>(new base::HashRefinementBoundaries());
  } else {
    // grid without boundaries
    abstractRefinement = std::unique_ptr<base::AbstractRefinement>(new base::HashRefinement());
  }

  base::GridStorage& gridStorage = grid.getStorage();
  // Set up linear system for hierarchization with the linear grid as
  // underlying grid, but set the internal grid storage to the one of the
  // B-spline/wavelet grid.
  // hierSLE is used for two purposes: 1. initial hierarchization
  // (with a linear solver),
  // 2. in the algorithm loop: evaluation of basis functions at grid points
  // (no linear solver due to forward substitution).
  // The right-hand side of the system does only matter for the
  // first purpose.
  base::HierarchisationSLE hierSLE(*linearGrid, gridStorage);

  // generate initial grid
  grid.getGenerator().regular(initialLevel);

  size_t currentN = gridStorage.getSize();
  // coeffs always has as much elements as there are grid points in
  // the grid, but fX has N elements (no resizing during the main loop)
  base::DataVector coeffs(currentN);
  // abbreviation (functionValues is a member variable of
  // IterativeGridGenerator)
  base::DataVector& fX = functionValues;

  fX.resize(std::max(N, currentN));
  fX.setAll(0.0);

  // parallel evaluation of f in the initial grid points
  evalFunction();

  for (size_t i = 0; i < currentN; i++) {
    coeffs[i] = fX[i];
  }

  // initial hierarchization
  if (currentN > 1) {
    base::DataVector fXCutoff(fX.getPointer(), currentN);
    base::sle_solver::BiCGStab sleSolver;

    // solve system
    base::Printer::getInstance().disableStatusPrinting();
    sleSolver.solve(hierSLE, fXCutoff, coeffs);
    base::Printer::getInstance().enableStatusPrinting();
  }

  // iteration counter
  size_t k = 0;

  // number of refinable points
  size_t refinablePtsCount;
  // number of points that should be refined
  size_t ptsToBeRefinedCount;
  // factor to be multiplied to pts_to_be_refined_count to get the
  // *real* number of points to be refined
  // (used like 1, 1/2, 1/4, 1/8, ... at the end of the algorithm to
  // fill the grid size up to N)
  double refineFactor = 1.0;
  // success or not?
  bool result = true;

  while (true) {
    // status printing
    {
      char str[10];
      snprintf(str, sizeof(str), "%.1f%%",
               static_cast<double>(currentN) / static_cast<double>(N) * 100.0);
      base::Printer::getInstance().printStatusUpdate(std::string(str) +
                                                     " (N = " + std::to_string(currentN) + ")");
    }

    // calculate number of points to be refined
    refinablePtsCount = abstractRefinement->getNumberOfRefinablePoints(gridStorage);
    ptsToBeRefinedCount =
        static_cast<int>(1.0 + refineFactor * gamma * static_cast<double>(refinablePtsCount));

    // refine
    base::SurplusRefinementFunctor refineFunc(coeffs, ptsToBeRefinedCount);
    abstractRefinement->free_refine(gridStorage, refineFunc);

    // new grid size
    size_t newN = gridStorage.getSize();

    if (newN == currentN) {
      // size unchanged ==> nothing refined (should not happen)
      base::Printer::getInstance().printStatusEnd(
          "error: size unchanged in IterativeGridGeneratorLinearSurplus");
      result = false;
      break;
    }

    if (newN > N) {
      // too many new points ==> undo refinement and try again with
      // refineFactor halved
      undoRefinement(currentN);

      if (ptsToBeRefinedCount == 1) {
        break;
      } else {
        refineFactor /= 2.0;
        k++;
        continue;
      }
    }

    coeffs.resize(newN);

    // evaluation of f in the new grid points
    evalFunction(currentN);

    // forward substitution
    // (hierSLE should always be a lower triangular matrix)
    for (size_t i = currentN; i < newN; i++) {
      coeffs[i] = fX[i];

      for (size_t j = 0; j < i; j++) {
        coeffs[i] -= hierSLE.getMatrixEntry(i, j) * coeffs[j];
      }
    }

    // next round
    currentN = newN;
    k++;
  }

  // delete superfluous entries in fX
  fX.resize(currentN);

  if (result) {
    base::Printer::getInstance().printStatusUpdate("100.0% (N = " + std::to_string(currentN) + ")");
    base::Printer::getInstance().printStatusEnd();
    return true;
  } else {
    return false;
  }
}
}  // namespace optimization
}  // namespace sgpp
