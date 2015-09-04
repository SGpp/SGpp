// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <cstring>
#include <stdexcept>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/gridgen/IterativeGridGeneratorLinearSurplus.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>

#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/sle/solver/BiCGStab.hpp>

namespace SGPP {
  namespace optimization {

    IterativeGridGeneratorLinearSurplus::IterativeGridGeneratorLinearSurplus(
      ScalarFunction& f, base::Grid& grid, size_t N, float_t adaptivity) :
      IterativeGridGenerator(f, grid, N),
      gamma(adaptivity) {
      if ((std::strcmp(grid.getType(),
                       "bspline") == 0) ||
          (std::strcmp(grid.getType(),
                       "wavelet") == 0) ||
          (std::strcmp(grid.getType(),
                       "linear") == 0) ||
          (std::strcmp(grid.getType(),
                       "fundamentalSpline") == 0)) {
        linearGrid = std::unique_ptr<base::Grid>(
                       new base::LinearGrid(f.getDimension()));
      } else if ((std::strcmp(grid.getType(),
                              "bsplineBoundary") == 0) ||
                 (std::strcmp(grid.getType(),
                              "waveletBoundary") == 0) ||
                 (std::strcmp(grid.getType(),
                              "linearBoundary") == 0)) {
        linearGrid = std::unique_ptr<base::Grid>(
                       new base::LinearBoundaryGrid(
                         f.getDimension()));
      } else if ((std::strcmp(grid.getType(),
                              "bsplineClenshawCurtis") == 0) ||
                 (std::strcmp(grid.getType(),
                              "linearClenshawCurtis") == 0)) {
        linearGrid = std::unique_ptr<base::Grid>(
                       new base::LinearClenshawCurtisGrid(f.getDimension()));
      } else if ((std::strcmp(grid.getType(),
                              "modBspline") == 0) ||
                 (std::strcmp(grid.getType(),
                              "modWavelet") == 0) ||
                 (std::strcmp(grid.getType(),
                              "modlinear") == 0) ||
                 (std::strcmp(grid.getType(),
                              "modFundamentalSpline") == 0)) {
        linearGrid = std::unique_ptr<base::Grid>(
                       new base::ModLinearGrid(f.getDimension()));
      } else if (std::strcmp(grid.getType(),
                             "modBsplineClenshawCurtis") == 0) {
        linearGrid = std::unique_ptr<base::Grid>(
                       new base::ModBsplineClenshawCurtisGrid(
                         f.getDimension(), 1));
      } else {
        throw std::invalid_argument("Grid type not supported.");
      }
    }

    float_t IterativeGridGeneratorLinearSurplus::getAdaptivity() const {
      return gamma;
    }

    void IterativeGridGeneratorLinearSurplus::setAdaptivity(float_t adaptivity) {
      this->gamma = adaptivity;
    }

    bool IterativeGridGeneratorLinearSurplus::generate() {
      printer.printStatusBegin("Adaptive grid generation (linear surplus)...");

      base::GridIndex::PointDistribution distr = base::GridIndex::Normal;

      if ((std::strcmp(grid.getType(), "bsplineClenshawCurtis") == 0) ||
          (std::strcmp(grid.getType(), "modBsplineClenshawCurtis") == 0) ||
          (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0)) {
        // Clenshaw-Curtis grid
        distr = base::GridIndex::ClenshawCurtis;
      }

      std::unique_ptr<base::AbstractRefinement> abstractRefinement;

      if ((std::strcmp(grid.getType(), "bsplineBoundary") == 0) ||
          (std::strcmp(grid.getType(), "waveletBoundary") == 0) ||
          (std::strcmp(grid.getType(), "linearBoundary") == 0) ||
          (std::strcmp(grid.getType(), "bsplineClenshawCurtis") == 0) ||
          (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0)) {
        // grid with boundaries
        abstractRefinement = std::unique_ptr<base::AbstractRefinement>(
                               new base::HashRefinementBoundaries());
      } else {
        // grid without boundaries
        abstractRefinement = std::unique_ptr<base::AbstractRefinement>(
                               new base::HashRefinement());
      }

      base::GridStorage& gridStorage = *grid.getStorage();
      // Set up linear system for hierarchization with the linear grid as
      // underlying grid, but set the internal grid storage to the one of the
      // B-spline/wavelet grid.
      // hierSLE is used for two purposes: 1. initial hierarchization
      // (with a linear solver),
      // 2. in the algorithm loop: evaluation of basis functions at grid points
      // (no linear solver due to forward substitution).
      // The right-hand side of the system does only matter for the
      // first purpose.
      HierarchisationSLE hierSLE(*linearGrid, gridStorage);

      // generate initial grid
      {
        std::unique_ptr<base::GridGenerator> gridGen(
          grid.createGridGenerator());
        gridGen->regular(3);
      }

      size_t currentN = gridStorage.size();
      // coeffs always has as much elements as there are grid points in
      // the grid, but fX has N elements (no resizing during the main loop)
      base::DataVector coeffs(currentN);
      // abbreviation (functionValues is a member variable of
      // IterativeGridGenerator)
      base::DataVector& fX = functionValues;

      fX.resize(std::max(N, currentN));
      fX.setAll(0.0);

      for (size_t i = 0; i < currentN; i++) {
        // set correct point distribution
        gridStorage[i]->setPointDistribution(distr);
      }

      // parallel evaluation of f in the initial grid points
      evalFunction();

      for (size_t i = 0; i < currentN; i++) {
        coeffs[i] = fX[i];
      }

      // initial hierarchization
      if (currentN > 1) {
        base::DataVector fXCutoff(fX.getPointer(), currentN);
        sle_solver::BiCGStab sleSolver;

        // solve system
        printer.disableStatusPrinting();
        sleSolver.solve(hierSLE, fXCutoff, coeffs);
        printer.enableStatusPrinting();
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
      float_t refineFactor = 1.0;
      // success or not?
      bool result = true;

      while (true) {
        // status printing
        {
          char str[10];
          snprintf(str, 10, "%.1f%%",
                   static_cast<float_t>(currentN) /
                   static_cast<float_t>(N) * 100.0);
          printer.printStatusUpdate(std::string(str) +
                                    " (N = " + std::to_string(currentN) + ")");
        }

        // calculate number of points to be refined
        refinablePtsCount =
          abstractRefinement->getNumberOfRefinablePoints(&gridStorage);
        ptsToBeRefinedCount =
          static_cast<int>(1.0 + refineFactor * gamma *
                           static_cast<float_t>(refinablePtsCount));

        // refine
        base::SurplusRefinementFunctor refineFunc(&coeffs,
            ptsToBeRefinedCount);
        abstractRefinement->free_refine(&gridStorage, &refineFunc);

        // new grid size
        size_t newN = gridStorage.size();

        if (newN == currentN) {
          // size unchanged ==> nothing refined (should not happen)
          printer.printStatusEnd(
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

        for (size_t i = currentN; i < newN; i++) {
          // set point distribution accordingly to
          // normal/Clenshaw-Curtis grids
          gridStorage[i]->setPointDistribution(distr);
        }

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
        printer.printStatusUpdate("100.0% (N = " +
                                  std::to_string(currentN) + ")");
        printer.printStatusEnd();
        return true;
      } else {
        return false;
      }
    }

  }
}
