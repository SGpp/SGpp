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
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearTruncatedBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>

#include <sgpp/optimization/sle/system/Hierarchisation.hpp>
#include <sgpp/optimization/sle/solver/BiCGStab.hpp>

namespace SGPP {
  namespace optimization {
    namespace gridgen {

      const float_t IterativeGridGeneratorLinearSurplus::DEFAULT_ALPHA = 0.2;

      IterativeGridGeneratorLinearSurplus::IterativeGridGeneratorLinearSurplus(
        function::Objective& f, base::Grid& grid, size_t N, float_t alpha) :
        IterativeGridGenerator(f, grid, N),
        alpha(alpha) {
        if ((std::strcmp(grid.getType(), "bspline") == 0) ||
            (std::strcmp(grid.getType(), "wavelet") == 0) ||
            (std::strcmp(grid.getType(), "linear") == 0)) {
          linearGrid = std::unique_ptr<base::Grid>(
                         new base::LinearGrid(f.getDimension()));
        } else if ((std::strcmp(grid.getType(), "bsplineTruncatedBoundary") == 0) ||
                   (std::strcmp(grid.getType(), "waveletTruncatedBoundary") == 0) ||
                   (std::strcmp(grid.getType(), "linearTruncatedBoundary") == 0)) {
          linearGrid = std::unique_ptr<base::Grid>(
                         new base::LinearTruncatedBoundaryGrid(f.getDimension()));
        } else if ((std::strcmp(grid.getType(), "bsplineClenshawCurtis") == 0) ||
                   (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0)) {
          linearGrid = std::unique_ptr<base::Grid>(
                         new base::LinearClenshawCurtisGrid(f.getDimension()));
        } else if ((std::strcmp(grid.getType(), "modBspline") == 0) ||
                   (std::strcmp(grid.getType(), "modWavelet") == 0) ||
                   (std::strcmp(grid.getType(), "modlinear") == 0)) {
          linearGrid = std::unique_ptr<base::Grid>(
                         new base::ModLinearGrid(f.getDimension()));
        } else {
          throw std::invalid_argument("Grid type not supported.");
        }
      }

      float_t IterativeGridGeneratorLinearSurplus::getAlpha() const {
        return alpha;
      }

      void IterativeGridGeneratorLinearSurplus::setAlpha(float_t alpha) {
        this->alpha = alpha;
      }

      bool IterativeGridGeneratorLinearSurplus::generate() {
        tools::printer.printStatusBegin("Adaptive grid generation (linear surplus)...");

        base::GridIndex::PointDistribution distr = base::GridIndex::Normal;

        if ((std::strcmp(grid.getType(), "bsplineClenshawCurtis") == 0) ||
            (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0)) {
          // Clenshaw-Curtis grid
          distr = base::GridIndex::ClenshawCurtis;
        }

        std::unique_ptr<base::AbstractRefinement> abstractRefinement;

        if ((std::strcmp(grid.getType(), "bsplineTruncatedBoundary") == 0) ||
            (std::strcmp(grid.getType(), "waveletTruncatedBoundary") == 0) ||
            (std::strcmp(grid.getType(), "linearTruncatedBoundary") == 0) ||
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

        base::GridStorage* gridStorage = grid.getStorage();
        // Set up linear system for hierarchisation with the linear grid as underlying grid,
        // but set the internal grid storage to the one of the B-spline/wavelet grid.
        // hier_system is used for two purposes: 1. initial hierarchisation (with a linear solver),
        // 2. in the algorithm loop: evaluation of basis functions at grid points
        // (no linear solver due to forward substitution).
        // The right-hand side of the system does only matter for the first purpose.
        sle::system::Hierarchisation hierSystem(*linearGrid);
        hierSystem.setGridStorage(gridStorage);

        // generate initial grid
        {
          std::unique_ptr<base::GridGenerator> gridGen(grid.createGridGenerator());
          gridGen->regular(3);
        }

        const size_t d = gridStorage->dim();
        size_t currentN = gridStorage->size();
        // coeffs always has as much elements as there are grid points in the grid,
        // but fX has N elements (no resizing during the main loop)
        base::DataVector coeffs(currentN);
        // abbreviation (functionValues is a member variable of IterativeGridGenerator)
        std::vector<float_t>& fX = functionValues;

        fX.assign(std::max(N, currentN), 0.0);

        for (size_t i = 0; i < currentN; i++) {
          // set correct point distribution
          gridStorage->get(i)->setPointDistribution(distr);
        }

        // parallel evaluation of f in the initial grid points
        #pragma omp parallel shared(fX, coeffs, currentN, gridStorage) default(none)
        {
          std::vector<float_t> x(d, 0.0);
          function::Objective* curFPtr = &f;
#ifdef _OPENMP
          std::unique_ptr<function::Objective> curF;

          if (omp_get_max_threads() > 1) {
            f.clone(curFPtr);
            curF = std::unique_ptr<function::Objective>(curFPtr);
          }

#endif

          #pragma omp for

          for (size_t i = 0; i < currentN; i++) {
            // convert grid point to coordinate vector
            #pragma omp critical
            {
              base::GridIndex* gp = gridStorage->get(i);

              for (size_t t = 0; t < d; t++) {
                x[t] = gp->getCoord(t);
              }
            }

            float_t fx = curFPtr->eval(x);

            #pragma omp critical
            {
              fX[i] = fx;
              coeffs[i] = fx;
            }
          }
        }

        // initial hierarchisation
        if (currentN > 1) {
          std::vector<float_t> fXCutoff(fX.begin(), fX.begin() + currentN);
          std::vector<float_t> coeffsVec;
          sle::solver::BiCGStab sleSolver;

          // solve system
          tools::printer.disableStatusPrinting();
          sleSolver.solve(hierSystem, fXCutoff, coeffsVec);
          tools::printer.enableStatusPrinting();

          // convert std::vector to base::DataVector
          coeffs = base::DataVector(&coeffsVec[0], currentN);
        }

        // iteration counter
        size_t k = 0;

        // number of refinable points
        size_t refinablePtsCount;
        // number of points that should be refined
        size_t ptsToBeRefinedCount;
        // factor to be multiplied to pts_to_be_refined_count to get the *real* number of
        // points to be refined
        // (used like 1, 1/2, 1/4, 1/8, ... at the end of the algorithm to fill the grid size up to N)
        float_t refineFactor = 1.0;
        // success or not?
        bool result = true;

        while (true) {
          // status printing
          {
            char str[10];
            snprintf(str, 10, "%.1f%%",
                     static_cast<float_t>(currentN) / static_cast<float_t>(N) * 100.0);
            tools::printer.printStatusUpdate(std::string(str) +
                                             " (N = " + toString(currentN) + ")");
          }

          // calculate number of points to be refined
          refinablePtsCount = abstractRefinement->getNumberOfRefinablePoints(gridStorage);
          ptsToBeRefinedCount = static_cast<int>(
                                  1.0 + refineFactor * alpha * static_cast<float_t>(refinablePtsCount));

          // refine
          base::SurplusRefinementFunctor refineFunc(&coeffs, ptsToBeRefinedCount);
          abstractRefinement->free_refine(gridStorage, &refineFunc);

          // new grid size
          size_t newN = gridStorage->size();

          if (newN == currentN) {
            // size unchanged ==> nothing refined (should not happen)
            tools::printer.printStatusEnd(
              "error: size unchanged in IterativeGridGeneratorLinearSurplus");
            result = false;
            break;
          }

          if (newN > N) {
            // too many new points ==> undo refinement and try again with refine_factor halved
            std::list<size_t> indicesToRemove;

            for (size_t i = currentN; i < newN; i++) {
              indicesToRemove.push_back(i);
            }

            gridStorage->deletePoints(indicesToRemove);

            if (ptsToBeRefinedCount == 1) {
              break;
            } else {
              refineFactor /= 2.0;
              k++;
              continue;
            }
          }

          coeffs.resize(newN);

          // parallel evaluation of f in the new grid points
          #pragma omp parallel shared(fX, currentN, newN, gridStorage, distr) default(none)
          {
            std::vector<float_t> x(d, 0.0);
            function::Objective* curFPtr = &f;
#ifdef _OPENMP
            std::unique_ptr<function::Objective> curF;

            if (omp_get_max_threads() > 1) {
              f.clone(curFPtr);
              curF = std::unique_ptr<function::Objective>(curFPtr);
            }

#endif

            #pragma omp for

            for (size_t i = currentN; i < newN; i++) {
              // convert grid point to coordinate vector
              #pragma omp critical
              {
                base::GridIndex* gp = gridStorage->get(i);
                // set point distribution accordingly to normal/Clenshaw-Curtis grids
                gp->setPointDistribution(distr);

                for (size_t t = 0; t < d; t++) {
                  x[t] = gp->getCoord(t);
                }
              }

              float_t fx = curFPtr->eval(x);

              #pragma omp critical
              fX[i] = fx;
            }
          }

          // forward substitution (hier_system should always be a lower triangular matrix)
          for (size_t i = currentN; i < newN; i++) {
            coeffs[i] = fX[i];

            for (size_t j = 0; j < i; j++) {
              coeffs[i] -= hierSystem.getMatrixEntry(i, j) * coeffs[j];
            }
          }

          // next round
          currentN = newN;
          k++;
        }

        // delete superfluous entries in fX
        fX.erase(fX.begin() + currentN, fX.end());

        if (result) {
          tools::printer.printStatusUpdate("100.0% (N = " + toString(currentN) + ")");
          tools::printer.printStatusEnd();
          return true;
        } else {
          return false;
        }
      }

    }
  }
}
