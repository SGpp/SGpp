/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include <cstring>
#include <stdexcept>

#include "opt/gridgen/IterativeGridGeneratorLinearSurplus.hpp"
#include "opt/tools/Printer.hpp"
#include "base/grid/generation/hashmap/HashRefinement.hpp"
#include "base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "base/grid/type/BsplineClenshawCurtisGrid.hpp"

#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/LinearClenshawCurtisGrid.hpp"
#include "base/grid/type/ModLinearGrid.hpp"

#include "opt/sle/system/Hierarchisation.hpp"
#include "opt/sle/solver/BiCGStab.hpp"

namespace sg {
  namespace opt {
    namespace gridgen {

      const double IterativeGridGeneratorLinearSurplus::DEFAULT_ALPHA = 0.2;

      IterativeGridGeneratorLinearSurplus::IterativeGridGeneratorLinearSurplus(
        function::Objective& f, base::Grid& grid, size_t N, double alpha) :
        IterativeGridGenerator(f, grid, N),
        alpha(alpha) {
        if ((std::strcmp(grid.getType(), "Bspline") == 0) ||
            (std::strcmp(grid.getType(), "Wavelet") == 0) ||
            (std::strcmp(grid.getType(), "linear") == 0)) {
          linear_grid = tools::SmartPointer<base::Grid>(
                          new base::LinearGrid(f.getDimension()));
        } else if ((std::strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0) ||
                   (std::strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0) ||
                   (std::strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)) {
          linear_grid = tools::SmartPointer<base::Grid>(
                          new base::LinearTrapezoidBoundaryGrid(f.getDimension()));
        } else if ((std::strcmp(grid.getType(), "BsplineClenshawCurtis") == 0) ||
                   (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0)) {
          linear_grid = tools::SmartPointer<base::Grid>(
                          new base::LinearClenshawCurtisGrid(f.getDimension()));
        } else if ((std::strcmp(grid.getType(), "modBspline") == 0) ||
                   (std::strcmp(grid.getType(), "modWavelet") == 0) ||
                   (std::strcmp(grid.getType(), "modlinear") == 0)) {
          linear_grid = tools::SmartPointer<base::Grid>(
                          new base::ModLinearGrid(f.getDimension()));
        } else {
          throw std::invalid_argument("Grid type not supported.");
        }
      }

      double IterativeGridGeneratorLinearSurplus::getAlpha() const {
        return alpha;
      }

      void IterativeGridGeneratorLinearSurplus::setAlpha(double alpha) {
        this->alpha = alpha;
      }

      bool IterativeGridGeneratorLinearSurplus::generate() {
        tools::printer.printStatusBegin("Adaptive grid generation (linear surplus)...");

        base::GridIndex::PointDistribution distr = base::GridIndex::Normal;

        if ((std::strcmp(grid.getType(), "BsplineClenshawCurtis") == 0) ||
            (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0)) {
          // Clenshaw-Curtis grid
          distr = base::GridIndex::ClenshawCurtis;
        }

        tools::SmartPointer<base::AbstractRefinement> abstract_refinement;

        if ((std::strcmp(grid.getType(), "BsplineTrapezoidBoundary") == 0) ||
            (std::strcmp(grid.getType(), "WaveletTrapezoidBoundary") == 0) ||
            (std::strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) ||
            (std::strcmp(grid.getType(), "BsplineClenshawCurtis") == 0) ||
            (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0)) {
          // grid with boundaries
          abstract_refinement = tools::SmartPointer<base::AbstractRefinement>(
                                  new base::HashRefinementBoundaries());
        } else {
          // grid without boundaries
          abstract_refinement = tools::SmartPointer<base::AbstractRefinement>(
                                  new base::HashRefinement());
        }

        base::GridStorage* grid_storage = grid.getStorage();
        // Set up linear system for hierarchisation with the linear grid as underlying grid,
        // but set the internal grid storage to the one of the B-spline/wavelet grid.
        // hier_system is used for two purposes: 1. initial hierarchisation (with a linear solver),
        // 2. in the algorithm loop: evaluation of basis functions at grid points
        // (no linear solver due to forward substitution).
        // The right-hand side of the system does only matter for the first purpose.
        sle::system::Hierarchisation hier_system(*linear_grid);
        hier_system.setGridStorage(grid_storage);

        // generate initial grid
        {
          tools::SmartPointer<base::GridGenerator> grid_gen(grid.createGridGenerator());
          grid_gen->regular(3);
        }

        const size_t d = grid_storage->dim();
        size_t current_N = grid_storage->size();
        // coeffs always has as much elements as there are grid points in the grid,
        // but fX has N elements (no resizing during the main loop)
        base::DataVector coeffs(current_N);
        // abbreviation (function_values is a member variable of IterativeGridGenerator)
        std::vector<double>& fX = function_values;

        fX.assign(std::max(N, current_N), 0.0);

        for (size_t i = 0; i < current_N; i++) {
          // set correct point distribution
          grid_storage->get(i)->setPointDistribution(distr);
        }

        // parallel evaluation of f in the initial grid points
        #pragma omp parallel shared(fX, coeffs, current_N, grid_storage) default(none)
        {
          std::vector<double> x(d, 0.0);
          tools::SmartPointer<function::Objective> cur_f(f.clone());

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

            double fx = cur_f->eval(x);

            #pragma omp critical
            {
              fX[i] = fx;
              coeffs[i] = fx;
            }
          }
        }

        // initial hierarchisation
        if (current_N > 1) {
          std::vector<double> fX_cutoff(fX.begin(), fX.begin() + current_N);
          std::vector<double> coeffs_vec;
          sle::solver::BiCGStab sle_solver;

          // solve system
          tools::printer.disableStatusPrinting();
          sle_solver.solve(hier_system, fX_cutoff, coeffs_vec);
          tools::printer.enableStatusPrinting();

          // convert std::vector to base::DataVector
          coeffs = base::DataVector(&coeffs_vec[0], current_N);
        }

        // iteration counter
        size_t k = 0;

        // number of refinable points
        size_t refinable_pts_count;
        // number of points that should be refined
        size_t pts_to_be_refined_count;
        // factor to be multiplied to pts_to_be_refined_count to get the *real* number of
        // points to be refined
        // (used like 1, 1/2, 1/4, 1/8, ... at the end of the algorithm to fill the grid size up to N)
        double refine_factor = 1.0;
        // success or not?
        bool result = true;

        while (true) {
          // status printing
          {
            char str[10];
            snprintf(str, 10, "%.1f%%",
                     static_cast<double>(current_N) / static_cast<double>(N) * 100.0);
            tools::printer.printStatusUpdate(std::string(str) +
                                             " (N = " + toString(current_N) + ")");
          }

          // calculate number of points to be refined
          refinable_pts_count = abstract_refinement->getNumberOfRefinablePoints(grid_storage);
          pts_to_be_refined_count = static_cast<int>(
                                      1.0 + refine_factor * alpha * static_cast<double>(refinable_pts_count));

          // refine
          base::SurplusRefinementFunctor refine_func(&coeffs, pts_to_be_refined_count);
          abstract_refinement->free_refine(grid_storage, &refine_func);

          // new grid size
          size_t new_N = grid_storage->size();

          if (new_N == current_N) {
            // size unchanged ==> nothing refined (should not happen)
            tools::printer.printStatusEnd(
              "error: size unchanged in IterativeGridGeneratorLinearSurplus");
            result = false;
            break;
          }

          if (new_N > N) {
            // too many new points ==> undo refinement and try again with refine_factor halved
            std::list<size_t> indices_to_remove;

            for (size_t i = current_N; i < new_N; i++) {
              indices_to_remove.push_back(i);
            }

            grid_storage->deletePoints(indices_to_remove);

            if (pts_to_be_refined_count == 1) {
              break;
            } else {
              refine_factor /= 2.0;
              k++;
              continue;
            }
          }

          coeffs.resize(new_N);

          // parallel evaluation of f in the new grid points
          #pragma omp parallel shared(fX, current_N, new_N, grid_storage, distr) default(none)
          {
            std::vector<double> x(d, 0.0);
            tools::SmartPointer<function::Objective> cur_f(f.clone());

            #pragma omp for
            for (size_t i = current_N; i < new_N; i++) {
              // convert grid point to coordinate vector
              #pragma omp critical
              {
                base::GridIndex* gp = grid_storage->get(i);
                // set point distribution accordingly to normal/Clenshaw-Curtis grids
                gp->setPointDistribution(distr);

                for (size_t t = 0; t < d; t++) {
                  x[t] = gp->abs(t);
                }
              }

              double fx = cur_f->eval(x);

              #pragma omp critical
              fX[i] = fx;
            }
          }

          // forward substitution (hier_system should always be a lower triangular matrix)
          for (size_t i = current_N; i < new_N; i++) {
            coeffs[i] = fX[i];

            for (size_t j = 0; j < i; j++) {
              coeffs[i] -= hier_system.getMatrixEntry(i, j) * coeffs[j];
            }
          }

          // next round
          current_N = new_N;
          k++;
        }

        // delete superfluous entries in fX
        fX.erase(fX.begin() + current_N, fX.end());

        if (result) {
          tools::printer.printStatusUpdate("100.0% (N = " + toString(current_N) + ")");
          tools::printer.printStatusEnd();
          return true;
        } else {
          return false;
        }
      }

    }
  }
}
