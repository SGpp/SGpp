// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
   * \page example_predictiveRefinement_cpp Spatially-Dimension-Adaptive Refinement in C++
   *
   *
   * We compute the sparse grid interpolant of the function \f$ f(x) =
   * \sin(\pi x).\f$ We perform spatially-dimension-adaptive
   * refinement of the sparse grid model, which means we refine a
   * particular grid point (locality) only in some dimensions
   * (dimensionality).
   *
   * For details on spatially-dimension-adaptive refinement see
   * \verbatim
   *  V. Khakhutskyy and M. Hegland: Spatially-Dimension-Adaptive Sparse Grids for Online Learning.
In D. Pflüger and J. Garcke (ed.), Sparse Grids and Applications - Stuttgart 2014, Volume 109 of
LNCSE, p. 133–162. Springer International Publishing, March 2016.
   * \endverbatim
   *
   *
   *
   * The example can be found in the file `predictiveRefinement.cpp`.
   */

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>

#include <cmath>
#include <iostream>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;
using sgpp::base::HashRefinement;
using sgpp::base::OperationEval;
using sgpp::base::PredictiveRefinement;
using sgpp::base::PredictiveRefinementIndicator;
using sgpp::base::OperationHierarchisation;

/**
   * We define the function \f$ f(x) =
   * \sin(\pi x).\f$ to interpolate.
   **/

double f(double x0, double x1) { return sin(x0 * M_PI); }

/**
   * Spatially-dimension-adaptive refinement uses squared prediction
   * error on a dataset to compute refinement indicators. Hence, here
   * we define a function to compute these squared errors.
   *
   */
DataVector& calculateError(const DataMatrix& dataSet, Grid& grid, const DataVector& alpha,
                           DataVector& error) {
  std::cout << "calculating error" << std::endl;

  // traverse dataSet
  DataVector vec(2);
  std::unique_ptr<OperationEval> opEval(sgpp::op_factory::createOperationEval(grid));

  for (unsigned int i = 0; i < dataSet.getNrows(); i++) {
    dataSet.getRow(i, vec);
    error[i] = pow(f(dataSet.get(i, 0), dataSet.get(i, 1)) - opEval->eval(alpha, vec), 2);
  }

  return error;
}

/**
   * Start with the main function
   *
   */

int main() {
  /**
     * create a two-dimensional piecewise bilinear grid
     */
  size_t dim = 2;
  std::unique_ptr<Grid> grid(Grid::createModLinearGrid(dim));
  GridStorage& gridStorage = grid->getStorage();
  std::cout << "dimensionality:                   " << gridStorage.getDimension() << std::endl;

  // create regular grid, level 3
  size_t level = 1;
  grid->getGenerator().regular(level);
  std::cout << "number of initial grid points:    " << gridStorage.getSize() << std::endl;

  /**
     * To create a dataset we use points on a regular 2d grid with a
     * step size of 1 / rows and 1 / cols.
     */

  int rows = 100;
  int cols = 100;

  DataMatrix dataSet(rows * cols, dim);
  DataVector vals(rows * cols);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      // xcoord
      dataSet.set(i * cols + j, 0, i * 1.0 / rows);
      // ycoord
      dataSet.set(i * cols + j, 1, j * 1.0 / cols);
      vals[i * cols + j] = f(i * 1.0 / rows, j * 1.0 / cols);
    }
  }

  /**
   * We refine adaptively 20 times. In every step we recompute the
   * vector of surpluses `alpha`, the vector with squared errors on
   * the dataset `errorVector`, and then call the refinement
   * routines.
   */

  // create coefficient vector
  DataVector alpha(gridStorage.getSize());
  alpha.setAll(0.0);
  std::cout << "length of alpha vector:           " << alpha.getSize() << std::endl;

  for (int step = 0; step < 20; step++) {
    /**
       * Step 1: calculate the surplus vector alpha. In data
       * mining with do it by solving a regression problem.
       * Here, the function can be evaluated at any point. Hence. we
       * simply evaluate it at the coordinates of the grid points to
       * obtain the nodal values. Then we use hierarchization to
       * obtain the surplus value.
       *
       */

    // set function values in alpha
    DataVector gridPointCoordinates(dim);
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
      gridStorage.getPoint(i).getStandardCoordinates(gridPointCoordinates);
      alpha[i] = f(gridPointCoordinates[0], gridPointCoordinates[1]);
    }

    // hierarchize
    std::unique_ptr<OperationHierarchisation>(
        sgpp::op_factory::createOperationHierarchisation(*grid))
        ->doHierarchisation(alpha);

    /**
       * Step 2: calculate squared errors.
       */
    DataVector errorVector(dataSet.getNrows());
    calculateError(dataSet, *grid, alpha, errorVector);

    /**
       * Step 3: call refinement routines. `PredictiveRefinement`
       * implements the decorator pattern and extends the
       * functionality of `HashRefinement`. `PredictiveRefinement`
       * requires a special kind of refinement functor --
       * `PredictiveRefinementIndicator` that can access the dataset
       * and the error vector. The refinement itself if performed by
       * calling `.free_refine()` same for normal refinement in
       * `HashRefinement`.
       *
       */
    // refinement  stuff
    HashRefinement refinement;
    PredictiveRefinement decorator(&refinement);

    // refine a single grid point each time
    std::cout << "Error over all = " << errorVector.sum() << std::endl;
    PredictiveRefinementIndicator indicator(*grid, dataSet, errorVector, 1);
    decorator.free_refine(gridStorage, indicator);

    std::cout << "Refinement step " << step + 1 << ", new grid size: " << gridStorage.getSize()
              << std::endl;

    // extend alpha vector (new entries uninitialized)
    alpha.resize(gridStorage.getSize());
  }
}

/**
   * The output of the program should look like this
   *
   * \verbatim
   * dimensionality:                   2
number of initial grid points:    1
length of alpha vector:           1
calculating error
Error over all = 2268.65
Refinement step 1, new grid size: 3
calculating error
Error over all = 264.09
Refinement step 2, new grid size: 5
calculating error
Error over all = 125.378
Refinement step 3, new grid size: 7
calculating error
Error over all = 3.48359
Refinement step 4, new grid size: 9
calculating error
Error over all = 1.99757
Refinement step 5, new grid size: 11
calculating error
Error over all = 0.845349
Refinement step 6, new grid size: 13
calculating error
Error over all = 0.464096
Refinement step 7, new grid size: 15
calculating error
Error over all = 0.0828432
Refinement step 8, new grid size: 17
calculating error
Error over all = 0.0828432
Refinement step 9, new grid size: 19
calculating error
Error over all = 0.068976
Refinement step 10, new grid size: 21
calculating error
Error over all = 0.0551672
Refinement step 11, new grid size: 23
calculating error
Error over all = 0.0413584
Refinement step 12, new grid size: 25
calculating error
Error over all = 0.0330229
Refinement step 13, new grid size: 27
calculating error
Error over all = 0.0230578
Refinement step 14, new grid size: 29
calculating error
Error over all = 0.0130926
Refinement step 15, new grid size: 31
calculating error
Error over all = 0.00856834
Refinement step 16, new grid size: 33
calculating error
Error over all = 0.00404405
Refinement step 17, new grid size: 35
calculating error
Error over all = 0.00404405
Refinement step 18, new grid size: 37
calculating error
Error over all = 0.00404405
Refinement step 19, new grid size: 41
calculating error
Error over all = 0.00404405
Refinement step 20, new grid size: 45

   * \endverbatim
   *
   */
