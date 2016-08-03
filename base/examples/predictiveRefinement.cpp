// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
   * \page example_predictiveRefinement_cpp Spatially-Dimension-Adaptive Refinement in C++
   *
   */

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>

#include <iostream>
#include <cmath>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;
using sgpp::base::HashRefinement;
using sgpp::base::OperationEval;
using sgpp::base::PredictiveRefinement;
using sgpp::base::PredictiveRefinementIndicator;

// function to interpolate
double f(double x0, double x1) {
  return sin(x0 * M_PI);
}

DataVector& calculateError(const DataMatrix& dataSet, Grid& grid,
                           const DataVector& alpha, DataVector& error) {
  std::cout << "calculating error" << std::endl;

  // traverse dataSet
  DataVector vec(2);
  std::unique_ptr<OperationEval> opEval(sgpp::op_factory::createOperationEval(grid));

  for (unsigned int i = 0; i < dataSet.getNrows(); i++) {
    dataSet.getRow(i, vec);
    error[i] = pow(f(dataSet.get(i, 0), dataSet.get(i, 1)) - opEval->eval(alpha,
                   vec), 2);
  }

  return error;
}

int main() {
  // create a two-dimensional piecewise bilinear grid
  size_t dim = 2;
  std::unique_ptr<Grid> grid(Grid::createModLinearGrid(dim));
  GridStorage& gridStorage = grid->getStorage();
  std::cout << "dimensionality:                   " << gridStorage.getDimension() << std::endl;

  // create regular grid, level 3
  size_t level = 1;
  grid->getGenerator().regular(level);
  std::cout << "number of initial grid points:    " << gridStorage.getSize() << std::endl;

  // create coefficient vector
  DataVector alpha(gridStorage.getSize());
  alpha.setAll(0.0);
  std::cout << "length of alpha vector:           " << alpha.getSize() << std::endl;

  int rows = 100;
  int cols = 100;

  DataMatrix dataSet(rows * cols, dim);
  DataVector vals(rows * cols);

  // Create a "List" of points where the error should be calculated.
  // This represents a regular 2d grid with a step size of 1 / rows and 1 / cols.
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
     * We refine adaptively 20 times
     */ 
  for (int step = 0; step < 20; step++) {

    /**
       * Now we need to calculate the surplus vector alpha. In data
       * mining with do it by solving a regression problem as shown in
       * example \ref example_classificationRefinementExample_cpp.
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
    sgpp::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(
      alpha);

    // calculate squared offset
    DataVector errorVector(dataSet.getNrows());
    calculateError(dataSet, *grid, alpha, errorVector);

    // refinement  stuff
    HashRefinement refinement;
    PredictiveRefinement decorator(&refinement);

    // refine a single grid point each time
    std::cout << "Error over all = "  << errorVector.sum() << std::endl;
    PredictiveRefinementIndicator indicator(*grid, dataSet, errorVector, 1);
    decorator.free_refine(gridStorage, indicator);

    std::cout << "Refinement step " << step + 1 << ", new grid size: " <<
         gridStorage.getSize() << std::endl;

    // plot grid

    // extend alpha vector (new entries uninitialized)
    alpha.resize(gridStorage.getSize());
  }
}
