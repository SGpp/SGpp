// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::GridGenerator;
using sgpp::base::GridPoint;
using sgpp::base::GridStorage;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::MakePositiveCandidateSearchAlgorithm;

double normal(DataVector& input) {
  double result = 1.;

  double mean = 0.5;
  double sigma = 0.02;
  size_t numDims = input.getSize();

  double norm = 1. / std::sqrt(std::pow(2 * M_PI, numDims) * std::pow(sigma, numDims));

  for (size_t i = 0; i < numDims; i++) {
    double x = (input[i] - mean);
    result *= std::exp(-0.5 * (x * x) / sigma);
  }
  return norm * result;
}

void testMakePositive(Grid& grid, size_t numDims, size_t level, size_t refnums,
                      MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm,
                      double (*f)(DataVector&), double tol = -1e-12) {
  // -------------------------------------------------------------------------------------------
  // interpolate the pdf
  // create a two-dimensional piecewise bilinear grid
  GridStorage& gridStorage = grid.getStorage();
  grid.getGenerator().regular(level);
  std::cout << "========================================================" << std::endl;
  std::cout << "dimensionality               : " << gridStorage.getDimension() << std::endl;
  std::cout << "level                        : " << level << std::endl;
  std::cout << "number of initial grid points: " << gridStorage.getSize() << std::endl;
  std::cout << "--------------------------------------------------------" << std::endl;

  // create coefficient vector
  DataVector alpha(gridStorage.getSize());
  alpha.setAll(0.0);

  DataVector x(gridStorage.getDimension());
  // set function values in alpha
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    gridStorage.getCoordinates(gridStorage.getPoint(i), x);
    alpha[i] = f(x);
  }

  // hierarchize
  sgpp::op_factory::createOperationHierarchisation(grid)->doHierarchisation(alpha);

  // refine adaptively 5 times
  for (size_t step = 0; step < refnums; step++) {
    // refine a single grid point each time
    SurplusRefinementFunctor functor(alpha, 1);
    grid.getGenerator().refine(functor);
    std::cout << "refinement step " << step + 1 << ", new grid size: " << alpha.getSize()
              << std::endl;

    alpha.resize(gridStorage.getSize());
    // extend alpha vector (new entries uninitialized)

    // set function values in alpha
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
      gridStorage.getCoordinates(gridStorage.getPoint(i), x);
      alpha[i] = f(x);
    }

    // hierarchize
    sgpp::op_factory::createOperationHierarchisation(grid)->doHierarchisation(alpha);
  }
  // -------------------------------------------------------------------------------------------
  // force the function to be positive
  Grid* positiveGrid = nullptr;
  DataVector positiveAlpha(alpha);
  sgpp::op_factory::createOperationMakePositive(grid, candidateSearchAlgorithm)
      ->makePositive(positiveGrid, positiveAlpha);

  std::cout << "(" << gridStorage.getDimension() << "," << level
            << ") : #grid points = " << grid.getSize() << " -> " << positiveGrid->getSize()
            << std::endl;

  // make sure that the sparse grid function is really positive
  auto fullGrid = sgpp::base::Grid::createLinearGrid(numDims);
  fullGrid->getGenerator().full(level);

  DataMatrix coordinates;
  DataVector nodalValues(fullGrid->getStorage().getSize());
  fullGrid->getStorage().getCoordinateArraysForEval(coordinates);
  sgpp::op_factory::createOperationMultipleEval(*positiveGrid, coordinates)
      ->mult(positiveAlpha, nodalValues);

  for (size_t i = 0; i < nodalValues.getSize(); ++i) {
    BOOST_CHECK_GE(nodalValues[i], tol);
  }

  delete positiveGrid;
}

BOOST_AUTO_TEST_SUITE(TestOperationMakePositive)

BOOST_AUTO_TEST_CASE(testOperationMakePositiveFullGridSearch) {
  // parameters
  size_t numDims = 4;
  size_t level = 4;
  size_t refnums = 0;

  // interpolate the normal pdf
  for (size_t idim = numDims; idim <= numDims; idim++) {
    for (size_t ilevel = level; ilevel <= level; ilevel++) {
      std::unique_ptr<Grid> grid = Grid::createLinearGrid(idim);
      testMakePositive(*grid, idim, ilevel, refnums, MakePositiveCandidateSearchAlgorithm::FullGrid,
                       &normal);
    }
  }
}

BOOST_AUTO_TEST_CASE(testOperationMakePositiveIntersections) {
  // parameters
  size_t numDims = 4;
  size_t level = 4;
  size_t refnums = 0;

  // interpolate the normal pdf
  for (size_t idim = numDims; idim <= numDims; idim++) {
    for (size_t ilevel = level; ilevel <= level; ilevel++) {
      std::unique_ptr<Grid> grid = Grid::createLinearGrid(idim);
      testMakePositive(*grid, idim, ilevel, refnums,
                       MakePositiveCandidateSearchAlgorithm::Intersections, &normal);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
