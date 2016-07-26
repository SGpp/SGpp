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

#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::GridGenerator;
using sgpp::base::GridPoint;
using sgpp::base::GridStorage;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::MakePositiveCandidateSearchAlgorithm;

BOOST_AUTO_TEST_SUITE(TestOperationLimitFunctionValueRange)

double sin(DataVector& x) {
  double result = 0;
  for (size_t i = 0; i < x.getSize(); i++) {
    result += std::sin(M_PI * (2 * x[i] - 1.0));
  }
  return result;
}

void testLimitFunctionValueRange(Grid& grid, size_t numDims, size_t level, size_t refnums,
                                 size_t side,
                                 MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm,
                                 double (*f)(DataVector&), double ylower, double yupper,
                                 double tol = 1e-12, bool verbose = true) {
  // -------------------------------------------------------------------------------------------
  // interpolate the pdf
  // create a two-dimensional piecewise bilinear grid
  GridStorage& gridStorage = grid.getStorage();
  grid.getGenerator().regular(level);
  if (verbose) {
    std::cout << "========================================================" << std::endl;
    std::cout << "dimensionality               : " << gridStorage.getDimension() << std::endl;
    std::cout << "level                        : " << level << std::endl;
    std::cout << "number of initial grid points: " << gridStorage.getSize() << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
  }

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

  // refine adaptively
  for (size_t step = 0; step < refnums; step++) {
    // refine a single grid point each time
    SurplusRefinementFunctor functor(alpha, 1);
    grid.getGenerator().refine(functor);
    alpha.resize(gridStorage.getSize());
    // extend alpha vector (new entries uninitialized)

    // set function values in alpha
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
      gridStorage.getCoordinates(gridStorage.getPoint(i), x);
      alpha[i] = f(x);
    }

    // hierarchize
    sgpp::op_factory::createOperationHierarchisation(grid)->doHierarchisation(alpha);
    if (verbose) {
      std::cout << "refinement step " << step + 1 << ", new grid size: " << alpha.getSize()
                << std::endl;
    }
  }

  size_t numFullGridPoints =
      static_cast<size_t>(std::pow(std::pow(2, gridStorage.getMaxLevel()) - 1, numDims));
  size_t maxLevel = gridStorage.getMaxLevel();
  if (verbose) {
    if (refnums > 0) {
      std::cout << "--------------------------------------------------------" << std::endl;
      std::cout << "level after refinement       : " << gridStorage.getMaxLevel() << std::endl;
      std::cout << "grid size after refinement   : " << gridStorage.getSize() << std::endl;
    }

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "number of full grid points   : " << numFullGridPoints
              << " (maxLevel = " << gridStorage.getMaxLevel() << ")" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
  }
  // -------------------------------------------------------------------------------------------
  // force the function to be positive
  Grid* positiveGrid = nullptr;
  DataVector positiveAlpha(alpha);
  auto opLimit =
      sgpp::op_factory::createOperationLimitFunctionValueRange(grid, candidateSearchAlgorithm);

  if (side == 0) {
    opLimit->doLowerLimitation(positiveGrid, positiveAlpha, ylower);
  } else if (side == 1) {
    opLimit->doUpperLimitation(positiveGrid, positiveAlpha, yupper);
  } else {  // both
    opLimit->doLimitation(positiveGrid, positiveAlpha, ylower, yupper);
  }

  if (verbose) {
    std::cout << "(" << gridStorage.getDimension() << "," << level
              << ") : #grid points = " << grid.getSize() << " -> " << positiveGrid->getSize()
              << " < " << numFullGridPoints << std::endl;
  }
  // make sure that the sparse grid function is really positive
  if (verbose) {
    std::cout << "check full grid for success: ";
  }
  auto fullGrid = sgpp::base::Grid::createLinearGrid(numDims);
  fullGrid->getGenerator().full(maxLevel);

  if (verbose) {
    std::cout << fullGrid->getSize() << std::endl;
    std::cout << "  evaluate at all grid points" << std::endl;
  }
  DataMatrix coordinates;
  DataVector nodalValues(fullGrid->getStorage().getSize());
  fullGrid->getStorage().getCoordinateArrays(coordinates);
  sgpp::op_factory::createOperationMultipleEval(*positiveGrid, coordinates)
      ->mult(positiveAlpha, nodalValues);
  if (verbose) {
    std::cout << "  check range: ";
  }
  for (size_t i = 0; i < nodalValues.getSize(); ++i) {
    if (side == 0 || side > 1) {
      BOOST_CHECK_GE(nodalValues[i], ylower - tol);
    }
    if (side == 1 || side > 1) {
      BOOST_CHECK_LE(nodalValues[i], yupper + tol);
    }
  }

  if (verbose) {
    std::cout << "Done" << std::endl;
  }
  delete positiveGrid;
}

BOOST_AUTO_TEST_CASE(testOperationLimitFunctionValueRangeLower) {
  // parameters
  size_t numDims = 2;
  size_t level = 2;
  size_t refnums = 0;

  for (size_t idim = numDims; idim <= numDims; idim++) {
    for (size_t ilevel = level; ilevel <= level; ilevel++) {
      std::unique_ptr<Grid> grid = Grid::createLinearGrid(idim);
      testLimitFunctionValueRange(*grid, idim, ilevel, refnums, 0,
                                  MakePositiveCandidateSearchAlgorithm::Intersections, &sin, -0.8,
                                  0.8);
    }
  }
}

BOOST_AUTO_TEST_CASE(testOperationLimitFunctionValueRangeUpper) {
  // parameters
  size_t numDims = 2;
  size_t level = 2;
  size_t refnums = 0;

  for (size_t idim = numDims; idim <= numDims; idim++) {
    for (size_t ilevel = level; ilevel <= level; ilevel++) {
      std::unique_ptr<Grid> grid = Grid::createLinearGrid(idim);
      testLimitFunctionValueRange(*grid, idim, ilevel, refnums, 1,
                                  MakePositiveCandidateSearchAlgorithm::Intersections, &sin, -0.8,
                                  0.8);
    }
  }
}

BOOST_AUTO_TEST_CASE(testOperationLimitFunctionValueRangeBothSides) {
  // parameters
  size_t numDims = 2;
  size_t level = 2;
  size_t refnums = 0;

  for (size_t idim = numDims; idim <= numDims; idim++) {
    for (size_t ilevel = level; ilevel <= level; ilevel++) {
      std::unique_ptr<Grid> grid = Grid::createLinearGrid(idim);
      testLimitFunctionValueRange(*grid, idim, ilevel, refnums, 2,
                                  MakePositiveCandidateSearchAlgorithm::Intersections, &sin, -0.8,
                                  0.8);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
