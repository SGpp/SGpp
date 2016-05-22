// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/common/DirichletGridConverter.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedGrid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <vector>

using sgpp::base::DataVector;
using sgpp::base::DirichletGridConverter;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridPoint;
using sgpp::base::GridStorage;
using sgpp::base::LinearGrid;

BOOST_AUTO_TEST_CASE(test_DirichletGridConverter) {
  /// number of the boundary grid's grid points
  size_t numTotalGridPoints;
  /// number of the inner grid's grid points
  size_t numInnerGridPoints;
  /**
  * array to store the position of i-th inner point in
  * the boundary grid's coefficients
  */
  size_t* conCoefArray;
  Grid* innerGridExact;

  int level = 3;
  int dimension = 2;
  std::unique_ptr<Grid> linearBoundaryGrid = Grid::createLinearBoundaryGrid(dimension);
  GridStorage& linearBoundaryGridStorageExact = linearBoundaryGrid->getStorage();
  linearBoundaryGrid->getGenerator().regular(level);

  // determine the number of grid points for both grids
  numTotalGridPoints = linearBoundaryGridStorageExact.getSize();
  numInnerGridPoints = linearBoundaryGridStorageExact.getNumberOfInnerPoints();

  DataVector boundaryGridCoeffs(numTotalGridPoints);
  for (size_t i = 0; i < numTotalGridPoints; ++i) {
    GridPoint& gp = linearBoundaryGridStorageExact.getGridPoint(i);
    boundaryGridCoeffs[i] =
        gp.getCoord(0) * (gp.getCoord(0)) + gp.getCoord(1) * (gp.getCoord(1));
  }

  sgpp::op_factory::createOperationHierarchisation(*linearBoundaryGrid)
      ->doHierarchisation(boundaryGridCoeffs);

  // allocate the translation array for the coefficients
  conCoefArray = new size_t[numInnerGridPoints];

  // Get the algorithmic dimensions
  std::vector<size_t> BSalgoDims = linearBoundaryGrid->getAlgorithmicDimensions();

  // create new inner Grid, with one grid point
  innerGridExact = new LinearGrid(linearBoundaryGrid->getBoundingBox());

  // Set algorithmic dimensions for inner Grid
  innerGridExact->setAlgorithmicDimensions(BSalgoDims);

  // create new DataVector for storing the inner grid's coefficients
  DataVector innerCoeffsExact(numInnerGridPoints);

  // Iterate through all grid points and filter inner points
  size_t numInner = 0;

  for (size_t i = 0; i < numTotalGridPoints; ++i) {
    GridPoint& curPoint = linearBoundaryGridStorageExact.getGridPoint(i);
    if (curPoint.isInnerPoint() == true) {
      // handle coefficients
      conCoefArray[numInner] = i;
      innerCoeffsExact.set(numInner, boundaryGridCoeffs.get(i));
      numInner++;
      // insert point into inner grid
      innerGridExact->getStorage().insert(curPoint);
    }
  }

  // Test BUILD
  DirichletGridConverter dirichGridConverter;
  Grid* innerGridActual;
  DataVector* innerCoeffsActual;
  dirichGridConverter.buildInnerGridWithCoefs(*linearBoundaryGrid, boundaryGridCoeffs,
                                              &innerGridActual, &innerCoeffsActual);

  BOOST_CHECK_EQUAL(innerCoeffsActual->getSize(), innerCoeffsExact.getSize());

  for (size_t i = 0; i < innerCoeffsActual->getSize(); ++i) {
    BOOST_CHECK_EQUAL(innerCoeffsActual->get(i), innerCoeffsExact.get(i));
  }

  int innerPointIndex = 0;
  {
    GridStorage& linearBoundaryGridStorageActual = innerGridActual->getStorage();
    for (size_t i = 0; i < numTotalGridPoints; ++i) {
      GridPoint& curPointExact = linearBoundaryGridStorageExact.getGridPoint(i);
      if (curPointExact.isInnerPoint() == true) {
        GridPoint& curPointActual = linearBoundaryGridStorageActual.getGridPoint(innerPointIndex);
        for (size_t curDim = 0; curDim < curPointExact.getDimension(); ++curDim) {
          BOOST_CHECK_EQUAL(curPointActual.getCoord(curDim), curPointExact.getCoord(curDim));
        }
        innerPointIndex++;
      }
    }
  }

  // TEST REBUILD
  dirichGridConverter.rebuildInnerGridWithCoefs(*linearBoundaryGrid, boundaryGridCoeffs,
                                                &innerGridActual, &innerCoeffsActual);

  BOOST_CHECK_EQUAL(innerCoeffsActual->getSize(), innerCoeffsExact.getSize());

  for (size_t i = 0; i < innerCoeffsActual->getSize(); ++i) {
    BOOST_CHECK_EQUAL(innerCoeffsActual->get(i), innerCoeffsExact.get(i));
  }

  innerPointIndex = 0;
  {
    GridStorage& linearBoundaryGridStorageActual = innerGridActual->getStorage();
    for (size_t i = 0; i < numTotalGridPoints; ++i) {
      GridPoint& curPointExact = linearBoundaryGridStorageExact.getGridPoint(i);
      if (curPointExact.isInnerPoint() == true) {
        GridPoint& curPointActual = linearBoundaryGridStorageActual.getGridPoint(innerPointIndex);
        for (size_t curDim = 0; curDim < curPointExact.getDimension(); ++curDim) {
          BOOST_CHECK_EQUAL(curPointActual.getCoord(curDim), curPointExact.getCoord(curDim));
        }
        innerPointIndex++;
      }
    }
  }

  // TEST calcInnerCoefs
  dirichGridConverter.calcInnerCoefs(boundaryGridCoeffs, *innerCoeffsActual);
  for (size_t i = 0; i < innerCoeffsActual->getSize(); ++i) {
    BOOST_CHECK_EQUAL(innerCoeffsActual->get(i), innerCoeffsExact.get(i));
  }

  // TEST updateBoundaryCoefs
  DataVector boundaryGridCoeffsActual(boundaryGridCoeffs);
  dirichGridConverter.updateBoundaryCoefs(boundaryGridCoeffsActual, *innerCoeffsActual);
  for (size_t i = 0; i < innerCoeffsActual->getSize(); ++i) {
    BOOST_CHECK_EQUAL(boundaryGridCoeffsActual.get(i), boundaryGridCoeffs.get(i));
  }

  delete conCoefArray;
  delete innerGridExact;
}
