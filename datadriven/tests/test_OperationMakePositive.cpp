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
#include <sgpp/datadriven/operation/hash/simple/OperationMakePositive.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <vector>
#include <list>
#include <utility>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::GridGenerator;
using sgpp::base::GridPoint;
using sgpp::base::GridStorage;
using sgpp::base::HashGridPoint;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::datadriven::MakePositiveCandidateSearchAlgorithm;
using sgpp::datadriven::MakePositiveInterpolationAlgorithm;

BOOST_AUTO_TEST_SUITE(TestOperationMakePositive)

double normal(DataVector& input) {
  double result = 1.;

  double mean = 0.5;
  double sigma = 0.02;
  size_t numDims = input.getSize();

  double norm =
      1. / std::sqrt(std::pow(2 * M_PI, numDims) * std::pow(sigma, numDims));

  for (size_t i = 0; i < numDims; i++) {
    double x = (input[i] - mean);
    result *= std::exp(-0.5 * (x * x) / sigma);
  }
  return norm * result;
}

double sin(DataVector& x) {
  double result = 0;
  for (size_t i = 0; i < x.getSize(); i++) {
    result += std::sin(M_PI * (2 * x[i] - 1.0)) + 0.8;
  }
  return result;
}

void testMakePositive(
    Grid& grid, size_t numDims, size_t level, size_t refnums,
    MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm,
    double (*f)(DataVector&), bool generateConsistentGrid = true,
    double tol = -1e-12, bool verbose = false) {
  // -------------------------------------------------------------------------------------------
  // interpolate the pdf
  // create a two-dimensional piecewise bilinear grid
  GridStorage& gridStorage = grid.getStorage();
  grid.getGenerator().regular(level);
  if (verbose) {
    std::cout << "========================================================"
              << std::endl;
    std::cout << "dimensionality               : " << gridStorage.getDimension()
              << std::endl;
    std::cout << "level                        : " << level << std::endl;
    std::cout << "number of initial grid points: " << gridStorage.getSize()
              << std::endl;
    std::cout << "--------------------------------------------------------"
              << std::endl;
  }

  // create coefficient vector
  DataVector alpha(gridStorage.getSize());

  DataVector x(gridStorage.getDimension());
  // set function values in alpha
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    gridStorage.getCoordinates(gridStorage.getPoint(i), x);
    alpha[i] = f(x);
  }

  // hierarchize
  std::unique_ptr<sgpp::base::OperationHierarchisation> opHier(
      sgpp::op_factory::createOperationHierarchisation(grid));
  opHier->doHierarchisation(alpha);

  // refine adaptively 5 times
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
    std::unique_ptr<sgpp::base::OperationHierarchisation> opHier(
        sgpp::op_factory::createOperationHierarchisation(grid));
    opHier->doHierarchisation(alpha);
    if (verbose) {
      std::cout << "refinement step " << step + 1
                << ", new grid size: " << alpha.getSize() << std::endl;
    }
  }

  size_t numFullGridPoints = static_cast<size_t>(
      std::pow(std::pow(2, gridStorage.getMaxLevel()) - 1, numDims));
  size_t maxLevel = gridStorage.getMaxLevel();
  if (verbose) {
    if (refnums > 0) {
      std::cout << "--------------------------------------------------------"
                << std::endl;
      std::cout << "level after refinement       : "
                << gridStorage.getMaxLevel() << std::endl;
      std::cout << "grid size after refinement   : " << gridStorage.getSize()
                << std::endl;
      std::cout << "--------------------------------------------------------"
                << std::endl;
    }

    std::cout << "number of full grid points   : " << numFullGridPoints
              << " (maxLevel = " << gridStorage.getMaxLevel() << ")"
              << std::endl;
    std::cout << "--------------------------------------------------------"
              << std::endl;
  }
  // -------------------------------------------------------------------------------------------
  // force the function to be positive
  std::unique_ptr<sgpp::datadriven::OperationMakePositive> opMakePositive(
      sgpp::op_factory::createOperationMakePositive(
          candidateSearchAlgorithm,
          MakePositiveInterpolationAlgorithm::SetToZero,
          generateConsistentGrid));
  std::unique_ptr<Grid> positiveGrid(grid.clone());
  DataVector positiveAlpha(alpha);
  opMakePositive->makePositive(*positiveGrid, positiveAlpha);

  if (verbose) {
    std::cout << "(" << gridStorage.getDimension() << ", " << level << ", "
              << refnums << ") : #grid points = " << grid.getSize() << " -> "
              << positiveGrid->getSize() << " < " << numFullGridPoints
              << std::endl;
  }

  // --------------------------------------------------------------------------------------------
  // make sure that the sparse grid function is really positive
  if (verbose) {
    std::cout << "check full grid for success: ";
  }
  std::unique_ptr<sgpp::base::Grid> fullGrid(
      sgpp::base::Grid::createLinearGrid(numDims));
  fullGrid->getGenerator().full(maxLevel);

  if (verbose) {
    std::cout << fullGrid->getSize() << std::endl;
    std::cout << "evaluate at all grid points" << std::endl;
  }
  DataMatrix coordinates;
  DataVector nodalValues(fullGrid->getStorage().getSize());
  fullGrid->getStorage().getCoordinateArrays(coordinates);

  // use streaming eval to support inconsistent grids
  sgpp::datadriven::OperationMultipleEvalConfiguration evalConfig(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT);
  std::unique_ptr<sgpp::base::OperationMultipleEval> opEval(
      sgpp::op_factory::createOperationMultipleEval(*positiveGrid, coordinates,
                                                    evalConfig));
  opEval->mult(positiveAlpha, nodalValues);

  if (verbose) {
    std::cout << "check for positivity: ";
  }
  for (size_t i = 0; i < nodalValues.getSize(); ++i) {
    BOOST_CHECK_GE(nodalValues[i], tol);
  }

  if (verbose) {
    std::cout << "Done" << std::endl;
  }

  // --------------------------------------------------------------------------------------------
  // make sure that the sparse grid is really minimal
  std::vector<size_t> addedGridPointIndices =
      opMakePositive->getAddedGridPointsForPositivity();
  std::vector<std::pair<std::shared_ptr<HashGridPoint>, double>>
      addedGridPoints;

  std::unique_ptr<Grid> coarsedGrid(positiveGrid->clone());
  GridStorage& coarsedGridStorage = coarsedGrid->getStorage();
  for (auto& ix : addedGridPointIndices) {
    std::shared_ptr<HashGridPoint> gp(
        new HashGridPoint(coarsedGridStorage.getPoint(ix)));
    addedGridPoints.push_back(std::make_pair(gp, positiveAlpha[ix]));
  }

  // remove each of the added grid points and check if the function value at
  // that point is negative
  std::list<size_t> toBeRemoved;
  DataVector coarsedAlpha(positiveAlpha);
  for (auto& values : addedGridPoints) {
    // load current point
    std::shared_ptr<HashGridPoint> gp = values.first;
    double coefficient = values.second;

    // delete grid point
    toBeRemoved.push_back(coarsedGridStorage.getSequenceNumber(*gp));
    auto remainingIndex = coarsedGridStorage.deletePoints(toBeRemoved);
    // restructure the coefficient vector
    coarsedAlpha.restructure(remainingIndex);

    // evaluate the function at this very point
    gp->getStandardCoordinates(x);
    std::unique_ptr<sgpp::base::OperationEval> opEvalNaive(
        sgpp::op_factory::createOperationEvalNaive(*coarsedGrid));
    double value = opEvalNaive->eval(coarsedAlpha, x);

    // make sure that it is smaller than zero
    BOOST_CHECK_LE(value, 0.0);

    // insert grid point again
    size_t ix = coarsedGridStorage.insert(*gp);

    // restructure alpha
    coarsedAlpha.insert(ix, coefficient);

    // clear helper list
    toBeRemoved.clear();
  }
}

BOOST_AUTO_TEST_CASE(testOperationMakePositiveConsistent) {
  // parameters
  size_t numDims = 3;
  size_t level = 3;
  size_t refIterations = 2;
  size_t refnums = 2;
  double tol = -1e-12;
  bool verbose = false;

  std::unique_ptr<Grid> gridIntersections;
  std::unique_ptr<Grid> gridIntersectionsJoin;
  std::unique_ptr<Grid> gridFull;

  for (size_t idim = 2; idim <= numDims; idim++) {
    for (size_t ilevel = 2; ilevel <= level; ilevel++) {
      for (size_t irefIteration = 0; irefIteration <= refIterations;
           irefIteration++) {
        // check whether both candidate search algorithms lead to the same
        // consistent grid
        gridIntersections.reset(Grid::createLinearGrid(idim));
        testMakePositive(*gridIntersections, idim, ilevel,
                         irefIteration * refnums,
                         MakePositiveCandidateSearchAlgorithm::Intersections,
                         &normal, true, tol, verbose);
        gridIntersectionsJoin.reset(Grid::createLinearGrid(idim));
        testMakePositive(
            *gridIntersectionsJoin, idim, ilevel, irefIteration * refnums,
            MakePositiveCandidateSearchAlgorithm::IntersectionsJoin, &normal,
            true, tol, verbose);
        if (irefIteration == 0) {
          gridFull.reset(Grid::createLinearGrid(idim));
          testMakePositive(*gridFull, idim, ilevel, irefIteration * refnums,
                           MakePositiveCandidateSearchAlgorithm::FullGrid,
                           &normal, true, tol, verbose);
          BOOST_CHECK_EQUAL(gridIntersections->getSize(), gridFull->getSize());
          BOOST_CHECK_EQUAL(gridIntersectionsJoin->getSize(),
                            gridFull->getSize());
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(testOperationMakePositiveInconsistent) {
  // parameters
  size_t numDims = 3;
  size_t level = 3;
  size_t refIterations = 2;
  size_t refnums = 2;
  double tol = -1e-12;
  bool verbose = false;

  std::unique_ptr<Grid> gridIntersections;
  std::unique_ptr<Grid> gridIntersectionsJoin;
  std::unique_ptr<Grid> gridFull;

  for (size_t idim = 2; idim <= numDims; idim++) {
    for (size_t ilevel = 2; ilevel <= level; ilevel++) {
      for (size_t irefIteration = 0; irefIteration <= refIterations;
           irefIteration++) {
        // check whether the candidate search algorithms lead to the same
        // inconsistent grid
        gridIntersections.reset(Grid::createLinearGrid(idim));
        testMakePositive(*gridIntersections, idim, ilevel,
                         irefIteration * refnums,
                         MakePositiveCandidateSearchAlgorithm::Intersections,
                         &normal, false, tol, verbose);
        gridIntersectionsJoin.reset(Grid::createLinearGrid(idim));
        testMakePositive(
            *gridIntersectionsJoin, idim, ilevel, irefIteration * refnums,
            MakePositiveCandidateSearchAlgorithm::IntersectionsJoin, &normal,
            true, tol, verbose);
        if (irefIteration == 0) {
          gridFull.reset(Grid::createLinearGrid(idim));
          testMakePositive(*gridFull, idim, ilevel, irefIteration * refnums,
                           MakePositiveCandidateSearchAlgorithm::FullGrid,
                           &normal, false, tol, verbose);
          BOOST_CHECK_EQUAL(gridIntersections->getSize(), gridFull->getSize());
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
