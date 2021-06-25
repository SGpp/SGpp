// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/application/SparseGridDensityEstimator.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/application/KernelDensityEstimator.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/optimization/operation/OptimizationOpFactory.hpp>

#include <vector>
#include <random>
#include <iostream>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;

using sgpp::datadriven::KernelDensityEstimator;
using sgpp::datadriven::KernelType;

double normal(DataVector& input, double mean = 0.5, double sigma = 0.1) {
  double result = 1.;
  double norm = 1. / (sigma * M_SQRT2PI);
  for (size_t i = 0; i < input.getSize(); i++) {
    result *= norm * std::exp((input[i] - mean) * (input[i] - mean) / (2. * sigma * sigma));
  }

  return result;
}

double parabola(DataVector& input) {
  double result = 1.;

  for (size_t i = 0; i < input.getSize(); i++) {
    result *= input[i] * (1. - input[i]) * 4.;
  }

  return result;
}

void hierarchize(Grid* grid, std::uint32_t level, DataVector& alpha, double (*func)(DataVector&)) {
  size_t dim = grid->getDimension();
  GridStorage& gs = grid->getStorage();
  gs.clear();
  grid->getGenerator().regular(level);
  alpha.resize(grid->getSize());

  DataVector coords(dim);
  for (size_t i = 0; i < grid->getSize(); i++) {
    gs.getPoint(i).getStandardCoordinates(coords);
    alpha[i] = func(coords);
  }
  if (grid->getType() == sgpp::base::GridType::Bspline ||
      grid->getType() == sgpp::base::GridType::ModBspline ||
      grid->getType() == sgpp::base::GridType::BsplineBoundary ||
      grid->getType() == sgpp::base::GridType::BsplineClenshawCurtis ||
      grid->getType() == sgpp::base::GridType::ModBsplineClenshawCurtis)
    sgpp::op_factory::createOperationMultipleHierarchisation(*grid)->doHierarchisation(alpha);
  else
    sgpp::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha);
}

void randu(DataVector& rvar, std::mt19937& generator) {
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (size_t j = 0; j < rvar.getSize(); ++j) {
    rvar[j] = distribution(generator);
  }
}

void randu(DataMatrix& rvar, std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t nsamples = rvar.getNrows(), ndim = rvar.getNcols();

  std::mt19937 generator(static_cast<std::mt19937::result_type>(seedValue));
  DataVector sample(ndim);
  for (size_t i = 0; i < nsamples; ++i) {
    randu(sample, generator);
    rvar.setRow(i, sample);
  }
}

void randn(DataVector& rvar, std::mt19937& generator) {
  std::normal_distribution<double> distribution(0.5, 0.1);
  for (size_t j = 0; j < rvar.getSize(); ++j) {
    rvar[j] = distribution(generator);
  }
}

void randn(DataMatrix& rvar, std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t nsamples = rvar.getNrows(), ndim = rvar.getNcols();

  std::mt19937 generator(static_cast<std::mt19937::result_type>(seedValue));
  DataVector sample(ndim);
  for (size_t i = 0; i < nsamples; ++i) {
    randn(sample, generator);
    rvar.setRow(i, sample);
  }
}

void testEqualityRosenblattInverseRosenblatt1D(
    Grid& grid, DataVector& alpha, size_t numSamples = 1000, double tolerance = 1e-12,
    std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t numDims = grid.getStorage().getDimension();
  DataVector u_vars(numSamples);
  std::mt19937 generator(static_cast<std::mt19937::result_type>(seedValue));

  // init samples to be transformed
  randu(u_vars, generator);

  std::unique_ptr<sgpp::datadriven::OperationTransformation1D> opInvRos(
      sgpp::op_factory::createOperationInverseRosenblattTransformation1D(grid));
  std::unique_ptr<sgpp::datadriven::OperationTransformation1D> opRos(
      sgpp::op_factory::createOperationRosenblattTransformation1D(grid));

  DataVector inversionErrors(numDims);
  for (size_t isample = 0; isample < numSamples; isample++) {
    // transform the u-space to x-space
    double x_var = opInvRos->doTransformation1D(&alpha, u_vars[isample]);
    // transform them back to the u-space
    double u_var_transformed = opRos->doTransformation1D(&alpha, x_var);

    // assert that x_vars and x_vars_transformed contain the same samples
    double inversionError = std::abs(u_vars[isample] - u_var_transformed) / u_vars[isample];
    BOOST_CHECK_SMALL(inversionError, tolerance);
  }
}

void testEqualityRosenblattInverseRosenblattDD(
    Grid& grid, DataVector& alpha, size_t numSamples = 1000, double tolerance = 1e-12,
    std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t numDims = grid.getStorage().getDimension();
  DataMatrix x_vars(numSamples, numDims);
  DataMatrix u_vars(numSamples, numDims);
  DataMatrix u_vars_transformed(numSamples, numDims);

  // init samples to be transformed
  randu(u_vars, seedValue);

  std::unique_ptr<sgpp::datadriven::OperationInverseRosenblattTransformation> opInvRos(
      sgpp::op_factory::createOperationInverseRosenblattTransformation(grid));
  std::unique_ptr<sgpp::datadriven::OperationRosenblattTransformation> opRos(
      sgpp::op_factory::createOperationRosenblattTransformation(grid));

  // transform the u-space to x-space
  opInvRos->doTransformation(&alpha, &u_vars, &x_vars);
  // transform them back to the u-space
  opRos->doTransformation(&alpha, &x_vars, &u_vars_transformed);

  DataVector u_sample(numDims);
  DataVector x_sample(numDims);
  DataVector u_sample_transformed(numDims);

  for (size_t isample = 0; isample < numSamples; isample++) {
    u_vars.getRow(isample, u_sample);
    x_vars.getRow(isample, x_sample);
    u_vars_transformed.getRow(isample, u_sample_transformed);
    for (size_t idim = 0; idim < numDims; idim++) {
      // assert that x_vars and x_vars_transformed contain the same samples
      double inversionError =
          std::abs(u_sample[idim] - u_sample_transformed[idim]) / u_sample[idim];
      // std::cout << "u_sample[idim]:" << u_sample[idim] << std::endl;
      // std::cout << "x_sample[idim]:" << x_sample[idim] << std::endl;
      // std::cout << "u_sample_transformed[idim]:" << u_sample_transformed[idim] << std::endl;
      BOOST_CHECK_SMALL(inversionError, tolerance);
      // assert that no negative values appear
      BOOST_CHECK(x_sample[idim] >= 0);
    }
  }
}

void testEqualityRosenblattInverseRosenblattKDE(
    KernelDensityEstimator& kde, size_t numSamples = 1000, double tolerance = 1e-12,
    std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t numDims = kde.getDim();
  DataMatrix x_vars(numSamples, numDims);
  DataMatrix u_vars(numSamples, numDims);
  DataMatrix u_vars_transformed(numSamples, numDims);

  // init samples to be transformed
  randu(u_vars, seedValue);

  auto opInvRos = sgpp::op_factory::createOperationInverseRosenblattTransformationKDE(kde);
  auto opRos = sgpp::op_factory::createOperationRosenblattTransformationKDE(kde);

  // transform the u-space to x-space
  opInvRos->doTransformation(u_vars, x_vars);
  // transform them back to the u-space
  opRos->doTransformation(x_vars, u_vars_transformed);

  DataVector u_sample(numDims);
  DataVector x_sample(numDims);
  DataVector u_sample_transformed(numDims);

  for (size_t isample = 0; isample < numSamples; isample++) {
    u_vars.getRow(isample, u_sample);
    x_vars.getRow(isample, x_sample);
    u_vars_transformed.getRow(isample, u_sample_transformed);
    for (size_t idim = 0; idim < numDims; idim++) {
      // assert that x_vars and x_vars_transformed contain the same samples
      double inversionError =
          std::abs(u_sample[idim] - u_sample_transformed[idim]) / u_sample[idim];
      BOOST_CHECK_SMALL(inversionError, tolerance);
    }
  }
}

BOOST_AUTO_TEST_SUITE(testRosenblattTransformation)

BOOST_AUTO_TEST_CASE(testRosenblattLinear1D) {
  Grid* grid = Grid::createLinearGrid(1);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattLinearDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createLinearGrid(dim);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattPoly1D) {
  Grid* grid = Grid::createPolyGrid(1, 3);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattPolyDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createPolyGrid(dim, 3);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattPolyBoundary1D) {
  Grid* grid = Grid::createPolyBoundaryGrid(1, 3);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattPolyBoundaryDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createPolyBoundaryGrid(dim, 3);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattModPoly1D) {
  Grid* grid = Grid::createModPolyGrid(1, 3);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattModPolyDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createModPolyGrid(dim, 3);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattPolyClenshawCurtis1D) {
  Grid* grid = Grid::createPolyClenshawCurtisGrid(1, 3);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattPolyClenshawCurtisDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createPolyClenshawCurtisGrid(dim, 3);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattPolyClenshawCurtisBoundary1D) {
  Grid* grid = Grid::createPolyClenshawCurtisBoundaryGrid(1, 3);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattPolyClenshawCurtisBoundaryDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createPolyClenshawCurtisBoundaryGrid(dim, 3);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattModPolyClenshawCurtis1D) {
  Grid* grid = Grid::createModPolyClenshawCurtisGrid(1, 3);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattModPolyClenshawCurtisDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createModPolyClenshawCurtisGrid(dim, 3);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattBspline1D) {
  Grid* grid = Grid::createBsplineGrid(1, 3);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattBsplineDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createBsplineGrid(dim, 3);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattBsplineBoundary1D) {
  Grid* grid = Grid::createBsplineBoundaryGrid(1, 3);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattBsplineBoundaryDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createBsplineBoundaryGrid(dim, 3);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattModBspline1D) {
  Grid* grid = Grid::createModBsplineGrid(1, 3);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattModBsplineDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createModBsplineGrid(dim, 3);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattBsplineClenshawCurtis1D) {
  Grid* grid = Grid::createBsplineClenshawCurtisGrid(1, 3);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattBsplineClenshawCurtisDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createBsplineClenshawCurtisGrid(dim, 3);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattModBsplineClenshawCurtis1D) {
  Grid* grid = Grid::createModBsplineClenshawCurtisGrid(1, 3);
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    hierarchize(grid, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
  }
  delete grid;
}

BOOST_AUTO_TEST_CASE(testRosenblattModBsplineClenshawCurtisDD) {
  DataVector alpha(20);
  std::uint32_t numSamples = 100;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    Grid* grid = Grid::createModBsplineClenshawCurtisGrid(dim, 3);
    for (std::uint32_t ilevel = 1; ilevel < 4; ilevel++) {
      hierarchize(grid, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
    }
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattKDE1D) {
  size_t numDims = 1;
  std::uint32_t numSamples = 100;
  while (numSamples <= 500) {
    // load samples
    DataMatrix samples(numSamples, numDims);
    randn(samples);

    // estimate density
    KernelDensityEstimator kde(samples);

    // do the test
    testEqualityRosenblattInverseRosenblattKDE(kde, numSamples);

    numSamples += 200;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattKDEDD) {
  std::uint32_t numSamples = 100;
  size_t numDims = 5;
  while (numSamples <= 500) {
    // load samples
    DataMatrix samples(numSamples, numDims);
    randn(samples);

    // estimate density
    KernelDensityEstimator kde(samples);

    // do the test
    testEqualityRosenblattInverseRosenblattKDE(kde, numSamples);

    numSamples += 200;
  }
}

BOOST_AUTO_TEST_SUITE_END()
