// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/application/GaussianKDE.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>
#include <random>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;

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

Grid* hierarchize(std::uint32_t dim, std::uint32_t level, DataVector& alpha,
                  double (*func)(DataVector&)) {
  std::unique_ptr<Grid> grid = Grid::createLinearGrid(dim);
  grid->getGenerator().regular(level);
  alpha.resize(grid->getSize());

  DataVector coords(dim);
  GridStorage& gs = grid->getStorage();
  for (size_t i = 0; i < grid->getSize(); i++) {
    gs.getGridPoint(i).getCoords(coords);
    alpha[i] = func(coords);
  }
  sgpp::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha);
  return grid.release();
}

void randu(DataVector& rvar, std::mt19937& generator) {
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (size_t j = 0; j < rvar.getSize(); ++j) {
    rvar[j] = distribution(generator);
  }
}

void randu(DataMatrix& rvar, std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t nsamples = rvar.getNrows(), ndim = rvar.getNcols();

  std::mt19937 generator(seedValue);
  DataVector sample(ndim);
  for (size_t i = 0; i < nsamples; ++i) {
    randu(sample, generator);
    rvar.setRow(i, sample);
  }
}

void testEqualityRosenblattInverseRosenblatt1D(
    Grid& grid, DataVector& alpha, size_t numSamples = 1000, double tolerance = 1e-14,
    std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t numDims = grid.getStorage().getDimension();
  DataVector u_vars(numSamples);
  std::mt19937 generator(seedValue);

  // init samples to be transformed
  randu(u_vars, generator);

  auto opInvRos = sgpp::op_factory::createOperationInverseRosenblattTransformation1D(grid);
  auto opRos = sgpp::op_factory::createOperationRosenblattTransformation1D(grid);

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
    Grid& grid, DataVector& alpha, size_t numSamples = 1000, double tolerance = 1e-14,
    std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t numDims = grid.getStorage().getDimension();
  DataMatrix x_vars(numSamples, numDims);
  DataMatrix u_vars(numSamples, numDims);
  DataMatrix u_vars_transformed(numSamples, numDims);

  // init samples to be transformed
  randu(u_vars, seedValue);

  auto opInvRos = sgpp::op_factory::createOperationInverseRosenblattTransformation(grid);
  auto opRos = sgpp::op_factory::createOperationRosenblattTransformation(grid);

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
      BOOST_CHECK_SMALL(inversionError, tolerance);
    }
  }
}

BOOST_AUTO_TEST_SUITE(testRosenblattTransformation)

BOOST_AUTO_TEST_CASE(testRosenblattLinear1D) {
  Grid* grid = NULL;
  DataVector alpha(20);
  std::uint32_t numSamples = 1000;
  for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
    grid = hierarchize(1, ilevel, alpha, &parabola);
    testEqualityRosenblattInverseRosenblatt1D(*grid, alpha, numSamples);
    delete grid;
  }
}

BOOST_AUTO_TEST_CASE(testRosenblattLinearDD) {
  Grid* grid = NULL;
  DataVector alpha(20);
  std::uint32_t numSamples = 1000;
  for (std::uint32_t dim = 2; dim < 4; dim++) {
    for (std::uint32_t ilevel = 1; ilevel < 5; ilevel++) {
      grid = hierarchize(dim, ilevel, alpha, &parabola);
      testEqualityRosenblattInverseRosenblattDD(*grid, alpha, numSamples);
      delete grid;
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
