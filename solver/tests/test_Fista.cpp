// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/solver/sle/fista/Fista.hpp>
#include <sgpp/solver/sle/fista/RidgeFunction.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <memory>
#include <random>
#include <vector>

using sgpp::base::DataVector;
using sgpp::base::DataMatrix;

BOOST_AUTO_TEST_SUITE(TestFista)

BOOST_AUTO_TEST_CASE(testFista) {
  const size_t dim = 2;
  const size_t level = 4;
  const auto numExamplesTrain = 1470;

  auto yTrain = DataVector(numExamplesTrain);
  auto datasetTrain = DataMatrix(numExamplesTrain, dim);
  for (auto i = 0; i < numExamplesTrain; ++i) {
    const double x1 = std::abs(std::sin(i));
    const double x2 = std::abs(std::sin(i * x1));
    const double comb = std::sinh(x1) + std::sinh(x2);
    const auto row = DataVector(std::vector<double>({x1, x2}));
    datasetTrain.setRow(i, row);
    yTrain.set(i, comb);
  }

  auto grid = sgpp::base::Grid::createModLinearGrid(dim);
  grid->getGenerator().regular(level);

  double lambda = 1e-3;
  const size_t maxIt = 2000;
  const double treshold = 1e-9;

  auto ridge = sgpp::solver::RidgeFunction(lambda);
  auto solver = sgpp::solver::Fista<decltype(ridge)>(ridge);
  auto op = sgpp::op_factory::createOperationMultipleEval(*grid, datasetTrain);

  auto weights = DataVector(grid->getSize(), 0.0);
  solver.solve(*op, weights, yTrain, maxIt, treshold);

  auto prediction = DataVector(numExamplesTrain);
  op->mult(weights, prediction);
  prediction.sub(yTrain);
  prediction.sqr();
  const double rmse = std::sqrt((1.0 / numExamplesTrain) * prediction.sum());
  BOOST_CHECK(rmse <= 0.001);

  // Check that solutions get smaller with increasing lambda!
  while (lambda < 10.0) {
    const double weightsNormBefore = weights.l2Norm();
    weights.setAll(0.0);
    lambda *= 10;
    auto ridge = sgpp::solver::RidgeFunction(lambda);
    auto solver = sgpp::solver::Fista<decltype(ridge)>(ridge);
    solver.solve(*op, weights, yTrain, maxIt, treshold);
    BOOST_CHECK(weights.l2Norm() <= weightsNormBefore);
  }
}

BOOST_AUTO_TEST_SUITE_END()
