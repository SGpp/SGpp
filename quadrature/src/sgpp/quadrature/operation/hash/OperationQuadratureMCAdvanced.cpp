// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/quadrature/operation/hash/OperationQuadratureMCAdvanced.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/Random.hpp>
#include <sgpp/quadrature/sampling/HaltonSampleGenerator.hpp>
#include <sgpp/quadrature/sampling/LatinHypercubeSampleGenerator.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>
#include <sgpp/quadrature/sampling/StratifiedSampleGenerator.hpp>

#include <cmath>
#include <iostream>
#include <vector>

namespace sgpp {
namespace quadrature {

OperationQuadratureMCAdvanced::OperationQuadratureMCAdvanced(sgpp::base::Grid& grid,
                                                             size_t numberOfSamples,
                                                             std::uint64_t seed)
    : grid(&grid), numberOfSamples(numberOfSamples), seed(seed) {
  dimensions = grid.getDimension();
  myGenerator = new sgpp::quadrature::NaiveSampleGenerator(dimensions, seed);
}

OperationQuadratureMCAdvanced::OperationQuadratureMCAdvanced(size_t dimensions,
                                                             size_t numberOfSamples,
                                                             std::uint64_t seed)
    : grid(NULL), numberOfSamples(numberOfSamples), dimensions(dimensions), seed(seed) {
  myGenerator = new sgpp::quadrature::NaiveSampleGenerator(dimensions, seed);
}

OperationQuadratureMCAdvanced::~OperationQuadratureMCAdvanced() {
  if (myGenerator != NULL) {
    delete myGenerator;
  }
}

void OperationQuadratureMCAdvanced::useNaiveMonteCarlo() {
  if (myGenerator != NULL) {
    delete myGenerator;
  }

  myGenerator = new sgpp::quadrature::NaiveSampleGenerator(dimensions, seed);
}

void OperationQuadratureMCAdvanced::useStratifiedMonteCarlo(
    std::vector<size_t>& strataPerDimension) {
  if (myGenerator != NULL) {
    delete myGenerator;
  }

  myGenerator = new sgpp::quadrature::StratifiedSampleGenerator(strataPerDimension, seed);
}

void OperationQuadratureMCAdvanced::useLatinHypercubeMonteCarlo() {
  if (myGenerator != NULL) {
    delete myGenerator;
  }

  myGenerator =
      new sgpp::quadrature::LatinHypercubeSampleGenerator(dimensions, numberOfSamples, seed);
}

void OperationQuadratureMCAdvanced::useQuasiMonteCarloWithHaltonSequences() {
  if (myGenerator != NULL) {
    delete myGenerator;
  }

  myGenerator = new sgpp::quadrature::HaltonSampleGenerator(dimensions);
}

double OperationQuadratureMCAdvanced::doQuadrature(sgpp::base::DataVector& alpha) {
  sgpp::base::DataMatrix dm(numberOfSamples, dimensions);

  myGenerator->getSamples(dm);

  sgpp::base::DataVector res = sgpp::base::DataVector(numberOfSamples);
  sgpp::op_factory::createOperationMultipleEval(*grid, dm)->mult(alpha, res);
  return res.sum() / static_cast<double>(numberOfSamples);
}

double OperationQuadratureMCAdvanced::doQuadratureFunc(FUNC func, void* clientdata) {
  // double* p = new double[dimensions];

  sgpp::base::DataMatrix dm(numberOfSamples, dimensions);
  myGenerator->getSamples(dm);

  // create number of paths (uniformly drawn from [0,1]^d)
  double res = 0;

  for (size_t i = 0; i < numberOfSamples; i++) {
    sgpp::base::DataVector dv(dimensions);
    dm.getRow(i, dv);
    res += func(static_cast<int>(dimensions), dv.getPointer(), clientdata);
  }

  // delete p;
  return res / static_cast<double>(numberOfSamples);
}

double OperationQuadratureMCAdvanced::doQuadratureL2Error(FUNC func, void* clientdata,
                                                          sgpp::base::DataVector& alpha) {
  // create number of paths (uniformly drawn from [0,1]^d)
  sgpp::base::DataMatrix dm(numberOfSamples, dimensions);
  myGenerator->getSamples(dm);

  double x;
  double* p = new double[dimensions];

  sgpp::base::DataVector point(dimensions);
  std::unique_ptr<sgpp::base::OperationEval> opEval(sgpp::op_factory::createOperationEval(*grid));
  double res = 0;

  for (size_t i = 0; i < numberOfSamples; i++) {
    for (size_t d = 0; d < dimensions; d++) {
      x = dm.get(i, d);
      p[d] = x;
      point[d] = x;
    }
    res += pow(
        func(static_cast<int>(dimensions), p, clientdata) - opEval->eval(alpha, point), 2);
  }

  delete[] p;
  return sqrt(res / static_cast<double>(numberOfSamples));
}

size_t OperationQuadratureMCAdvanced::getDimensions() { return dimensions; }

}  // namespace quadrature
}  // namespace sgpp
