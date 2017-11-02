// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <random>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::GridType;

void randu(DataVector& rvar, std::mt19937& generator) {
  std::normal_distribution<double> distribution(0.5, 0.1);
  for (size_t j = 0; j < rvar.getSize(); ++j) {
    rvar[j] = distribution(generator);
    while (rvar[j] < 0 || rvar[j] > 1) // make sure the sample is in the unit cube
      rvar[j] = distribution(generator);
  }
}

void createSamples(DataMatrix& rvar, std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t nsamples = rvar.getNrows(), ndim = rvar.getNcols();

  std::mt19937 generator(seedValue);
  DataVector sample(ndim);
  for (size_t i = 0; i < nsamples; ++i) {
    randu(sample, generator);
    rvar.setRow(i, sample);
  }
}

void solve(DataMatrix samples, sgpp::base::RegularGridConfiguration gridConfig, double lambda) {
}


int main(int argc, char** argv) {
  size_t d = 2;
  int level = 3;
  GridType gridType = sgpp::base::GridType::Linear;
  double lambda = 0.0;

  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdpativityConfiguration adaptivityConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  sgpp::datadriven::CrossvalidationForRegularizationConfiguration crossvalidationConfig;
  sgpp::datadriven::SGDEConfiguration sgdeConfig;
  gridConfig.dim_ = d;
  gridConfig.level_ = level;
  gridConfig.type_ = gridType;
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.noPoints_ = 15;
  solverConfig.threshold_ = 1e-15;
  solverConfig.verbose_ = false;
  regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Laplace;
  crossvalidationConfig.enable_ = false;
  sgdeConfig.makePositive_ = false;

  sgpp::datadriven::LearnerSGDE learnerSGDE(gridConfig, adaptivityConfig, solverConfig, regularizationConfig,
                                          crossvalidationConfig, sgdeConfig);
  DataMatrix trainSamples;
  createSamples(trainSamples);
  learnerSGDE.initialize(trainSamples);
  sgpp::base::Grid grid = learnerSGDE.getGrid();
  DataVector alpha(learnerSGDE.getAlpha());

  return 0;
}
