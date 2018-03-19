/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * RosenblattTransformation.cpp
 *
 *  Created on: 22.01.2018
 *      Author: Lars Wolfsteller
 */

#include <sgpp/datadriven/datamining/modules/dataSource/RosenblattTransformation.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation.hpp>

#include <random>

namespace sgpp {
namespace datadriven {

RosenblattTransformation::RosenblattTransformation(Dataset* dataset, size_t numSamples)
  : grid(nullptr),
    alpha(nullptr),
    datasetTransformed(nullptr),
    datasetInvTransformed(nullptr),
    numSamples(numSamples) {
  // Sample #numSamples random samples from dataset
  DataMatrix samples(numSamples, dataset->getDimension());
  DataVector currSample(dataset->getDimension());

  std::mt19937 generator;
  std::uniform_real_distribution<double> distr(
      0, static_cast<double>(dataset->getNumberInstances()-1));

  for (unsigned int i = 0; i < numSamples; i++) {
      dataset->getData().getRow(static_cast<size_t>(distr(generator)), currSample);
      samples.setRow(i, currSample);
    }

  // Approximate probability density function (PDF)
  double lambda = 1e-5;
  size_t dim = dataset->getDimension();
  size_t level = 3;
  sgpp::datadriven::LearnerSGDE learner = createSGDELearner(dim, level, lambda);
  learner.initialize(samples);
  learner.train();

  // Get grid and alpha
  grid = learner.getGrid();
  alpha = learner.getSurpluses();
}


Dataset* RosenblattTransformation::doTransformation(Dataset* dataset) {
  OperationRosenblattTransformation* opRos(
      sgpp::op_factory::createOperationRosenblattTransformation(*this->grid));
  datasetTransformed = new Dataset{dataset->getNumberInstances(), dataset->getDimension()};

  opRos->doTransformation(this->alpha.get(), &dataset->getData(), &datasetTransformed->getData());

  return datasetTransformed;
}


Dataset* RosenblattTransformation::doInverseTransformation(Dataset* dataset) {
  std::unique_ptr<sgpp::datadriven::OperationInverseRosenblattTransformation> opInvRos(
      sgpp::op_factory::createOperationInverseRosenblattTransformation(*this->grid));
  datasetInvTransformed = new Dataset{dataset->getNumberInstances(), dataset->getDimension()};

  opInvRos->doTransformation(
      this->alpha.get(), &dataset->getData(), &datasetInvTransformed->getData());

  return datasetInvTransformed;
}


sgpp::datadriven::LearnerSGDE RosenblattTransformation::createSGDELearner(size_t dim, size_t level,
                                                double lambda) {
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = dim;
  gridConfig.level_ = static_cast<int>(level);
  gridConfig.type_ = sgpp::base::GridType::Linear;

  // configure adaptive refinement
  sgpp::base::AdpativityConfiguration adaptConfig;
  adaptConfig.numRefinements_ = 0;

  // configure solver
  sgpp::solver::SLESolverConfiguration solverConfig;
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-10;
  solverConfig.threshold_ = 1e-10;

  // configure regularization
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.regType_ =
    sgpp::datadriven::RegularizationType::Laplace;

  // configure learner
  sgpp::datadriven::CrossvalidationForRegularizationConfiguration
    crossvalidationConfig;
  crossvalidationConfig.enable_ = false;

  sgpp::datadriven::LearnerSGDE learner(gridConfig,
                                        adaptConfig,
                                        solverConfig,
                                        regularizationConfig,
                                        crossvalidationConfig);
  return learner;
}
} /* namespace datadriven */
} /* namespace sgpp */
