// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/dataSource/RosenblattTransformation.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>

#include <random>

namespace sgpp {
namespace datadriven {

RosenblattTransformation::RosenblattTransformation()
    : grid(nullptr), alpha(nullptr), datasetTransformed(nullptr), datasetInvTransformed(nullptr) {}

void RosenblattTransformation::initialize(Dataset *dataset, DataTransformationConfig config) {
  RosenblattTransformationConfig rbConfig = config.rosenblattConfig_;

  // Sample #numSamples random samples from dataset
  DataMatrix samples(rbConfig.numSamples_, dataset->getDimension());
  DataVector currSample(dataset->getDimension());

  std::mt19937 generator;
  std::uniform_real_distribution<double> distr(
      0, static_cast<double>(dataset->getNumberInstances() - 1));

  for (unsigned int i = 0; i < rbConfig.numSamples_; i++) {
    dataset->getData().getRow(static_cast<size_t>(distr(generator)), currSample);
    samples.setRow(i, currSample);
  }

  // Approximate probability density function (PDF)
  size_t dim = dataset->getDimension();
  LearnerSGDE learner = createSGDELearner(dim, rbConfig);

  learner.initialize(samples);
  learner.train();

  // Get grid and alpha
  grid = learner.getSharedGrid();
  alpha = learner.getSharedSurpluses();

  std::cout << "Rosenblatt transformation initialized" << std::endl;
}

Dataset *RosenblattTransformation::doTransformation(Dataset *dataset) {
  std::cout << "Performing Rosenblatt transformation" << std::endl;
  OperationRosenblattTransformation *opRos(
      sgpp::op_factory::createOperationRosenblattTransformation(*this->grid));
  datasetTransformed = new Dataset{dataset->getNumberInstances(), dataset->getDimension()};

  opRos->doTransformation(this->alpha.get(), &dataset->getData(), &datasetTransformed->getData());

  return datasetTransformed;
}

Dataset *RosenblattTransformation::doInverseTransformation(Dataset *dataset) {
  std::cout << "Performing Rosenblatt inverse transformation" << std::endl;
  std::unique_ptr<sgpp::datadriven::OperationInverseRosenblattTransformation> opInvRos(
      sgpp::op_factory::createOperationInverseRosenblattTransformation(*this->grid));
  datasetInvTransformed = new Dataset{dataset->getNumberInstances(), dataset->getDimension()};

  opInvRos->doTransformation(this->alpha.get(), &dataset->getData(),
                             &datasetInvTransformed->getData());

  return datasetInvTransformed;
}

sgpp::datadriven::LearnerSGDE RosenblattTransformation::createSGDELearner(
    size_t dim, RosenblattTransformationConfig config) {
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = dim;
  gridConfig.level_ = static_cast<int>(config.gridLevel_);
  gridConfig.type_ = sgpp::base::GridType::Linear;

  // configure adaptive refinement
  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  adaptivityConfig.numRefinements_ = 0;

  // configure solver
  sgpp::solver::SLESolverConfiguration solverConfig;
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = config.solverMaxIterations_;
  solverConfig.eps_ = config.solverEps_;
  solverConfig.threshold_ = config.solverThreshold_;

  // configure regularization
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Laplace;

  // configure learner
  sgpp::datadriven::CrossvalidationConfiguration crossValidationConfig;
  crossValidationConfig.enable_ = false;

  sgpp::datadriven::LearnerSGDE learner(gridConfig, adaptivityConfig, solverConfig,
                                        regularizationConfig, crossValidationConfig);
  return learner;
}
} /* namespace datadriven */
} /* namespace sgpp */
