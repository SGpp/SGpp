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

#include <ctime>

namespace sgpp {
namespace datadriven {

int RandomNumber () { return (std::rand()%100); }

RosenblattTransformation::RosenblattTransformation(Dataset* dataset, size_t numSamples)
  : DataTransformation{},
    grid(nullptr),
    alpha(nullptr),
    datasetTransformed(nullptr),
    datasetInvTransformed(nullptr),
    numSamples(numSamples) {
  // Sample #numSamples random samples from dataset
  DataMatrix samples(numSamples, dataset->getDimension());
  DataVector currSample(dataset->getDimension());
  DataVector randVector(numSamples);

  // Fill vector with random numbers
  std::srand(std::time(0));
  std::generate(randVector.begin(), randVector.end(), RandomNumber);

  for (unsigned int i = 0; i < numSamples; i++) {
    dataset->getData().getRow(static_cast<size_t>(randVector[i]), currSample);
    samples.setRow(i, currSample);
  }

  // Approximate probability density function (PDF)
  double lambda = 1e-5;
  size_t dim = dataset->getDimension();
  size_t level = 4;
  sgpp::datadriven::LearnerSGDE learner = createSGDELearner(dim, level, lambda);
  learner.initialize(samples);

  // Get grid and alpha
  grid = learner.getGrid().get();
  alpha = learner.getSurpluses().get();
}


Dataset* RosenblattTransformation::doTransformation(Dataset* dataset) {
  std::unique_ptr<sgpp::datadriven::OperationRosenblattTransformation> opRos(
      sgpp::op_factory::createOperationRosenblattTransformation(*grid));

  datasetTransformed = new Dataset{dataset->getNumberInstances(), dataset->getDimension()};

  opRos->doTransformation(alpha, &dataset->getData(), &datasetTransformed->getData());
  return datasetTransformed;
}


Dataset* RosenblattTransformation::doInverseTransformation(Dataset* dataset) {
  std::unique_ptr<sgpp::datadriven::OperationInverseRosenblattTransformation> opInvRos(
      sgpp::op_factory::createOperationInverseRosenblattTransformation(*grid));

  datasetInvTransformed = new Dataset{dataset->getNumberInstances(), dataset->getDimension()};

  opInvRos->doTransformation(alpha, &dataset->getData(), &datasetInvTransformed->getData());
  return datasetInvTransformed;
}


/**
 * Helper function
 * It configures and creates a SGDE learner with meaningful parameters
 */
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
