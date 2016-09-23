// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/application/ClassificationLearner.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <string>
#include <vector>
#include <exception>
#include <limits>
#include <ostream>

/**
 * @brief getLearner
 * @param dimension the number of new dimensions
 * @return a classification for dimension dimension
 */
sgpp::datadriven::ClassificationLearner getLearner(size_t dimension) {
  auto gridConfig = sgpp::base::RegularGridConfiguration();
  gridConfig.dim_ = dimension;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::ModLinear;

  auto adaptivityConfig = sgpp::base::AdpativityConfiguration();
  adaptivityConfig.noPoints_ = 0;
  adaptivityConfig.numRefinements_ = 0;

  auto solverConfig = sgpp::solver::SLESolverConfiguration();
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 500;
  solverConfig.eps_ = 1e-8;

  auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
  regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Diagonal;
  regularizationConfig.lambda_ = 0.00001;
  regularizationConfig.exponentBase_ = 0.25;

  return sgpp::datadriven::ClassificationLearner(gridConfig, adaptivityConfig, solverConfig,
                                                 solverConfig, regularizationConfig);
}

/**
 * @brief main
   Creates a sparse grid classification learner and returns the training accuracy for the ripley dataset.
 * @return
 */
int main(int argc, char** argv) {
  const auto filenameTrain =
      std::string("../../datasets/ripley/ripleyGarcke.train.arff");

  auto dataTrain = sgpp::datadriven::ARFFTools::readARFF(filenameTrain);
  std::cout << "Read file " << filenameTrain << "." << std::endl;
  auto xTrain = dataTrain.getData();
  auto yTrain = dataTrain.getTargets();
  const auto dimensions = dataTrain.getDimension();

  auto learner = getLearner(dimensions);
  learner.train(xTrain, yTrain);
  const auto accuracy = learner.getAccuracy(xTrain, yTrain);
  std::cout << "Best config got a training acc of " << accuracy << "!" << std::endl;
}
