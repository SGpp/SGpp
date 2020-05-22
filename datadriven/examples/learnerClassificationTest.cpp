// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_learnerClassificationTest_cpp Learner Classification Test
 *
 * This represents a small example how to use sparse grids for classification
 * problems. It uses the artificial Ripley dataset.
 */

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/ClassificationLearner.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <exception>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

/**
 * @brief getLearner creates a sparse grid classification learner.
 * @param dimension the number of dimensions.
 * @return a classification for dimension dimension
 */
sgpp::datadriven::ClassificationLearner getLearner(size_t dimension) {
  auto gridConfig = sgpp::base::RegularGridConfiguration();
  gridConfig.dim_ = dimension;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::ModLinear;

  auto adaptivityConfig = sgpp::base::AdaptivityConfiguration();
  adaptivityConfig.numRefinementPoints_ = 0;
  adaptivityConfig.numCoarseningPoints_ = 0;
  adaptivityConfig.numRefinements_ = 0;

  auto solverConfig = sgpp::solver::SLESolverConfiguration();
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 500;
  solverConfig.eps_ = 1e-8;

  auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Diagonal;
  regularizationConfig.lambda_ = 0.00001;
  regularizationConfig.exponentBase_ = 0.25;

  return sgpp::datadriven::ClassificationLearner(gridConfig, adaptivityConfig, solverConfig,
                                                 solverConfig, regularizationConfig);
}

/**
 * @brief main
   Creates a sparse grid classification learner and prints the training accuracy
 for the ripley
 dataset.
 * @return
 */
int main(int argc, char** argv) {
  const auto filenameTrain = std::string("../datasets/ripley/ripleyGarcke.train.arff");

  auto dataTrain = sgpp::datadriven::ARFFTools::readARFFFromFile(filenameTrain);
  std::cout << "Read file " << filenameTrain << "." << std::endl;
  auto xTrain = dataTrain.getData();
  auto yTrain = dataTrain.getTargets();
  const auto dimensions = dataTrain.getDimension();

  auto learner = getLearner(dimensions);
  learner.train(xTrain, yTrain);
  const auto accuracy = learner.getAccuracy(xTrain, yTrain);
  std::cout << "Best config got a training acc of " << accuracy << "!" << std::endl;
}
