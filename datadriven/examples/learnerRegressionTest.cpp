// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/application/RegressionLearner.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <string>
#include <exception>

sgpp::datadriven::RegressionLearner getLearner(sgpp::datadriven::Dataset& data) {
  auto gridConfig = sgpp::base::RegularGridConfiguration();
  gridConfig.dim_ = data.getDimension();
  gridConfig.level_ = 2;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  auto adaptivityConfig = sgpp::base::AdpativityConfiguration();
  adaptivityConfig.noPoints_ = 0;
  adaptivityConfig.numRefinements_ = 0;

  auto regularizationType = sgpp::datadriven::RegularizationType::Diagonal;
  auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
  regularizationConfig.regType_ = regularizationType;
  regularizationConfig.lambda = 0.1;

  auto solverConfig = sgpp::solver::SLESolverConfiguration();
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 5000;
  solverConfig.eps_ = 1e-6;

  std::cout << "Initializing the learner." << std::endl;
  return sgpp::datadriven::RegressionLearner(gridConfig, adaptivityConfig, solverConfig,
                                             regularizationConfig);
}

int main(int argc, char** argv) {
  const auto filenameTrain = std::string("../tests/data/friedman3_10k_train.arff");
  const auto filenameValidation = std::string("../tests/data/friedman3_10k_validation.arff");

  auto dataTrain = sgpp::datadriven::ARFFTools::readARFF(filenameTrain);
  std::cout << "Read file " << filenameTrain << "." << std::endl;
  auto xTrain = dataTrain.getData();
  auto yTrain = dataTrain.getTargets();

  auto learner = getLearner(dataTrain);

  std::cout << "Training the learner." << std::endl;
  learner.train(xTrain, yTrain);
  std::cout << "Finished training." << std::endl;

  const double MSETrain = learner.getMSE(xTrain, yTrain);
  std::cout << "Training MSE = " << MSETrain << std::endl;

  auto dataValidation = sgpp::datadriven::ARFFTools::readARFF(filenameValidation);
  std::cout << "Read file " << filenameValidation << "." << std::endl;
  auto xValidation = dataValidation.getData();
  auto yValidation = dataValidation.getTargets();

  const double MSEValidation = learner.getMSE(xValidation, yValidation);
  std::cout << "Validation MSE = " << MSEValidation << std::endl;
}
