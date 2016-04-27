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
#include <vector>
#include <exception>
#include <limits>
#include <ostream>

sgpp::datadriven::RegressionLearner getLearner(
    size_t dimension, sgpp::datadriven::RegularizationConfiguration regularizationConfig) {
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

  return sgpp::datadriven::RegressionLearner(gridConfig, adaptivityConfig, solverConfig,
                                             regularizationConfig);
}

std::string showRegularizationConfiguration(
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig) {
  std::ostringstream ss;
  switch (regularizationConfig.regType_) {
    case sgpp::datadriven::RegularizationType::Diagonal:
      ss << "type: DiagonalMatrix\t";
      break;
    case sgpp::datadriven::RegularizationType::Identity:
      ss << "type: IdentityMatrix\t";
      break;
    case sgpp::datadriven::RegularizationType::Laplace:
      ss << "type: Laplace\t";
      break;
    default:
      ss << "type: unknown\t";
  }

  ss << "lambda: " << regularizationConfig.lambda
     << "\tmultiplicationFactor: " << regularizationConfig.multiplicationFactor;
  return ss.str();
}

sgpp::datadriven::RegularizationConfiguration gridSearch(
    std::vector<sgpp::datadriven::RegularizationConfiguration> configs, size_t dimension,
    sgpp::base::DataMatrix& xTrain, sgpp::base::DataVector& yTrain,
    sgpp::base::DataMatrix& xValidation, sgpp::base::DataVector& yValidation) {
  double bestMSE = std::numeric_limits<double>::max();
  sgpp::datadriven::RegularizationConfiguration bestConfig;
  for (const auto& config : configs) {
    // Step 1: Create a learner
    auto learner = getLearner(dimension, config);
    // Step 2: Train it with the hyperparameter
    learner.train(xTrain, yTrain);
    // Step 3: Evaluate accuracy
    const double curMSE = learner.getMSE(xValidation, yValidation);
    std::cout << "Tested parameters are\n" << showRegularizationConfiguration(config) << ".\n";
    if (curMSE < bestMSE) {
      std::cout << "Better! MSE is now " << curMSE << std::endl;
      bestConfig = config;
      bestMSE = curMSE;
    } else {
      std::cout << "Worse!  MSE is now " << curMSE << std::endl;
    }
  }
  std::cout << "gridSearch finished with parameters " << showRegularizationConfiguration(bestConfig)
            << std::endl;
  return bestConfig;
}

std::vector<sgpp::datadriven::RegularizationConfiguration> getConfigs() {
  decltype(getConfigs()) result;
  std::vector<double> lambdas = {0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125};

  std::vector<double> multiplicationFactors = {0.5, 0.25, 0.125};
  for (const auto lambda : lambdas) {
    // Identity
    const auto regularizationType = sgpp::datadriven::RegularizationType::Identity;
    auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
    regularizationConfig.regType_ = regularizationType;
    regularizationConfig.lambda = lambda;
    regularizationConfig.multiplicationFactor = 0.25;
    result.push_back(regularizationConfig);

    // Diagonal
    for (const auto multiplicationFactor : multiplicationFactors) {
      const auto regularizationType = sgpp::datadriven::RegularizationType::Diagonal;
      auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
      regularizationConfig.regType_ = regularizationType;
      regularizationConfig.lambda = lambda;
      regularizationConfig.multiplicationFactor = multiplicationFactor;
      result.push_back(regularizationConfig);
    }
  }

  return result;
}

int main(int argc, char** argv) {
  const auto filenameTrain = std::string("../tests/data/friedman3_10k_train.arff");
  const auto filenameValidation = std::string("../tests/data/friedman3_10k_validation.arff");
  const auto filenameTest = std::string("../tests/data/friedman3_10k_testing.arff");

  auto dataTrain = sgpp::datadriven::ARFFTools::readARFF(filenameTrain);
  std::cout << "Read file " << filenameTrain << "." << std::endl;
  auto xTrain = dataTrain.getData();
  auto yTrain = dataTrain.getTargets();
  const auto dimensions = dataTrain.getDimension();

  auto dataValidation = sgpp::datadriven::ARFFTools::readARFF(filenameValidation);
  std::cout << "Read file " << filenameValidation << "." << std::endl;
  auto xValidation = dataValidation.getData();
  auto yValidation = dataValidation.getTargets();

  const auto configs = getConfigs();
  const auto bestConfig = gridSearch(configs, dimensions, xTrain, yTrain, xValidation, yValidation);

  auto dataTest = sgpp::datadriven::ARFFTools::readARFF(filenameTest);
  std::cout << "Read file " << filenameTest << "." << std::endl;
  auto xTest = dataTest.getData();
  auto yTest = dataTest.getTargets();

  auto learner = getLearner(dimensions, bestConfig);
  learner.train(xTrain, yTrain);
  const auto MSETest = learner.getMSE(xTest, yTest);
  std::cout << "Best config got a testing MSE of " << MSETest << "!" << std::endl;
}
