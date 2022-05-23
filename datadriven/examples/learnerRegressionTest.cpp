// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_learnerRegressionTest_cpp Regression Learner
 * This example demonstrates sparse grid regression learning.
 */

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegressionLearner.hpp>
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
 * @brief getLearner
 * @param dimension is the number of dimensions
 * @param regularizationConfig
 * @return a sparse grid regression learner
 */
sgpp::datadriven::RegressionLearner getLearner(
    size_t dimension, sgpp::datadriven::RegularizationConfiguration regularizationConfig) {
  auto gridConfig = sgpp::base::RegularGridConfiguration();
  gridConfig.dim_ = dimension;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::ModLinear;

  // gridConfig.type_ = sgpp::base::GridType::ModNakBspline;
  gridConfig.maxDegree_ = 3;

  auto adaptivityConfig = sgpp::base::AdaptivityConfiguration();
  adaptivityConfig.numRefinementPoints_ = 0;
  adaptivityConfig.numRefinements_ = 0;

  auto solverConfig = sgpp::solver::SLESolverConfiguration();
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-8;
  solverConfig.threshold_ = 1e-5;

  return sgpp::datadriven::RegressionLearner(gridConfig, adaptivityConfig, solverConfig,
                                             solverConfig, regularizationConfig);
}

/**
 * @brief showRegularizationConfiguration
 * @param regularizationConfig
 * @return type of the regularization method as string
 */
std::string showRegularizationConfiguration(
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig) {
  std::ostringstream ss;
  const auto regType = regularizationConfig.type_;
  if (regType == sgpp::datadriven::RegularizationType::Diagonal) {
    ss << "type: DiagonalMatrix\t";
  } else if (regType == sgpp::datadriven::RegularizationType::Identity) {
    ss << "type: IdentityMatrix\t";
  } else if (regType == sgpp::datadriven::RegularizationType::Laplace) {
    ss << "type: Laplace\t";
  } else {
    ss << "type: unknown\t";
  }

  ss << "lambda: " << regularizationConfig.lambda_
     << "\tmultiplicationFactor: " << regularizationConfig.exponentBase_;
  return ss.str();
}

/**
 * @brief gridSearch performs a hyper-parameter grid search over configs using a
 * holdout validation set.
 * @param configs are the regularization configs that will be tried
 * @param dimension is the number of dimensions
 * @param xTrain are the training predictors
 * @param yTrain is the training target
 * @param xValidation are the validation predictors
 * @param yValidation is the validation target
 * @return best found regularization configuration
 */
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
      std::cout << "Better! RMSE is now " << std::sqrt(curMSE) << std::endl;
      bestConfig = config;
      bestMSE = curMSE;
    } else {
      std::cout << "Worse!  RMSE is now " << std::sqrt(curMSE) << std::endl;
    }
  }
  std::cout << "gridSearch finished with parameters " << showRegularizationConfiguration(bestConfig)
            << std::endl;
  return bestConfig;
}

/**
 * @brief getConfigs
 * @return some regularization configurations for seven lambdas between 1 and 0.000001
 * and for exponent bases 1.0, 0.5, 0.25, 0.125
 */
std::vector<sgpp::datadriven::RegularizationConfiguration> getConfigs() {
  decltype(getConfigs()) result;
  std::vector<double> lambdas = {1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};

  std::vector<double> exponentBases = {1.0, 0.5, 0.25, 0.125};
  for (const auto lambda : lambdas) {
    // Identity
    const auto regularizationType = sgpp::datadriven::RegularizationType::Identity;
    auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
    regularizationConfig.type_ = regularizationType;
    regularizationConfig.lambda_ = lambda;
    regularizationConfig.exponentBase_ = 0.25;
    result.push_back(regularizationConfig);
    {
      // Laplace
      const auto regularizationType = sgpp::datadriven::RegularizationType::Laplace;
      auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
      regularizationConfig.type_ = regularizationType;
      regularizationConfig.lambda_ = lambda;
      regularizationConfig.exponentBase_ = 0.25;
      result.push_back(regularizationConfig);
    }
    // Diagonal
    for (const auto exponentBase : exponentBases) {
      const auto regularizationType = sgpp::datadriven::RegularizationType::Diagonal;
      auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
      regularizationConfig.type_ = regularizationType;
      regularizationConfig.lambda_ = lambda;
      regularizationConfig.exponentBase_ = exponentBase;
      result.push_back(regularizationConfig);
    }
  }

  return result;
}

/**
 * @brief main is an example for the RegressionLearner. It performs a grid search for the best
 * hyper-parameter for the Friedman3 dataset using the diagonal Tikhonov regularization method.
 */
int main(int argc, char** argv) {
  const auto filenameTrain = std::string("../datasets/friedman/friedman3_10k_train.arff");
  const auto filenameValidation = std::string("../datasets/friedman/friedman3_10k_validation.arff");
  const auto filenameTest = std::string("../datasets/friedman/friedman3_10k_test.arff");

  auto dataTrain = sgpp::datadriven::ARFFTools::readARFFFromFile(filenameTrain);
  std::cout << "Read file " << filenameTrain << "." << std::endl;
  auto xTrain = dataTrain.getData();
  auto yTrain = dataTrain.getTargets();
  const auto dimensions = dataTrain.getDimension();

  auto dataValidation = sgpp::datadriven::ARFFTools::readARFFFromFile(filenameValidation);
  std::cout << "Read file " << filenameValidation << "." << std::endl;
  auto xValidation = dataValidation.getData();
  auto yValidation = dataValidation.getTargets();

  const auto configs = getConfigs();
  const auto bestConfig = gridSearch(configs, dimensions, xTrain, yTrain, xValidation, yValidation);

  auto dataTest = sgpp::datadriven::ARFFTools::readARFFFromFile(filenameTest);
  std::cout << "Read file " << filenameTest << "." << std::endl;
  auto xTest = dataTest.getData();
  auto yTest = dataTest.getTargets();

  auto learner = getLearner(dimensions, bestConfig);
  learner.train(xTrain, yTrain);
  const auto MSETest = learner.getMSE(xTest, yTest);
  std::cout << "Best config got a testing MSE of " << MSETest << "!" << std::endl;
}
