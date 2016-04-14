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

int main(int argc, char** argv) {
  const auto filename = std::string("../tests/data/friedman_4d_2000.arff");
  auto dataset = sgpp::datadriven::ARFFTools::readARFF(filename);
  std::cout << "Read file " << filename << " ." << std::endl;

  auto data = dataset.getData();
  auto classes = dataset.getTargets();

  auto gridConfig = sgpp::base::RegularGridConfiguration();
  gridConfig.dim_ = dataset.getDimension();
  gridConfig.level_ = 5;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  // TODO(ARGH): Validate!
  auto adaptivityConfig = sgpp::base::AdpativityConfiguration();
  adaptivityConfig.noPoints_ = 0;
  adaptivityConfig.numRefinements_ = 0;

  auto regularizationType = sgpp::datadriven::RegularizationType::Identity;
  auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
  regularizationConfig.regType_ = regularizationType;

  auto solverConfig = sgpp::solver::SLESolverConfiguration();
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-5;

  std::cout << "Initializing the learner." << std::endl;
  auto learner = sgpp::datadriven::RegressionLearner(gridConfig, adaptivityConfig, solverConfig,
                                                     regularizationConfig);

  std::cout << "Training the learner." << std::endl;
  try {
    learner.train(data, classes);
  } catch (const std::exception& e) {
    std::cout << "std exception " << e.what() << std::endl;
  } catch (...) {
    std::cout << "unknown exception" << std::endl;
  }

  std::cout << "Finished training." << std::endl;
}
