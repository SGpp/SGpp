// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunctionGradient.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/datadriven/application/RegressionLearner.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>

using namespace sgpp::base;

// Friedmann function. We should replace this with the Friedman dataset
double objectiveFunc(DataVector x) {
  return 10 * sin(M_PI * x[0] * x[1]) + 20 * (x[2] - 0.5) * (x[2] - 0.5) + 10 * x[3] + 5 * x[4];
}

int main() {
  int numPoints = 1000;
  int dim = 5;
  DataMatrix evaluationPoints(numPoints, dim);
  DataVector functionValues(numPoints);

  // create random dataset
  RandomNumberGenerator::getInstance().setSeed(42);
  DataVector v(dim);
  for (int i = 0; i < numPoints; i++) {
    RandomNumberGenerator::getInstance().getUniformRV(v, 0.0, 1.0);
    evaluationPoints.setRow(i, v);
    functionValues[i] = objectiveFunc(v);
  }

  //  const auto regularizationType = sgpp::datadriven::RegularizationType::Diagonal;
  const auto regularizationType = sgpp::datadriven::RegularizationType::Identity;
  double lambda = 0.01;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig{};
  regularizationConfig.type_ = regularizationType;
  regularizationConfig.lambda_ = lambda;

  int level = 3;
  auto gridConfig = RegularGridConfiguration();
  gridConfig.dim_ = dim;
  gridConfig.level_ = level;
  // gridConfig.type_ = GridType::Linear;
  // gridConfig.type_ = GridType::NakBsplineExtended;
  gridConfig.type_ = GridType::NakPBspline;
  gridConfig.maxDegree_ = 3;

  auto adaptivityConfig = AdaptivityConfiguration();  // no adaptivity
  adaptivityConfig.numRefinements_ = 0;

  auto solverConfig = sgpp::solver::SLESolverConfiguration();
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-8;
  solverConfig.threshold_ = 1e-5;
  auto learner = sgpp::datadriven::RegressionLearner(gridConfig, adaptivityConfig, solverConfig,
                                                     solverConfig, regularizationConfig);
  learner.train(evaluationPoints, functionValues);
  std::shared_ptr<Grid> grid = learner.getGridPtr();
  DataVector coefficients = learner.getWeights();

  auto interpolant = std::make_shared<InterpolantScalarFunction>(*grid, coefficients);
  //   auto interpolantGradient =
  //       std::make_shared<InterpolantScalarFunctionGradient>(*grid, coefficients);

  // L2 error
  int numMCPoints = 1000;
  double err = 0;
  double trueVal;
  double approxVal;
  for (int i = 0; i < numMCPoints; i++) {
    RandomNumberGenerator::getInstance().getUniformRV(v, 0.0, 1.0);
    trueVal = objectiveFunc(v);
    approxVal = interpolant->eval(v);
    err += (trueVal - approxVal) * (trueVal - approxVal);
  }
  err = sqrt(err / numMCPoints);
  std::cout << "L2 error: " << err << "\n";
}