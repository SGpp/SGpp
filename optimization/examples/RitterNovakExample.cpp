// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// Example file demonstrating the use of IterativeGridGeneratorRitterNovak
// and IterativeGridGeneratorFullAdaptiveRitterNovak from the SGpp optimization
// module.

// To fully enable this example, ensure that the following dependencies are installed (otherwise you
// will get a slightly limited functionality):
// 1. [libgp](https://github.com/mblum/libgp#building-and-testing)
// 2. [BayesOpt](https://rmcantin.github.io/bayesopt/html/install.html)
// Additionally, compile SGpp with the `USE_LIBGP` and `USE_BAYESOPT` options enabled.
// Note: Before compiling SGpp, verify that libgp and BayesOpt are installed system-wide.

#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorFullAdaptiveRitterNovak.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Rastrigin.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <array>  // For std::array
#include <cassert>
#include <cmath>
#include <iostream>   // For std::cout, std::endl
#include <memory>     // For std::unique_ptr
#include <random>     // For std::mt19937, std::uniform_real_distribution
#include <stdexcept>  // For std::runtime_error
#include <utility>    // For std::pair, std::make_pair
#include <vector>

/**
 * @brief Demonstrates the standard IterativeGridGeneratorRitterNovak.
 *
 * This function creates a B-spline grid, initializes the Ritter-Novak
 * grid generator, runs the generation, and then solves the resulting
 * system of linear equations to get the hierarchical coefficients.
 * Finally, it prints all grid points and their corresponding coefficients.
 */
void gridGenerationRitterNovak() {
  std::cout << "--- Running gridGenerationRitterNovak ---" << std::endl;
  const size_t dimension = 2;
  const size_t bSplineDegree = 3;
  const size_t maxGridPoints = 50;
  const double adaptivity = 0.15;

  std::unique_ptr<sgpp::base::ScalarFunction> function =
      std::make_unique<sgpp::optimization::test_problems::RastriginObjective>(dimension);

  std::unique_ptr<sgpp::base::Grid> grid(
      sgpp::base::Grid::createBsplineGrid(dimension, bSplineDegree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  sgpp::optimization::IterativeGridGeneratorRitterNovak gridGen(*function, *grid, maxGridPoints,
                                                                adaptivity);

  if (!gridGen.generate()) {
    throw std::runtime_error("Grid generation failed, exiting");
  }

  sgpp::base::HierarchisationSLE hierSLE(*grid);
  const size_t N = hierSLE.getDimension();

  sgpp::base::sle_solver::Auto sleSolver;
  sgpp::base::DataVector alpha(N);
  sgpp::base::DataVector b = gridGen.getFunctionValues();

  if (!sleSolver.solve(hierSLE, b, alpha)) {
    throw std::runtime_error("Solving failed");
  }

  // if success print all point on the grid with their coefficients
  std::cout << "Grid generation successful. Points: " << gridStorage.getSize() << std::endl;
  for (size_t i = 0; i < gridStorage.getSize(); ++i) {
    sgpp::base::DataVector point(dimension);
    gridStorage[i].getStandardCoordinates(point);
    std::cout << "Point " << i << ": " << point.toString() << ", Coefficient: " << alpha[i]
              << std::endl;
  }
  std::cout << "------------------------------------------" << std::endl;
}

/**
 * @brief Demonstrates IterativeGridGeneratorFullAdaptiveRitterNovak
 * with minimal setup.
 *
 * This function uses the fully adaptive Ritter-Novak generator.
 * It does not require solving an SLE externally, as the generator
 * computes the hierarchical coefficients (alpha) internally.
 */
void gridGenerationFullAdaptiveRitterNovakMinimal() {
  std::cout << "--- Running gridGenerationFullAdaptiveRitterNovakMinimal ---" << std::endl;
  const size_t dimension = 2;
  const size_t bSplineDegree = 3;
  const size_t maxGridPoints = 50;

  std::unique_ptr<sgpp::base::ScalarFunction> function =
      std::make_unique<sgpp::optimization::test_problems::RastriginObjective>(dimension);

  std::unique_ptr<sgpp::base::Grid> grid(
      sgpp::base::Grid::createBsplineGrid(dimension, bSplineDegree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  std::vector<std::pair<double, double>> domain(dimension);
  for (size_t t = 0; t < dimension; ++t) {
    domain[t] = std::make_pair(0.0, 1.0);
  }

  // Example hyperparameters, the last hyperparameter needs to be 0 (because the gaussian processes
  // are not available)
  const sgpp::optimization::IGGFARNHelper::Hyperparameters hyperparameters(
      0.49825808, 0.00042716, 0.00003296, 0.14564321, 0.11427497, 0);

  sgpp::optimization::IterativeGridGeneratorFullAdaptiveRitterNovak gridGen(
      *function, *grid, maxGridPoints, hyperparameters, domain);

  if (!gridGen.generate()) {
    throw std::runtime_error("Grid generation failed, exiting");
  }

  const sgpp::base::DataVector alpha = gridGen.getAlpha();

  // if success print all point on the grid with their coefficients
  std::cout << "Grid generation successful. Points: " << gridStorage.getSize() << std::endl;
  for (size_t i = 0; i < gridStorage.getSize(); ++i) {
    sgpp::base::DataVector point(dimension);
    gridStorage[i].getStandardCoordinates(point);
    std::cout << "Point " << i << ": " << point.toString() << ", Coefficient: " << alpha[i]
              << std::endl;
  }
  std::cout << "------------------------------------------------------" << std::endl;
}

/**
 * @brief Demonstrates IterativeGridGeneratorFullAdaptiveRitterNovak
 * with a Gaussian Process surrogate model and a stopping criterion.
 *
 * This is a more advanced example showing how to inject a surrogate model
 * (from libgp) into the grid generation process. It also demonstrates
 * how to set up a Monte Carlo-based stopping criterion.
 */
void gridGenerationFullAdaptiveRitterNovak() {
  std::cout << "--- Running gridGenerationFullAdaptiveRitterNovak (GP) ---" << std::endl;
  const size_t dimension = 2;
  const size_t bSplineDegree = 3;
  const size_t maxGridPoints = 50;
  const size_t stoppingCriterionMonteCarloPointCount = 50000;
  const bool useStoppingCriterion = true;
  std::unique_ptr<sgpp::base::ScalarFunction> function =
      std::make_unique<sgpp::optimization::test_problems::RastriginObjective>(dimension);

  std::unique_ptr<sgpp::base::Grid> grid(
      sgpp::base::Grid::createBsplineGrid(dimension, bSplineDegree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  std::vector<std::pair<double, double>> domain(dimension);
  for (size_t t = 0; t < dimension; ++t) {
    domain[t] = std::make_pair(0.0, 1.0);
  }

  const sgpp::optimization::IGGFARNHelper::Hyperparameters hyperparameters(
      0.49825808, 0.00042716, 0.00003296, 0.14564321, 0.11427497, 0.24136362);

  sgpp::optimization::IterativeGridGeneratorFullAdaptiveRitterNovak gridGen(
      *function, *grid, maxGridPoints, hyperparameters, domain);

  if (useStoppingCriterion) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    std::vector<sgpp::base::DataVector> validationPoints;
    for (size_t i = 0; i < stoppingCriterionMonteCarloPointCount; ++i) {
      sgpp::base::DataVector point(dimension);
      for (size_t t = 0; t < dimension; ++t) {
        const double rand01 = dis(gen);
        point[t] = domain[t].first + rand01 * (domain[t].second - domain[t].first);
      }
      validationPoints.push_back(point);
    }
    gridGen.activateStoppingCriterion(validationPoints);
  }

  if (!gridGen.generate()) {
    throw std::runtime_error("Grid generation failed, exiting");
  }

  const sgpp::base::DataVector alpha = gridGen.getAlpha();

  // if success print all point on the grid with their coefficients
  std::cout << "Grid generation successful. Points: " << gridStorage.getSize() << std::endl;
  for (size_t i = 0; i < gridStorage.getSize(); ++i) {
    sgpp::base::DataVector point(dimension);
    gridStorage[i].getStandardCoordinates(point);
    std::cout << "Point " << i << ": " << point.toString() << ", Coefficient: " << alpha[i]
              << std::endl;
  }
  std::cout << "------------------------------------------------------" << std::endl;
}

/**
 * @brief Main function to run the grid generation examples.
 */
int main(int argc, const char* argv[]) {
  try {
    gridGenerationRitterNovak();
    gridGenerationFullAdaptiveRitterNovakMinimal();
    if (sgpp::optimization::IterativeGridGeneratorFullAdaptiveRitterNovak::LIBGP_AVAILABLE) {
      gridGenerationFullAdaptiveRitterNovak();
    } else {
      std::cout << "Skipping gridGenerationFullAdaptiveRitterNovak (GP) example since SGpp was not "
                   "compiled with USE_LIBGP."
                << std::endl;
    }
  } catch (const std::exception& e) {
    std::cerr << "An exception occurred: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
