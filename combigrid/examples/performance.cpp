// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridOptimizedPCESummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridPCESummationStrategy.hpp>

#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>

#include <sgpp/combigrid/pce/PolynomialChaosExpansion.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModelFactory.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>

#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>

#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <cmath>

#include <iostream>
#include <vector>

double f(sgpp::base::DataVector const &x) {
  double prod = 1.0;
  for (size_t dim = 0; dim < x.getSize(); ++dim) {
    prod *= exp(-x[dim] * x[dim]);
  }
  return prod;
}

auto func = sgpp::combigrid::MultiFunction(f);

void testEfficientPCE() {
  size_t d = 5;

  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  config.polyParameters.alpha_ = 10.0;
  config.polyParameters.beta_ = 5.0;

  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  for (size_t q = 2; q < 8; ++q) {
    // interpolate on adaptively refined grid
    // auto op = sgpp::combigrid::CombigridOperation::createLinearLejaPolynomialInterpolation(d,
    // func);
    auto op = sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(
        d, func);
    sgpp::combigrid::Stopwatch stopwatch;
    stopwatch.start();
    op->getLevelManager()->addRegularLevels(q);
    //    op->getLevelManager()->addLevelsAdaptiveParallel(1000, 4);
    stopwatch.log();
    // compute the variance

    // initialize the surrogate model
    sgpp::combigrid::CombigridSurrogateModelConfiguration config;
    config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_CHAOS_EXPANSION;
    config.loadFromCombigridOperation(op);
    config.basisFunction = functionBasis;

    stopwatch.start();
    auto pce = sgpp::combigrid::createCombigridSurrogateModel(config);
    double mean = pce->mean();
    double variance = pce->variance();
    sgpp::base::DataVector sobol_indices;
    sgpp::base::DataVector total_sobol_indices;
    // pce.getComponentSobolIndices(sobol_indices);
    // pce.getTotalSobolIndices(total_sobol_indices);
    std::cout << "Time: " << stopwatch.elapsedSeconds() / static_cast<double>(op->numGridPoints())
              << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "#gp = " << op->getLevelManager()->numGridPoints() << std::endl;
    std::cout << "E(u) = " << mean << std::endl;
    std::cout << "Var(u) = " << variance << std::endl;
    // std::cout << "Sobol indices = " << sobol_indices.toString() << std::endl;
    // std::cout << "Sum Sobol indices = " << sobol_indices.sum() << std::endl;
    // std::cout << "Total Sobol indices = " << total_sobol_indices.toString() << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
  }
}

void testFullGridPCE(sgpp::combigrid::FullGridSummationStrategyType summationType) {
  size_t d = 5;

  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  config.polyParameters.alpha_ = 10.0;
  config.polyParameters.beta_ = 5.0;

  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);
  std::vector<std::shared_ptr<sgpp::combigrid::AbstractInfiniteFunctionBasis1D>> functionBases(
      d, functionBasis);

  for (size_t q = 2; q < 8; ++q) {
    // interpolate on adaptively refined grid
    // std::vector<std::shared_ptr<sgpp::combigrid::AbstractPointHierarchy>> pointHierarchies(
    //     d, sgpp::combigrid::CombiHierarchies::linearLeja());
    std::vector<std::shared_ptr<sgpp::combigrid::AbstractPointHierarchy>> pointHierarchies(
        d, sgpp::combigrid::CombiHierarchies::expClenshawCurtis());
    std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> storage =
        std::make_shared<sgpp::combigrid::CombigridTreeStorage>(pointHierarchies, true, func);

    auto levelManager = std::make_shared<sgpp::combigrid::AveragingLevelManager>();

    sgpp::combigrid::CombiEvaluators::TensorCollection tensorEvaluators(
        d, sgpp::combigrid::CombiEvaluators::tensorInterpolation(functionBasis));

    auto op = std::make_shared<sgpp::combigrid::CombigridTensorOperation>(
        pointHierarchies, tensorEvaluators, levelManager, storage, summationType);

    sgpp::combigrid::Stopwatch stopwatch;
    stopwatch.start();
    levelManager->addRegularLevels(q);
    auto result = op->getResult();
    //    levelManager->addLevelsAdaptiveParallel(1000, 4);
    std::cout << "Time: "
              << stopwatch.elapsedSeconds() / static_cast<double>(levelManager->numGridPoints())
              << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "#gp = " << levelManager->numGridPoints() << std::endl;
    std::cout << "E(u) = " << result.get(sgpp::combigrid::MultiIndex(d, 0)) << "\n";
    std::cout << "Var(u) = " << std::pow(result.norm(), 2) << "\n";
    std::cout << "---------------------------------------------------------" << std::endl;
  }
}

void oldSpeedtest() {
  size_t dim = 6;
  auto op =
      sgpp::combigrid::CombigridOperation::createLinearL2LejaPolynomialInterpolation(dim, func);
  sgpp::base::DataVector parameter(std::vector<double>{0.1, 0.2, 0.3, 0.4, 0.5, 0.6});

  {
    sgpp::combigrid::Stopwatch sw;
    op->evaluate(10, parameter);
    sw.log();
    std::cout << "Number of grid points: " << op->numGridPoints() << "\n";
  }

  {
    sgpp::combigrid::Stopwatch sw;
    for (size_t i = 0; i < 10; ++i) {
      op->evaluate(10, parameter);
    }
    sw.log();
    std::cout << "Number of grid points: " << op->numGridPoints() << "\n";
  }
}

int main() {
  try {
    // testEfficientPCE();
    testFullGridPCE(sgpp::combigrid::FullGridSummationStrategyType::ONEDSUBSPACEPCE);
    testFullGridPCE(sgpp::combigrid::FullGridSummationStrategyType::FULLSUBSPACEDPCE);
  }
  catch (sgpp::base::generation_exception& exc)  {
    std::cout << "Exception: " << exc.what() << std::endl;
    std::cout << "Skipping example..." << std::endl;
  }

  return 0;
}
