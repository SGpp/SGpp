// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialScalarProductEvaluator.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModelFactory.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <cmath>
#include <iostream>
#include <vector>

int main() {
  sgpp::combigrid::AtanBeta model;
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = model.bounds[0];
  config.polyParameters.upperBound_ = model.bounds[1];
  config.polyParameters.alpha_ = model.alpha1;
  config.polyParameters.beta_ = model.beta1;
  auto functionBasis1 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::JACOBI;
  config.polyParameters.lowerBound_ = model.bounds[0];
  config.polyParameters.upperBound_ = model.bounds[1];
  config.polyParameters.alpha_ = model.alpha2;
  config.polyParameters.beta_ = model.beta2;
  auto functionBasis2 = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  sgpp::combigrid::OrthogonalBasisFunctionsCollection basisFunctions;
  basisFunctions.push_back(functionBasis1);
  basisFunctions.push_back(functionBasis2);

  sgpp::combigrid::MultiFunction func(model.eval);

  sgpp::combigrid::CombiHierarchies::Collection grids{
      model.numDims, sgpp::combigrid::CombiHierarchies::expClenshawCurtis()};

  sgpp::combigrid::CombiEvaluators::TensorCollection evaluators;
  evaluators.push_back(sgpp::combigrid::CombiEvaluators::tensorInterpolation(functionBasis1));
  evaluators.push_back(sgpp::combigrid::CombiEvaluators::tensorInterpolation(functionBasis2));

  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager =
      std::make_shared<sgpp::combigrid::AveragingLevelManager>();

  auto tensor_op = std::make_shared<sgpp::combigrid::CombigridTensorOperation>(
      grids, evaluators, levelManager, func, false);

  auto tensor_levelManager = tensor_op->getLevelManager();

  // initialize the surrogate model
  sgpp::combigrid::CombigridSurrogateModelConfiguration sc_config;
  sc_config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_STOCHASTIC_COLLOCATION;
  sc_config.loadFromCombigridOperation(tensor_op);
  sc_config.basisFunctions = basisFunctions;
  auto sc = sgpp::combigrid::createCombigridSurrogateModel(sc_config);

  sgpp::combigrid::CombigridSurrogateModelConfiguration update_config;

  sgpp::combigrid::Stopwatch stopwatch;
  for (size_t q = 0; q < 6; ++q) {
    tensor_levelManager->addRegularLevels(q);
    //    tensor_levelManager->addLevelsAdaptive(100);
    // compute the variance
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "compute mean and variance of stochastic collocation" << std::endl;
    std::cout << "#gp = " << tensor_levelManager->numGridPoints() << std::endl;
    stopwatch.start();
    update_config.levelStructure = tensor_levelManager->getLevelStructure();
    sc->updateConfig(update_config);
    double mean = sc->mean();
    double variance = sc->variance();
    stopwatch.log();
    std::cout << "|mu - E(u)|        = " << std::abs(model.mean - mean) << std::endl;
    std::cout << "|sigma^2 - Var(u)| = " << std::abs(model.variance - variance) << std::endl;
  }
}
