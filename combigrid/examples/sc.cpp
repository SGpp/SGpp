// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
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

  sgpp::combigrid::ProbabilityDensityFunction1DConfiguration config;
  config.pdfParameters.type_ = sgpp::combigrid::ProbabilityDensityFunctionType::BETA;
  config.pdfParameters.lowerBound_ = model.bounds[0];
  config.pdfParameters.upperBound_ = model.bounds[1];
  config.pdfParameters.alpha_ = model.alpha1;
  config.pdfParameters.beta_ = model.beta1;
  auto pdf1 = std::make_shared<sgpp::combigrid::ProbabilityDensityFunction1D>(config);

  config.pdfParameters.type_ = sgpp::combigrid::ProbabilityDensityFunctionType::BETA;
  config.pdfParameters.lowerBound_ = model.bounds[0];
  config.pdfParameters.upperBound_ = model.bounds[1];
  config.pdfParameters.alpha_ = model.alpha2;
  config.pdfParameters.beta_ = model.beta2;
  auto pdf2 = std::make_shared<sgpp::combigrid::ProbabilityDensityFunction1D>(config);

  sgpp::combigrid::WeightFunctionsCollection weightFunctions;
  weightFunctions.push_back(pdf1->getWeightFunction());
  weightFunctions.push_back(pdf2->getWeightFunction());

  sgpp::combigrid::MultiFunction func(model.eval);

  sgpp::combigrid::CombiHierarchies::Collection grids{
      model.numDims, sgpp::combigrid::CombiHierarchies::expClenshawCurtis()};

  sgpp::combigrid::CombiEvaluators::Collection evaluators;
  evaluators.push_back(sgpp::combigrid::CombiEvaluators::polynomialInterpolation());
  evaluators.push_back(sgpp::combigrid::CombiEvaluators::polynomialInterpolation());

  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager =
      std::make_shared<sgpp::combigrid::AveragingLevelManager>();

  auto op = std::make_shared<sgpp::combigrid::CombigridOperation>(grids, evaluators, levelManager,
                                                                  func, false);

  auto op_levelManager = op->getLevelManager();

  // initialize the surrogate model
  std::vector<double> bounds{pdf1->lowerBound(), pdf1->upperBound(), pdf2->lowerBound(),
                             pdf2->upperBound()};
  sgpp::combigrid::CombigridSurrogateModelConfiguration sc_config;
  sc_config.type = sgpp::combigrid::CombigridSurrogateModelsType::POLYNOMIAL_STOCHASTIC_COLLOCATION;
  sc_config.loadFromCombigridOperation(op);
  sc_config.weightFunctions = weightFunctions;
  sc_config.bounds = sgpp::base::DataVector(bounds);
  auto sc = sgpp::combigrid::createCombigridSurrogateModel(sc_config);

  sgpp::combigrid::CombigridSurrogateModelConfiguration update_config;

  sgpp::combigrid::Stopwatch stopwatch;
  for (size_t q = 0; q < 8; ++q) {
    op_levelManager->addRegularLevels(q);
    //    tensor_levelManager->addLevelsAdaptive(100);
    // compute the variance
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "compute mean and variance of stochastic collocation" << std::endl;
    std::cout << "#gp = " << op_levelManager->numGridPoints() << std::endl;
    stopwatch.start();
    update_config.levelStructure = op_levelManager->getLevelStructure();
    sc->updateConfig(update_config);
    double mean = sc->mean();
    double variance = sc->variance();
    stopwatch.log();
    std::cout << "|mu - E(u)|        = " << std::abs(model.mean - mean) << std::endl;
    std::cout << "|sigma^2 - Var(u)| = " << std::abs(model.variance - variance) << std::endl;
  }
}
