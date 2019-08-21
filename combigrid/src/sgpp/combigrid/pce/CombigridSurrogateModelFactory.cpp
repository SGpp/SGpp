// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/pce/CombigridSurrogateModelFactory.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/pce/PolynomialChaosExpansion.hpp>
#include <sgpp/combigrid/pce/PolynomialStochasticCollocation.hpp>

#include <sgpp/base/exception/application_exception.hpp>

namespace sgpp {
namespace combigrid {

std::shared_ptr<CombigridSurrogateModel> createCombigridSurrogateModel(
    CombigridSurrogateModelConfiguration& config) {
  switch (config.type) {
    case CombigridSurrogateModelsType::POLYNOMIAL_CHAOS_EXPANSION:
      return std::make_shared<PolynomialChaosExpansion>(config);
    case CombigridSurrogateModelsType::POLYNOMIAL_STOCHASTIC_COLLOCATION:
      return std::make_shared<PolynomialStochasticCollocation>(config);
    case CombigridSurrogateModelsType::BSPLINE_STOCHASTIC_COLLOCATION:
      return std::make_shared<PolynomialChaosExpansion>(config);
  }

  throw sgpp::base::application_exception(
      "createCombigridSurrogateModel: surrogate model type is unknown.");
}

} /* namespace combigrid */
} /* namespace sgpp */
