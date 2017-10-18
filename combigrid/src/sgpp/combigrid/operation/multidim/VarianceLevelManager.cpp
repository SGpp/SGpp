// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/VarianceLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>

#include <algorithm>
#include <vector>

namespace sgpp {
namespace combigrid {

VarianceLevelManager::VarianceLevelManager(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
    std::shared_ptr<AbstractCombigridStorage> storage,
    std::vector<std::shared_ptr<OrthogonalPolynomialBasis1D>> basisFunctions)
    : LevelManager(), pointHierarchies(pointHierarchies), storage(storage) {
  std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>> evaluatorPrototypes(
      basisFunctions.size());
  size_t i = 0;
  for (auto& basisFunction : basisFunctions) {
    evaluatorPrototypes[i++] = CombiEvaluators::tensorInterpolation(basisFunction);
  }
  fullGridEval = std::make_shared<FullGridLinearCallbackEvaluator<FloatTensorVector>>(
      storage, evaluatorPrototypes, pointHierarchies);
}

VarianceLevelManager::~VarianceLevelManager() {}

double VarianceLevelManager::computePriority(const MultiIndex& level) {
  auto predecessors = getPredecessors(level);

  if (predecessors.empty()) {
    return 0.0;
  }

  double sum = 0.0;

  for (auto& predLevel : predecessors) {
    // compute the local variance
    double local_variance = std::pow(fullGridEval->eval(predLevel).norm(), 2);
    sum += local_variance;
  }

  return sum / static_cast<double>(predecessors.size());
}

std::shared_ptr<LevelManager> VarianceLevelManager::clone() {
  return std::make_shared<VarianceLevelManager>(*this);
}

} /* namespace combigrid */
} /* namespace sgpp */
