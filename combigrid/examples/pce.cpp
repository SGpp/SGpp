// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/pce/PolynomialChaosExpansion.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>

#include <cmath>

#include <iostream>
#include <vector>

int main() {
  sgpp::combigrid::Ishigami ishigamiModel;
  sgpp::combigrid::MultiFunction func(ishigamiModel.eval);

  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  for (size_t q = 6; q < 7; ++q) {
    // interpolate on adaptively refined grid
    auto op = sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(
        ishigamiModel.numDims, func);
    sgpp::combigrid::Stopwatch stopwatch;
    stopwatch.start();
    op->getLevelManager()->addRegularLevels(q);
    //    op->getLevelManager()->addLevelsAdaptiveParallel(1000, 4);
    stopwatch.log();
    // compute the variance
    stopwatch.start();
    sgpp::combigrid::PolynomialChaosExpansion pce(op, functionBasis);

    pce.getCombigridTensorOperation()->getLevelManager()->addLevelsAdaptiveByNumLevels(10);
    pce.getCombigridTensorOperation()->getLevelManager()->addLevelsAdaptiveByNumLevels(10);
    //    pce.getCombigridTensorOperation()->getLevelManager()->addLevelsAdaptiveByNumLevels(10);

    stopwatch.log();
    // compute the variance
    stopwatch.start();

    double mean = pce.mean();
    double variance = pce.variance();
    sgpp::base::DataVector sobol_indices;
    sgpp::base::DataVector total_sobol_indices;
    pce.getComponentSobolIndices(sobol_indices);
    pce.getTotalSobolIndices(total_sobol_indices);
    std::cout << "Time: " << stopwatch.elapsedSeconds() / static_cast<double>(op->numGridPoints())
              << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "#gp = " << op->getLevelManager()->numGridPoints() << std::endl;
    std::cout << "E(u) = " << mean << std::endl;
    std::cout << "Var(u) = " << variance << std::endl;
    std::cout << "Sobol indices = " << sobol_indices.toString() << std::endl;
    std::cout << "Sum Sobol indices = " << sobol_indices.sum() << std::endl;
    std::cout << "Total Sobol indices = " << total_sobol_indices.toString() << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
  }
}
