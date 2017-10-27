// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <cmath>

#include <iostream>
#include <vector>

double f(sgpp::base::DataVector const &v) {
  // return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  // return v[0] * v[0] + v[1] * v[1];
  //  double ans = 1.0;
  //  for (size_t idim = 0; idim < v.getSize(); idim++) {
  //    ans *= 4 * v[idim] * (1.0 - v[idim]);
  //  }
  //  return ans;
  return exp(3 * v[0] * v[0] + v[1]) * atan(10 * v[2]) + sin(3 * v[1] + v[2]);
}

int main() {
  auto func = sgpp::combigrid::MultiFunction(f);
  size_t d = 3;

  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  config.polyParameters.alpha_ = 10.0;
  config.polyParameters.beta_ = 5.0;

  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  for (size_t q = 0; q <= 8; ++q) {
    // interpolate on adaptively refined grid
    auto op = sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(
        d, func);
    sgpp::combigrid::Stopwatch stopwatch;
    stopwatch.start();
    op->getLevelManager()->addRegularLevels(q);
    op->getLevelManager()->addLevelsAdaptive(1000);
    stopwatch.log();
    // compute the variance
    stopwatch.start();
    auto tensor_op =
        sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
            op->getPointHierarchies(), op->getStorage(), op->getLevelManager(), functionBasis);
    auto tensor_result = tensor_op->getResult();
    std::cout << "Time " << stopwatch.elapsedSeconds() / static_cast<double>(op->numGridPoints())
              << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "#gp = " << op->getLevelManager()->numGridPoints() << std::endl;
    std::cout << "E(u) = " << tensor_result.get(sgpp::combigrid::MultiIndex(d, 0)) << std::endl;
    std::cout << "Var(u) = " << std::pow(tensor_result.norm(), 2) << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
  }
}
