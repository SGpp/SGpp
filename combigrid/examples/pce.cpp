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
  //  return exp(3 * v[0] * v[0] + v[1]) * atan(10 * v[2]) + sin(3 * v[1] + v[2]);

  // Ishigami function in (-pi, pi)
  double a = 7.0, b = 0.1;
  sgpp::base::DataVector x(v);
  sgpp::base::DataVector pi(v.getSize());
  pi.setAll(M_PI);
  x.mult(2 * M_PI);
  x.sub(pi);
  return std::sin(x[0]) + a * std::sin(x[1]) * std::sin(x[1]) +
         b * x[2] * x[2] * x[2] * x[2] * std::sin(x[0]);
}

double f_variance(double a = 7., double b = 0.1) {
  double pi_4 = M_PI * M_PI * M_PI * M_PI;
  return a * a / 8. + b * pi_4 / 5 + b * b * pi_4 * pi_4 / 18. + 0.5;
}

int main() {
  auto func = sgpp::combigrid::MultiFunction(f);
  size_t d = 3;

  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration config;
  config.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  config.polyParameters.alpha_ = 10.0;
  config.polyParameters.beta_ = 5.0;

  auto functionBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(config);

  for (size_t q = 5; q < 6; ++q) {
    // interpolate on adaptively refined grid
    auto op = sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(
        d, func);
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
