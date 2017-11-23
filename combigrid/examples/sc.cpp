// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/pce/PolynomialStochasticCollocation.hpp>
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
  sgpp::base::DataVector pi(v.size());
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

  auto op =
      sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(d, func);
  auto op_levelManager = op->getLevelManager();
  sgpp::combigrid::PolynomialStochasticCollocation sc(op, functionBasis);
  auto tensor_levelManager = sc.getCombigridTensorOperation()->getLevelManager();
  for (size_t q = 0; q < 8; ++q) {
    sgpp::combigrid::Stopwatch stopwatch;
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "add regular levels " << q << " to interpolation operation" << std::endl;
    stopwatch.start();
    op_levelManager->addRegularLevels(q);
    stopwatch.log();
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "add regular levels " << q << " to tensor operation" << std::endl;
    stopwatch.start();
    tensor_levelManager->addRegularLevels(q);
    stopwatch.log();
    std::cout << "#gp = " << op_levelManager->numGridPoints() << std::endl;
    // compute the variance
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "compute mean and variance of stochastic collocation" << std::endl;
    stopwatch.start();
    double mean = sc.mean();
    double variance = sc.variance();
    stopwatch.log();
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "#gp = " << op_levelManager->numGridPoints() << std::endl;
    std::cout << "|mu - E(u)|        = " << std::abs(3.5 - mean) << std::endl;
    std::cout << "|sigma^2 - Var(u)| = " << std::abs(f_variance() - variance) << std::endl;
  }
}
