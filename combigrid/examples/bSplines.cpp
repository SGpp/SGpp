// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/QuadratureEvaluator.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

double f(sgpp::base::DataVector const& v) {
  //  return v[0];  // * v[0] * v[0];
  //  return sin(1. / (1. + v[0] * v[0])) * v[1];
  return v[0] * v[0];
  //  return std::atan(50 * (v[0] - .35)) + M_PI / 2 + 4 * std::pow(v[1], 3) +
  //         std::exp(v[0] * v[1] - 1);
}

void interpolate(size_t maxlevel, double& max_err, double& L2_err) {
  size_t numDimensions = 2;
  size_t degree = 5;
  sgpp::combigrid::MultiFunction func(f);
  auto operation =
      sgpp::combigrid::CombigridOperation::createExpUniformBoundaryBsplineInterpolation(
          numDimensions, func, degree);

  //  auto operation =
  //      sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(
  //          numDimensions, func);

  //  auto operation =
  //  sgpp::combigrid::CombigridOperation::createExpUniformBoundaryLinearInterpolation(
  //      numDimensions, func);

  double diff = 0.0;
  // generator generates num_points random points in [0,1]^numDimensions
  size_t num_points = 1000;
  sgpp::quadrature::NaiveSampleGenerator generator(numDimensions);
  sgpp::base::DataVector p(numDimensions, 0);

  for (size_t i = 0; i < num_points; i++) {
    generator.getSample(p);
    diff = fabs(operation->evaluate(maxlevel, p) - f(p));
    max_err = (diff > max_err) ? diff : max_err;
    L2_err += pow(diff, 2);
  }
  L2_err = sqrt(L2_err / static_cast<double>(num_points));

  std::cout << "# grid points: " << operation->numGridPoints() << " ";
}

/**
 * @param level level of the underlying 1D subspace
 * @return vector containing the integrals of all basisfunctions
 */
std::vector<double> integrate(size_t level) {
  size_t numDimensions = 1;
  size_t degree = 3;
  sgpp::combigrid::CombiHierarchies::Collection grids(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  bool needsOrdered = true;

  auto evaluator = sgpp::combigrid::CombiEvaluators::BSplineQuadrature(degree);
  evaluator->setGridPoints(grids[0]->getPoints(level, needsOrdered));
  std::vector<sgpp::combigrid::FloatScalarVector> weights = evaluator->getBasisValues();
  std::vector<double> integrals(weights.size());
  for (size_t i = 0; i < weights.size(); i++) {
    integrals[i] = weights[i].value();
  }
  return integrals;
}

double integrinterpolate(size_t level) {
  size_t numDimensions = 1;
  size_t degree = 3;
  sgpp::combigrid::MultiFunction func(f);
  auto operation = sgpp::combigrid::CombigridOperation::createExpUniformBoundaryBsplineIntegration(
      numDimensions, func, degree);
  return operation->evaluate(level);
}

int main() {
  // Interpolation
  //  sgpp::base::SGppStopwatch watch;
  //  watch.start();
  //  size_t minLevel = 0;
  //  size_t maxLevel = 8;
  //  std::vector<double> maxErr(maxLevel + 1, 0);
  //  std::vector<double> L2Err(maxLevel + 1, 0);
  //  for (size_t l = minLevel; l < maxLevel + 1; l++) {
  //    interpolate(l, maxErr[l], L2Err[l]);
  //    std::cout << "level: " << l << " max err " << maxErr[l] << " L2 err " << L2Err[l] <<
  //    std::endl;
  //  }
  //  std::cout << " Total Runtime: " << watch.stop() << " s" << std::endl;

  // Integrate B-Splines
  //  size_t level = 4;
  //  std::vector<double> integrals = integrate(level);
  //  for (size_t i = 0; i < integrals.size(); i++) {
  //    std::cout << integrals[i] << " ";
  //  }
  //  std::cout << "\n";

  // Integrate Interpolants
  size_t level = 6;
  double integral = integrinterpolate(level);
  std::cout << integral << std::endl;

  return 0;
}
