// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/integration/MCIntegrator.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridQuadraticSummationStrategy.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineRoutines.hpp>
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
  //  return 1;
  return v[0] * sin(v[0] + v[1]) * exp(v[1] * v[2]);
  //  return v[0] * sin(v[0]);
  //  return std::atan(50 * (v[0] - .35)) + M_PI / 2 + 4 * std::pow(v[1], 3) +
  //         std::exp(v[0] * v[1] - 1);
}

void interpolate(size_t maxlevel, size_t numDimensions, size_t degree, double& max_err,
                 double& L2_err) {
  sgpp::combigrid::MultiFunction func(f);
  max_err = 0.0;
  L2_err = 0.0;
  auto operation =
      sgpp::combigrid::CombigridOperation::createExpUniformBoundaryBsplineInterpolation(
          numDimensions, func, degree);

  //  auto operation =
  //  sgpp::combigrid::CombigridOperation::createLinearL2LejaPolynomialInterpolation(
  //      numDimensions, func, 2);

  double diff = 0.0;
  // generator generates num_points random points in [0,1]^numDimensions
  size_t num_points = 1000;
  sgpp::quadrature::NaiveSampleGenerator generator(numDimensions);
  sgpp::base::DataVector p(numDimensions, 0);
  std::vector<sgpp::base::DataVector> params;

  for (size_t i = 0; i < num_points; i++) {
    generator.getSample(p);
    params.push_back(p);
    diff = func(p) - operation->evaluate(maxlevel, p);
    max_err = (fabs(diff) > max_err) ? fabs(diff) : max_err;
    L2_err += diff * diff;
  }

  L2_err = sqrt(L2_err / static_cast<double>(num_points));

  //  std::cout << "# grid points: " << operation->numGridPoints() << " ";

  auto multiOperation =
      sgpp::combigrid::CombigridMultiOperation::createExpUniformBoundaryBsplineInterpolation(
          numDimensions, func, degree);

  //  auto multiOperation =
  //      sgpp::combigrid::CombigridMultiOperation::createLinearL2LejaPolynomialInterpolation(
  //          numDimensions, func, 2);

  double MultiL2_err = 0.0;
  auto result = multiOperation->evaluate(maxlevel, params);
  for (size_t i = 0; i < params.size(); ++i) {
    //    std::cout << params[i][0] << " " << result[i] << std::endl;
    diff = func(params[i]) - result[i];
    MultiL2_err += diff * diff;
  }

  MultiL2_err = sqrt(MultiL2_err / static_cast<double>(num_points));

  std::cout << std::scientific << maxlevel << " " << L2_err << " " << MultiL2_err << std::endl;

  //  std::string plotstr = "/home/rehmemk/SGS_Sync/Plotting/combigrid_bsplines/interpolant.dat";
  //  remove(plotstr.c_str());
  //  std::ofstream plotfile;
  //  plotfile.open(plotstr.c_str(), std::ios::app);
  //  plotfile << "#Basis functions  \n";
  //  std::vector<sgpp::base::DataVector> plotparams;
  //  for (double k = 0; k < num_points + 1; ++k) {
  //    p[0] = k / (double)num_points;
  //    plotparams.push_back(p);
  //  }
  //  auto multiresult = multiOperation->evaluate(maxlevel, plotparams);
  //  for (double k = 0; k < num_points + 1; k++) {
  //    plotfile << plotparams[k][0] << ", " << operation->evaluate(maxlevel, plotparams[k]) << ", "
  //             << multiresult[k] << "\n";
  //  }
  //
  //  plotfile.close();
}
/**
 * @param level level of the underlying 1D subspace
 * @return vector containing the integrals of all basisfunctions
 */
std::vector<double> integrateBasisFunctions(size_t level, size_t numDimensions, size_t degree) {
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

double integrate(size_t level, size_t numDimensions, size_t degree) {
  sgpp::combigrid::MultiFunction func(f);
  auto operation = sgpp::combigrid::CombigridOperation::createExpUniformBoundaryBsplineQuadrature(
      numDimensions, func, degree);
  return operation->evaluate(level);
}

double interpolate_and_integrate(size_t level, size_t numDimensions, size_t degree) {
  sgpp::combigrid::MultiFunction func(f);
  sgpp::combigrid::CombiHierarchies::Collection grids(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::Collection evaluators(numDimensions);
  evaluators[0] = sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree);
  evaluators[1] = sgpp::combigrid::CombiEvaluators::BSplineQuadrature(degree);

  auto operation = sgpp::combigrid::CombigridOperation::auxiliaryBsplineFunction(
      numDimensions, func, grids, evaluators, degree);

  sgpp::base::DataVector p(1, 0.3);
  double res = operation->evaluate(level, p);

  return res;
}

double variance(size_t level, size_t numDimensions, size_t degree) {
  sgpp::combigrid::MultiFunction func(f);
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineMixedQuadrature(degree));
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::WeightedRatioLevelManager());
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);
  bool exploitNesting = false;
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::QUADRATIC;
  auto varianceOperation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, evaluators, levelManager, gf, exploitNesting, summationStrategyType);

  sgpp::base::DataVector var = varianceOperation->evaluate(level);
  return var[0];
}

int main() {
  size_t numDimensions = 3;
  size_t degree = 3;
  size_t level = 4;

  // Interpolation
  //  sgpp::base::SGppStopwatch watch;
  //  watch.start();
  //  size_t minLevel = 2;
  //  size_t maxLevel = 2;
  //
  //  std::vector<double> maxErr(maxLevel + 1, 0);
  //  std::vector<double> L2Err(maxLevel + 1, 0);
  //  for (size_t l = minLevel; l < maxLevel + 1; l++) {
  //    interpolate(l, numDimensions, degree, maxErr[l], L2Err[l]);
  //    //    std::cout << "level: " << l << " max err " << maxErr[l] << " L2 err " << L2Err[l] <<
  //    //    std::endl;
  //  }

  //  std::cout << " Total Runtime: " << watch.stop() << " s" << std::endl;

  // Integration
  //  double integral = integrate(level, numDimensions, degree);
  //  std::cout << "integral:  " << integral << std::endl;

  // Integrate basis functions
  //  std::vector<double> integrals = integrateBasisFunctions(level, numDimensions, degree);
  //  std::cout << "------------------------------------" << std::endl;
  //  for (size_t i = 0; i < integrals.size(); i++) {
  //    std::cout << integrals[i] << " ";
  //  }

  //  Interpolate in one direction and integrate in the other
  //  double res = interpolate_and_integrate(level, numDimensions, degree);
  //  std::cout << res << std::endl;

  // Calculate variance (integral of func^2 respectively)
  sgpp::base::SGppStopwatch watch;
  double var = variance(level, numDimensions, degree);
  std::cout << "variance : " << var << "  (currently variance is simply integral(f^2)" << std::endl;
  std::cout << "Total Runtime: " << watch.stop() << " s" << std::endl;

  return 0;
}
