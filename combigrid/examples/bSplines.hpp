// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryCombigridBasis.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/common/GridConversion.hpp>
#include <sgpp/combigrid/integration/MCIntegrator.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/OperationConfiguration.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/RegularLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridGridBasedEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridLinearSummationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/FullGridQuadraticSummationStrategy.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineScalarProductEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/PolynomialQuadratureEvaluator.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "../src/sgpp/combigrid/utils/BSplineRoutines.hpp"

size_t numDimensions = 2;
double f(sgpp::base::DataVector const& v) {
  //  return 1;
  //  return v[0] * sin(v[0] + v[1]) * exp(v[1] * v[2]);
  //  return std::atan(50 * (v[0] - .35));
  return std::atan(50 * (v[0] - .35)) + M_PI / 2 + 4 * std::pow(v[1], 3) +
         std::exp(v[0] * v[1] - 1);
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

double integrateSquare(size_t level, size_t numDimensions, size_t degree) {
  sgpp::combigrid::MultiFunction func(f);
  auto quadratureSquareOperation =
      sgpp::combigrid::CombigridMultiOperation::createExpUniformBoundaryBsplineSquareQuadrature(
          numDimensions, func, degree);

  sgpp::base::DataVector res = quadratureSquareOperation->evaluate(level);
  return res[0];
}

double variance(size_t level, size_t numDimensions, size_t degree) {
  sgpp::combigrid::MultiFunction func(f);
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());

  sgpp::combigrid::EvaluatorConfiguration evalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineScalarProduct, degree);

  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(evalConfig));
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::AveragingLevelManager());
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);
  bool exploitNesting = false;
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::VARIANCE;
  auto operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, evaluators, levelManager, gf, exploitNesting, summationStrategyType);
  sgpp::base::DataVector res = operation->evaluate(level);
  return res[0];
}

double interpolateOneLevel(sgpp::combigrid::MultiIndex level, size_t degree) {
  size_t numDimensions = level.size();
  sgpp::combigrid::MultiFunction func(f);
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::Collection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);
  bool exploitNesting = false;
  auto summationStrategyType = sgpp::combigrid::FullGridSummationStrategyType::LINEAR;

  auto storage = std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage>(
      new sgpp::combigrid::CombigridTreeStorage(pointHierarchies, exploitNesting));

  std::shared_ptr<sgpp::combigrid::AbstractFullGridEvaluator<sgpp::combigrid::FloatScalarVector>>
      fullGridEval = std::make_shared<
          sgpp::combigrid::FullGridGridBasedEvaluator<sgpp::combigrid::FloatScalarVector>>(
          storage, evaluators, pointHierarchies, gf, summationStrategyType);

  std::vector<sgpp::combigrid::FloatScalarVector> params;
  params.push_back(sgpp::combigrid::FloatScalarVector(0.2));
  params.push_back(sgpp::combigrid::FloatScalarVector(0.7));

  fullGridEval->setParameters(params);

  auto result = fullGridEval->eval(level);

  double res = result.value();
  return res;
}

double integrateOneLevel(sgpp::combigrid::MultiIndex level, size_t degree) {
  size_t numDimensions = level.size();
  sgpp::combigrid::MultiFunction func(f);
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::Collection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineQuadrature(degree));
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);
  bool exploitNesting = false;
  auto summationStrategyType = sgpp::combigrid::FullGridSummationStrategyType::LINEAR;

  auto storage = std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage>(
      new sgpp::combigrid::CombigridTreeStorage(pointHierarchies, exploitNesting));

  std::shared_ptr<sgpp::combigrid::AbstractFullGridEvaluator<sgpp::combigrid::FloatScalarVector>>
      fullGridEval = std::make_shared<
          sgpp::combigrid::FullGridGridBasedEvaluator<sgpp::combigrid::FloatScalarVector>>(
          storage, evaluators, pointHierarchies, gf, summationStrategyType);

  auto result = fullGridEval->eval(level);

  double res = result.value();
  return res;
}

void printLevelstructure(std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> levelstructure) {
  auto it = levelstructure->getStoredDataIterator();
  while (it->isValid()) {
    sgpp::combigrid::MultiIndex index = it->getMultiIndex();
    for (auto& i : index) {
      std::cout << i << " ";
    }
    std::cout << "\n";
    it->moveToNext();
  }
}

std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> createVarianceLevelStructure(
    size_t numlevels, size_t degree, sgpp::combigrid::CombiHierarchies::Collection pointHierarchies,
    sgpp::combigrid::GridFunction gf, bool exploitNesting) {
  // Operation creating a variance adaptive level structure
  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineScalarProduct, degree);

  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
      new sgpp::combigrid::AveragingLevelManager());

  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::VARIANCE;
  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, levelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);

  // create level structure
  // First step is to guarantee the existence of level (1,..,1). Otherwise a conversion to an SG
  // grid wouldn't be possible
  Operation->getLevelManager()->addRegularLevels(numDimensions);
  Operation->getLevelManager()->addLevelsAdaptiveByNumLevels(numlevels);
  auto levelStructure = Operation->getLevelManager()->getLevelStructure();
  return levelStructure;
}

double interpolateVarianceAdaptively(size_t numlevels, size_t degree) {
  sgpp::combigrid::MultiFunction func(f);
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  bool exploitNesting = false;
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);

  auto levelStructure =
      createVarianceLevelStructure(numlevels, degree, pointHierarchies, gf, exploitNesting);
  printLevelstructure(levelStructure);

  // create interpolation Operation
  sgpp::combigrid::EvaluatorConfiguration evalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, degree);
  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(evalConfig));
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  std::shared_ptr<sgpp::combigrid::LevelManager> dummylevelManager(
      new sgpp::combigrid::AveragingLevelManager());

  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, evaluators, dummylevelManager, gf, exploitNesting, summationStrategyType);

  // calculate error
  size_t num_MCpoints = 1000;
  sgpp::quadrature::NaiveSampleGenerator generator(numDimensions);
  sgpp::base::DataVector p(numDimensions, 0);
  sgpp::base::DataMatrix params(numDimensions, 0);
  sgpp::base::DataVector funceval;
  for (size_t i = 0; i < num_MCpoints; i++) {
    generator.getSample(p);
    params.appendCol(p);
    funceval.push_back(func(p));
  }
  //  diff = func(p) - Operation->evaluate(maxlevel, p);
  //  max_err = (fabs(diff) > max_err) ? fabs(diff) : max_err;
  //  L2_err += diff * diff;

  Operation->setParameters(params);
  Operation->getLevelManager()->addLevelsFromStructure(levelStructure);
  sgpp::base::DataVector interpoleval = Operation->getResult();
  // std::cout << params.size() <<  " " << interpoleval.size() << " " << funceval.size() <<
  // std::endl;
  interpoleval.sub(funceval);
  return interpoleval.l2Norm();
}

void BsplineTimeComparison() {
  // deg = 5, level = 27:
  // SNak:  1.30 s
  // Combi: 1.36s
  // Kein nennenswerter Vorteil SNak erkennbar
  size_t degree = 5;
  std::vector<double> SNakeval;
  std::vector<double> Combieval;
  sgpp::base::SNakBsplineBoundaryCombigridBase basis(degree);
  unsigned int level = 27;
  sgpp::combigrid::Stopwatch watch;
  watch.start();
  double x = 0.714;
  for (unsigned int i = 1; i <= std::pow(2, level) - 1; i++) {
    SNakeval.push_back(basis.eval(level, i, x));
  }
  std::cout << watch.elapsedSeconds() << " " << std::endl;

  std::vector<double> xValues;
  for (double i = 1; i <= std::pow(2, level) - 1; i = i + 2) {
    xValues.push_back(i / std::pow(2, level));
  }
  watch.start();
  std::vector<double> xi = createNakKnots(xValues, degree);
  for (unsigned int i = 1; i <= std::pow(2, level) - 1; i++) {
    Combieval.push_back(sgpp::combigrid::nonUniformBSpline(x, degree, i, xi));
  }
  std::cout << watch.elapsedSeconds() << " " << std::endl;

  //  for (size_t i = 0; i < Combieval.size(); i++) {
  //    std::cout << fabs(SNakeval[i] - Combieval[i]) << std::endl;
  //  }
}
