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
  return v[0];
  //    return v[0] * sin(v[1]) ;
  //  return std::atan(50 * (v[0] - .35)) + M_PI / 2 + 4 * std::pow(v[1], 3) +
  //         std::exp(v[0] * v[1] - 1);

  //  double prod = 1.0;
  //
  //  for (size_t i = 0; i < v.getSize(); ++i) {
  //    double x = v[i];
  //    prod *= std::cos(-x * x / static_cast<double>((i + 1) * (i + 1)));
  //  }
  //  return prod;
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
  size_t num_points = 10;
  sgpp::quadrature::NaiveSampleGenerator generator(numDimensions);
  sgpp::base::DataVector p(numDimensions, 0);
  std::vector<sgpp::base::DataVector> params;

  for (size_t i = 0; i < num_points; i++) {
    generator.getSample(p);
    params.push_back(p);
    diff = func(p) - operation->evaluate(maxlevel, p);
    std::cout << diff << " ";
    max_err = (fabs(diff) > max_err) ? fabs(diff) : max_err;
    L2_err += diff * diff;
  }
  std::cout << "\n";
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
    diff = func(params[i]) - result[i];
    std::cout << diff << " ";
    MultiL2_err += diff * diff;
  }
  std::cout << "\n \n";

  MultiL2_err = sqrt(MultiL2_err / static_cast<double>(num_points));

  //  std::cout << std::scientific << maxlevel << " " << L2_err << " " << MultiL2_err << std::endl;
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
//
//// This example is not ready yet. FullGridQuadratureSummationstrategy must be finished
// double variance(size_t level, size_t numDimensions, size_t degree) {
//  sgpp::combigrid::MultiFunction func(f);
//  //  std::vector<std::shared_ptr<sgpp::combigrid::AbstractPointHierarchy>> pointHierarchies(
//  //      numDimensions);
//  //  std::vector<
//  // std::shared_ptr<sgpp::combigrid::AbstractLinearEvaluator<sgpp::combigrid::FloatArrayVector>>>
//  //      evaluators(numDimensions);
//  //  evaluators[0] =
//  //      std::make_shared<ArrayEvaluator<sgpp::combigrid::BSplineQuadratureMixedEvaluator>>;
//
//  sgpp::combigrid::CombiHierarchies::Collection grids(
//      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
//  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
//      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineMixedQuadrature(degree));
//  auto storage = std::make_shared<sgpp::combigrid::CombigridTreeStorage>(grids, func);
//
//  // stores the values of the objective function
//  auto funcStorage = std::make_shared<sgpp::combigrid::CombigridTreeStorage>(grids, func);
//
//  sgpp::combigrid::GridFunction gf([=](std::shared_ptr<sgpp::combigrid::TensorGrid> grid) {
//    sgpp::combigrid::CombiEvaluators::Collection interpolEvaluators(
//        numDimensions, sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));
//    size_t numDimensions = grid->getDimension();
//    auto coefficientTree = std::make_shared<sgpp::combigrid::TreeStorage<double>>(numDimensions);
//    auto level = grid->getLevel();
//    std::vector<size_t> numGridPointsVec = grid->numPoints();
//    size_t numGridPoints = 1;
//    for (size_t i = 0; i < numGridPointsVec.size(); i++) {
//      numGridPoints *= numGridPointsVec[i];
//    }
//
//    sgpp::combigrid::CombiEvaluators::Collection evalCopy(numDimensions);
//    for (size_t dim = 0; dim < numDimensions; ++dim) {
//      evalCopy[dim] = interpolEvaluators[dim]->cloneLinear();
//      bool needsSorted = evalCopy[dim]->needsOrderedPoints();
//      auto gridPoints = grids[dim]->getPoints(level[dim], needsSorted);
//      evalCopy[dim]->setGridPoints(gridPoints);
//    }
//    sgpp::base::DataMatrix A(numGridPoints, numGridPoints);
//    sgpp::base::DataVector coefficients_sle(numGridPoints);
//    sgpp::base::DataVector functionValues(numGridPoints);
//
//    // Creates an iterator that yields the multi-indices of all grid points in the grid.
//    sgpp::combigrid::MultiIndexIterator it(grid->numPoints());
//    auto funcIter =
//        funcStorage->getGuidedIterator(level, it, std::vector<bool>(numDimensions, true));
//
//    for (size_t ixEvalPoints = 0; funcIter->isValid(); ++ixEvalPoints, funcIter->moveToNext()) {
//      auto gridPoint = grid->getGridPoint(funcIter->getMultiIndex());
//      functionValues[ixEvalPoints] = funcIter->value();
//
//      std::vector<std::vector<double>> basisValues;
//      for (size_t dim = 0; dim < numDimensions; ++dim) {
//        evalCopy[dim]->setParameter(sgpp::combigrid::FloatScalarVector(gridPoint[dim]));
//        auto basisValues1D = evalCopy[dim]->getBasisValues();
//        // basis values at gridPoint
//        std::vector<double> basisValues1D_vec(basisValues1D.size());
//        for (size_t i = 0; i < basisValues1D.size(); i++) {
//          basisValues1D_vec[i] = basisValues1D[i].value();
//        }
//        basisValues.push_back(basisValues1D_vec);
//      }
//
//      sgpp::combigrid::MultiIndexIterator innerIter(grid->numPoints());
//      for (size_t ixBasisFunctions = 0; innerIter.isValid();
//           ++ixBasisFunctions, innerIter.moveToNext()) {
//        double splineValue = 1.0;
//        auto innerIndex = innerIter.getMultiIndex();
//        for (size_t dim = 0; dim < numDimensions; ++dim) {
//          splineValue *= basisValues[dim][innerIndex[dim]];
//        }
//        A.set(ixEvalPoints, ixBasisFunctions, splineValue);
//      }
//    }
//
//    sgpp::optimization::FullSLE sle(A);
//    sgpp::optimization::sle_solver::Auto solver;
//    sgpp::optimization::Printer::getInstance().setVerbosity(-1);
//    bool solved = solver.solve(sle, functionValues, coefficients_sle);
//
//    /*std::cout << A.toString() << std::endl;
//    std::cout << "fct: ";
//    for (size_t i = 0; i < functionValues.size(); i++) {
//      std::cout << functionValues[i] << " ";
//    }
//    std::cout << "\ncoeff: ";
//    for (size_t i = 0; i < coefficients_sle.size(); i++) {
//      std::cout << coefficients_sle[i] << " ";
//    }
//    std::cout << "\n";
//    std::cout << "--------" << std::endl;
//    */
//    if (!solved) {
//      exit(-1);
//    }
//
//    it.reset();
//    for (size_t vecIndex = 0; it.isValid(); ++vecIndex, it.moveToNext()) {
//      coefficientTree->set(it.getMultiIndex(), coefficients_sle[vecIndex]);
//    }
//
//    return coefficientTree;
//  });
//
//  // Kann keine normale CombigridOperation für MixedQuadrature nehmen, da CombigridOperation als
//  // Templateparameter FloatScalarVector benutzt. MixedQuadrature benötigt jedoch
//  FloarArrayVector.
//  // Ist es sinnvoll CombigridOperation zu templatisieren? Es soll ja ein einfaches Interface für
//  // die Operationen bieten. Gibt es noch andere Anwendungsfälle wo eine templatisierte Version
//  von
//  // Nutzen wäre?
//  // Um die MixedQuadrature hier durchzuführen überspringe ich den Operation-Hilfsschritt und
//  nutze
//  // direkt addRegularLevels (setParameters in diesem Fall nicht nötig, da Quadratur)
//
//  // QUATSCH. Nutze MultiOperation!
//
//  auto summationStrategy = std::make_shared<
//      sgpp::combigrid::FullGridQuadraticSummationStrategy<sgpp::combigrid::FloatArrayVector>>(
//      storage, evaluators, grids, gf);
//  auto combiGridEval =
//      std::make_shared<sgpp::combigrid::CombigridEvaluator<sgpp::combigrid::FloatArrayVector>>(
//          numDimensions, summationStrategy);
//  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
//      new sgpp::combigrid::WeightedRatioLevelManager(combiGridEval));
//
//  levelManager->addRegularLevels(level);
//  sgpp::combigrid::FloatArrayVector innerSums = combiGridEval->getValue();
//
//  double var = 0.0;
//  for (size_t i = 0; i < innerSums.size(); i++) {
//    //    std::cout << innerSums[i].value() << " ";
//    var += innerSums[i].value();
//  }
//  //  std::cout << "\n";
//
//  return var;
//}

int main() {
  size_t numDimensions = 1;
  size_t degree = 3;
  size_t level = 3;

  // Interpolation
  sgpp::base::SGppStopwatch watch;
  watch.start();
  size_t minLevel = 0;
  size_t maxLevel = 7;

  std::vector<double> maxErr(maxLevel + 1, 0);
  std::vector<double> L2Err(maxLevel + 1, 0);
  for (size_t l = minLevel; l < maxLevel + 1; l++) {
    interpolate(l, numDimensions, degree, maxErr[l], L2Err[l]);
    //    std::cout << "level: " << l << " max err " << maxErr[l] << " L2 err " << L2Err[l] <<
    //    std::endl;
  }

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
  //  double var = variance(level, numDimensions, degree);
  //  std::cout << "variance : " << var << "  (currently variance is simply integral(f^2)" <<
  //  std::endl;

  return 0;
}
