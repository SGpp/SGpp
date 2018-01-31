// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModelFactory.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineQuadratureEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineInterpolationCoefficientEvaluator.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>

#include <sgpp/optimization/sle/solver/UMFPACK.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>

#include <cmath>

#include <iostream>
#include <vector>

double u(double x) { return std::log(std::exp(-x) * std::cos(4. * x * (1. - x))); }

int main() {
  size_t q = 5;
  size_t n = static_cast<size_t>(std::pow(2, q)) + 1;
  double h = 1. / static_cast<double>(n - 1);

  sgpp::base::DataVector rhs(n);
  std::vector<double> gridPoints(n);
  for (size_t i = 0; i < n; i++) {
    gridPoints[i] = static_cast<double>(i) * h;
    rhs[i] = u(gridPoints[i]);
  }

  // compute interpolation matrix
  sgpp::combigrid::BSplineInterpolationEvaluator evaluator;
  evaluator.setGridPoints(gridPoints);

  // compute the interpolation matrix
  sgpp::base::DataMatrix V(n, n);
  for (size_t i = 0; i < n; ++i) {
    evaluator.setParameter(sgpp::combigrid::FloatScalarVector(gridPoints[i]));
    auto result = evaluator.getBasisValues();
    for (size_t j = 0; j < n; ++j) {
      V.set(i, j, result[j].getValue());
    }
  }

  // compute coeffs for standard bspline basis
  sgpp::base::DataVector coeffs(n);
  sgpp::optimization::FullSLE sle(V);
  sgpp::optimization::sle_solver::UMFPACK solver;
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  solver.solve(sle, rhs, coeffs);

  // compute Tensor for orthogonal Bspline basis
  sgpp::combigrid::BSplineInterpolationCoefficientEvaluator eval;
  eval.setGridPoints(gridPoints);
  std::vector<sgpp::combigrid::FloatTensorVector> basisValues = eval.getBasisValues();

  // solve interpolation problem -> multiply the rhs with the basisValues of the tensor
  sgpp::base::DataVector ocoeffs(n);
  for (size_t i = 0; i < basisValues.size(); i++) {
    // vector vector multiplication
    for (size_t j = 0; j < rhs.size(); j++) {
      ocoeffs[i] += basisValues[j].get(sgpp::combigrid::MultiIndex{i}).getValue() * rhs[j];
    }
  }

  // compute the mean
  sgpp::combigrid::BSplineQuadratureEvaluator quadEval;
  quadEval.setGridPoints(gridPoints);
  std::vector<sgpp::combigrid::FloatScalarVector> quadResult = quadEval.getBasisValues();

  double mean = 0.0;
  for (size_t i = 0; i < quadResult.size(); i++) {
    mean += coeffs[i] * quadResult[i].getValue();
  }

  double variance = 0.0;
  for (size_t i = 0; i < quadResult.size(); i++) {
    variance += ocoeffs[i] * ocoeffs[i];
  }
  variance -= mean * mean;

  std::cout << "---------------------------------------------------------" << std::endl;
  std::cout << "E(u)   = " << mean << std::endl;
  std::cout << "Var(u) = " << variance << std::endl;
  std::cout << "---------------------------------------------------------" << std::endl;
}
