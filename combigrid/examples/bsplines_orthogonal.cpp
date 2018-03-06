// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineInterpolationCoefficientEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineInterpolationEvaluator.hpp>
#include <sgpp/combigrid/operation/onedim/BSplineQuadratureEvaluator.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModelFactory.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <cmath>

#include <iomanip>
#include <iostream>
#include <vector>

double u(double x) {
  return 4. * x * (1. - x);
  //  return std::log(std::exp(-x) * std::cos(4. * x * (1. - x)));
}

int main() {
  size_t q = 3;
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
  sgpp::optimization::sle_solver::Auto solver;
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  solver.solve(sle, rhs, coeffs);

  std::cout << "----------------------------------" << std::endl;
  std::cout << "coeffs" << std::endl;
  std::cout << "[";
  for (size_t i = 0; i < n - 1; i++) {
    std::cout << std::setw(15) << std::setprecision(10) << coeffs[i] << ", ";
  }
  std::cout << coeffs[coeffs.size() - 1] << "]" << std::endl;

  // compute Tensor for orthogonal Bspline basis
  sgpp::combigrid::BSplineInterpolationCoefficientEvaluator eval;
  eval.setGridPoints(gridPoints);
  std::vector<sgpp::combigrid::FloatTensorVector> basisValues = eval.getBasisValues();

  // solve interpolation problem -> multiply the rhs with the basisValues of the tensor
  size_t offset = getUniqueIndex(q, 0);
  sgpp::base::DataMatrix V_ocoeffs(n, n);
  for (size_t i = 0; i < basisValues.size(); i++) {
    // vector vector multiplication
    for (size_t j = 0; j < rhs.size(); j++) {
      double value = basisValues[j].get(sgpp::combigrid::MultiIndex{offset + i}).getValue();
      V_ocoeffs.set(i, j, value);
    }
  }

  std::cout << "----------------------------------" << std::endl;
  std::cout << "V" << std::endl;
  std::cout << "[";
  for (size_t i = 0; i < n; i++) {
    std::cout << "[";
    for (size_t j = 0; j < n; j++) {
      std::cout << std::setw(15) << std::setprecision(10) << V_ocoeffs(i, j) << ", ";
    }
    std::cout << "]" << std::endl << " ";
  }
  std::cout << "]" << std::endl;
  std::cout << "----------------------------------" << std::endl;
  std::cout << "rhs" << std::endl;
  std::cout << "[";
  for (size_t i = 0; i < n - 1; i++) {
    std::cout << std::setw(15) << std::setprecision(10) << rhs[i] << ", ";
  }
  std::cout << rhs[rhs.size() - 1] << "]" << std::endl;
  std::cout << "----------------------------------" << std::endl;

  // compute coeffs for standard bspline basis
  sgpp::base::DataVector ocoeffs(n);
  sgpp::optimization::FullSLE sle_ocoeffs(V_ocoeffs);
  solver.solve(sle_ocoeffs, rhs, ocoeffs);

  // compute the mean
  sgpp::combigrid::BSplineQuadratureEvaluator quadEval;
  quadEval.setGridPoints(gridPoints);
  std::vector<sgpp::combigrid::FloatScalarVector> quadResult = quadEval.getBasisValues();

  double mean = 0.0;
  for (size_t i = 0; i < quadResult.size(); i++) {
    mean += coeffs[i] * quadResult[i].getValue();
  }

  double variance = 0.0;
  std::cout << " ----------------------------------" << std::endl;
  for (size_t i = 0; i < ocoeffs.size(); i++) {
    std::cout << i << ": (" << (offset + i) << ") -> " << ocoeffs[i] << std::endl;
    variance += ocoeffs[i] * ocoeffs[i];
  }
  variance -= mean * mean;

  std::cout << "---------------------------------------------------------" << std::endl;
  std::cout << "E(u)   = " << mean << std::endl;
  std::cout << "Var(u) = " << variance << std::endl;
  std::cout << "---------------------------------------------------------" << std::endl;
}
