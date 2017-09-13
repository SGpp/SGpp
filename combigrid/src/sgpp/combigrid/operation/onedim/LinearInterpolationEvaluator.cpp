// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <vector>

namespace sgpp {
namespace combigrid {

LinearInterpolationEvaluator::LinearInterpolationEvaluator()
    : evaluationPoint(0.0), basisValues(), xValues() {}

LinearInterpolationEvaluator::~LinearInterpolationEvaluator() {}

LinearInterpolationEvaluator::LinearInterpolationEvaluator(
    const LinearInterpolationEvaluator& other)
    : evaluationPoint(other.evaluationPoint),
      basisValues(other.basisValues),
      xValues(other.xValues) {}

void LinearInterpolationEvaluator::computeBasisValues() {
  size_t numPoints = xValues.size();

  basisValues = std::vector<FloatScalarVector>(numPoints, FloatScalarVector(0.0));

  if (numPoints == 0) {
    return;
  }

  // returns linear interpolation between function values in evaluationPoint
  if (basistype == BasisCoefficientsComputationType::FUNCTION_VALUES) {
    if (evaluationPoint <= xValues[0]) {
      if (std::abs(xValues[0]) > 1e-14) {
        basisValues[0] = FloatScalarVector(evaluationPoint / xValues[0]);
      }
      return;
    }

    // TODO(holzmudd): could be optimized by using binary search
    for (size_t i = 1; i < numPoints; ++i) {
      if (evaluationPoint <= xValues[i]) {
        double x0 = xValues[i - 1];
        double x1 = xValues[i];

        double weight = (x1 - evaluationPoint) / (x1 - x0);
        basisValues[i - 1] = FloatScalarVector(weight);
        basisValues[i] = FloatScalarVector(1 - weight);

        return;
      }
    }

    // if we did not return in the loop, then evaluationPoint > all xValues...
    if (std::abs(xValues[numPoints - 1] - 1.0) > 1e-14) {
      basisValues[numPoints - 1] =
          FloatScalarVector((evaluationPoint - 1.0) / (xValues[numPoints - 1] - 1.0));
    }

    // returns evaluation of corresponding B-Spline in evaluationPoint
  } else if (basistype == BasisCoefficientsComputationType::SLE) {
    //    std::cout << "current level: " << level << "  evaluation Point: " << evaluationPoint
    //              << " grid: [";

    //    for (size_t i = 0; i < numPoints; ++i) {
    //      std::cout << xValues[i] << " ";
    //    }
    //    std::cout << "]" << std::endl;

    // ToDo (rehmemk) Bspline degree should be chooseable somewhere outside when setting
    // BasisCoeffCompType to SLE
    size_t degree = 3;

    std::unique_ptr<base::BsplineBasis<int, int>> bsplineBasis;
    //    bsplineBasis = std::unique_ptr<base::SBsplineBase>(new sgpp::base::SBsplineBase(degree));

    //    std::cout << "basis values: ";
    // hInv = 2^level
    const double hInv = static_cast<double>(1 << (level + 1));
    for (int i = 0; i < (int)numPoints; i++) {
      // cardinal B-Spline:  b(x/h - k)  with x the evaluation point, h the gridwidth and k the
      // index of the B-Spline.
      // the index of the B-Spline that has its midpoint in 0 is -2, the one that has it in h has
      // index -1, 2h has 0, 3h has 1 and so on...

      // Choosing the index like this sets the first Bspline the point h
      double index = i - 1;
      basisValues[i] = bsplineBasis->uniformBSpline(evaluationPoint * hInv - index, degree);
      //      std::cout << bsplineBasis->uniformBSpline(evaluationPoint * hInv - index, degree) <<
      //      ", ";
    }
    //    std::cout << "\n\n";
  }
}

void LinearInterpolationEvaluator::setGridPoints(const std::vector<double>& newXValues) {
  xValues = newXValues;
  computeBasisValues();
}

void LinearInterpolationEvaluator::setFunctionValuesAtGridPoints(
    std::vector<double>& functionValues) {
  basisCoefficients = functionValues;
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>
LinearInterpolationEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>(
      new LinearInterpolationEvaluator(*this));
}

bool LinearInterpolationEvaluator::needsOrderedPoints() { return true; }

bool LinearInterpolationEvaluator::needsParameter() { return true; }

void LinearInterpolationEvaluator::setParameter(const FloatScalarVector& param) {
  evaluationPoint = param.value();
  computeBasisValues();
}

} /* namespace combigrid */
} /* namespace sgpp*/
