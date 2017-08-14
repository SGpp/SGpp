// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/onedim/ZetaHermiteInterpolationEvaluator.hpp>
#include <vector>

namespace sgpp {
namespace combigrid {

ZetaHermiteInterpolationEvaluator::ZetaHermiteInterpolationEvaluator()
    : evaluationPoint(0.0), basisCoefficients(), xValues() {}

ZetaHermiteInterpolationEvaluator::~ZetaHermiteInterpolationEvaluator() {}

ZetaHermiteInterpolationEvaluator::ZetaHermiteInterpolationEvaluator(
    const ZetaHermiteInterpolationEvaluator& other)
    : evaluationPoint(other.evaluationPoint),
      basisCoefficients(other.basisCoefficients),
      xValues(other.xValues) {}

void ZetaHermiteInterpolationEvaluator::computeBasisCoefficients() {
  size_t numPoints = xValues.size();

  basisCoefficients = std::vector<FloatScalarVector>(numPoints, FloatScalarVector(0.0));

  if (numPoints == 0) {
    return;
  }


  // without this evaluation always zero on this point
  /*if(evaluationPoint== 0 && xValues[0]==0){
    basisCoefficients[0] = FloatScalarVector(evalUniform(evaluationPoint  ));

    return;
  }
  */

  if (evaluationPoint <= xValues[0]) {
    if (std::abs(xValues[0]) > 1e-14) {
      basisCoefficients[0] = FloatScalarVector(evalUniform((evaluationPoint/ xValues[0]) -1 )*xValues[0]);
    }
    return;
  }

  // TODO(holzmudd): could be optimized by using binary search
  for (size_t i = 1; i < numPoints; ++i) {
    if (evaluationPoint <= xValues[i]) {
      double x0 = xValues[i - 1];
      double x1 = xValues[i];

      // translate to uniformpoint
      double weight = (evaluationPoint - x0) / (x1 - x0);
      basisCoefficients[i - 1] = FloatScalarVector(evalUniform(weight)*(x1 - x0));
      basisCoefficients[i] = FloatScalarVector(evalUniform(weight - 1)*(x1 - x0));

      return;
    }
  }

  // if we did not return in the loop, then evaluationPoint > all xValues...
  if (std::abs(xValues[numPoints - 1] - 1.0) > 1e-14) {
    double weight =(evaluationPoint - xValues[numPoints - 1]) / (1.0- xValues[numPoints - 1]);
    basisCoefficients[numPoints - 1] =
        FloatScalarVector(evalUniform(weight)*(1.0- xValues[numPoints - 1]));
  }

}

void ZetaHermiteInterpolationEvaluator::setGridPoints(const std::vector<double>& newXValues) {
  xValues = newXValues;
  computeBasisCoefficients();
}

std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >
ZetaHermiteInterpolationEvaluator::cloneLinear() {
  return std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector> >(
      new ZetaHermiteInterpolationEvaluator(*this));
}

bool ZetaHermiteInterpolationEvaluator::needsOrderedPoints() { return true; }

bool ZetaHermiteInterpolationEvaluator::needsParameter() { return true; }

void ZetaHermiteInterpolationEvaluator::setParameter(const FloatScalarVector& param) {
  evaluationPoint = param.value();
  computeBasisCoefficients();
}

double ZetaHermiteInterpolationEvaluator::evalUniform(double x) {

  if (x > 1 || x < -1)
      return 0;

    else if (x  < 0) {
      return pow(x, 3) + 2 * pow(x, 2) + x;
      
      
    } else {
      return pow(x, 3) - 2 * pow(x, 2) + x;
      
    }
}

} /* namespace combigrid */
} /* namespace sgpp*/
