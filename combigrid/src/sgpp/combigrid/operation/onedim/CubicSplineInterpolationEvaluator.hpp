// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class CubicSplineInterpolationEvaluator : public AbstractLinearEvaluator<FloatScalarVector> {
 public:
  CubicSplineInterpolationEvaluator();
  virtual ~CubicSplineInterpolationEvaluator();
  CubicSplineInterpolationEvaluator(CubicSplineInterpolationEvaluator const &other);

  virtual std::vector<FloatScalarVector> getBasisCoefficients() { return basisCoefficients; }

  virtual void setGridPoints(std::vector<double> const &newXValues);
  virtual std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear();
  virtual bool needsOrderedPoints();
  virtual bool needsParameter();
  virtual void setParameter(FloatScalarVector const &param);

 private:
  void computeBasisCoefficients();

  double evaluationPoint;
  std::vector<FloatScalarVector> basisCoefficients;
  std::vector<double> xValues;

  struct SplineSet {
    double a;
    double b;
    double c;
    double d;
  };

  std::vector<std::vector<SplineSet>> gridCoefficients;
};

} /* namespace combigrid */
} /* namespace sgpp */
