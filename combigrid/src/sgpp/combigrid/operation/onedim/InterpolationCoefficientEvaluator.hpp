// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class InterpolationCoefficientEvaluator : public AbstractLinearEvaluator<FloatTensorVector> {
  std::vector<FloatTensorVector> basisValues;
  std::vector<double> basisCoefficients;
  std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis;

  std::vector<double> xValues;

  void computeBasisValues();

 public:
  InterpolationCoefficientEvaluator();
  explicit InterpolationCoefficientEvaluator(
      std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis);
  ~InterpolationCoefficientEvaluator() override;
  InterpolationCoefficientEvaluator(InterpolationCoefficientEvaluator const &other);

  std::vector<FloatTensorVector> getBasisValues() override { return basisValues; }
  std::vector<double> getBasisCoefficients() override { return basisCoefficients; }

  void setGridPoints(std::vector<double> const &newXValues) override;
  void setBasisCoefficientsAtGridPoints(std::vector<double> &functionValues) override;
  std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>> cloneLinear() override;
  bool needsOrderedPoints() override;
  bool needsParameter() override;
  void setParameter(FloatTensorVector const &param) override;
};

} /* namespace combigrid */
} /* namespace sgpp */
