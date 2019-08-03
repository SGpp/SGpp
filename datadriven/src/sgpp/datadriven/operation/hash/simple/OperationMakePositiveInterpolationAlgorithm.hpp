// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>
#include "../../../../../../../base/src/sgpp/base/function/scalar/ScalarFunction.hpp"

namespace sgpp {
namespace datadriven {

class OperationMakePositiveInterpolationAlgorithm {
 public:
  OperationMakePositiveInterpolationAlgorithm();
  virtual ~OperationMakePositiveInterpolationAlgorithm();

  virtual void initialize(base::Grid& grid);
  virtual void computeHierarchicalCoefficients(base::Grid& grid, base::DataVector& alpha,
                                               std::vector<size_t>& addedGridPoints,
                                               double tol = -1e-14) = 0;

 protected:
  double offset;
};

// -------------------------------------------------------------------------------------------
class OperationMakePositiveSetToZero : public OperationMakePositiveInterpolationAlgorithm {
 public:
  OperationMakePositiveSetToZero();
  virtual ~OperationMakePositiveSetToZero();

  void computeHierarchicalCoefficients(base::Grid& grid, base::DataVector& alpha,
                                       std::vector<size_t>& addedGridPoints,
                                       double tol = -1e-14) override;
};

// -------------------------------------------------------------------------------------------
class OperationMakePositiveInterpolateExp : public OperationMakePositiveInterpolationAlgorithm {
 public:
  void computeHierarchicalCoefficients(base::Grid& grid, base::DataVector& alpha,
                                       std::vector<size_t>& addedGridPoints,
                                       double tol = -1e-14) override;
};
// -------------------------------------------------------------------------------------------
class OperationMakePositiveInterpolateBoundaryOfSupport
    : public OperationMakePositiveInterpolationAlgorithm {
 public:
  void computeHierarchicalCoefficients(base::Grid& grid, base::DataVector& alpha,
                                       std::vector<size_t>& addedGridPoints,
                                       double tol = -1e-14) override;

 private:
  double computeMinimum(base::Grid& grid, base::DataVector& alpha, base::HashGridPoint& gp);
};
// -------------------------------------------------------------------------------------------
class OperationMakePositiveInterpolateFunction
    : public OperationMakePositiveInterpolationAlgorithm {
 public:
  explicit OperationMakePositiveInterpolateFunction(base::ScalarFunction* f);
  virtual ~OperationMakePositiveInterpolateFunction();

  void computeHierarchicalCoefficients(base::Grid& grid, base::DataVector& alpha,
                                       std::vector<size_t>& addedGridPoints,
                                       double tol = -1e-14) override;

 private:
  base::ScalarFunction* f;
};

} /* namespace datadriven */
} /* namespace sgpp */
