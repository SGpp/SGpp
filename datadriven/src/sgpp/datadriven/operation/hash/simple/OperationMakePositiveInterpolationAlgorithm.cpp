// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationMakePositiveInterpolationAlgorithm.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>
#include <algorithm>
#include <limits>

namespace sgpp {
namespace datadriven {

OperationMakePositiveInterpolationAlgorithm::OperationMakePositiveInterpolationAlgorithm()
    : offset(0.0) {}

OperationMakePositiveInterpolationAlgorithm::~OperationMakePositiveInterpolationAlgorithm() {}

void OperationMakePositiveInterpolationAlgorithm::initialize(base::Grid& grid) { offset = 0.0; }

// -------------------------------------------------------------------------------------------

OperationMakePositiveSetToZero::OperationMakePositiveSetToZero() {}
OperationMakePositiveSetToZero::~OperationMakePositiveSetToZero() {}

void OperationMakePositiveSetToZero::computeHierarchicalCoefficients(
    base::Grid& grid, base::DataVector& alpha, std::vector<size_t>& addedGridPoints, double tol) {
  // compute the nodal values of the newly added grid points and subtract the
  // nodal value from the hierarchical coefficient. This sets the function value at that point
  // to zero
  base::HashGridStorage& gridStorage = grid.getStorage();

  auto opEval = std::unique_ptr<base::OperationEval>(op_factory::createOperationEvalNaive(grid));
  base::DataVector x(gridStorage.getDimension());

  for (auto& i : addedGridPoints) {
    gridStorage.getPoint(i).getStandardCoordinates(x);
    double yi = opEval->eval(alpha, x);
    if (yi < tol) {
      alpha[i] = alpha[i] - yi + offset;
    } else {
      alpha[i] = 0.0;
    }
  }
}

// -------------------------------------------------------------------------------------------
void OperationMakePositiveInterpolateExp::computeHierarchicalCoefficients(
    base::Grid& grid, base::DataVector& alpha, std::vector<size_t>& addedGridPoints, double tol) {
  // compute the nodal values of the newly added grid points and subtract the
  // nodal value from the hierarchical coefficient. This sets the function value at that point
  // to zero. We additionaly add the value from the log interpolant.
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto opEval = std::unique_ptr<base::OperationEval>(op_factory::createOperationEvalNaive(grid));
  base::DataVector x(gridStorage.getDimension());

  for (auto& i : addedGridPoints) {
    gridStorage.getPoint(i).getStandardCoordinates(x);
    double yi = opEval->eval(alpha, x);
    if (yi < tol) {
      // eval exp of yi
      double expyi = 0.0;
      if (yi < -0.1) {
        expyi = std::exp(yi);
      }

      alpha[i] = alpha[i] - yi + expyi;
    } else {
      alpha[i] = 0.0;
    }
  }
}

// -------------------------------------------------------------------------------------------
double OperationMakePositiveInterpolateBoundaryOfSupport::computeMinimum(base::Grid& grid,
                                                                         base::DataVector& alpha,
                                                                         base::HashGridPoint& gp) {
  double minimumInterpolatedValue = std::numeric_limits<double>::max();

  base::HashGridStorage& gs = grid.getStorage();
  size_t numDims = gs.getDimension();

  std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEvalNaive(grid));
  base::DataVector x(numDims);

  double leftValue = 0.0;
  double rightValue = 0.0;

  for (size_t idim = 0; idim < numDims; idim++) {
    // get right boundary point in the current dimension
    auto level = gp.getLevel(idim);
    auto index = gp.getIndex(idim);

    // check grid point at the left boundary
    gp.getLeftBoundaryPoint(idim);
    if (gs.isContaining(gp)) {
      gp.getStandardCoordinates(x);
      leftValue = opEval->eval(alpha, x);
    } else {
      leftValue = 0.0;
    }

    // reset grid point
    gp.set(idim, level, index);

    // check grid point at the right boundary
    gp.getRightBoundaryPoint(idim);
    if (gs.isContaining(gp)) {
      gp.getStandardCoordinates(x);
      rightValue = opEval->eval(alpha, x);
    } else {
      rightValue = 0.0;
    }

    double interpolatedValue = std::abs(leftValue + rightValue) / 2.0;

    if (idim == 0 || (interpolatedValue > 0 && interpolatedValue < minimumInterpolatedValue)) {
      minimumInterpolatedValue = interpolatedValue;
    }

    // reset grid point
    gp.set(idim, level, index);
  }

  return minimumInterpolatedValue;
}

void OperationMakePositiveInterpolateBoundaryOfSupport::computeHierarchicalCoefficients(
    base::Grid& grid, base::DataVector& alpha, std::vector<size_t>& addedGridPoints, double tol) {
  // compute the nodal values of the newly added grid points and subtract the
  // nodal value from the hierarchical coefficient. This sets the function value at that point
  // to zero. We additionaly interpolate the boundary points linearly in each dimension and
  // add the minimum to the current function value.
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto opEval = std::unique_ptr<base::OperationEval>(op_factory::createOperationEvalNaive(grid));

  base::DataVector x(gridStorage.getDimension());

  for (auto& i : addedGridPoints) {
    auto gp = gridStorage.getPoint(i);
    gp.getStandardCoordinates(x);
    double yi = opEval->eval(alpha, x);
    if (yi < tol) {
      alpha[i] = alpha[i] - yi + computeMinimum(grid, alpha, gp);
    } else {
      alpha[i] = 0.0;
    }
  }
}
// -------------------------------------------------------------------------------------------

OperationMakePositiveInterpolateFunction::OperationMakePositiveInterpolateFunction(
    base::ScalarFunction* f)
    : f(f) {}
OperationMakePositiveInterpolateFunction::~OperationMakePositiveInterpolateFunction() {}

void OperationMakePositiveInterpolateFunction::computeHierarchicalCoefficients(
    base::Grid& grid, base::DataVector& alpha, std::vector<size_t>& addedGridPoints, double tol) {
  base::HashGridStorage& gridStorage = grid.getStorage();
  auto opEval = std::unique_ptr<base::OperationEval>(op_factory::createOperationEvalNaive(grid));
  base::DataVector x(gridStorage.getDimension());

  for (auto& i : addedGridPoints) {
    gridStorage.getPoint(i).getStandardCoordinates(x);
    double yi = opEval->eval(alpha, x);
    if (yi < tol) {
      alpha[i] = alpha[i] - yi + std::abs(f->eval(x));
    } else {
      alpha[i] = 0.0;
    }
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
