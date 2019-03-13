// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/parallel/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>

#include <sgpp/globaldef.hpp>
#include <cstring>

namespace sgpp {

namespace parallel {

LearnerVectorizedPerformance LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(
    sgpp::base::Grid& grid, size_t numInstances, sgpp::solver::SLESolverType solver,
    size_t numIterations, size_t sizeDatatype) {
  LearnerVectorizedPerformance result;

  result.GByte_ = 0.0;
  result.GFlop_ = 0.0;

  size_t nDim = grid.getStorage().getDimension();
  size_t nGridsize = grid.getSize();

  if (grid.getType() == base::GridType::ModLinear) {
    for (size_t g = 0; g < grid.getSize(); g++) {
      sgpp::base::GridPoint& curPoint = grid.getStorage().getPoint(g);

      for (size_t h = 0; h < nDim; h++) {
        sgpp::base::level_t level;
        sgpp::base::index_t index;

        curPoint.get(h, level, index);

        if (level == 1) {
        } else if (index == 1) {
          result.GFlop_ +=
              1e-9 * 8.0 * static_cast<double>(numIterations) * static_cast<double>(numInstances);
          result.GByte_ += 1e-9 * 4.0 * static_cast<double>(numIterations) *
                           static_cast<double>(sizeDatatype) * static_cast<double>(numInstances);
        } else if (index ==
                   static_cast<sgpp::base::index_t>(
                       (1 << static_cast<sgpp::base::index_t>(level)) -
                       1)) {
          result.GFlop_ +=
              1e-9 * 10.0 * static_cast<double>(numIterations) * static_cast<double>(numInstances);
          result.GByte_ += 1e-9 * 6.0 * static_cast<double>(numIterations) *
                           static_cast<double>(sizeDatatype) * static_cast<double>(numInstances);
        } else {
          result.GFlop_ +=
              1e-9 * 12.0 * static_cast<double>(numIterations) * static_cast<double>(numInstances);
          result.GByte_ += 1e-9 * 6.0 * static_cast<double>(numIterations) *
                           static_cast<double>(sizeDatatype) * static_cast<double>(numInstances);
        }
      }
    }

    // GBytes for EvalTrans (coefficients)
    result.GByte_ += 1e-9 * static_cast<double>(numIterations) *
                     ((static_cast<double>(nGridsize) * static_cast<double>(numInstances + 1) *
                       static_cast<double>(sizeDatatype)));

    // GBytes for Eval (coefficients)
    result.GByte_ += 1e-9 * static_cast<double>(numIterations) *
                     ((static_cast<double>(nGridsize + 1) * static_cast<double>(numInstances) *
                       static_cast<double>(sizeDatatype)));
  } else {
    // GFlops
    result.GFlop_ += 2.0 * 1e-9 * static_cast<double>(nGridsize) *
                     static_cast<double>(numInstances) * static_cast<double>(nDim) * 6.0 *
                     static_cast<double>(numIterations);

    // GBytes
    result.GByte_ += 2.0 * 1e-9 * static_cast<double>(nGridsize) *
                     static_cast<double>(numInstances) * static_cast<double>(nDim) * 3.0 *
                     static_cast<double>(numIterations) * static_cast<double>(sizeDatatype);

    // GBytes for EvalTrans (coefficients)
    result.GByte_ += 1e-9 * static_cast<double>(numIterations) *
                     ((static_cast<double>(nGridsize) * static_cast<double>(numInstances) *
                       static_cast<double>(sizeDatatype)));

    // GBytes for Eval (coefficients)
    result.GByte_ += 1e-9 * static_cast<double>(numIterations) *
                     ((static_cast<double>(nGridsize) * static_cast<double>(numInstances) *
                       static_cast<double>(sizeDatatype)));
  }

  if (solver == sgpp::solver::SLESolverType::BiCGSTAB) {
    result.GFlop_ = result.GFlop_ * 2.0;
    result.GByte_ = result.GByte_ * 2.0;
  }

  return result;
}
}  // namespace parallel
}  // namespace sgpp