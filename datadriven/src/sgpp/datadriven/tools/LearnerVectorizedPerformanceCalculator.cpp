// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>

#include <sgpp/globaldef.hpp>

#include <cstring>

namespace sgpp {
namespace datadriven {

LearnerVectorizedPerformance LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(
    base::Grid& grid, size_t numInstances, solver::SLESolverType solver, size_t numIterations,
    size_t sizeDatatype, bool reuseAlpha, bool verbose) {
  LearnerVectorizedPerformance result;

  result.GByte_ = 0.0;
  result.GFlop_ = 0.0;

  size_t nDim = grid.getDimension();
  size_t nGridsize = grid.getSize();

  // iterations + right hand side + initial residual
  double actualIterations = static_cast<double>(numIterations) + 0.5 + 1.0;
  if (reuseAlpha) {
    // if the last alpha is reused another mult() operation is required
    actualIterations += 1.0;
  }

  // each additional recalculation of the residual after 50 iterations has the cost of another
  // iteration
  size_t additionalResidualRecalculations = numIterations / 50;
  actualIterations += static_cast<double>(additionalResidualRecalculations);

  if (verbose) {
    if (reuseAlpha) {
      std::cout << "note: performance calculation: coefficients are not reused" << std::endl;
    } else {
      std::cout << "note: performance calculation: coefficients are reused" << std::endl;
    }
    std::cout << "note: performance calculation: iterations for calculation: " << actualIterations
              << std::endl;
    std::cout << "note: performance calculation: additional iteration for residual recalculation: "
              << additionalResidualRecalculations << std::endl;
  }

  if (grid.getType() == base::GridType::ModLinear) {
    for (size_t g = 0; g < grid.getSize(); g++) {
      base::GridPoint& curPoint = grid.getStorage().getPoint(g);

      for (size_t h = 0; h < nDim; h++) {
        base::level_t level;
        base::index_t index;

        curPoint.get(h, level, index);

        if (level == 1) {
        } else if (index == 1) {
          result.GFlop_ += 1e-9 * 8.0 * actualIterations * static_cast<double>(numInstances);
          result.GByte_ += 1e-9 * 4.0 * actualIterations * static_cast<double>(sizeDatatype) *
                           static_cast<double>(numInstances);
        } else if (index ==
                   static_cast<base::index_t>(
                       (1 << static_cast<base::index_t>(level)) - 1)) {
          result.GFlop_ += 1e-9 * 10.0 * actualIterations * static_cast<double>(numInstances);
          result.GByte_ += 1e-9 * 6.0 * actualIterations * static_cast<double>(sizeDatatype) *
                           static_cast<double>(numInstances);
        } else {
          result.GFlop_ += 1e-9 * 12.0 * actualIterations * static_cast<double>(numInstances);
          result.GByte_ += 1e-9 * 6.0 * actualIterations * static_cast<double>(sizeDatatype) *
                           static_cast<double>(numInstances);
        }
      }
    }

    // GBytes for EvalTrans (coefficients)
    result.GByte_ += 1e-9 * actualIterations *
                     ((static_cast<double>(nGridsize) * static_cast<double>(numInstances + 1) *
                       static_cast<double>(sizeDatatype)));

    // GBytes for Eval (coefficients)
    result.GByte_ += 1e-9 * actualIterations *
                     ((static_cast<double>(nGridsize + 1) * static_cast<double>(numInstances) *
                       static_cast<double>(sizeDatatype)));
  } else {
    // GFlops
    result.GFlop_ += 2.0 * 1e-9 * static_cast<double>(nGridsize) *
                     static_cast<double>(numInstances) * static_cast<double>(nDim) * 6.0 *
                     actualIterations;

    // GBytes
    result.GByte_ += 2.0 * 1e-9 * static_cast<double>(nGridsize) *
                     static_cast<double>(numInstances) * static_cast<double>(nDim) * 3.0 *
                     actualIterations * static_cast<double>(sizeDatatype);

    // GBytes for EvalTrans (coefficients)
    result.GByte_ += 1e-9 * actualIterations *
                     ((static_cast<double>(nGridsize) * static_cast<double>(numInstances) *
                       static_cast<double>(sizeDatatype)));

    // GBytes for Eval (coefficients)
    result.GByte_ += 1e-9 * actualIterations *
                     ((static_cast<double>(nGridsize) * static_cast<double>(numInstances) *
                       static_cast<double>(sizeDatatype)));
  }

  if (solver == solver::SLESolverType::BiCGSTAB) {
    result.GFlop_ = result.GFlop_ * 2.0;
    result.GByte_ = result.GByte_ * 2.0;
  }

  return result;
}

}  // namespace datadriven
}  // namespace sgpp
