/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/grid/GridStorage.hpp"
#include "parallel/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp"
#include <cstring>

namespace sg {

  namespace parallel {

    LearnerVectorizedPerformance LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(sg::base::Grid& Grid,
        size_t numInstances, sg::solver::SLESolverType solver, size_t numIterations, size_t sizeDatatype) {
      LearnerVectorizedPerformance result;

      result.GByte_ = 0.0;
      result.GFlop_ = 0.0;

      size_t nDim = Grid.getStorage()->dim();
      size_t nGridsize = Grid.getSize();

      if (strcmp(Grid.getType(), "modlinear") == 0) {
        for (size_t g = 0; g < Grid.getSize(); g++) {
          sg::base::GridIndex* curPoint = Grid.getStorage()->get(g);

          for (size_t h = 0; h < nDim; h++) {
            sg::base::GridStorage::index_type::level_type level;
            sg::base::GridStorage::index_type::index_type index;

            curPoint->get(h, level, index);

            if (level == 1) {
            } else if (index == 1) {
              result.GFlop_ += 1e-9 * 8.0 * static_cast<double>(numIterations) * static_cast<double>(numInstances);
              result.GByte_ += 1e-9 * 4.0 * static_cast<double>(numIterations) * static_cast<double>(sizeDatatype) * static_cast<double>(numInstances);
            } else if (index == static_cast<sg::base::GridStorage::index_type::index_type>((1 << static_cast<sg::base::GridStorage::index_type::index_type>(level)) - 1)) {
              result.GFlop_ += 1e-9 * 10.0 * static_cast<double>(numIterations) * static_cast<double>(numInstances);
              result.GByte_ += 1e-9 * 6.0 * static_cast<double>(numIterations) * static_cast<double>(sizeDatatype) * static_cast<double>(numInstances);
            } else {
              result.GFlop_ += 1e-9 * 12.0 * static_cast<double>(numIterations) * static_cast<double>(numInstances);
              result.GByte_ += 1e-9 * 6.0 * static_cast<double>(numIterations) * static_cast<double>(sizeDatatype) * static_cast<double>(numInstances);
            }
          }
        }

        // GBytes for EvalTrans (coefficients)
        result.GByte_ += 1e-9 * static_cast<double>(numIterations)
                         * ((static_cast<double>(nGridsize) * static_cast<double>(numInstances + 1) * static_cast<double>(sizeDatatype)));

        // GBytes for Eval (coefficients)
        result.GByte_ += 1e-9 * static_cast<double>(numIterations)
                         * ((static_cast<double>(nGridsize + 1) * static_cast<double>(numInstances) * static_cast<double>(sizeDatatype)));
      } else {
        // GFlops
        result.GFlop_ += 2.0 * 1e-9 * static_cast<double>(nGridsize) * static_cast<double>(numInstances) * static_cast<double>(nDim) * 6.0 * static_cast<double>(numIterations);

        // GBytes
        result.GByte_ += 2.0 * 1e-9 * static_cast<double>(nGridsize) * static_cast<double>(numInstances) * static_cast<double>(nDim) * 3.0 * static_cast<double>(numIterations) * static_cast<double>(sizeDatatype);

        // GBytes for EvalTrans (coefficients)
        result.GByte_ += 1e-9 * static_cast<double>(numIterations)
                         * ((static_cast<double>(nGridsize) * static_cast<double>(numInstances) * static_cast<double>(sizeDatatype)));

        // GBytes for Eval (coefficients)
        result.GByte_ += 1e-9 * static_cast<double>(numIterations)
                         * ((static_cast<double>(nGridsize) * static_cast<double>(numInstances) * static_cast<double>(sizeDatatype)));
      }

      if (solver == sg::solver::BiCGSTAB) {
        result.GFlop_ = result.GFlop_ * 2.0;
        result.GByte_ = result.GByte_ * 2.0;
      }

      return result;
    }

  }

}
