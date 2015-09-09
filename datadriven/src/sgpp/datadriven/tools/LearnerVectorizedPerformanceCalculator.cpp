// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>
#include <cstring>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    LearnerVectorizedPerformance LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(SGPP::base::Grid& grid,
        size_t numInstances, SGPP::solver::SLESolverType solver, size_t numIterations, size_t sizeDatatype) {
      LearnerVectorizedPerformance result;

      result.GByte_ = 0.0;
      result.GFlop_ = 0.0;

      size_t nDim = grid.getStorage()->dim();
      size_t nGridsize = grid.getSize();

      if (grid.getType() == base::GridType::ModLinear) {
        for (size_t g = 0; g < grid.getSize(); g++) {
          SGPP::base::GridIndex* curPoint = grid.getStorage()->get(g);

          for (size_t h = 0; h < nDim; h++) {
            SGPP::base::GridStorage::index_type::level_type level;
            SGPP::base::GridStorage::index_type::index_type index;

            curPoint->get(h, level, index);

            if (level == 1) {
            } else if (index == 1) {
              result.GFlop_ += 1e-9 * 8.0 * static_cast<float_t>(numIterations)
                               * static_cast<float_t>(numInstances);
              result.GByte_ += 1e-9 * 4.0 * static_cast<float_t>(numIterations) * static_cast<float_t>(sizeDatatype)
                               * static_cast<float_t>(numInstances);
            } else if (index
                       == static_cast<SGPP::base::GridStorage::index_type::index_type>((1
                           << static_cast<SGPP::base::GridStorage::index_type::index_type>(level)) - 1)) {
              result.GFlop_ += 1e-9 * 10.0 * static_cast<float_t>(numIterations)
                               * static_cast<float_t>(numInstances);
              result.GByte_ += 1e-9 * 6.0 * static_cast<float_t>(numIterations) * static_cast<float_t>(sizeDatatype)
                               * static_cast<float_t>(numInstances);
            } else {
              result.GFlop_ += 1e-9 * 12.0 * static_cast<float_t>(numIterations)
                               * static_cast<float_t>(numInstances);
              result.GByte_ += 1e-9 * 6.0 * static_cast<float_t>(numIterations) * static_cast<float_t>(sizeDatatype)
                               * static_cast<float_t>(numInstances);
            }
          }
        }

        // GBytes for EvalTrans (coefficients)
        result.GByte_ += 1e-9 * static_cast<float_t>(numIterations)
                         * ((static_cast<float_t>(nGridsize) * static_cast<float_t>(numInstances + 1)
                             * static_cast<float_t>(sizeDatatype)));

        // GBytes for Eval (coefficients)
        result.GByte_ += 1e-9 * static_cast<float_t>(numIterations)
                         * ((static_cast<float_t>(nGridsize + 1) * static_cast<float_t>(numInstances)
                             * static_cast<float_t>(sizeDatatype)));
      } else {
        // GFlops
        result.GFlop_ += 2.0 * 1e-9 * static_cast<float_t>(nGridsize) * static_cast<float_t>(numInstances)
                         * static_cast<float_t>(nDim) * 6.0 * static_cast<float_t>(numIterations);

        // GBytes
        result.GByte_ += 2.0 * 1e-9 * static_cast<float_t>(nGridsize) * static_cast<float_t>(numInstances)
                         * static_cast<float_t>(nDim) * 3.0 * static_cast<float_t>(numIterations)
                         * static_cast<float_t>(sizeDatatype);

        // GBytes for EvalTrans (coefficients)
        result.GByte_ += 1e-9 * static_cast<float_t>(numIterations)
                         * ((static_cast<float_t>(nGridsize) * static_cast<float_t>(numInstances)
                             * static_cast<float_t>(sizeDatatype)));

        // GBytes for Eval (coefficients)
        result.GByte_ += 1e-9 * static_cast<float_t>(numIterations)
                         * ((static_cast<float_t>(nGridsize) * static_cast<float_t>(numInstances)
                             * static_cast<float_t>(sizeDatatype)));
      }

      if (solver == SGPP::solver::BiCGSTAB) {
        result.GFlop_ = result.GFlop_ * 2.0;
        result.GByte_ = result.GByte_ * 2.0;
      }

      return result;
    }

  }
}
