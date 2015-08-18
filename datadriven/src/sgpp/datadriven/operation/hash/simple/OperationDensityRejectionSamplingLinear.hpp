// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYREJECTIONSAMPLINGLINEAR_HPP
#define OPERATIONDENSITYREJECTIONSAMPLINGLINEAR_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityRejectionSampling.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * Sampling with rejection sampling method
     */

    class OperationDensityRejectionSamplingLinear : public OperationDensityRejectionSampling {
      public:
        OperationDensityRejectionSamplingLinear(base::Grid* grid) : grid(grid) {}
        virtual ~OperationDensityRejectionSamplingLinear() {}

        /**
         * Rejection sampling
         *
         * @param alpha Coefficient vector for current grid
         * @param samples Output DataMatrix (rows: # of samples, columns: # of dims)
         * @param num_samples # of samples to draw
         * @param trial_max maximum # of trials for drawing a sample (exceeding will cause operation to stop)
         */
        void doSampling(base::DataVector* alpha, base::DataMatrix*& samples, size_t num_samples, size_t trial_max);

      protected:
        base::Grid* grid;
    };

  }
}
#endif /* OPERATIONDENSITYREJECTIONSAMPLINGLINEAR_HPP */





