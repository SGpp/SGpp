// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYREJECTIONSAMPLING_HPP_
#define OPERATIONDENSITYREJECTIONSAMPLING_HPP_

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * Sampling on all dimensions
     */

    class OperationDensityRejectionSampling {
      public:
        OperationDensityRejectionSampling() {}
        virtual ~OperationDensityRejectionSampling() {}

        /**
         * Rejection sampling
         *
         * @param alpha Coefficient vector for current grid
         * @param samples Output DataMatrix (rows: # of samples, columns: # of dims)
         * @param num_samples # of samples to draw
         * @param trial_max maximum # of trials for drawing a sample (exceeding will cause operation to stop)
         */
        virtual void doSampling(base::DataVector* alpha, base::DataMatrix*& samples, size_t num_samples, size_t trial_max) = 0;
    };

  }
}
#endif /* OPERATIONDENSITYREJECTIONSAMPLING_HPP_ */