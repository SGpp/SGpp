// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYSAMPLING1DLINEAR_HPP_
#define OPERATIONDENSITYSAMPLING1DLINEAR_HPP_

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/OperationDensitySampling1D.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {
    class OperationDensitySampling1DLinear : public OperationDensitySampling1D {
      protected:
        base::Grid* grid;
      public:
        OperationDensitySampling1DLinear(base::Grid* grid);
        virtual ~OperationDensitySampling1DLinear();

        /**
         * Sampling on 1D grid
         *
         * @param alpha Coefficient vector for current grid (1D grid)
         * @param num_samples # of samples to draw
         * @param samples Output DataVector
        * @param seedp seed
         */
        void doSampling1D(base::DataVector* alpha, size_t num_samples, base::DataVector*& samples, unsigned int* seedp);
    };

  }
}

#endif /* OPERATIONDENSITYSAMPLING1D_HPP_ */