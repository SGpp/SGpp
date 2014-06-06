/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Benjamin Peherstorfer (pehersto@in.tum.de)

#ifndef OPERATIONDENSITYSAMPLING1D_HPP
#define OPERATIONDENSITYSAMPLING1D_HPP

#include "base/grid/Grid.hpp"
#include <cstring>

namespace sg {
  namespace datadriven {

    /**
     * Sample 1D Probability Density Function
     */

    class OperationDensitySampling1D {
      public:
        OperationDensitySampling1D() {}
        virtual ~OperationDensitySampling1D() {}

        /**
         * Sampling on 1D grid
         *
         * @param alpha Coefficient vector for current grid (1D grid)
         * @param num_samples # of samples to draw
         * @param samples Output DataVector
        * @param seedp seed
         */
        virtual void doSampling1D(base::DataVector* alpha, size_t num_samples, base::DataVector*& samples, unsigned int* seedp) = 0;
    };

  }
}
#endif /* OPERATIONDENSITYSAMPLING1D_HPP */
