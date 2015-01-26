/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author A. Mo-Hellenbrand

#ifndef OPERATIONDENSITYSAMPLING_HPP_
#define OPERATIONDENSITYSAMPLING_HPP_

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * Sampling on all dimensions
     */

    class OperationDensitySampling {
      public:
        OperationDensitySampling() {}
        virtual ~OperationDensitySampling() {}

        /**
         * Sampling with mixed starting dimensions
         *
         * @param alpha Coefficient vector for current grid
         * @param samples Output DataMatrix (rows: # of samples, columns: # of dims)
         * @param num_samples # of samples to draw
         */
        virtual void doSampling(base::DataVector* alpha, base::DataMatrix*& samples, size_t num_samples) = 0;

        /**
         * Sampling with specified starting dimension
         *
         * @param alpha Coefficient vector for current grid
         * @param samples Output DataMatrix (rows: # of samples, columns: # of dims)
         * @param num_samples # of samples to draw
         * @param dim_x Starting dimension
         */
        virtual void doSampling(base::DataVector* alpha, base::DataMatrix*& samples, size_t num_samples, size_t dim_x) = 0;
    };

  }
}
#endif /* OPERATIONDENSITYSAMPLING_HPP_ */
