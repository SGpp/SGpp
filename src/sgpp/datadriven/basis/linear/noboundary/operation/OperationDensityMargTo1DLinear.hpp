/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author A. Mo-Hellenbrand

#ifndef OPERATIONDENSITYMARGTO1DLINEAR_HPP
#define OPERATIONDENSITYMARGTO1DLINEAR_HPP

#include "base/grid/Grid.hpp"
#include "datadriven/operation/OperationDensityMargTo1D.hpp"
#include <cstring>

namespace sg {
  namespace datadriven {

    /**
     * keep applying marginalize to function until it's reduced to only 1 dimension
     */

    class OperationDensityMargTo1DLinear : public OperationDensityMargTo1D {
      public:
        OperationDensityMargTo1DLinear(base::Grid* grid) : grid(grid) {}
        virtual ~OperationDensityMargTo1DLinear() {}

        /**
         * Keep applying marginalizes to (Density) Functions, until it's reduced to 1 dimension (dim_x)
         *
         * @param alpha Coefficient vector for current grid
         * @param grid_x output 1D grid pointer
         * @param alpha_x Coefficient vector for new grid (grid_x). Will be initialized.
         * @param dim_x Target dimension, all other dimensions will be marginalized
         */
        void margToDimX(base::DataVector* alpha, base::Grid*& grid_x, base::DataVector*& alpha_x, size_t dim_x);

      protected:
        base::Grid* grid;
        void marg_next_dim(base::Grid* g_in, base::DataVector* a_in, base::Grid*& g_out, base::DataVector*& a_out, size_t dims, size_t dim_x, size_t& count);
    };

  }
}
#endif /* OPERATIONDENSITYMARGTO1DLINEAR_HPP */





