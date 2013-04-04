/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef PHIPHIUPMODLINEAR_HPP
#define PHIPHIUPMODLINEAR_HPP

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace pde {



    /**
     * Implementation of sweep operator (): 1D Up for
     * Bilinearform \f$\int_{x} \phi(x) \phi(x) dx\f$
     * on mod-linear grids
     */
    class PhiPhiUpModLinear {
      protected:
        typedef sg::base::GridStorage::grid_iterator grid_iterator;
        /// Pointer to the sg::base::GridStorage object
        sg::base::GridStorage* storage;

      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        PhiPhiUpModLinear(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        ~PhiPhiUpModLinear();

        /**
         * This operations performs the calculation of up in the direction of dimension <i>dim</i>
         *
         * @param source sg::base::DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
         * @param result sg::base::DataVector that contains the result of the up operation
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         */
        void operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim);

      protected:

        /**
         * recursive function for the calculation of Up
         *
         * @param source sg::base::DataVector that contains the coefficients of the ansatzfunction
         * @param result sg::base::DataVector in which the result of the operation is stored
         * @param index reference to a griditerator object that is used navigate through the grid
         * @param dim the dimension in which the operation is executed
         * @param fl function value on the left boundary
         * @param fr function value on the right boundary
         */
        void rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr);
    };



  }
}

#endif /* PHIPHIUPMODLINEAR_HPP */
