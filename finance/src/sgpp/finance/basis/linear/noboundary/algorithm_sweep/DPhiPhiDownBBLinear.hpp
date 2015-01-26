/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#ifndef DPHIPHIDOWNBBLINEAR_HPP
#define DPHIPHIDOWNBBLINEAR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sg {
  namespace finance {



    /**
     * Implementation of sweep operator (): 1D Down for
     * Bilinearform \f$\int_{x} \frac{\partial \phi(x)}{x} \phi(x) dx\f$
     */
    class DPhiPhiDownBBLinear {
      protected:
        typedef sg::base::GridStorage::grid_iterator grid_iterator;

        /// Pointer to the sg::base::GridStorage Object
        sg::base::GridStorage* storage;
        /// Pointer to the bounding box Obejct
        sg::base::BoundingBox* boundingBox;


      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        DPhiPhiDownBBLinear(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~DPhiPhiDownBBLinear();

        /**
         * This operations performs the calculation of down in the direction of dimension <i>dim</i>
         *
         * @param source sg::base::DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
         * @param result sg::base::DataVector that contains the result of the down operation
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         */
        virtual void operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim);

      protected:

        /**
         * recursive function for the calculation of Down without and with Bounding Box support
         * (calculations are independent from bounding box)
         *
         * @param source sg::base::DataVector that contains the coefficients of the ansatzfunction
         * @param result sg::base::DataVector in which the result of the operation is stored
         * @param index reference to a griditerator object that is used navigate through the grid
         * @param dim the dimension in which the operation is executed
         * @param fl function value on the left boundary
         * @param fr function value on the right boundary
         */
        void rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double fl, double fr);
    };

    // namespace detail

  } // namespace sg
}

#endif /* PHIDPHIDOWNBBLINEAR_HPP */
