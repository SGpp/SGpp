/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#ifndef SQXDPHIDPHIUPBBLINEARSTRETCHED_HPP
#define SQXDPHIDPHIUPBBLINEARSTRETCHED_HPP

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace finance {



    /**
     * Implementation of sweep operator (): 1D Up for
     * Bilinearform \f$\int_{x} x^{2} \frac{\partial \phi(x)}{x} \frac{\partial \phi(x)}{x} dx\f$
     */
    class SqXdPhidPhiUpBBLinearStretched {
      protected:
        typedef sg::base::GridStorage::grid_iterator grid_iterator;

        /// Pointer to sg::base::GridStorage object
        sg::base::GridStorage* storage;
        /// Pointer to the stretching Obejct
        sg::base::Stretching* stretching;

      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        SqXdPhidPhiUpBBLinearStretched(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~SqXdPhidPhiUpBBLinearStretched();


        /**
         * This operations performs the calculation of up in the direction of dimension <i>dim</i>
         * on a grid with fix Dirichlet 0 boundary conditions
         *
         * @param source sg::base::DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
         * @param result sg::base::DataVector that contains the result of the up operation
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         */
        virtual void operator()(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim);

      protected:

        /**
         * recursive function for the calculation of Up without Bounding Box support
         *
         * @param source sg::base::DataVector that contains the coefficients of the ansatzfunction
         * @param result sg::base::DataVector in which the result of the operation is stored
         * @param index reference to a griditerator object that is used navigate through the grid
         * @param dim the dimension in which the operation is executed
         * @param fl function value on the left boundary, reference parameter
         * @param fr function value on the right boundary, reference parameter
         */
        void rec(sg::base::DataVector& source, sg::base::DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr);
    };

    // namespace detail

  } // namespace sg
}

#endif /* SQXDPHIDPHIUPBBLINEARSTRETCHEDBOUNDARY_HPP */
