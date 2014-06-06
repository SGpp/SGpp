/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (dirk.pflueger@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#ifndef PHIPHIUPBBLINEARSTRETCHED_HPP
#define PHIPHIUPBBLINEARSTRETCHED_HPP

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace pde {



    /**
     * Implementation of sweep operator (): 1D Up for
     * Bilinearform \f$\int_{x} \phi(x) \phi(x) dx\f$
     */
    class PhiPhiUpBBLinearStretched {
      protected:
        typedef sg::base::GridStorage::grid_iterator grid_iterator;

        /// Pointer to sg::base::GridStorage object
        sg::base::GridStorage* storage;
        /// Pointer to the sg::base::Stretching Object
        sg::base::Stretching* stretching;

      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        PhiPhiUpBBLinearStretched(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~PhiPhiUpBBLinearStretched();

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
         * recursive function for the calculation of Up without bounding Box support
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

#endif /* PHIPHIUPBBLINEARSTRETCHED_HPP */
