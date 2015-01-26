/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */

#ifndef UPDOWNONEOPDIMWITHSHADOW_HPP
#define UPDOWNONEOPDIMWITHSHADOW_HPP

#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/operation/OperationMatrix.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>


namespace sg {
  namespace pde {

    /**
     * Implements the Up/Down scheme with one dimension with a special operation. Before the actual
     * operation starts, all grid points from the shadow storage are copied into the actual grid.
     * After the calculation, all shadow grid points are deleted, the result vector is adapted
     * accordingly.
     *
     * @version $HEAD$
     */
    class UpDownOneOpDimWithShadow: public sg::base::OperationMatrix {
      public:

        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         * @param shadowStorage shadow storage
         */
        UpDownOneOpDimWithShadow(sg::base::GridStorage* storage, sg::base::GridStorage* shadowStorage);

        /**
         * Destructor
         */
        virtual ~UpDownOneOpDimWithShadow();


        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

      protected:
        typedef sg::base::GridStorage::grid_iterator grid_iterator;

        /// Pointer to the grid's storage object
        sg::base::GridStorage* storage;
        sg::base::GridStorage* shadowStorage;

        /**
         * Recursive procedure for updown().
         *
         * @param dim the current dimension
         * @param op_dim the dimension in which a special operation is applied
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        void updown(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim);

        /**
         * This functions adds all grid points of the shadow storage into the actual grid.
         */
        void expandGrid();

        /**
         * Removes the previously added shadow grid points from the actual grid. Thus, the
         * grid is the same shape as before the call of mult.
         */
        void shrinkGrid();

        /**
         * All calculations for op_dim.
         *
         * @param alpha the coefficients of the grid points
         * @param result the result of the operations
         * @param dim the current dimension in the recursion
         * @param op_dim the dimension in that a special operation is applied
         */
        virtual void specialOP(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim);


        /**
         * std 1D up operation
         *
         * @param dim dimension in which to apply the up-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * std 1D down operation
         *
         * @param dim dimension in which to apply the down-part
         * @param alpha vector of coefficients
         * @param result vector to store the results in
         */
        virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * special 1D down operation that is only executed in one direction
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that down-Gradient is applied
         */
        virtual void downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;

        /**
         * special 1D up operation that is only executed in one direction
         *
         * @param alpha the coefficients of the gridpoints
         * @param result vector with the result of this operation
         * @param dim the dimension in that up-Gradient is applied
         */
        virtual void upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
    };

  }
}

#endif /* UPDOWNONEOPDIMWITHSHADOW_HPP */
