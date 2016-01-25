// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DEHIERARCHISATIONPOLYBOUNDARY_HPP
#define DEHIERARCHISATIONPOLYBOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace base {

    /**
     * Class that implements the dehierarchisation on a polynomial sparse grid. Therefore
     * the ()operator has to be implement in order to use the sweep algorithm for
     * the grid traversal
     */
    class DehierarchisationPolyBoundary {
      protected:
        typedef GridStorage::grid_iterator grid_iterator;
        typedef GridStorage::index_type::level_type level_type;
        typedef GridStorage::index_type::index_type index_type;

        /// the grid object
        GridStorage* storage;

        /// the base
        SPolyBoundaryBase* base;

      public:
        /**
         * Constructor, must be bind to a grid
         *
         * @param storage the grid storage object of the the grid, on which the dehierarchisation should be executed
         * @param base Polynomial basis
         */
        DehierarchisationPolyBoundary(GridStorage* storage, SPolyBoundaryBase* base);
        /**
         * Destructor
         */
        virtual ~DehierarchisationPolyBoundary();

        /**
         * Implements operator() needed by the sweep class during the grid traversal. This function
         * is applied to the whole grid.
         *
         * @param source this DataVector holds the linear base coefficients of the sparse grid's ansatz-functions
         * @param result this DataVector holds the node base coefficients of the function that should be applied to the sparse grid
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         */
        virtual void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim);

      protected:

        /**
         * Recursive dehierarchisation algorithm, this algorithms works in-place -> source should be equal to result
         *
         * @param source this DataVector holds the linear base coefficients of the sparse grid's ansatz-functions
         * @param result this DataVector holds the node base coefficients of the function that should be applied to the sparse grid
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         * @param coeffs nodal coefficients computed so far
         */
        void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, DataVector& coeffs);
    };

    // namespace detail

  } // namespace sg
}

#endif /* DEHIERARCHISATIONPOLYBOUNDARY_HPP */
