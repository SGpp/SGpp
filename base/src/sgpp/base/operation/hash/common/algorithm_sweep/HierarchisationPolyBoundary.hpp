// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HIERARCHISATIONPOLYBOUNDARY_HPP
#define HIERARCHISATIONPOLYBOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {

  namespace base {

    /**
     * Class that implements the hierarchisation on a polynomial sparse grid. Therefore
     * the ()operator has to be implement in order to use the sweep algorithm for
     * the grid traversal
     */
    class HierarchisationPolyBoundary {
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
         * @param storage the grid storage object of the the grid, on which the hierarchisation should be executed
         * @param base The polynomial basis functions
         */
        HierarchisationPolyBoundary(GridStorage* storage, SPolyBoundaryBase* base);

        /**
         * Destructor
         */
        virtual ~HierarchisationPolyBoundary();

        /**
         * Implements operator() needed by the sweep class during the grid traversal. This function
         * is applied to the whole grid.
         *
         * @param source this DataVector holds the node base coefficients of the function that should be applied to the sparse grid
         * @param result this DataVector holds the linear base coefficients of the sparse grid's ansatz-functions
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         */
        virtual void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim);

      protected:

        /**
         * Recursive hierarchisaton algorithm, this algorithms works in-place -> source should be equal to result
         *
         * @todo add graphical explanation here
         *
         * @param source this DataVector holds the node base coefficients of the function that should be applied to the sparse grid
         * @param result this DataVector holds the linear base coefficients of the sparse grid's ansatz-functions
         * @param index a iterator object of the grid
         * @param dim current fixed dimension of the 'execution direction'
         * @param coeffs coefficients of the basis functions as calculated so far
         */
        void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, DataVector& coeffs);

    };

  } // namespace base

} // namespace sg

#endif /* HIERARCHISATIONPOLY_HPP */
