// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONMODWAVELET_HPP
#define OPERATIONHIERARCHISATIONMODWAVELET_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Hierarchisation on sparse grid, mod wavelet case
     */
    class OperationHierarchisationModWavelet : public OperationHierarchisation {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         */
        OperationHierarchisationModWavelet(GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationHierarchisationModWavelet() {}

        /**
         * Implements the hierarchisation on a sprase grid with mod wavelets base functions
         *
         * @param node_values the functions values in the node base
         *
         */
        virtual void doHierarchisation(DataVector& node_values);

        /**
         * Implements the dehierarchisation on a sprase grid with mod wavelets base functions
         *
         * @param alpha the coefficients of the sparse grid's base functions
         *
         */
        virtual void doDehierarchisation(DataVector& alpha);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONHIERARCHISATIONMODWAVELET_HPP */
