// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPERATION_HASH_OPERATIONMULTIPLEHIERARCHISATIONWAVELET_HPP
#define SGPP_OPTIMIZATION_OPERATION_HASH_OPERATIONMULTIPLEHIERARCHISATIONWAVELET_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisation.hpp>
#include <sgpp/base/grid/type/WaveletGrid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Hierarchisation operation for wavelet basis functions on Noboundary grids.
     */
    class OperationMultipleHierarchisationWavelet :
      public OperationMultipleHierarchisation {
      public:
        /**
         * Constructor.
         *
         * @param grid      grid
         */
        OperationMultipleHierarchisationWavelet(base::WaveletGrid& grid) : grid(grid) {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationMultipleHierarchisationWavelet() {
        }

        /**
         * @param[in,out] nodeValues    before: vector of function values at the grid points,
         *                              after: vector of hierarchical coefficients
         */
        virtual void doHierarchisation(std::vector<base::DataVector*> nodeValues);

        /**
         * @param[in,out] alpha         before: vector of hierarchical coefficients,
         *                              after: vector of function values at the grid points
         */
        virtual void doDehierarchisation(std::vector<base::DataVector*> alpha);

      protected:
        /// storage of the sparse grid
        base::WaveletGrid& grid;
    };

  }
}

#endif
