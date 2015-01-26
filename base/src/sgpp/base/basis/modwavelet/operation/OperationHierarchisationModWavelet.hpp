/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONHIERARCHISATIONMODWAVELET_HPP
#define OPERATIONHIERARCHISATIONMODWAVELET_HPP

#include <sgpp/base/operation/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

namespace sg {
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
         * @todo (heinecke, nice) Implement the hierarchisation on the sparse grid with mod wavelets base functions
         */
        virtual void doHierarchisation(DataVector& node_values);

        /**
         * Implements the dehierarchisation on a sprase grid with mod wavelets base functions
         *
         * @param alpha the coefficients of the sparse grid's base functions
         *
         * @todo (heinecke, nice) Implement the dehierarchisation on the sparse grid with mod wavelets base functions
         */
        virtual void doDehierarchisation(DataVector& alpha);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONHIERARCHISATIONMODWAVELET_HPP */
