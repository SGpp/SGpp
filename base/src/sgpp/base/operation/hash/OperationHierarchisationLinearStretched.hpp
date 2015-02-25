// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONLINEARSTRETCHED_HPP
#define OPERATIONHIERARCHISATIONLINEARSTRETCHED_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Hierarchisation on sparse grid, linear grid without boundaries
     */
    class OperationHierarchisationLinearStretched : public OperationHierarchisation {
      public:
        /**
         * Constructor of OperationHierarchisationLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationHierarchisationLinearStretched(GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationHierarchisationLinearStretched() {}

        virtual void doHierarchisation(DataVector& node_values);
        virtual void doDehierarchisation(DataVector& alpha);

      protected:
        /// Pointer to the grid's GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONHIERARCHISATIONLINEARSTRETCHED_HPP */