// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONLINEARBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONLINEARBOUNDARY_HPP

#include <sgpp/base/operation/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Hierarchisation on sparse grid, linear case with boundaries
     *
     * @version $HEAD$
     */
    class OperationHierarchisationLinearBoundary : public OperationHierarchisation {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         */
        OperationHierarchisationLinearBoundary(GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationHierarchisationLinearBoundary() {}

        virtual void doHierarchisation(DataVector& node_values);
        virtual void doDehierarchisation(DataVector& alpha);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONHIERARCHISATIONLINEARBOUNDARY_HPP */