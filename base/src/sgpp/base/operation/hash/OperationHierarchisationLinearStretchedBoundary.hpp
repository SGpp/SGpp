// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONLINEARSTRETCHEDBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONLINEARSTRETCHEDBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Hierarchisation on sparse grid, linear stretched case with boundaries
     *
     * @version $HEAD$
     */
    class OperationHierarchisationLinearStretchedBoundary : public OperationHierarchisation {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         */
        OperationHierarchisationLinearStretchedBoundary(GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationHierarchisationLinearStretchedBoundary() {}

        virtual void doHierarchisation(DataVector& node_values);
        virtual void doDehierarchisation(DataVector& alpha);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONHIERARCHISATIONLINEARSTRETCHEDBOUNDARY_HPP */