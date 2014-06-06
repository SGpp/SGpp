/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONHIERARCHISATIONLINEAR_HPP
#define OPERATIONHIERARCHISATIONLINEAR_HPP

#include "base/operation/OperationHierarchisation.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg {
  namespace base {

    /**
     * Hierarchisation on sparse grid, linear grid without boundaries
     */
    class OperationHierarchisationLinear : public OperationHierarchisation {
      public:
        /**
         * Constructor of OperationHierarchisationLinear
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationHierarchisationLinear(GridStorage* storage) : storage(storage) {}

        /**
         * Destructor
         */
        virtual ~OperationHierarchisationLinear() {}

        virtual void doHierarchisation(DataVector& node_values);
        virtual void doDehierarchisation(DataVector& alpha);

      protected:
        /// Pointer to the grid's GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONHIERARCHISATION_HPP */
