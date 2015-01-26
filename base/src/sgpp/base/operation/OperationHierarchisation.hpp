// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATION_HPP
#define OPERATIONHIERARCHISATION_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class implements the hierarchisation and dehierarchisation on the sparse grid
     */
    class OperationHierarchisation {
      public:
        /**
         * Constructor
         */
        OperationHierarchisation() {}

        /**
         * Destructor
         */
        virtual ~OperationHierarchisation() {}

        /**
         * Implements the hierarchisation on a sparse grid
         *
         * @param node_values the function's values in the nodal basis
         */
        virtual void doHierarchisation(DataVector& node_values) = 0;

        /**
         * Implements the dehierarchisation on a sparse grid
         *
         * @param alpha the coefficients of the sparse grid's basis functions
         */
        virtual void doDehierarchisation(DataVector& alpha) = 0;
    };

  }
}

#endif /* OPERATIONHIERARCHISATION_HPP */