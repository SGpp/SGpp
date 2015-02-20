// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPERATION_HASH_OPERATIONMULTIPLEHIERARCHISATION_HPP
#define SGPP_OPTIMIZATION_OPERATION_HASH_OPERATIONMULTIPLEHIERARCHISATION_HPP

#include <vector>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Abstract operation for hierarchisation and dehierarchisation for
     * multiple sets of function values at the grid nodes.
     */
    class OperationMultipleHierarchisation :
      public base::OperationHierarchisation {
      public:
        /**
         * Constructor.
         */
        OperationMultipleHierarchisation() {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationMultipleHierarchisation() {
        }

        /**
         * Virtual method for hierarchising for one set of function values.
         * Defaults to calling the other doHierarchisation() by doing two
         * copy operations.
         *
         * @param[in,out] nodeValues before: vector of function values at
         *                           the grid points,
         *                           after: vector of hierarchical coefficients
         */
        virtual void doHierarchisation(base::DataVector& nodeValues) {
          std::vector<base::DataVector*> nodeValuesVec;
          nodeValuesVec.push_back(&nodeValues);
          doHierarchisation(nodeValuesVec);
        }

        /**
         * Virtual method for dehierarchising for one set of function values.
         * Defaults to calling the other doDehierarchisation() by doing two
         * copy operations.
         *
         * @param[in,out] alpha before: vector of hierarchical coefficients,
         *                      after: vector of function values at
         *                      the grid points
         */
        virtual void doDehierarchisation(base::DataVector& alpha) {
          std::vector<base::DataVector*> alphaVec;
          alphaVec.push_back(&alpha);
          doDehierarchisation(alphaVec);
        }

        /**
         * Pure virtual method for hierarchising for multiple sets of
         * function values.
         *
         * @param[in,out] nodeValues before: vector of function values at
         *                           the grid points,
         *                           after: vector of hierarchical coefficients
         */
        virtual void doHierarchisation(
          std::vector<base::DataVector*> nodeValues) = 0;

        /**
         * Pure virtual method for dehierarchising for multiple sets of
         * coefficients.
         *
         * @param[in,out] alpha before: vector of hierarchical coefficients,
         *                      after: vector of function values at
         *                      the grid points
         */
        virtual void doDehierarchisation(
          std::vector<base::DataVector*> alpha) = 0;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_OPERATION_HASH_OPERATIONMULTIPLEHIERARCHISATION_HPP */
