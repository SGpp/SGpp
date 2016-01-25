// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONMODPOLY_HPP
#define OPERATIONHIERARCHISATIONMODPOLY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Hierarchisation on sparse grid, mod poly case
     *
     */
    class OperationHierarchisationModPoly : public OperationHierarchisation {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         * @param degree the polynom's max. degree
         */
        OperationHierarchisationModPoly(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

        /**
         * Destructor
         */
        virtual ~OperationHierarchisationModPoly() override {}

        /**
         * Implements the hierarchisation on a sprase grid with mod poly base functions
         *
         * @param node_values the functions values in the node base
         *
         */
        virtual void doHierarchisation(DataVector& node_values) override;

        /**
         * Implements the dehierarchisation on a sprase grid with mod poly base functions
         *
         * @param alpha the coefficients of the sparse grid's base functions
         *
         */
        virtual void doDehierarchisation(DataVector& alpha) override;

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
        /// Mod Poly Basis object
        SPolyModifiedBase base;
    };

  }
}

#endif /* OPERATIONHIERARCHISATIONMODPOLY_HPP */
