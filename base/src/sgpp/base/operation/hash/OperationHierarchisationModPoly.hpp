// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONMODPOLY_HPP
#define OPERATIONHIERARCHISATIONMODPOLY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/basis/modpoly/ModifiedPolyBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Hierarchisation on sparse grid, mod poly case
     *
     * @version $HEAD$
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
        virtual ~OperationHierarchisationModPoly() {}

        /**
         * Implements the hierarchisation on a sprase grid with mod poly base functions
         *
         * @param node_values the functions values in the node base
         *
         * @todo Implement the hierarchisation on the sparse grid with mod poly base functions
         */
        virtual void doHierarchisation(DataVector& node_values);

        /**
         * Implements the dehierarchisation on a sprase grid with mod poly base functions
         *
         * @param alpha the coefficients of the sparse grid's base functions
         *
         * @todo Implement the dehierarchisation on the sparse grid with mod poly base functions
         */
        virtual void doDehierarchisation(DataVector& alpha);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
        /// Mod Poly Basis object
        SModPolyBase base;
    };

  }
}

#endif /* OPERATIONHIERARCHISATIONMODPOLY_HPP */