/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONHIERARCHISATIONMODPOLY_HPP
#define OPERATIONHIERARCHISATIONMODPOLY_HPP

#include "base/operation/OperationHierarchisation.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/modpoly/ModifiedPolyBasis.hpp"
#include "base/datatypes/DataVector.hpp"


namespace sg {
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
