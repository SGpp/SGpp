/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONHIERARCHISATIONMODBSPLINE_HPP
#define OPERATIONHIERARCHISATIONMODBSPLINE_HPP

#include "base/operation/OperationHierarchisation.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/basis/modbspline/ModifiedBsplineBasis.hpp"
#include "base/datatypes/DataVector.hpp"


namespace sg {
  namespace base {

    /**
     * Hierarchisation on sparse grid, mod bspline case
     *
     * @version $HEAD$
     */
    class OperationHierarchisationModBspline : public OperationHierarchisation {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's GridStorage object
         * @param degree the bsplinenom's max. degree
         */
        OperationHierarchisationModBspline(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

        /**
         * Destructor
         */
        virtual ~OperationHierarchisationModBspline() {}

        /**
         * Implements the hierarchisation on a sprase grid with mod bspline base functions
         *
         * @param node_values the functions values in the node base
         *
         * @todo Implement the hierarchisation on the sparse grid with mod bspline base functions
         */
        virtual void doHierarchisation(DataVector& node_values);

        /**
         * Implements the dehierarchisation on a sprase grid with mod bspline base functions
         *
         * @param alpha the coefficients of the sparse grid's base functions
         *
         * @todo Implement the dehierarchisation on the sparse grid with mod bspline base functions
         */
        virtual void doDehierarchisation(DataVector& alpha);

      protected:
        /// Pointer to GridStorage object
        GridStorage* storage;
        /// Mod Bspline Basis object
        SModBsplineBase base;
    };

  }
}

#endif /* OPERATIONHIERARCHISATIONMODBSPLINE_HPP */
