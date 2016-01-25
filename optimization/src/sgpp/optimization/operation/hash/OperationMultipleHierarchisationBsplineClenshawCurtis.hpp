// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPERATION_HASH_OPERATIONMULTIPLEHIERARCHISATIONBSPLINECLENSHAWCURTIS_HPP
#define SGPP_OPTIMIZATION_OPERATION_HASH_OPERATIONMULTIPLEHIERARCHISATIONBSPLINECLENSHAWCURTIS_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisation.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Hierarchisation operation for B-spline basis functions on
     * Clenshaw-Curtis grids.
     */
    class OperationMultipleHierarchisationBsplineClenshawCurtis :
      public OperationMultipleHierarchisation {
      public:
        /**
         * Constructor.
         *
         * @param grid      grid
         */
        OperationMultipleHierarchisationBsplineClenshawCurtis(
          base::BsplineClenshawCurtisGrid& grid);

        /**
         * Virtual destructor.
         */
        virtual ~OperationMultipleHierarchisationBsplineClenshawCurtis();

        /**
         * @param[in,out] nodeValues before: vector of function values at
         *                           the grid points,
         *                           after: vector of hierarchical coefficients
         * @return                   whether hierarchisation was successful
         */
        virtual bool doHierarchisation(base::DataVector& nodeValues);

        /**
         * @param[in,out] alpha before: vector of hierarchical coefficients,
         *                      after: vector of function values at
         *                      the grid points
         */
        virtual void doDehierarchisation(base::DataVector& alpha);

        /**
         * @param[in,out] nodeValues before: matrix of function values at
         *                           the grid points,
         *                           after: matrix of hierarchical coefficients
         * @return                   whether hierarchisation was successful
         */
        virtual bool doHierarchisation(base::DataMatrix& nodeValues);

        /**
         * @param[in,out] alpha before: matrix of hierarchical coefficients,
         *                      after: matrix of function values at
         *                      the grid points
         */
        virtual void doDehierarchisation(base::DataMatrix& alpha);

      protected:
        /// storage of the sparse grid
        base::BsplineClenshawCurtisGrid& grid;
    };

  }
}

#endif /* SGPP_OPTIMIZATION_OPERATION_HASH_OPERATIONMULTIPLEHIERARCHISATIONBSPLINECLENSHAWCURTIS_HPP */
