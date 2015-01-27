// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef BASE_OP_FACTORY_HPP
#define BASE_OP_FACTORY_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/operation/hash/OperationSecondMoment.hpp>
#include <sgpp/base/operation/hash/OperationConvert.hpp>
#include <sgpp/base/operation/hash/OperationIdentity.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/hash/OperationStencilHierarchisation.hpp>

/*
 * This file contains factory methods for operations.
 */

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace op_factory {

    /**
     * Factory method, returning an OperationHierarchisation for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used for hierarchisation
     * @return Pointer to the new OperationHierarchisation object for the Grid grid
     */
    base::OperationHierarchisation* createOperationHierarchisation(base::Grid& grid);
    /**
     * Factory method, returning an OperationQuadrature for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used for quadrature
     * @return Pointer to the new OperationQuadrature for the Grid grid
     */
    base::OperationQuadrature* createOperationQuadrature(base::Grid& grid);
    /**
     * Factory method, returning an OperationFirstMoment for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used for quadrature
     * @return Pointer to the new OperationFirstMoment for the Grid grid
     */
    base::OperationFirstMoment* createOperationFirstMoment(base::Grid& grid);
    /**
     * Factory method, returning an OperationSecondMoment for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used for quadrature
     * @return Pointer to the new OperationSecondMoment for the Grid grid
     */
    base::OperationSecondMoment* createOperationSecondMoment(base::Grid& grid);
    /**
     * Factory method, returning an OperationConvert for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used for conversion
     * @return Pointer to the new OperationConvert object for the Grid grid
     */
    base::OperationConvert* createOperationConvert(base::Grid& grid);
    /**
     * Factory method, returning an OperationIdentity for the grid at hand.
     * Note: object has to be freed after use.
     * Just calls OperationIdentity() independent of grid; factory method
     * provided for uniform use.
     *
     * @return Pointer to the new OperationIdentity object
     */
    base::OperationMatrix* createOperationIdentity(base::Grid& grid);
    /**
     * Factory method, returning an OperationEval for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used
     * @return Pointer to the new OperationEval object for the Grid grid
     */
    base::OperationEval* createOperationEval(base::Grid& grid);
    /**
     * Factory method, returning an OperationMultipleEval for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used
     * @param dataset The dataset (DataMatrix, one datapoint per row) that is to be evaluated for
     * the sparse grid function
     * @return Pointer to the new OperationMultipleEval object for the Grid grid
     */
    base::OperationMultipleEval* createOperationMultipleEval(base::Grid& grid, base::DataMatrix &dataset);
  }

}

#endif /*BASE_OP_FACTORY_HPP*/