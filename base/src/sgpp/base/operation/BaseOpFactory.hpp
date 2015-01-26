/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef BASE_OP_FACTORY_HPP
#define BASE_OP_FACTORY_HPP

#include "base/grid/Grid.hpp"

#include "base/operation/OperationHierarchisation.hpp"
#include "base/operation/OperationQuadrature.hpp"
#include "base/operation/OperationFirstMoment.hpp"
#include "base/operation/OperationSecondMoment.hpp"
#include "base/operation/OperationConvert.hpp"
#include "base/operation/OperationIdentity.hpp"
#include "base/operation/OperationMatrix.hpp"
#include "base/operation/OperationEval.hpp"
#include "base/operation/OperationMultipleEval.hpp"
#include "base/operation/OperationStencilHierarchisation.hpp"

/*
 * This file contains factory methods for operations.
 */

namespace sg {

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
