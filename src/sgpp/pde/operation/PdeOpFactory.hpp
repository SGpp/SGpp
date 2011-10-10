/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef PDE_OP_FACTORY_HPP
#define PDE_OP_FACTORY_HPP

#include "grid/Grid.hpp"

#include "operation/common/OperationMatrix.hpp"

/*
 * This file contains factory methods for operations.
 */
 
namespace sg
{

namespace op_factory
{
  /**
   * Factory method, returning an OperationLaplace (OperationMatrix) for the grid at hand.
   * Note: object has to bee freed after use.
   *
   * @param grid Grid which is to be used
   * @return Pointer to the new OperationMatrix object for the Grid grid
   */
  base::OperationMatrix* createOperationLaplace(base::Grid& grid);
  /**
   * Factory method, returning an OperationLaplace (OperationMatrix) for the grid at hand.
   * Note: object has to bee freed after use.
   *
   * @param grid Grid which is to be used
   * @param coef Coefficient vector for OperationLaplace
   * @return Pointer to the new OperationMatrix object for the Grid grid
   */
  base::OperationMatrix* createOperationLaplace(base::Grid& grid, sg::base::DataVector& coef);
  /**
   * Factory method, returning an OperationLTwoDotProduct (OperationMatrix) for the grid at hand.
   * Note: object has to bee freed after use.
   *
   * @param grid Grid which is to be used
   * @return Pointer to the new OperationMatrix object for the Grid grid
   */
  base::OperationMatrix* createOperationLTwoDotProduct(base::Grid& grid);
  /**
   * Factory method, returning an OperationUpDownTest (OperationMatrix) for the grid at hand.
   * Note: object has to bee freed after use.
   *
   * @param grid Grid which is to be used
   * @return Pointer to the new OperationMatrix object for the Grid grid
   */
  base::OperationMatrix* createOperationUpDownTest(base::Grid& grid);

}

}

#endif /*PDE_OP_FACTORY_HPP*/
