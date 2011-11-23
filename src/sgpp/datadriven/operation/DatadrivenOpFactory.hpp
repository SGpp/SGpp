/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef DATADRIVEN_OP_FACTORY_HPP
#define DATADRIVEN_OP_FACTORY_HPP

#include "base/grid/Grid.hpp"

#include "datadriven/operation/OperationTest.hpp"
#include "base/operation/OperationMatrix.hpp"

/*
 * This file contains factory methods for operations.
 */
 
namespace sg
{
namespace op_factory
{
  /**
   * Factory method, returning an OperationTest for the grid at hand.
   * Note: object has to be freed after use.
   *
   * @param grid Grid which is to be used for OperationTest
   * @return Pointer to the new OperationTest object for the Grid grid
   */
  datadriven::OperationTest* createOperationTest(base::Grid& grid);

  /**
   * Factory method, returning an OperationTest for the grid at hand.
   * Note: object has to be freed after use.
   *
   * @param grid Grid which is to be used for OperationRegularizationDiagonal
   * @return Pointer to the new OperationTest object for the Grid grid
   */
  base::OperationMatrix* createOperationRegularizationDiagonal(base::Grid& grid, int mode, double k);

}
}

#endif /*DATADRIVEN_OP_FACTORY_HPP*/
