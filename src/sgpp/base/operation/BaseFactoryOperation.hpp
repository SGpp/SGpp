/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef BASE_FACTORY_OPERATION_HPP
#define BASE_FACTORY_OPERATION_HPP

#include "grid/Grid.hpp"

#include "operation/common/OperationHierarchisation.hpp"
#include "operation/common/OperationQuadrature.hpp"
#include "operation/common/OperationConvert.hpp"
#include "operation/common/OperationIdentity.hpp"
#include "operation/common/OperationMatrix.hpp"
#include "operation/common/OperationEval.hpp"
#include "operation/datadriven/OperationMultipleEval.hpp"

namespace sg
{

namespace GOperationFactory
{

  /**
   * Factory method, returning an OperationHierarchisation for the grid at hand.
   * Note: object has to bee freed after use.
   *
   * @param grid Grid which is to be used for hierarchisation
   * @return Pointer to the new OperationHierarchisation object for the Grid grid
   */
  base::OperationHierarchisation* createOperationHierarchisation(base::Grid& grid);
  /**
   * Factory method, returning an OperationQuadrature for the grid at hand.
   * Note: object has to bee freed after use.
   *
   * @param grid Grid which is to be used for quadrature
   * @return Pointer to the new OperationQuadrature for the Grid grid
   */
  base::OperationQuadrature* createOperationQuadrature(base::Grid& grid);
  /**
   * Factory method, returning an OperationConvert for the grid at hand.
   * Note: object has to bee freed after use.
   *
   * @param grid Grid which is to be used for conversion
   * @return Pointer to the new OperationConvert object for the Grid grid
   */
  base::OperationConvert* createOperationConvert(base::Grid& grid);
  /**
   * Factory method, returning an OperationIdentity for the grid at hand.
   * Note: object has to bee freed after use.
   * Just calls OperationIdentity() independent of grid; factory method 
   * provided for uniform use.
   *
   * @return Pointer to the new OperationIdentity object
   */
  base::OperationMatrix* createOperationIdentity(base::Grid& grid);
  /**
   * Factory method, returning an OperationEval for the grid at hand.
   * Note: object has to bee freed after use.
   *
   * @param grid Grid which is to be used
   * @return Pointer to the new OperationEval object for the Grid grid
   */
  base::OperationEval* createOperationEval(base::Grid& grid);
  /**
   * Factory method, returning an OperationMultipleEval for the grid at hand.
   * Note: object has to bee freed after use.
   *
   * @param grid Grid which is to be used
   * @param dataset The dataset (DataMatrix, one datapoint per row) that is to be evaluated for
   * the sparse grid function
   * @return Pointer to the new OperationMultipleEval object for the Grid grid
   */
  base::OperationMultipleEval* createOperationMultipleEval(base::Grid& grid, base::DataMatrix* dataset);
}

}

#endif /*BASE_FACTORY_OPERATION_HPP*/
