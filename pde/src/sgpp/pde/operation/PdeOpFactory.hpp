// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PDE_OP_FACTORY_HPP
#define PDE_OP_FACTORY_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/base/operation/hash/OperationMatrix.hpp>

/*
 * This file contains factory methods for operations.
 */

#include <sgpp/globaldef.hpp>


namespace SGPP {

namespace op_factory {
/**
 * Factory method, returning an OperationLaplace (OperationMatrix) for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used
 * @return Pointer to the new OperationMatrix object for the Grid grid
 */
base::OperationMatrix* createOperationLaplace(base::Grid& grid);
/**
 * Factory method, returning an OperationLaplace (OperationMatrix) for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used
 * @param coef Coefficient vector for OperationLaplace
 * @return Pointer to the new OperationMatrix object for the Grid grid
 */
base::OperationMatrix* createOperationLaplace(base::Grid& grid,
    SGPP::base::DataVector& coef);
/**
 * Factory method, returning an OperationLTwoDotProduct (OperationMatrix) for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used
 * @return Pointer to the new OperationMatrix object for the Grid grid
 */
base::OperationMatrix* createOperationLTwoDotProduct(base::Grid& grid);

/**
   * Factory method, returning an OperationLTwoDotExplicit (OperationMatrix) for the grid at hand.
   * Note: object has to be freed after use.
   *
   * @param grid Grid which is to be used
   * @return Pointer to the new OperationMatrix object for the Grid grid
   */
base::OperationMatrix* createOperationLTwoDotExplicit(base::Grid& grid);

/**
   * Factory method, returning an OperationLTwoDotExplicit (OperationMatrix) for the grid at hand.
   * Note: object has to be freed after use. The DataMatrix m is used to store the the matrix and
   * is not destroyed if the OperationMatrix is destroyed.
   *
   * @param grid Grid which is to be used
   * @param m DataMatrix in which the data is stored
   * @return Pointer to the new OperationMatrix object for the Grid grid
   */
base::OperationMatrix* createOperationLTwoDotExplicit(base::DataMatrix* m,
    base::Grid& grid);

/**
 * Factory method, returning an OperationLaplace (OperationMatrix) for the grid at hand.
 * Note: object has to be freed after use.
 *
 * This Laplacian is implemented by a multi-dimensional sweep operation
 *
 * @param grid Grid which is to be used
 * @return Pointer to the new OperationMatrix object for the Grid grid
 */
base::OperationMatrix* createOperationLaplaceEnhanced(base::Grid& grid);

/**
 * Factory method, returning an OperationLaplace (OperationMatrix) for the grid at hand.
 * Note: object has to be freed after use.
 *
 * This Laplacian is implemented by a multi-dimensional sweep operation
 *
 * @param grid Grid which is to be used
 * @param coef Coefficient vector for OperationLaplace
 * @return Pointer to the new OperationMatrix object for the Grid grid
 */
base::OperationMatrix* createOperationLaplaceEnhanced(base::Grid& grid,
    SGPP::base::DataVector& coef);

}

}

#endif /*PDE_OP_FACTORY_HPP*/
