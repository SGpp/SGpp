/* *****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                          *
* This file is part of the SG++ project. For conditions of distribution and    *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp           *
***************************************************************************** */
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#ifndef MISC_OP_FACTORY_HPP
#define MISC_OP_FACTORY_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/OperationMatrix.hpp>

/*
 * This file contains factory methods for operations.
 */
namespace sg {

  namespace op_factory {
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
    base::OperationMatrix* createOperationLaplaceEnhanced(base::Grid& grid, sg::base::DataVector& coef);
  }

}

#endif /* EXPERIMENTAL_OP_FACTORY_HPP */
