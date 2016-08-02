// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// SGPPpp.sparsegrids.org

#ifndef COMBIGRID_OP_FACTORY_HPP
#define COMBIGRID_OP_FACTORY_HPP

/*
 * This file contains factory methods for operations.
 */

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace op_factory {

    /**
     * Factory method, returning an OperationQuadrature for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used for quadrature
     * @return Pointer to the new OperationQuadrature for the Grid grid
     */
    // base::OperationQuadrature* createOperationQuadrature(base::Grid& grid);

    /**
     * Factory method, returning an OperationEval for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used
     * @return Pointer to the new OperationEval object for the Grid grid
     */
    // base::OperationEval* createOperationEval(base::Grid& grid);
  }

}

#endif /*COMBIGRID_OP_FACTORY_HPP*/
