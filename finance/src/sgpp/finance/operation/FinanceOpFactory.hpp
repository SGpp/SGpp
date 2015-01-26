// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef FINANCE_OP_FACTORY_HPP
#define FINANCE_OP_FACTORY_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/base/operation/OperationMatrix.hpp>

/*
 * This file contains factory methods for operations.
 */

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace op_factory {
    /**
     * Factory method, returning an OperationGamma (OperationMatrix) for the grid at hand.
     * Note: object has to be freed after use.
     *
     * This operation calculates the following bilinear form which is
     * needed to solve the Black Scholes equation:
     * \f$ \int_{\Omega} S_i S_j \frac{\partial u(\vec{s}}{\partial S_i} \frac{\partial v(\vec{s}}{\partial S_j} d \vec{s}\f$
     *
     * @param grid Grid which is to be used
     * @param coef Reference to a DataMatrix object that contains the constant coeffecients of this bilinear form
     * @return Pointer to the new OperationMatrix object for the Grid grid
     */
    base::OperationMatrix* createOperationGamma(base::Grid& grid, base::DataMatrix& coef);
    /**
     * Factory method, returning an OperationGammaLog (OperationMatrix) for the grid at hand.
     * Note: object has to be freed after use.
     *
     * This operation calculates the following bilinear form which is
     * needed to solve the Black Scholes equation:
     * \f$ \int_{\Omega} \frac{\partial u(\vec{s}}{\partial S_i} \frac{\partial v(\vec{s}}{\partial S_j} d \vec{s}\f$
     *
     * @param grid Grid which is to be used
     * @param coef Reference to a DataMatrix object that contains the constant coeffecients of this bilinear form
     * @return Pointer to the new OperationMatrix object for the Grid grid
     */
    base::OperationMatrix* createOperationGammaLog(base::Grid& grid, base::DataMatrix& coef);
    /**
     * Factory method, returning an OperationLB (OperationMatrix) for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used
     * @return Pointer to the new OperationMatrix object for the Grid grid
     */
    base::OperationMatrix* createOperationLB(base::Grid& grid);
    /**
     * Factory method, returning an OperationLE (OperationMatrix) for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used
     * @return Pointer to the new OperationMatrix object for the Grid grid
     */
    base::OperationMatrix* createOperationLE(base::Grid& grid);
    /**
     * Factory method, returning an OperationLD (OperationMatrix) for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used
     * @return Pointer to the new OperationMatrix object for the Grid grid
     */
    base::OperationMatrix* createOperationLD(base::Grid& grid);
    /**
     * Factory method, returning an OperationLF (OperationMatrix) for the grid at hand.
     * Note: object has to be freed after use.
     *
     * @param grid Grid which is to be used
     * @return Pointer to the new OperationMatrix object for the Grid grid
     */
    base::OperationMatrix* createOperationLF(base::Grid& grid);
    /**
     * Factory method, returning an OperationDelta (OperationMatrix) for the grid at hand.
     * Note: object has to be freed after use.
     *
     * This operation calculates the following bilinear form which is
     * needed to solve the Black Scholes equation:
     * \f$ \int_{\Omega} S_i v(\vec{s}) \frac{\partial u(\vec{s}}{\partial S_i} d \vec{s}\f$
     *
     * @param grid Grid which is to be used
     * @param coef Reference to a DataMatrix object that contains the constant coeffecients of this bilinear form
     * @return Pointer to the new OperationMatrix object for the Grid grid
     */
    base::OperationMatrix* createOperationDelta(base::Grid& grid, base::DataVector& coef);
    /**
     * Factory method, returning an OperationDeltaLog (OperationMatrix) for the grid at hand.
     * Note: object has to be freed after use.
     *
     * This operation calculates the following bilinear form which is
     * needed to solve the Black Scholes equation:
     * \f$ \int_{\Omega} \frac{\partial u(\vec{s}}{\partial S_i} v(\vec{s}) d \vec{s}\f$
     *
     * @param grid Grid which is to be used
     * @param coef Reference to a DataMatrix object that contains the constant coeffecients of this bilinear form
     * @return Pointer to the new OperationMatrix object for the Grid grid
     */
    base::OperationMatrix* createOperationDeltaLog(base::Grid& grid, base::DataVector& coef);

    base::OperationMatrix* createOperationHestonBLog(base::Grid& grid, base::DataMatrix& coef);
    base::OperationMatrix* createOperationHestonCLog(base::Grid& grid, base::DataMatrix& coef);
    base::OperationMatrix* createOperationHestonDLog(base::Grid& grid, base::DataVector& coef);
    base::OperationMatrix* createOperationHestonELog(base::Grid& grid, base::DataVector& coef);
    base::OperationMatrix* createOperationHestonFLog(base::Grid& grid, base::DataVector& coef);
    base::OperationMatrix* createOperationHestonGLog(base::Grid& grid, base::DataVector& coef);
    base::OperationMatrix* createOperationHestonHLog(base::Grid& grid, base::DataMatrix& coef);
    base::OperationMatrix* createOperationHestonKLog(base::Grid& grid, double**** * coef);
    base::OperationMatrix* createOperationHestonX(base::Grid& grid, base::DataMatrix& coef);
    base::OperationMatrix* createOperationHestonY(base::Grid& grid, base::DataMatrix& coef);
    base::OperationMatrix* createOperationHestonW(base::Grid& grid, base::DataMatrix& coef);
    base::OperationMatrix* createOperationHestonZ(base::Grid& grid, base::DataVector& coef);

  }

}

#endif /*FINANCE_OP_FACTORY_HPP*/