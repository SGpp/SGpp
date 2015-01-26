/* ****************************************************************************
* Copyright (C) 2010-2014 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONELLIPTICPDESOLVERSYSTEMFREEBOUNDARIES_HPP
#define OPERATIONELLITPICPDESOLVERSYSTEMFREEBOUNDARIES_HPP

#include "pde/operation/OperationEllipticPDESolverSystem.hpp"

namespace sg {
  namespace pde {

    /**
     * Defines a System that is used to solve elliptic partial
     * differential equations. So an instance of this class has to pass to
     * any SLE Solver used in SGpp, here degrees of freedom exists on
     * the boundaries!
     *
     * \f$L \vec{u} = rhs\f$
     *
     * L: space discretization (L-Operator)
     * rhs: right hand sider)
     *
     */
    class OperationEllipticPDESolverSystemFreeBoundaries : public OperationEllipticPDESolverSystem {
      protected:

        /**
         * applies the PDE's system matrix, on complete grid - with boundaries
         *
         * @param alpha the coefficients of the sparse grid's ansatzfunctions
         * @param result reference to the sg::base::DataVector into which the result is written
         */
        virtual void applyLOperator(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;

      public:
        /**
         * Constructor
         *
         * @param SparseGrid the grid, for which the system should be solved
         * @param rhs the right hand side of the corresponding system
         */
        OperationEllipticPDESolverSystemFreeBoundaries(sg::base::Grid& SparseGrid, sg::base::DataVector& rhs);

        /**
         * Destructor
         */
        virtual ~OperationEllipticPDESolverSystemFreeBoundaries();

        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual sg::base::DataVector* generateRHS();
    };

  }
}

#endif /* OPERATIONELLITPTICPDESOLVERMATRIXFREEBOUNDARIES_HPP */
