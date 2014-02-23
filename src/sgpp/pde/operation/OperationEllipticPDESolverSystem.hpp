/* ****************************************************************************
* Copyright (C) 2010-2014 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONELLITPICPDESOLVERSYSTEM_HPP
#define OPERATIONELLIPTICPDESOLVERSYSTEM_HPP

#include "base/grid/Grid.hpp"
#include "base/operation/OperationMatrix.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace pde {

    /**
     * Abstract definition of a System that is used to solve elliptic partial
     * differential equations. So an instance of this class has to pass to
     * any SLE Solver used in SGpp.
     *
     * \f$L \vec{u} = rhs\f$
     *
     * L: space discretization (L-Operator)
     * rhs: right hand side
     */
    class OperationEllipticPDESolverSystem : public sg::base::OperationMatrix {
      protected:
        /// Pointer to the grid object
        sg::base::Grid* BoundGrid;
        /// the right hand side of the system
        sg::base::DataVector* rhs;
        /// Stores number of gridpoints, inner grid
        size_t numGridpointsInner;
        /// Stores number of gridpoints, complete grid
        size_t numGridpointsComplete;

      public:
        /**
         * Constructor
         *
         * @param SparseGrid the grid, for which the system should be solved
         * @param rhs the right hand side of the corresponding system
         */
        OperationEllipticPDESolverSystem(sg::base::Grid& SparseGrid, sg::base::DataVector& rhs);

        /**
         * Destructor
         */
        virtual ~OperationEllipticPDESolverSystem();

        /**
         * Multiplicates a vector with the matrix \f$ L \f$
         *
         * @param alpha sg::base::DataVector that contains the ansatzfunctions' coefficients
         * @param result sg::base::DataVector into which the result of the space discretization operation is stored
         */
        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;

        /**
         * generates the right hand side of the system
         *
         * @return returns the rhs
         */
        virtual sg::base::DataVector* generateRHS() = 0;

        /**
         * Returns the number of grid points for the complete grid
         *
         * @return the number of grid points for the complete grid
         */
        size_t getNumGridPointsComplete();

        /**
         * Returns the number of grid points for the inner grid
         *
         * @return the number of grid points for the inner grid
         */
        size_t getNumGridPointsInner();
    };

  }
}

#endif /* OPERATIONELLITPICPDESOLVERMATRIX_HPP */
