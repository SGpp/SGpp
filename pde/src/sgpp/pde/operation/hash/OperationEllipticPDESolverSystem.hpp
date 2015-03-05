// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONELLITPICPDESOLVERSYSTEM_HPP
#define OPERATIONELLIPTICPDESOLVERSYSTEM_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
    class OperationEllipticPDESolverSystem : public SGPP::base::OperationMatrix {
      protected:
        /// Pointer to the grid object
        SGPP::base::Grid* BoundGrid;
        /// the right hand side of the system
        SGPP::base::DataVector* rhs;
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
        OperationEllipticPDESolverSystem(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& rhs);

        /**
         * Destructor
         */
        virtual ~OperationEllipticPDESolverSystem();

        /**
         * Multiplicates a vector with the matrix \f$ L \f$
         *
         * @param alpha SGPP::base::DataVector that contains the ansatzfunctions' coefficients
         * @param result SGPP::base::DataVector into which the result of the space discretization operation is stored
         */
        virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;

        /**
         * generates the right hand side of the system
         *
         * @return returns the rhs
         */
        virtual SGPP::base::DataVector* generateRHS() = 0;

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