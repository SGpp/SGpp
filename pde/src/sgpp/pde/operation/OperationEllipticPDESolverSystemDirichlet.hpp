/* ****************************************************************************
* Copyright (C) 2010-2014 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP
#define OPERATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP

#include <sgpp/pde/operation/OperationEllipticPDESolverSystem.hpp>
#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/base/grid/common/DirichletGridConverter.hpp>

#include <string>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    /**
     * Defines a System that is used to solve elliptic partial
     * differential equations. So an instance of this class has to pass to
     * any SLE Solver used in SGpp.
     *
     * \f$L \vec{u} = rhs\f$
     *
     * L: space discretization (L-Operator)
     * rhs: right hand side
     *
     * This class is a specialized version of OperationEllipticPDESolverSystem which
     * exploits Dirichlet boundary conditions. Since there are no degrees of freedom
     * on on the boundaries the iterative solver (CG or BiCGSTAB) has only to take
     * inner grid points into account.
     *
     * The inner grid is constructed during the constructor call!
     */
    class OperationEllipticPDESolverSystemDirichlet : public OperationEllipticPDESolverSystem {
      protected:
        /// Pointer to the alphas (ansatzfunctions' coefficients; inner points only)
        SGPP::base::DataVector* alpha_inner;
        /// Routine to modify the boundaries/inner points of the grid
        SGPP::base::DirichletUpdateVector* BoundaryUpdate;
        /// Class that allows a simple conversion between a grid with and a without boundary points
        SGPP::base::DirichletGridConverter* GridConverter;
        /// Pointer to the inner grid object
        SGPP::base::Grid* InnerGrid;
        /// rhs for the inner grid
        SGPP::base::DataVector* rhs_inner;

        /**
         * applies the PDE's system matrix, on complete grid - with boundaries
         *
         * @param alpha the coefficients of the sparse grid's ansatzfunctions
         * @param result reference to the SGPP::base::DataVector into which the result is written
         */
        virtual void applyLOperatorComplete(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;

        /**
         * applies the PDE's system matrix, on inner grid only
         *
         * @param alpha the coefficients of the sparse grid's ansatzfunctions
         * @param result reference to the SGPP::base::DataVector into which the result is written
         */
        virtual void applyLOperatorInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;

      public:
        /**
         * Constructor
         *
         * @param SparseGrid the grid, for which the system should be solved
         * @param rhs the right hand side of the corresponding system
         */
        OperationEllipticPDESolverSystemDirichlet(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& rhs);

        /**
         * Destructor
         */
        virtual ~OperationEllipticPDESolverSystemDirichlet();

        virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

        virtual SGPP::base::DataVector* generateRHS();

        /**
         * gets a pointer to the sparse grids coefficients used in the CG method to solve
         * one timestep. This is useful because (direchlet) boundaries can be skipped when
         * solving the system.
         *
         * @return alpha vector for CG method
         */
        virtual SGPP::base::DataVector* getGridCoefficientsForCG();

        /**
         * Gets the solution for the complete grid
         *
         * @param Solution SGPP::base::DataVector that must have a dimension equal to the bound's grid dimension, the result is written to Solution
         * @param SolutionInner Solution on the inner grid
         */
        virtual void getSolutionBoundGrid(SGPP::base::DataVector& Solution, SGPP::base::DataVector& SolutionInner);
    };
  }
}

#endif /* OPERATIONELLIPTICPDESOLVERMATRIXDIRICHLET_HPP */
