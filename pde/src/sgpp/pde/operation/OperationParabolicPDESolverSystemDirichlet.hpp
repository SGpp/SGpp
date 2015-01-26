/* ****************************************************************************
* Copyright (C) 2010-2014 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONPARABOLICPDESOLVERSYSTEMDIRICHLET_HPP
#define OPERATIONPARABOLICPDESOLVERSYSTEMDIRICHLET_HPP

#include <sgpp/pde/operation/OperationParabolicPDESolverSystem.hpp>
#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/base/grid/common/DirichletGridConverter.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    /**
     * Defines a System that is used to solve parabolic partial
     * differential equations. So an instance of this class has to pass to
     * any ODE Solver used in SGpp.
     *
     * \f$A \dot{u} = L \vec{u}\f$
     *
     * A: mass matrix
     * L: space discretization (L-Operator)
     *
     * This class defines an elliptic problem in every timestep which is solved
     * using an iterative SLE solver, that solving step is integrated in the
     * ODE Solver.
     *
     * This class is a specialized version of OperationParabolicPDESolverSystem which
     * exploit Dirichlet boundary conditions. Hence there are no degrees of freedom
     * on on the boundaries the iterative solver (CG or BiCGSTAB) has only to take
     * inner grid points into account.
     */
    class OperationParabolicPDESolverSystemDirichlet : public OperationParabolicPDESolverSystem {
      protected:
        /// Pointer to the alphas (ansatzfunctions' coefficients; inner points only)
        SGPP::base::DataVector* alpha_inner;
        /// Routine to modify the boundaries/inner points of the grid
        SGPP::base::DirichletUpdateVector* BoundaryUpdate;
        /// Class that allows a simple conversion between a grid with and a without boundary points
        SGPP::base::DirichletGridConverter* GridConverter;
        /// Pointer to the inner grid object
        SGPP::base::Grid* InnerGrid;

        /**
         * applies the PDE's mass matrix, on complete grid - with boundaries
         *
         * @param alpha the coefficients of the sparse grid's ansatzfunctions
         * @param result reference to the SGPP::base::DataVector into which the result is written
         */
        virtual void applyMassMatrixComplete(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;

        /**
         * applies the PDE's system matrix, on complete grid - with boundaries
         *
         * @param alpha the coefficients of the sparse grid's ansatzfunctions
         * @param result reference to the SGPP::base::DataVector into which the result is written
         */
        virtual void applyLOperatorComplete(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;

        /**
         * applies the PDE's mass matrix, on inner grid only
         *
         * @param alpha the coefficients of the sparse grid's ansatzfunctions
         * @param result reference to the SGPP::base::DataVector into which the result is written
         */
        virtual void applyMassMatrixInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;

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
         */
        OperationParabolicPDESolverSystemDirichlet();

        /**
         * Destructor
         */
        virtual ~OperationParabolicPDESolverSystemDirichlet();

        /**
         * Multiplicates a vector with the matrix
         *
         * @param alpha SGPP::base::DataVector that contains the ansatzfunctions' coefficients
         * @param result SGPP::base::DataVector into which the result of the space discretization operation is stored
         */
        virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

        /**
         * generates the right hand side of the system
         *
         * @return returns the rhs
         */
        virtual SGPP::base::DataVector* generateRHS();

        virtual SGPP::base::DataVector* getGridCoefficientsForCG();
    };

  }
}

#endif /* OPERATIONPARABOLICPDESOLVERSYSTEMDIRICHLET_HPP */
