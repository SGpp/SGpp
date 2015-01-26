// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONPARABOLICPDESOLVERSYSTEMFREEBOUNDARIES_HPP
#define OPERATIONPARABOLICPDESOLVERSYSTEMFREEBOUNDARIES_HPP

#include <sgpp/pde/operation/OperationParabolicPDESolverSystem.hpp>

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
     */
    class OperationParabolicPDESolverSystemFreeBoundaries : public OperationParabolicPDESolverSystem {
      protected:
        /**
         * applies the PDE's mass matrix, on complete grid - with boundaries
         *
         * @param alpha the coefficients of the sparse grid's ansatzfunctions
         * @param result reference to the SGPP::base::DataVector into which the result is written
         */
        virtual void applyMassMatrix(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;

        /**
         * applies the PDE's system matrix, on complete grid - with boundaries
         *
         * @param alpha the coefficients of the sparse grid's ansatzfunctions
         * @param result reference to the SGPP::base::DataVector into which the result is written
         */
        virtual void applyLOperator(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;

      public:
        /**
         * Constructor
         */
        OperationParabolicPDESolverSystemFreeBoundaries();

        /**
         * Destructor
         */
        virtual ~OperationParabolicPDESolverSystemFreeBoundaries();

        virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

        virtual SGPP::base::DataVector* generateRHS();

        virtual SGPP::base::DataVector* getGridCoefficientsForCG();
    };

  }
}

#endif /* OPERATIONPARABOLICPDESOLVERSYSTEMFREEBOUNDARIES_HPP */