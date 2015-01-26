/* ****************************************************************************
* Copyright (C) 2010-2014 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONPARABOLICPDESOLVERSYSTEMFREEBOUNDARIES_HPP
#define OPERATIONPARABOLICPDESOLVERSYSTEMFREEBOUNDARIES_HPP

#include <sgpp/pde/operation/OperationParabolicPDESolverSystem.hpp>

namespace sg {
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
         * @param result reference to the sg::base::DataVector into which the result is written
         */
        virtual void applyMassMatrix(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;

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
         */
        OperationParabolicPDESolverSystemFreeBoundaries();

        /**
         * Destructor
         */
        virtual ~OperationParabolicPDESolverSystemFreeBoundaries();

        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual sg::base::DataVector* generateRHS();

        virtual sg::base::DataVector* getGridCoefficientsForCG();
    };

  }
}

#endif /* OPERATIONPARABOLICPDESOLVERSYSTEMFREEBOUNDARIES_HPP */
