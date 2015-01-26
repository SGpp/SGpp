/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETVECTORIZEDMPI_HPP
#define POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETVECTORIZEDMPI_HPP

#include <sgpp/pde/operation/OperationEllipticPDESolverSystemDirichlet.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    /**
     * This class uses OperationEllipticPDESolverSystemDirichlet
     * to define a solver system for the Poission Equation.
     *
     * For the mult-routine only the Laplace-Operator is required
     */
    class PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI : public SGPP::pde::OperationEllipticPDESolverSystemDirichlet {
      protected:
        SGPP::base::OperationMatrix* Laplace_Inner;
        SGPP::base::OperationMatrix* Laplace_Complete;

        void applyLOperatorComplete(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

        void applyLOperatorInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

      public:
        /**
         * Constructor
         *
         * @param SparseGrid reference to a sparse grid on which the Poisson Equation should be solved
         * @param rhs the right hand side for solving the elliptic PDE
         */
        PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& rhs);

        /**
         * Destructor
         */
        virtual ~PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI();
    };

  }
}

#endif /* POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETVECTORIZEDMPI_HPP */
