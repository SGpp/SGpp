/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETPARALLELMPI_HPP
#define POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETPARALLELMPI_HPP

#include "operation/pde/OperationEllipticPDESolverSystemDirichlet.hpp"

namespace sg
{

/**
 * This class uses OperationEllipticPDESolverSystemDirichlet
 * to define a solver system for the Poission Equation.
 *
 * For the mult-routine only the Laplace-Operator is required
 *
 * There is a parallelization over all operators, required
 * to solve the poisson equation.
 */
class PoissonEquationEllipticPDESolverSystemDirichletParallelMPI : public OperationEllipticPDESolverSystemDirichlet
{
protected:
	OperationMatrix* Laplace_Inner;
	OperationMatrix* Laplace_Complete;

	void applyLOperatorComplete(DataVector& alpha, DataVector& result);

	void applyLOperatorInner(DataVector& alpha, DataVector& result);

public:
	/**
	 * Constructor
	 *
	 * @param SparseGrid reference to a sparse grid on which the Poisson Equation should be solved
	 * @param rhs the right hand side for solving the elliptic PDE
	 */
	PoissonEquationEllipticPDESolverSystemDirichletParallelMPI(Grid& SparseGrid, DataVector& rhs);

	/**
	 * Destructor
	 */
	virtual ~PoissonEquationEllipticPDESolverSystemDirichletParallelMPI();

	void mult(DataVector& alpha, DataVector& result);

	DataVector* generateRHS();
};

}

#endif /* POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLETPARALLELMPI_HPP */
