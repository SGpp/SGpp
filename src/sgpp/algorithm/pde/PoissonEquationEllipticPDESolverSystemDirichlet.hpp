/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP
#define POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP

#include "operation/pde/OperationEllipticPDESolverSystemDirichlet.hpp"

namespace sg
{

/**
 * This class uses OperationEllipticPDESolverSystemDirichlet
 * to define a solver system for the Poission Equation.
 *
 * For the mult-routine only the Laplace-Operator is required
 */
class PoissonEquationEllipticPDESolverSystemDirichlet : public OperationEllipticPDESolverSystemDirichlet
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
	 * @param rhs the right hand side for solving the elliptic PDE
	 */
	PoissonEquationEllipticPDESolverSystemDirichlet(Grid& SparseGrid, DataVector& rhs);

	/**
	 * Destructor
	 */
	virtual ~PoissonEquationEllipticPDESolverSystemDirichlet();
};

}

#endif /* POISSONEQUATIONELLIPTICPDESOLVERSYSTEMDIRICHLET_HPP */
