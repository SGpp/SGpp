/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONELLITPICPDESOLVERSYSTEM_HPP
#define OPERATIONELLIPTICPDESOLVERSYSTEM_HPP

#include "grid/Grid.hpp"
#include "operation/common/OperationMatrix.hpp"
#include "data/DataVector.hpp"
using namespace sg::base;

namespace sg
{

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
class OperationEllipticPDESolverSystem : public OperationMatrix
{
protected:
	/// the right hand side of the system
	DataVector* rhs;
	/// Pointer to the grid object
	Grid* BoundGrid;
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
	OperationEllipticPDESolverSystem(Grid& SparseGrid, DataVector& rhs);

	/**
	 * Destructor
	 */
	virtual ~OperationEllipticPDESolverSystem();

	/**
	 * Multiplicates a vector with the matrix \f$ L \f$
	 *
	 * @param alpha DataVector that contains the ansatzfunctions' coefficients
	 * @param result DataVector into which the result of the space discretization operation is stored
	 */
	virtual void mult(DataVector& alpha, DataVector& result) = 0;

	/**
	 * generates the right hand side of the system
	 *
	 * @return returns the rhs
	 */
	virtual DataVector* generateRHS() = 0;

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

#endif /* OPERATIONELLITPICPDESOLVERMATRIX_HPP */
