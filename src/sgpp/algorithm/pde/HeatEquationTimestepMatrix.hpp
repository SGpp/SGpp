/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HEATEQUATIONTIMESTEPMATRIX_HPP
#define HEATEQUATIONTIMESTEPMATRIX_HPP

#include "data/DataVector.hpp"
#include "grid/Grid.hpp"
#include "operation/common/OperationODESolverMatrix.hpp"

#include <string>

namespace sg
{

/**
 * @todo (heinecke) add description here
 */
class HeatEquationTimestepMatrix : public OperationODESolverMatrix
{
private:
	/// the heat coefficient
	double a;
	/// the Laplace Operation (Stiffness Matrix)
	OperationMatrix* OpLaplace;
	/// the LTwoDotProduct Operation (Mass Matrix)
	OperationMatrix* OpMass;
	/**
	 *  specifies in which solver this matrix is used, valid values are:
	 *  ExEul for explicit Euler
	 *  ImEul for implicit Euler
	 *  CrNic for Crank Nicolson solver
	 */
	std::string tOperationMode;
	/// the size of one timestep used in the ODE Solver
	double TimestepSize;
	/// Pointer to the grid object
	Grid* myGrid;

	/**
	 * @todo (heinecke) add description
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyMassMatrix(DataVector& alpha, DataVector& result);

	/**
	 * @todo (heinecke) add description
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyLOperator(DataVector& alpha, DataVector& result);

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param a the heat coefficient
	 * @param TimestepSize the size of one timestep used in the ODE Solver
	 * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
	 *  							ImEul for implicit Euler, CrNic for Crank Nicolson solver
	 */
	HeatEquationTimestepMatrix(Grid& SparseGrid, double a, double TimestepSize, std::string OperationMode = "ExEul");

	/**
	 * Std-Destructor
	 */
	virtual ~HeatEquationTimestepMatrix();

	virtual void mult(DataVector& alpha, DataVector& result);

	virtual void generateRHS(DataVector& data, DataVector& rhs);

	virtual void finishTimestep(DataVector& alpha);

	virtual void startTimestep(DataVector& alpha);

	virtual Grid* getGrid();
};

}

#endif /* HEATEQUATIONTIMESTEPMATRIX_HPP */
