/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi (qic@in.tum.de)

#ifndef HULLWHITEODESOLVERSYSTEM_HPP
#define HULLWHITEODESOLVERSYSTEM_HPP

#include "grid/Grid.hpp"
#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"
#include "grid/common/DirichletUpdateVector.hpp"
#include "operation/pde/OperationODESolverSystemNeumann.hpp"

namespace sg
{

/**
 * This class implements the ODESolverSystem for the HullWhite
 * Equation.
 */
class HullWhiteODESolverSystem : public OperationODESolverSystemNeumann
{
protected:
	double theta;
	double sigma;
	double a;
	/// the B matrix Operation, on boundary grid
	OperationMatrix* OpBBound;
	/// the D matrix Operation, on boundary grid
	OperationMatrix* OpDBound;
	/// the E matrix Operation, on boundary grid
	OperationMatrix* OpEBound;
	/// the F matrix Operation, on boundary grid
	OperationMatrix* OpFBound;
	/// the LTwoDotProduct Operation (Mass Matrix A), on boundary grid
	OperationMatrix* OpLTwoBound;
	/// use coarsening between timesteps in order to reduce gridsize
	bool useCoarsen;
	/// Threshold used to decide if a grid point should be deleted
	double coarsenThreshold;
	/// Percent how many of the removable points should be tested for deletion
	double coarsenPercent;
	/// denotes the number of complete coarsen procedures per timestep
	size_t numExecCoarsen;

	std::vector<size_t> HWalgoDims;
	/// Routine to modify the boundaries/inner points of the grid
	DirichletUpdateVector* BoundaryUpdate;

	virtual void applyLOperator(DataVector& alpha, DataVector& result);

	virtual void applyMassMatrix(DataVector& alpha, DataVector& result);

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param alpha the ansatzfunctions' coefficients
	 * @param theta reference to the theta
	 * @param sigma reference to the sigma
	 * @param a reference to the a
	 * @param TimestepSize the size of one timestep used in the ODE Solver
	 * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
	 *  							ImEul for implicit Euler, CrNic for Crank Nicolson solver
	 * @param useCoarsen specifies if the grid should be coarsened between timesteps
	 * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
	 * @param coarsenPercent Number of removable grid points that should be tested for deletion
	 * @param numExecCoarsen denotes the number of complete coarsen procedures per timestep
	 */
	HullWhiteODESolverSystem(Grid& SparseGrid, DataVector& alpha, double sigma, double theta,
		    double a, double TimestepSize, std::string OperationMode = "ExEul",
			bool useCoarsen = false, double coarsenThreshold = 0.0, double coarsenPercent = 0.0,
			size_t numExecCoarsen = 0);

	/**
	 * Std-Destructor
	 */
	virtual ~HullWhiteODESolverSystem();

	void finishTimestep(bool isLastTimestep = false);

	void startTimestep();
};

}

#endif /* HULLWHITEODESOLVERSYSTEM_HPP */
