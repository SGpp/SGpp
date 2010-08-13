/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HULLWHITEODESOLVERSYSTEM_HPP
#define HULLWHITEODESOLVERSYSTEM_HPP

#include "grid/Grid.hpp"
#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"
#include "operation/pde/OperationODESolverSystem.hpp"

namespace sg
{

/**
 * This class implements the ODESolverSystem for the HullWhite
 * Equation.
 */
class HullWhiteODESolverSystem : public OperationODESolverSystem
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
	/// the B matrix Operation, on Inner grid
	OperationMatrix* OpBInner;
	/// the D matrix Operation, on Inner grid
	OperationMatrix* OpDInner;
	/// the E matrix Operation, on Inner grid
	OperationMatrix* OpEInner;
	/// the F matrix Operation, on Inner grid
	OperationMatrix* OpFInner;
	/// the LTwoDotProduct Operation (Mass Matrix), on Inner grid
	OperationMatrix* OpLTwoInner;
	/// Pointer to the coefficients of operation Delta
	//DataVector* deltaCoef;
	/// Pointer to the coefficients ot operation Gamma
	//DataMatrix* gammaCoef;
	/// use coarsening between timesteps in order to reduce gridsize
	bool useCoarsen;
	/// Threshold used to decide if a grid point should be deleted
	double coarsenThreshold;
	/// Percent how many of the removable points should be tested for deletion
	double coarsenPercent;
	/// denotes the number of complete coarsen procedures per timestep
	size_t numExecCoarsen;

	virtual void applyLOperatorInner(DataVector& alpha, DataVector& result);

	virtual void applyLOperatorComplete(DataVector& alpha, DataVector& result);

	virtual void applyMassMatrixInner(DataVector& alpha, DataVector& result);

	virtual void applyMassMatrixComplete(DataVector& alpha, DataVector& result);

	/**
	 * Build the coefficients for the Gamma Operation, which
	 * are the assets' covariance matrix multiplied by 0.5
	 *
	 * this routine handles also the symmtrie of the
	 * gamma operation
	 */
	//void buildGammaCoefficients();

	/**
	 * Build the coefficients for the combined Delta Operation
	 */
	//void buildDeltaCoefficients();

	/**
	 * Build the coefficients for the Gamma Operation, which
	 * are the assets' covariance matrix multiplied by 0.5
	 *
	 * this routine handles also the symmtrie of the
	 * gamma operation
	 *
	 * This function builds the coefficients for the Log Transformed Black Scholes Equation
	 */
	//void buildGammaCoefficientsLogTransform();

	/**
	 * Build the coefficients for the combined Delta Operation
	 *
	 * This function builds the coefficients for the Log Transformed Black Scholes Equation
	 */
	//void buildDeltaCoefficientsLogTransform();

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param alpha the ansatzfunctions' coefficients
	 * @param mu reference to the mus
	 * @param sigma reference to the sigmas
	 * @param rho reference to the rhos
	 * @param r the riskfree interest rate
	 * @param TimestepSize the size of one timestep used in the ODE Solver
	 * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
	 *  							ImEul for implicit Euler, CrNic for Crank Nicolson solver
	 * @param bLogTransform indicates that this system belongs to a log-transformed Black Scholes Equation
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

#endif /* BLACKSCHOLESODESOLVERSYSTEM_HPP */
