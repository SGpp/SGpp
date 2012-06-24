/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Sam Maurus (MA thesis)

#ifndef HESTONPARABOLICPDESOLVERSYSTEMEUROAMER_HPP
#define HESTONPARABOLICPDESOLVERSYSTEMEUROAMER_HPP

#include "base/grid/Grid.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/GridStorage.hpp"
#include "pde/operation/OperationParabolicPDESolverSystemDirichlet.hpp"
#include "finance/tools/Hedging.hpp"

namespace sg
{
namespace finance
{

/**
 * This class implements the ParabolicPDESolverSystem for the BlackScholes
 * Equation.
 *
 * Here European or American Options with fix Dirichlet boundaries are solved.
 */
class HestonParabolicPDESolverSystemEuroAmer : public sg::pde::OperationParabolicPDESolverSystemDirichlet
{
protected:

	/// the riskfree interest rate
	double r;

	// The various Heston operators
	sg::base::OperationMatrix* OpABound;
	sg::base::OperationMatrix* OpAInner;
	sg::base::OperationMatrix* OpBBound;
	sg::base::OperationMatrix* OpBInner;
	sg::base::OperationMatrix* OpCBound;
	sg::base::OperationMatrix* OpCInner;
	sg::base::OperationMatrix* OpDBound;
	sg::base::OperationMatrix* OpDInner;
	sg::base::OperationMatrix* OpEBound;
	sg::base::OperationMatrix* OpEInner;
	sg::base::OperationMatrix* OpFBound;
	sg::base::OperationMatrix* OpFInner;
	sg::base::OperationMatrix* OpGBound;
	sg::base::OperationMatrix* OpGInner;
	sg::base::OperationMatrix* OpHBound;
	sg::base::OperationMatrix* OpHInner;
	sg::base::OperationMatrix* OpKBound;
	sg::base::OperationMatrix* OpKInner;
	sg::base::OperationMatrix* OpXBound;
	sg::base::OperationMatrix* OpXInner;
	sg::base::OperationMatrix* OpYBound;
	sg::base::OperationMatrix* OpYInner;
	sg::base::OperationMatrix* OpWBound;
	sg::base::OperationMatrix* OpWInner;
	sg::base::OperationMatrix* OpZBound;
	sg::base::OperationMatrix* OpZInner;
	sg::base::OperationMatrix* OpLTwoBound;
	sg::base::OperationMatrix* OpLTwoInner;

	// Pointer to the vector containing the volatility of volatility values
	sg::base::DataVector* volvols;
	sg::base::DataVector* kappas;
	sg::base::DataVector* thetas;

	/// Pointer to the rhos
	sg::base::DataMatrix* hMatrix;

	//	Coefficient collections for the operators
	// Up/down one op-dim
	sg::base::DataVector* dCoeff;
	sg::base::DataVector* eCoeff;
	sg::base::DataVector* fCoeff;
	sg::base::DataVector* gCoeff;
	sg::base::DataVector* zCoeff;

	// Up/down two op-dims
	sg::base::DataMatrix* bCoeff;
	sg::base::DataMatrix* cCoeff;
	sg::base::DataMatrix* hCoeff;
	sg::base::DataMatrix* xCoeff;
	sg::base::DataMatrix* yCoeff;
	sg::base::DataMatrix* wCoeff;

	// Up/down four op dims
	double**** kCoeff;

	// Use coarsening between timesteps in order to reduce gridsize
	bool useCoarsen;
	// Adaptive mode during solving Heston Equation: coarsen, refine, coarsenNrefine
	std::string adaptSolveMode;
	// Number of points the are coarsened in each coarsening-step
	int numCoarsenPoints;
	// Threshold used to decide if a grid point should be deleted
	double coarsenThreshold;
	// Threshold used to decide if a grid point should be refined
	double refineThreshold;
	// Refine mode during solving the Heston Equation: classic or maxLevel
	std::string refineMode;
	// MaxLevel max. Level of refinement
	size_t refineMaxLevel;
	// The algorithmic dimensions used in this system
	std::vector<size_t> HestonAlgoDims;
	// The number of assets (half the number of problem dimensions)
	size_t nAssets;
	// Stores the number of executed timesteps
	size_t nExecTimesteps;
	// The strike of the current option
	double dStrike;
	// The type of the current option
	std::string option_type;
	// store whether log coordinates are used
	bool b_log_transform;

	//
	virtual void applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result);
	virtual void applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result);
	virtual void applyMassMatrixInner(sg::base::DataVector& alpha, sg::base::DataVector& result);
	virtual void applyMassMatrixComplete(sg::base::DataVector& alpha, sg::base::DataVector& result);

	//
	void buildDCoefficients();
	void buildFCoefficients();
	void buildGCoefficients();
	void buildXCoefficients();
	void buildYCoefficients();
	void buildWCoefficients();
	void buildZCoefficients();

	void buildACoefficientsLogTransform();
	void buildBCoefficientsLogTransform();
	void buildCCoefficientsLogTransform();
	void buildDCoefficientsLogTransform();
	void buildECoefficientsLogTransform();
	void buildFCoefficientsLogTransform();
	void buildGCoefficientsLogTransform();
	void buildHCoefficientsLogTransform();
	void buildKCoefficientsLogTransform();

	void create4dEqualDimSizeArray(size_t dimSize, double***** array);
	void delete4dEqualDimSizeArray(size_t dimSize, double***** array);
	void setAll4dEqualDimSizeArray(size_t dimSize, double***** array, double value);

public:
	/**
	 * Std-Constructor
	 *
	 * Todo: comment
	 */
	HestonParabolicPDESolverSystemEuroAmer(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, sg::base::DataVector& thetas, sg::base::DataVector& volvols,
			sg::base::DataVector& kappas,
			sg::base::DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
			double dStrike, std::string option_type,
			bool bLogTransform = false, bool useCoarsen = false, double coarsenThreshold = 0.0, std::string adaptSolveMode ="none",
			int numCoarsenPoints = -1, double refineThreshold = 0.0, std::string refineMode = "classic", size_t refineMaxLevel = 0);

	/**
	 * Std-Destructor
	 */
	virtual ~HestonParabolicPDESolverSystemEuroAmer();

	virtual void finishTimestep();

	virtual void coarsenAndRefine(bool isLastTimestep = false);

	void startTimestep();
};

}

}

#endif /* HESTONPARABOLICPDESOLVERSYSTEMEUROAMER_HPP */
