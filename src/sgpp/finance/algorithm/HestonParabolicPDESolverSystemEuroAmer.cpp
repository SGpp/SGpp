/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "finance/algorithm/HestonParabolicPDESolverSystemEuroAmer.hpp"
#include "base/exception/algorithm_exception.hpp"
#include "base/grid/generation/functors/SurplusCoarseningFunctor.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "pde/operation/PdeOpFactory.hpp"
#include "finance/operation/FinanceOpFactory.hpp"
#include <cmath>

//#define NOBOUNDARYDISCOUNT

namespace sg
{
namespace finance
{

/*
 * Constructor. Checks for the sanity of the supplied parameters and initialises the required members.
 */
HestonParabolicPDESolverSystemEuroAmer::HestonParabolicPDESolverSystemEuroAmer(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, sg::base::DataVector& thetas, sg::base::DataVector& volvols,
		sg::base::DataVector& kappas,
		sg::base::DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
		double dStrike, std::string option_type,
		bool bLogTransform, bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
		int numCoarsenPoints, double refineThreshold, std::string refineMode, size_t refineMaxLevel)
{
	this->BoundGrid = &SparseGrid;
	this->alpha_complete = &alpha;

	this->alpha_complete_old = new sg::base::DataVector(*this->alpha_complete);
	this->alpha_complete_tmp = new sg::base::DataVector(*this->alpha_complete);
	this->oldGridStorage = new sg::base::GridStorage(*(this->BoundGrid)->getStorage());
	this->secondGridStorage = new sg::base::GridStorage(*(this->BoundGrid)->getStorage());

	this->InnerGrid = NULL;
	this->alpha_inner = NULL;
	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->TimestepSize_old = TimestepSize;
	this->BoundaryUpdate = new sg::base::DirichletUpdateVector(SparseGrid.getStorage());
	this->GridConverter = new sg::base::DirichletGridConverter();
	this->r = r;
	this->thetas = &thetas;
	this->volvols = &volvols;
	this->kappas = &kappas;
	this->hMatrix = &rho;
	this->HestonAlgoDims = this->BoundGrid->getAlgorithmicDimensions();
	this->nAssets = this->HestonAlgoDims.size() / 2;
	this->nExecTimesteps = 0;

	// throw exception if algorithmic dimensions isn't even
	if ((this->HestonAlgoDims.size()%2) != 0 )
	{
		throw sg::base::algorithm_exception("HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : Number of algorithmic dimensions is not even.");
	}

	// throw exception if grid dimensions not equal algorithmic dimensions
	if (this->HestonAlgoDims.size() != this->BoundGrid->getStorage()->dim())
	{
		throw sg::base::algorithm_exception("HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : Number of algorithmic dimensions is not equal to the number of grid's dimensions.");
	}

	// Test if 2*dimC = dimG, where dimC is the number of dimensions in the coefficient vectors and dimG is the number of grid dimensions
	if (this->BoundGrid->getStorage()->dim() != (2*this->thetas->getSize()) || this->BoundGrid->getStorage()->dim() != (2*this->kappas->getSize()) || this->BoundGrid->getStorage()->dim() != (2*this->volvols->getSize()))
	{
		throw sg::base::algorithm_exception("HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : Dimension of theta/volvol/kappa parameters != half of grid's dimensions!");
	}

	// Test if number of dimensions in the coefficients match the numbers of grid dimensions (hmatrix)
	if (this->BoundGrid->getStorage()->dim() != this->hMatrix->getNrows() || this->BoundGrid->getStorage()->dim() != this->hMatrix->getNcols())
	{
		throw sg::base::algorithm_exception("HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : Row or col of hmatrix parameter don't match the grid's dimensions!");
	}

	// Test if all algorithmic dimensions are inside the grid's dimensions
	for (size_t i = 0; i < this->HestonAlgoDims.size(); i++)
	{
		if (this->HestonAlgoDims[i] >= this->BoundGrid->getStorage()->dim())
		{
			throw sg::base::algorithm_exception("HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : Minimum one algorithmic dimension is not inside the grid's dimensions!");
		}
	}

	// Test if there are double algorithmic dimensions
	std::vector<size_t> tempAlgoDims(this->HestonAlgoDims);
	for (size_t i = 0; i < this->HestonAlgoDims.size(); i++)
	{
		size_t dimCount = 0;
		for (size_t j = 0; j < tempAlgoDims.size(); j++)
		{
			if (this->HestonAlgoDims[i] == tempAlgoDims[j])
			{
				dimCount++;
			}
		}

		if (dimCount > 1)
		{
			throw sg::base::algorithm_exception("HestonParabolicPDESolverSystemEuropean::HestonParabolicPDESolverSystemEuropean : There is minimum one doubled algorithmic dimension!");
		}
	}

	// Build the coefficient collections for the operations
	// Some of the coefficient collections are vectors, some are matrices and one is a tensor (for the Up/Down Four Op Dim case)
	size_t coefficientVectorSize = this->HestonAlgoDims.size();
	this->dCoeff = new sg::base::DataVector(coefficientVectorSize);
	this->eCoeff = new sg::base::DataVector(coefficientVectorSize);
	this->fCoeff = new sg::base::DataVector(coefficientVectorSize);
	this->gCoeff = new sg::base::DataVector(coefficientVectorSize);
	this->zCoeff = new sg::base::DataVector(coefficientVectorSize);
	this->bCoeff = new sg::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
	this->cCoeff = new sg::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
	this->hCoeff = new sg::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
	this->xCoeff = new sg::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
	this->yCoeff = new sg::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
	this->wCoeff = new sg::base::DataMatrix(coefficientVectorSize, coefficientVectorSize);
	create4dEqualDimSizeArray(coefficientVectorSize, &(this->kCoeff));

	// create the inner grid
	// Here we have to make sure that the alpha_inner doesn't come out with zeros at the other end.
	this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

	if (bLogTransform == false)
	{
		// Build coefficients
		buildXCoefficients();
		buildYCoefficients();
		buildWCoefficients();
		buildZCoefficients();
		buildGCoefficients();
		buildFCoefficients();
		buildDCoefficients();

		// Create operators
		this->OpXBound = sg::op_factory::createOperationHestonX(*this->BoundGrid, *this->xCoeff);
		this->OpXInner = sg::op_factory::createOperationHestonX(*this->InnerGrid, *this->xCoeff);
		this->OpYBound = sg::op_factory::createOperationHestonY(*this->BoundGrid, *this->yCoeff);
		this->OpYInner = sg::op_factory::createOperationHestonY(*this->InnerGrid, *this->yCoeff);
		this->OpWBound = sg::op_factory::createOperationHestonW(*this->BoundGrid, *this->wCoeff);
		this->OpWInner = sg::op_factory::createOperationHestonW(*this->InnerGrid, *this->wCoeff);
		this->OpZBound = sg::op_factory::createOperationHestonZ(*this->BoundGrid, *this->zCoeff);
		this->OpZInner = sg::op_factory::createOperationHestonZ(*this->InnerGrid, *this->zCoeff);
		this->OpGBound = sg::op_factory::createOperationHestonGLog(*this->BoundGrid, *this->gCoeff);
		this->OpGInner = sg::op_factory::createOperationHestonGLog(*this->InnerGrid, *this->gCoeff);
		this->OpDBound = sg::op_factory::createOperationHestonDLog(*this->BoundGrid, *this->dCoeff);
		this->OpDInner = sg::op_factory::createOperationHestonDLog(*this->InnerGrid, *this->dCoeff);
		this->OpFBound = sg::op_factory::createOperationHestonFLog(*this->BoundGrid, *this->fCoeff);
		this->OpFInner = sg::op_factory::createOperationHestonFLog(*this->InnerGrid, *this->fCoeff);
	}
	else
	{
		// Log-transformed case
		// coefficients
		buildBCoefficientsLogTransform();
		buildCCoefficientsLogTransform();
		buildDCoefficientsLogTransform();
		buildECoefficientsLogTransform();
		buildFCoefficientsLogTransform();
		buildGCoefficientsLogTransform();
		buildHCoefficientsLogTransform();
		buildKCoefficientsLogTransform();

		// operations
		this->OpBBound = sg::op_factory::createOperationHestonBLog(*this->BoundGrid, *this->bCoeff);
		this->OpBInner = sg::op_factory::createOperationHestonBLog(*this->InnerGrid, *this->bCoeff);
		this->OpCBound = sg::op_factory::createOperationHestonCLog(*this->BoundGrid, *this->cCoeff);
		this->OpCInner = sg::op_factory::createOperationHestonCLog(*this->InnerGrid, *this->cCoeff);
		this->OpDBound = sg::op_factory::createOperationHestonDLog(*this->BoundGrid, *this->dCoeff);
		this->OpDInner = sg::op_factory::createOperationHestonDLog(*this->InnerGrid, *this->dCoeff);
		this->OpEBound = sg::op_factory::createOperationHestonELog(*this->BoundGrid, *this->eCoeff);
		this->OpEInner = sg::op_factory::createOperationHestonELog(*this->InnerGrid, *this->eCoeff);
		this->OpFBound = sg::op_factory::createOperationHestonFLog(*this->BoundGrid, *this->fCoeff);
		this->OpFInner = sg::op_factory::createOperationHestonFLog(*this->InnerGrid, *this->fCoeff);
		this->OpGBound = sg::op_factory::createOperationHestonGLog(*this->BoundGrid, *this->gCoeff);
		this->OpGInner = sg::op_factory::createOperationHestonGLog(*this->InnerGrid, *this->gCoeff);
		this->OpHBound = sg::op_factory::createOperationHestonHLog(*this->BoundGrid, *this->hCoeff);
		this->OpHInner = sg::op_factory::createOperationHestonHLog(*this->InnerGrid, *this->hCoeff);
		this->OpKBound = sg::op_factory::createOperationHestonKLog(*this->BoundGrid, &(this->kCoeff));
		this->OpKInner = sg::op_factory::createOperationHestonKLog(*this->InnerGrid, &(this->kCoeff));
	}

	// Create operations that are needed by both the log-transformed and Cartesian case
	this->OpAInner = sg::op_factory::createOperationLTwoDotProduct(*this->InnerGrid);
	this->OpABound = sg::op_factory::createOperationLTwoDotProduct(*this->BoundGrid);

	// right hand side if System
	this->rhs = NULL;

	// set coarsen settings
	this->useCoarsen = useCoarsen;
	this->coarsenThreshold = coarsenThreshold;
	this->refineThreshold = refineThreshold;
	this->adaptSolveMode = adaptSolveMode;
	this->numCoarsenPoints = numCoarsenPoints;
	this->refineMode = refineMode;
	this->refineMaxLevel = refineMaxLevel;

	// init Number of AverageGridPoins
	this->numSumGridpointsInner = 0;
	this->numSumGridpointsComplete = 0;

	// init option type and strike
	this->dStrike = dStrike;
	this->option_type = option_type;

	// save coordinate transformations
	this->b_log_transform = bLogTransform;
}

HestonParabolicPDESolverSystemEuroAmer::~HestonParabolicPDESolverSystemEuroAmer()
{
	// Todo: update

	//	delete this->OpDeltaBound;
	//	delete this->OpGammaBound;
	//	delete this->OpLTwoBound;
	//	delete this->OpDeltaInner;
	//	delete this->OpGammaInner;
	//	delete this->OpLTwoInner;
	//	delete this->gammaCoef;
	//	delete this->deltaCoef;
	//	delete this->BoundaryUpdate;
	//	delete this->GridConverter;
	//	if (this->InnerGrid != NULL)
	//	{
	//		delete this->InnerGrid;
	//	}
	//	if (this->alpha_inner != NULL)
	//	{
	//		delete this->alpha_inner;
	//	}
	//	if (this->rhs != NULL)
	//	{
	//		delete this->rhs;
	//	}
	//	delete this->alpha_complete_old;
	//	delete this->alpha_complete_tmp;
	delete4dEqualDimSizeArray(this->HestonAlgoDims.size(), &(this->kCoeff));

	//	std::cout << "test1" << std::endl;
	//
	//	for(size_t i=0;i<this->HestonAlgoDims.size();i++)
	//	{
	//		for(size_t j=0;j<this->HestonAlgoDims.size();j++)
	//		{
	//			for(size_t k=0;k<this->HestonAlgoDims.size();k++)
	//			{
	//				std::cout << "test2" << std::endl;
	//				delete[] this->kCoeff[i][j][k];
	//			}
	//			std::cout << "test3" << std::endl;
	//			delete[] this->kCoeff[i][j];
	//			std::cout << "test4" << std::endl;
	//		}
	//		std::cout << "test5" << std::endl;
	//		delete[] this->kCoeff[i];
	//		std::cout << "test6" << std::endl;
	//	}
	//	std::cout << "test7" << std::endl;
	//	delete[] this->kCoeff;
	//	std::cout << "test8" << std::endl;

}

void HestonParabolicPDESolverSystemEuroAmer::applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the riskfree rate
	if (this->r != 0.0)
	{
		this->OpABound->mult(alpha, temp);
		result.axpy((-1.0)*this->r, temp);
	}

	if(this->b_log_transform)
	{
		// Log-transformed coordinates

		// Apply the B operator
		this->OpBBound->mult(alpha, temp);
		result.sub(temp);

		// Apply the C operator
		this->OpCBound->mult(alpha, temp);
		result.sub(temp);

		// Apply the D operator
		this->OpDBound->mult(alpha, temp);
		result.sub(temp);

		// Apply the E operator
		this->OpEBound->mult(alpha, temp);
		result.add(temp);

		// Apply the F operator
		this->OpFBound->mult(alpha, temp);
		result.add(temp);

		// Apply the G operator
		this->OpGBound->mult(alpha, temp);
		result.sub(temp);

		// Apply the H operator
		this->OpHBound->mult(alpha, temp);
		result.sub(temp);

		// Apply the K operator
		this->OpKBound->mult(alpha, temp);
		result.sub(temp);
	}
	else
	{
		// Cartesian coordinates

		// Apply the X method
		this->OpXBound->mult(alpha, temp);
		result.sub(temp);

		// Apply the Y method
		this->OpYBound->mult(alpha, temp);
		result.sub(temp);

		this->OpGBound->mult(alpha, temp);
		result.sub(temp);

		this->OpDBound->mult(alpha, temp);
		result.sub(temp);

		this->OpFBound->mult(alpha, temp);
		result.add(temp);

		this->OpWBound->mult(alpha, temp);
		result.sub(temp);

		this->OpZBound->mult(alpha, temp);
		result.add(temp);
	}
}

void HestonParabolicPDESolverSystemEuroAmer::applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the riskfree rate
	if (this->r != 0.0)
	{
		this->OpAInner->mult(alpha, temp);
		result.axpy((-1.0)*this->r, temp);
	}

	if(this->b_log_transform)
	{
		// Log-transformed coordinates

		// Apply the B method
		this->OpBInner->mult(alpha, temp);
		result.sub(temp);

		// Apply the C method
		this->OpCInner->mult(alpha, temp);
		result.sub(temp);

		// Apply the D method
		this->OpDInner->mult(alpha, temp);
		result.sub(temp);

		// Apply the E method
		this->OpEInner->mult(alpha, temp);
		result.add(temp);

		// Apply the F method
		this->OpFInner->mult(alpha, temp);
		result.add(temp);

		// Apply the G method
		this->OpGInner->mult(alpha, temp);
		result.sub(temp);

		// Apply the H method
		this->OpHInner->mult(alpha, temp);
		result.sub(temp);

		// Apply the K method
		this->OpKInner->mult(alpha, temp);
		result.sub(temp);
	}
	else
	{
		// Cartesian coordinates

		// Apply the X method
		this->OpXInner->mult(alpha, temp);
		result.sub(temp);

		// Apply the Y method
		this->OpYInner->mult(alpha, temp);
		result.sub(temp);

		this->OpGInner->mult(alpha, temp);
		result.sub(temp);

		this->OpDInner->mult(alpha, temp);
		result.sub(temp);

		this->OpFInner->mult(alpha, temp);
		result.add(temp);

		this->OpWInner->mult(alpha, temp);
		result.sub(temp);

		this->OpZInner->mult(alpha, temp);
		result.add(temp);
	}
}

void HestonParabolicPDESolverSystemEuroAmer::applyMassMatrixComplete(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpABound->mult(alpha, temp);

	result.add(temp);
}

void HestonParabolicPDESolverSystemEuroAmer::applyMassMatrixInner(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpAInner->mult(alpha, temp);

	result.add(temp);
}

/*
 * Post-processing following the main timestep logic.
 * Updates the boundary coefficients and performs discounting.
 */
void HestonParabolicPDESolverSystemEuroAmer::finishTimestep()
{
	this->nExecTimesteps++;

	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

#ifndef NOBOUNDARYDISCOUNT
	// Adjust the boundaries with the riskfree rate
	if (this->r != 0.0)
	{
		if (this->tOperationMode == "ExEul" || this->tOperationMode == "AdBas")
		{
			this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete, exp(((-1.0)*(this->r*this->TimestepSize))));
		}
	}
#endif
}

/*
 * Coarsens and refines the sparse grid based on the relevant thresholds.
 * @param isLastTimestep: True if the current timestep is the last. If this is true, integrated coarsening is not performed.
 */
void HestonParabolicPDESolverSystemEuroAmer::coarsenAndRefine(bool isLastTimestep)
{
	// add number of Gridpoints
	this->numSumGridpointsInner += this->InnerGrid->getSize();
	this->numSumGridpointsComplete += this->BoundGrid->getSize();

	if (this->useCoarsen == true && isLastTimestep == false)
	{
		///////////////////////////////////////////////////
		// Start integrated refinement & coarsening
		///////////////////////////////////////////////////

		size_t originalGridSize = this->BoundGrid->getStorage()->size();

		// Coarsen the grid
		sg::base::GridGenerator* myGenerator = this->BoundGrid->createGridGenerator();

		if (this->adaptSolveMode == "refine" || this->adaptSolveMode == "coarsenNrefine")
		{
			size_t numRefines = myGenerator->getNumberOfRefinablePoints();
			sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(this->alpha_complete, numRefines, this->refineThreshold);
			if (this->refineMode == "maxLevel")
			{
				myGenerator->refineMaxLevel(myRefineFunc, this->refineMaxLevel);
				this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
			}
			if (this->refineMode == "classic")
			{
				myGenerator->refine(myRefineFunc);
				this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
			}
			delete myRefineFunc;
		}

		if (this->adaptSolveMode == "coarsen" || this->adaptSolveMode == "coarsenNrefine")
		{
			size_t numCoarsen = myGenerator->getNumberOfRemovablePoints();
			sg::base::SurplusCoarseningFunctor* myCoarsenFunctor = new sg::base::SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, this->coarsenThreshold);
			myGenerator->coarsenNFirstOnly(myCoarsenFunctor, this->alpha_complete, originalGridSize);
			delete myCoarsenFunctor;
		}

		delete myGenerator;

		///////////////////////////////////////////////////
		// End integrated refinement & coarsening
		///////////////////////////////////////////////////

		// rebuild the inner grid + coefficients
		this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);
	}
}

/*
 * Pre-processing operation before the main timestep operations begin. This method adjusts the boundaries by discounting them.
 */
void HestonParabolicPDESolverSystemEuroAmer::startTimestep()
{
#ifndef NOBOUNDARYDISCOUNT
	// Adjust the boundaries with the riskfree rate
	if (this->r != 0.0)
	{
		if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul")
		{

			sg::base::GridStorage* storage = this->BoundGrid->getStorage();

			for (size_t i = 0; i < storage->size(); i++)
			{
				sg::base::GridIndex* curPoint = (*storage)[i];

				sg::base::DataVector pointCoords(storage->dim());
				curPoint->getCoords(pointCoords);

				// Discount the boundary points
				if(!curPoint->isInnerPoint())
				{
					this->alpha_complete->set(i, this->alpha_complete->get(i)*exp(((-1.0)*(this->r*this->TimestepSize))));
				}
			}
		}
	}
#endif
}

/*
 * Builds the coefficient object for the G operator.
 * Note that only vanilla options are supported for this (non-transformed) case.
 */
void HestonParabolicPDESolverSystemEuroAmer::buildGCoefficients()
{
	double volvol = volvols->get(0);
	double rho = this->hMatrix->get(0,1);
	double kappa = this->kappas->get(0);

	this->gCoeff->setAll(0.0);
	this->gCoeff->set(1, (kappa) + rho*volvol);
}

/*
 * Builds the coefficient object for the (non-log-transformed) F operator.
 * Note that only vanilla options are supported for this (non-transformed) case.
 */
void HestonParabolicPDESolverSystemEuroAmer::buildFCoefficients()
{
	double volvol = volvols->get(0);
	double kappa = this->kappas->get(0);
	double theta = this->thetas->get(0);

	this->fCoeff->setAll(0.0);
	this->fCoeff->set(1, kappa*theta - 0.5*pow(volvol,2.0));
}

/*
 * Builds the coefficient object for the (non-log-transformed) D operator.
 * Note that only vanilla options are supported for this (non-transformed) case.
 */
void HestonParabolicPDESolverSystemEuroAmer::buildDCoefficients()
{
	double volvol = volvols->get(0);

	this->dCoeff->setAll(0.0);
	this->dCoeff->set(1, (0.5)*pow(volvol,2.0));
}

/*
 * Builds the coefficient object for the X operator.
 * Note that only vanilla options are supported for this (non-transformed) case.
 */
void HestonParabolicPDESolverSystemEuroAmer::buildXCoefficients()
{
	this->xCoeff->setAll(0.0);
	this->xCoeff->set(0,1,1.0);
}

/*
 * Builds the coefficient object for the W operator.
 * Note that only vanilla options are supported for this (non-transformed) case.
 */
void HestonParabolicPDESolverSystemEuroAmer::buildWCoefficients()
{
	double volvol = volvols->get(0);
	double rho = this->hMatrix->get(0,1);

	this->wCoeff->setAll(0.0);
	this->wCoeff->set(0,1, rho*volvol);
}

/*
 * Builds the coefficient object for the Z operator.
 * Note that only vanilla options are supported for this (non-transformed) case.
 */
void HestonParabolicPDESolverSystemEuroAmer::buildZCoefficients()
{
	this->zCoeff->setAll(0.0);
	this->zCoeff->set(0,this->r);
}

/*
 * Builds the coefficient object for the Y operator.
 * Note that only vanilla options are supported for this (non-transformed) case.
 */
void HestonParabolicPDESolverSystemEuroAmer::buildYCoefficients()
{
	this->yCoeff->setAll(0.0);
	this->yCoeff->set(0,1,0.5);
}

/*
 * Builds the coefficients object for the B operator.
 * Operator B has two custom 1D operators, so the coefficients object is a matrix.
 * The operator has a non-zero coefficient for each pairing of a stock and its variance, e.g. elements (1,2), (3,4), (5,6) etc. (zero-based array elements (0,1), (2,3), (4,5) etc.)
 */
void HestonParabolicPDESolverSystemEuroAmer::buildBCoefficientsLogTransform()
{
	this->bCoeff->setAll(0.0);
	for (size_t i = 0; i < this->nAssets; i++)
	{
		this->bCoeff->set(2*i,2*i+1, 0.5);	// Coefficient is a constant of 0.5
	}
}
/*
 * Builds the coefficients object for the C operator.
 * Operator C has two custom 1D operators, so the coefficients object is a matrix.
 * The operator has a non-zero coefficient for each pairing of a stock and its variance, e.g. elements (1,2), (3,4), (5,6) etc. (zero-based array elements (0,1), (2,3), (4,5) etc.)
 */
void HestonParabolicPDESolverSystemEuroAmer::buildCCoefficientsLogTransform()
{
	double volvol, rho;
	this->cCoeff->setAll(0.0);
	for (size_t i = 0; i < this->nAssets; i++)
	{
		// Determine and set the coefficient for this asset
		volvol = volvols->get(i);
		rho = this->hMatrix->get(2*i,2*i+1);
		this->cCoeff->set(2*i, 2*i+1,volvol*rho);
	}
}

/*
 * Builds the coefficients object for the D operator.
 * The D operator is oneOpDim, so there is only a vector of coefficients to set.
 * Only the variance dimensions are active in this operator, i.e. elements 2, 4, 6, 8 etc. (zero-based array elements 1, 3, 5, 7 etc.)
 */
void HestonParabolicPDESolverSystemEuroAmer::buildDCoefficientsLogTransform()
{
	this->dCoeff->setAll(0.0);
	for(size_t i = 0; i < this->nAssets;i++)
	{
		// Determine and set the coefficient for this asset
		double volvol = volvols->get(i);
		this->dCoeff->set(2*i+1, (0.5)*pow(volvol,2.0));
	}
}

/*
 * Builds the coefficients object for the E operator.
 * The E operator is oneOpDim, so there is only a vector of coefficients to set.
 * Only the stock-price dimensions are active in this operator, i.e. elements 1, 3, 5, 7 etc. (zero-based array elements 0, 2, 4, 6 etc.)
 */
void HestonParabolicPDESolverSystemEuroAmer::buildECoefficientsLogTransform()
{
	this->eCoeff->setAll(0.0);
	for (size_t i = 0; i < this->nAssets; i++)
	{
		this->eCoeff->set(2*i, this->r);	// Coefficient is independent of which asset it is.
	}
}

/*
 * Builds the coefficients object for the F operator.
 * The F operator is oneOpDim, so there is only a vector of coefficients to set.
 * Only the variance elements are active in this operator, i.e. elements 2, 4, 6, 8 etc. (zero-based array elements 1, 3, 5, 7 etc.)
 */
void HestonParabolicPDESolverSystemEuroAmer::buildFCoefficientsLogTransform()
{
	double theta, kappa, volvol;
	this->fCoeff->setAll(0.0);
	for (size_t i = 0; i < nAssets; i++)
	{
		theta = this->thetas->get(i);
		kappa = this->kappas->get(i);
		volvol = this->volvols->get(i);
		this->fCoeff->set(2*i+1, kappa*theta - 0.5*pow(volvol,2.0));
	}
}

/*
 * Builds the coefficients object for the G operator.
 * The G operator is oneOpDim, so there is only a vector of coefficients to set.
 * Only the variance elements are active in this operator, i.e. elements 2, 4, 6, 8 etc. (zero-based vector elements 1, 3, 5, 7 etc.)
 */
void HestonParabolicPDESolverSystemEuroAmer::buildGCoefficientsLogTransform()
{
	this->gCoeff->setAll(0.0);
	double kappa;
	for (size_t i = 0; i < nAssets; i++)
	{
		kappa = this->kappas->get(i);
		this->gCoeff->set(2*i+1, kappa);
	}
}

/*
 * Builds the coefficients object for the H operator.
 * Operator H has two custom 1D operators, so the coefficients object is a matrix.
 * The operator has a non-zero coefficient for each pairing of a stock and its variance, e.g. elements (1,2), (3,4), (5,6) etc. (zero-based array elements (0,1), (2,3), (4,5) etc.)
 */
void HestonParabolicPDESolverSystemEuroAmer::buildHCoefficientsLogTransform()
{
	this->hCoeff->setAll(0.0);
	for (size_t i = 0; i < nAssets; i++)
	{
		this->hCoeff->set(2*i, 2*i+1, 0.5);
	}
}

/*
 * Builds the coefficients object for the K operator.
 * This operator is more exciting. It has four custom 1D operators, so the coefficients object is a 4d array.
 * The operator has a non-zero coefficient for each quad-pairing of a stock and its variance, and a DIFFERENT stock and its variance, e.g. for three assets: (1,2,3,4), (1,2,5,6), (3,4,5,6). (zero-based array elements (0,1,2,3), (0,1,4,5), (2,3,4,5))
 */
void HestonParabolicPDESolverSystemEuroAmer::buildKCoefficientsLogTransform()
{
	double rho;
	setAll4dEqualDimSizeArray(this->HestonAlgoDims.size(), &(this->kCoeff), 0.0);
	for (size_t i = 0; i < nAssets; i++)
	{
		for (size_t j = i+1; j < nAssets; j++)
		{
			rho = this->hMatrix->get(2*i,2*j);
			this->kCoeff[2*i][2*i+1][2*j][2*j+1] = rho;
		}
	}
}

/*
 * Utility method for creating a 4d array of equal size in each dimension. The array is allocated dynamically based on the provided size.
 */
void HestonParabolicPDESolverSystemEuroAmer::create4dEqualDimSizeArray(size_t dimSize, double***** array)
{
	(*array) = (double****) calloc(dimSize, sizeof(double***));
	for(size_t i=0;i<dimSize;i++)
	{
		(*array)[i] = (double***) calloc(dimSize, sizeof(double**));

		for(size_t j=0;j<dimSize;j++)
		{
			(*array)[i][j] = (double**) calloc(dimSize, sizeof(double*));

			for(size_t k=0;k<dimSize;k++)
			{
				(*array)[i][j][k] = (double*) calloc(dimSize, sizeof(double));

				for (size_t m=0;m<dimSize;m++)
				{
					(*array)[i][j][k][m] = 0.0;
				}
			}
		}
	}
}

/*
 * Utility method to free the memory allocated to a 4d array of equal size in each dimension.
 */
void HestonParabolicPDESolverSystemEuroAmer::delete4dEqualDimSizeArray(size_t dimSize, double***** array)
{
	for(size_t i=0;i<dimSize;i++)
	{
		for(size_t j=0;j<dimSize;j++)
		{
			for(size_t k=0;k<dimSize;k++)
			{
				delete[] (*array)[i][j][k];
			}
			delete[] (*array)[i][j];
		}
		delete[] (*array)[i];
	}
	delete[] (*array);
}

/*
 * Utility method to set all values in a 4d array of equal size to one value (e.g. this method can be used to set all the values to zero).
 */
void HestonParabolicPDESolverSystemEuroAmer::setAll4dEqualDimSizeArray(size_t dimSize, double***** array, double value)
{
	for(size_t i=0;i<dimSize;i++)
	{
		for(size_t j=0;j<dimSize;j++)
		{
			for(size_t k=0;k<dimSize;k++)
			{
				for (size_t m=0;m<dimSize;m++)
				{
					(*array)[i][j][k][m] = value;
				}
			}
		}
	}
}

}
}
