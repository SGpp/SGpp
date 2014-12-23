/******************************************************************************
 * Copyright (C) 2012 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/ModLinearGrid.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "base/operation/OperationMultipleEval.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/exception/application_exception.hpp"
#include "base/tools/PrecisionConverter.hpp"
#include "base/tools/GridPrinter.hpp"

#include "solver/sle/ConjugateGradientsSP.hpp"
#include "solver/sle/BiCGStabSP.hpp"

#include "datadriven/application/LearnerBaseSP.hpp"

#include <iostream>

#include "parallel/tools/MPI/SGppMPITools.hpp"

namespace sg {

namespace datadriven {

LearnerBaseSP::LearnerBaseSP(const bool isRegression, const bool isVerbose) :
		alpha_(NULL), grid_(NULL), isVerbose_(isVerbose), isRegression_(isRegression), isTrained_(false), execTime_(
				0.0), GFlop_(0.0), GByte_(0.0) {
#ifdef USE_MPI

	// suppress output from all process but proc0,
	// output is (in the normal, correctly working
	// case) the same for all MPI processes, so no
	// need to see output more than once
	if (sg::parallel::myGlobalMPIComm->getMyRank() != 0) {
		this->isVerbose_ = false;
	}

#endif
}

LearnerBaseSP::LearnerBaseSP(const std::string tGridFilename, const std::string tAlphaFilename, const bool isRegression,
		const bool isVerbose) :
		alpha_(NULL), grid_(NULL), isVerbose_(isVerbose), isRegression_(isRegression), isTrained_(false), execTime_(
				0.0), GFlop_(0.0), GByte_(0.0) {
	// @TODO (heinecke) implement
	throw base::application_exception("LearnerBaseSP::LearnerBaseSP: This construct isn't implemented, yet!");
}

bool isVerbose_;
/// is regression selected
bool isRegression_;
/// is the grid trained
bool isTrained_;
/// execution time
double execTime_;
/// number of executed Floating Point operations
double GFlop_;
/// number of transferred Gbytes
double GByte_;

LearnerBaseSP::LearnerBaseSP(const LearnerBaseSP& copyMe) :
		isVerbose_(copyMe.isVerbose_), isRegression_(copyMe.isRegression_), isTrained_(false), execTime_(0.0), GFlop_(
				0.0), GByte_(0.0) {
	this->isRegression_ = copyMe.isRegression_;
	this->isTrained_ = false;
	this->GFlop_ = 0.0;
	this->GByte_ = 0.0;
	this->execTime_ = 0.0;

	// safety, should not happen
	if (alpha_ != NULL)
		delete alpha_;

	if (grid_ != NULL)
		delete grid_;

	// @TODO (heinecke) grid copy constructor
	grid_ = sg::base::Grid::unserialize(copyMe.grid_->serialize());
	alpha_ = new sg::base::DataVectorSP(*(copyMe.alpha_));
}

LearnerBaseSP::~LearnerBaseSP() {
	// if user does no cleaning
	if (alpha_ != NULL)
		delete alpha_;

	if (grid_ != NULL)
		delete grid_;
}

void LearnerBaseSP::InitializeGrid(const sg::base::RegularGridConfiguration& GridConfig) {
	if (GridConfig.type_ == sg::base::LinearTrapezoidBoundary) {
		grid_ = new sg::base::LinearTrapezoidBoundaryGrid(GridConfig.dim_);
	} else if (GridConfig.type_ == sg::base::ModLinear) {
		grid_ = new sg::base::ModLinearGrid(GridConfig.dim_);
	} else if (GridConfig.type_ == sg::base::Linear) {
		grid_ = new sg::base::LinearGrid(GridConfig.dim_);
	} else {
		grid_ = NULL;
		throw base::application_exception("LearnerBaseSP::InitializeGrid: An unsupported grid type was chosen!");
	}

	// Generate regular Grid with LEVELS Levels
	sg::base::GridGenerator* myGenerator = grid_->createGridGenerator();
	myGenerator->regular(GridConfig.level_);
	delete myGenerator;

	// Create alpha
	alpha_ = new sg::base::DataVectorSP(grid_->getSize());
	alpha_->setAll(0.0);
}

void LearnerBaseSP::preProcessing() {
}

void LearnerBaseSP::postProcessing(const sg::base::DataMatrixSP& trainDataset, const sg::solver::SLESolverType& solver,
		const size_t numNeededIterations) {
	if (this->isVerbose_) {
		std::cout << std::endl;
		std::cout << "Current Execution Time: " << execTime_ << std::endl;
		std::cout << std::endl;
	}
}

LearnerTiming LearnerBaseSP::train(sg::base::DataMatrixSP& trainDataset, sg::base::DataVectorSP& classes,
		const sg::base::RegularGridConfiguration& GridConfig,
		const sg::solver::SLESolverSPConfiguration& SolverConfigRefine,
		const sg::solver::SLESolverSPConfiguration& SolverConfigFinal,
		const sg::base::AdpativityConfiguration& AdaptConfig, const bool testAccDuringAdapt, const float lambda) {
	LearnerTiming result;

	if (trainDataset.getNrows() != classes.getSize()) {
		throw base::application_exception("LearnerBaseSP::train: length of classes vector does not match to dataset!");
	}

	result.timeComplete_ = 0.0;
	result.timeMultComplete_ = 0.0;
	result.timeMultCompute_ = 0.0;
	result.timeMultTransComplete_ = 0.0;
	result.timeMultTransCompute_ = 0.0;
	result.timeRegularization_ = 0.0;
	result.GFlop_ = 0.0;
	result.GByte_ = 0.0;

	execTime_ = 0.0;
	GFlop_ = 0.0;
	GByte_ = 0.0;

	double oldAcc = 0.0;

	// Construct Grid
	if (alpha_ != NULL)
		delete alpha_;

	if (grid_ != NULL)
		delete grid_;

	if (isTrained_ == true)
		isTrained_ = false;

	InitializeGrid(GridConfig);

	// check if grid was created
	if (grid_ == NULL)
		return result;

	// create DMSystem
	sg::datadriven::DMSystemMatrixBaseSP* DMSystem = createDMSystem(trainDataset, lambda);

	// check if System was created
	if (DMSystem == NULL)
		return result;

	sg::solver::SLESolverSP* myCG;

	if (SolverConfigRefine.type_ == sg::solver::CG) {
		myCG = new sg::solver::ConjugateGradientsSP(SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
	} else if (SolverConfigRefine.type_ == sg::solver::BiCGSTAB) {
		myCG = new sg::solver::BiCGStabSP(SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
	} else {
		throw base::application_exception("LearnerBaseSP::train: An unsupported SLE solver type was chosen!");
	}

	// Pre-Procession
	preProcessing();

	if (isVerbose_)
		std::cout << "Starting Learning...." << std::endl;

	// execute adaptsteps
	sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();

	for (size_t i = 0; i < AdaptConfig.numRefinements_ + 1; i++) {
		if (isVerbose_)
			std::cout << std::endl << "Doing refinement: " << i << std::endl;

#ifdef USE_MPI
		// This barrier is needed since just the time measurement
		// of process 0 is printed
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		myStopwatch->start();

		// Do Refinements
		if (i > 0) {
#ifdef USE_MPI

			if (parallel::myGlobalMPIComm->getMyRank() == 0) {
#endif
			sg::base::DataVector alphaDP(alpha_->getSize());
			sg::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_, alphaDP);
			sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(&alphaDP,
					AdaptConfig.noPoints_, AdaptConfig.threshold_);
			grid_->createGridGenerator()->refine(myRefineFunc);
			delete myRefineFunc;
#ifdef USE_MPI
			std::string serialized_grid = grid_->getStorage()->serialize();

			parallel::myGlobalMPIComm->broadcastGridStorage(serialized_grid);
		} else {
			std::string serialized_grid = "";

			parallel::myGlobalMPIComm->receiveGridStorage(serialized_grid);

			grid_->getStorage()->emptyStorage();
			grid_->getStorage()->unserialize_noAlgoDims(serialized_grid);
		}

#endif
			DMSystem->rebuildLevelAndIndex();

			if (isVerbose_)
				std::cout << "New Grid Size: " << grid_->getSize() << std::endl;

			alpha_->resizeZero(grid_->getSize());
		} else {
			if (isVerbose_)
				std::cout << "Grid Size: " << grid_->getSize() << std::endl;
		}

		sg::base::DataVectorSP b(alpha_->getSize());
		DMSystem->generateb(classes, b);

		if (i == AdaptConfig.numRefinements_) {
			myCG->setMaxIterations(SolverConfigFinal.maxIterations_);
			myCG->setEpsilon(SolverConfigFinal.eps_);
		}

		myCG->solve(*DMSystem, *alpha_, b, true, false, 0.0);

#ifdef USE_MPI
		// This barrier is needed since just the time measurement
		// of process 0 is printed
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		execTime_ += myStopwatch->stop();

		if (isVerbose_) {
			std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
			std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;
		}

		// use post-processing to determine Flops and time
		if (i < AdaptConfig.numRefinements_) {
			postProcessing(trainDataset, SolverConfigRefine.type_, myCG->getNumberIterations());
		} else {
			postProcessing(trainDataset, SolverConfigFinal.type_, myCG->getNumberIterations());
		}

		double tmp1, tmp2, tmp3, tmp4;
		DMSystem->getTimers(tmp1, tmp2, tmp3, tmp4);
		result.timeComplete_ = execTime_;
		result.timeMultComplete_ = tmp1;
		result.timeMultCompute_ = tmp2;
		result.timeMultTransComplete_ = tmp3;
		result.timeMultTransCompute_ = tmp4;
		// @TODO fix regularization timings, if needed
		result.timeRegularization_ = 0.0;
		result.GFlop_ = GFlop_;
		result.GByte_ = GByte_;

		if (testAccDuringAdapt) {
			double acc = getAccuracy(trainDataset, classes);

			if (isVerbose_) {
				if (isRegression_) {
					if (isVerbose_)
						std::cout << "MSE (train): " << acc << std::endl;
				} else {
					if (isVerbose_)
						std::cout << "Acc (train): " << acc << std::endl;
				}
			}

			if (isRegression_) {
				if ((i > 0) && (oldAcc <= acc)) {
					if (isVerbose_)
						std::cout << "The grid is becoming worse --> stop learning" << std::endl;

					break;
				}
			} else {
				if ((i > 0) && (oldAcc >= acc)) {
					if (isVerbose_)
						std::cout << "The grid is becoming worse --> stop learning" << std::endl;

					break;
				}
			}

			oldAcc = acc;
		}
	}

	if (isVerbose_) {
		std::cout << "Finished Training!" << std::endl << std::endl;
		std::cout << "Training took: " << execTime_ << " seconds" << std::endl << std::endl;
	}

	isTrained_ = true;

	delete myStopwatch;
	delete myCG;
	delete DMSystem;

	return result;
}

LearnerTiming LearnerBaseSP::train(sg::base::DataMatrixSP& trainDataset, sg::base::DataVectorSP& classes,
		const sg::base::RegularGridConfiguration& GridConfig, const sg::solver::SLESolverSPConfiguration& SolverConfig,
		const float lambda) {
	sg::base::AdpativityConfiguration AdaptConfig;

	AdaptConfig.maxLevelType_ = false;
	AdaptConfig.noPoints_ = 0;
	AdaptConfig.numRefinements_ = 0;
	AdaptConfig.percent_ = 0.0;
	AdaptConfig.threshold_ = 0.0;

	return train(trainDataset, classes, GridConfig, SolverConfig, SolverConfig, AdaptConfig, false, lambda);
}

sg::base::DataVectorSP LearnerBaseSP::predict(sg::base::DataMatrixSP& testDataset) {
	sg::base::DataVectorSP classesComputed(testDataset.getNrows());

	sg::base::DataVector classesComputedDP(testDataset.getNrows());
	sg::base::DataVector alphaDP(grid_->getSize());
	sg::base::DataMatrix testDatasetDP(testDataset.getNrows(), testDataset.getNcols());

	sg::base::PrecisionConverter::convertDataMatrixSPToDataMatrix(testDataset, testDatasetDP);
	sg::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_, alphaDP);

	sg::base::OperationMultipleEval* MultEval = sg::op_factory::createOperationMultipleEval(*grid_, testDatasetDP);
	MultEval->mult(alphaDP, classesComputedDP);
	delete MultEval;

	sg::base::PrecisionConverter::convertDataVectorToDataVectorSP(classesComputedDP, classesComputed);

	return classesComputed;
}

void LearnerBaseSP::store(std::string tGridFilename, std::string tAlphaFilename) {
	// @TODO (heinecke) implement
	throw base::application_exception("LearnerBaseSP::store: This method isn't implemented, yet!");
}

double LearnerBaseSP::getAccuracy(sg::base::DataMatrixSP& testDataset, const sg::base::DataVectorSP& classesReference,
		const float threshold) {
	// evaluate test dataset
	sg::base::DataVectorSP classesComputed = predict(testDataset);

	return getAccuracy(classesComputed, classesReference, threshold);
}

double LearnerBaseSP::getAccuracy(const sg::base::DataVectorSP& classesComputed,
		const sg::base::DataVectorSP& classesReference, const float threshold) {
	double result = -1.0;

	if (classesComputed.getSize() != classesReference.getSize()) {
		throw base::application_exception("LearnerBaseSP::getAccuracy: lengths of classes vectors do not match!");
	}

	if (isRegression_) {
		sg::base::DataVectorSP tmp(classesComputed);
		tmp.sub(classesReference);
		tmp.sqr();
		result = tmp.sum();
		result /= static_cast<double>(tmp.getSize());
	} else {
		size_t correct = 0;

		for (size_t i = 0; i < classesComputed.getSize(); i++) {
			if ((classesComputed.get(i) >= threshold && classesReference.get(i) >= 0.0f)
					|| (classesComputed.get(i) < threshold && classesReference.get(i) < 0.0f)) {
				correct++;
			}
		}

		result = static_cast<double>(correct) / static_cast<double>(classesComputed.getSize());
	}

	return result;
}

ClassificatorQuality LearnerBaseSP::getCassificatorQuality(sg::base::DataMatrixSP& testDataset,
		const sg::base::DataVectorSP& classesReference, const float threshold) {
	// evaluate test dataset
	sg::base::DataVectorSP classesComputed = predict(testDataset);

	return getCassificatorQuality(classesComputed, classesReference, threshold);
}

ClassificatorQuality LearnerBaseSP::getCassificatorQuality(const sg::base::DataVectorSP& classesComputed,
		const sg::base::DataVectorSP& classesReference, const float threshold) {
	ClassificatorQuality result;

	if (isRegression_) {
		throw base::application_exception(
				"LearnerBaseSP::getCassificatorQuality: this method is not valid for regression problems!");
	}

	if (classesComputed.getSize() != classesReference.getSize()) {
		throw base::application_exception(
				"LearnerBaseSP::getCassificatorQuality: lengths of classes vectors do not match!");
	}

	result.truePositive_ = 0;
	result.trueNegative_ = 0;
	result.falsePositive_ = 0;
	result.falseNegative_ = 0;

	for (size_t i = 0; i < classesComputed.getSize(); i++) {
		if ((classesComputed.get(i) >= threshold && classesReference.get(i) >= 0.0f)) {
			result.truePositive_++;
		} else if ((classesComputed.get(i) < threshold && classesReference.get(i) < 0.0f)) {
			result.trueNegative_++;
		} else if ((classesComputed.get(i) >= threshold && classesReference.get(i) < 0.0f)) {
			result.falsePositive_++;
		} else { // ( (classesComputed.get(i) < threshold && classesReference.get(i) >= 0) )
			result.falseNegative_++;
		}
	}

	return result;
}

void LearnerBaseSP::dumpGrid(std::string tFilename) {
	if (isTrained_) {
		sg::base::GridPrinter myPlotter(*grid_);
		sg::base::DataVector tmp_alpha(alpha_->getSize());
		sg::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_, tmp_alpha);
		myPlotter.printSparseGrid(tmp_alpha, tFilename, false);
	}
}

void LearnerBaseSP::dumpFunction(std::string tFilename, size_t resolution) {
	if (isTrained_ && grid_->getStorage()->dim() <= 2) {
		sg::base::GridPrinter myPlotter(*grid_);
		sg::base::DataVector tmp_alpha(alpha_->getSize());
		sg::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_, tmp_alpha);
		myPlotter.printGrid(tmp_alpha, tFilename, resolution);
	}
}

bool LearnerBaseSP::getIsRegression() const {
	return isRegression_;
}

bool LearnerBaseSP::getIsVerbose() const {
	return isVerbose_;
}

void LearnerBaseSP::setIsVerbose(const bool isVerbose) {
	isVerbose_ = isVerbose;
}

}

}
