// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "ModelFittingLeastSquares.hpp"

#include <sgpp/datadriven/algorithm/SystemMatrixLeastSquaresIdentity.hpp>
#include <sgpp/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>

//TODO: use the system matrix with flexible regularization

namespace SGPP {
namespace datadriven {

ModelFittingLeastSquares::ModelFittingLeastSquares(
SGPP::datadriven::DataMiningConfigurationLeastSquares config) :
		datadriven::ModelFittingBase(), configuration(config) {
}

ModelFittingLeastSquares::~ModelFittingLeastSquares() {
}

datadriven::DMSystemMatrixBase*
ModelFittingLeastSquares::createDMSystem(base::DataMatrix& trainDataset,
		float_t lambda) {
	if (this->grid) {
		//TODO: improve
		throw base::application_exception("grid is null");
	}

	datadriven::SystemMatrixLeastSquaresIdentity* systemMatrix =
			new datadriven::SystemMatrixLeastSquaresIdentity(*(this->grid),
					trainDataset, lambda);
	systemMatrix->setImplementation(this->implementationConfiguration);
	return systemMatrix;
}

void ModelFittingLeastSquares::fit(datadriven::Dataset& dataset) {
	/*
	 LearnerTiming result;

	 if (trainDataset.getNrows() != classes.getSize()) {
	 throw base::application_exception(
	 "LearnerBase::train: length of classes vector does not match "
	 "to dataset!");
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

	 float_t oldAcc = 0.0;

	 // Construct Grid
	 if (alpha != NULL)
	 delete alpha;

	 if (grid != NULL)
	 delete grid;

	 if (isTrained_ == true)
	 isTrained_ = false;

	 InitializeGrid(GridConfig);

	 // check if grid was created
	 if (grid == NULL)
	 return result;

	 // create DMSystem
	 SGPP::datadriven::DMSystemMatrixBase* DMSystem = createDMSystem(
	 trainDataset, lambdaRegularization);

	 // check if System was created
	 if (DMSystem == NULL)
	 return result;

	 SGPP::solver::SLESolver* myCG;

	 if (SolverConfigRefine.type_ == SGPP::solver::SLESolverType::CG) {
	 myCG = new SGPP::solver::ConjugateGradients(
	 SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
	 } else if (SolverConfigRefine.type_
	 == SGPP::solver::SLESolverType::BiCGSTAB) {
	 myCG = new SGPP::solver::BiCGStab(SolverConfigRefine.maxIterations_,
	 SolverConfigRefine.eps_);
	 } else {
	 throw base::application_exception(
	 "LearnerBase::train: An unsupported SLE solver type was "
	 "chosen!");
	 }

	 // Pre-Procession
	 preProcessing();

	 if (isVerbose_)
	 std::cout << "Starting Learning...." << std::endl;

	 // execute adaptsteps
	 SGPP::base::SGppStopwatch* myStopwatch = new SGPP::base::SGppStopwatch();
	 SGPP::base::SGppStopwatch* myStopwatch2 = new SGPP::base::SGppStopwatch();

	 for (size_t i = 0; i < AdaptConfig.numRefinements_ + 1; i++) {
	 if (isVerbose_)
	 std::cout << std::endl << "Doing refinement: " << i << std::endl;

	 this->currentRefinementStep = i;

	 myStopwatch->start();

	 // Do Refinements
	 if (i > 0) {
	 myStopwatch2->start();

	 // disable refinement here!
	 SGPP::base::SurplusRefinementFunctor* myRefineFunc = new
	 SGPP::base::SurplusRefinementFunctor(alpha_, AdaptConfig.noPoints_,
	 AdaptConfig.threshold_);
	 grid_->createGridGenerator()->refine(myRefineFunc);
	 delete myRefineFunc;

	 // tell the SLE manager that the grid changed (for interal data structures)
	 DMSystem->prepareGrid();

	 alpha_->resizeZero(grid_->getSize());
	 float_t refineTime = myStopwatch2->stop();

	 if (isVerbose_)
	 std::cout << "New Grid Size: " << grid_->getSize()
	 << " (Refinement took " << refineTime << " secs)"
	 << std::endl;
	 } else {
	 if (isVerbose_)
	 std::cout << "Grid Size: " << grid_->getSize() << std::endl;
	 }

	 SGPP::base::DataVector b(alpha_->getSize());
	 DMSystem->generateb(classes, b);

	 if (i == AdaptConfig.numRefinements_) {
	 myCG->setMaxIterations(SolverConfigFinal.maxIterations_);
	 myCG->setEpsilon(SolverConfigFinal.eps_);
	 }

	 myCG->solve(*DMSystem, *alpha_, b, true, false, 0.0);

	 float_t stopTime = myStopwatch->stop();
	 this->execTime_ += stopTime;
	 this->stepExecTime_ = stopTime;

	 if (isVerbose_) {
	 std::cout << std::endl;
	 std::cout << "Needed Iterations: " << myCG->getNumberIterations()
	 << std::endl;
	 std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;
	 }

	 // use post-processing to determine Flops and time
	 if (i < AdaptConfig.numRefinements_) {
	 postProcessing(trainDataset, SolverConfigRefine.type_,
	 myCG->getNumberIterations());
	 } else {
	 postProcessing(trainDataset, SolverConfigFinal.type_,
	 myCG->getNumberIterations());
	 }

	 float_t tmp1, tmp2, tmp3, tmp4;
	 DMSystem->getTimers(tmp1, tmp2, tmp3, tmp4);
	 result.timeComplete_ = execTime_;
	 result.timeMultComplete_ = tmp1;
	 result.timeMultCompute_ = tmp2;
	 result.timeMultTransComplete_ = tmp3;
	 result.timeMultTransCompute_ = tmp4;
	 result.timeRegularization_ = 0.0;
	 result.GFlop_ = GFlop_;
	 result.GByte_ = GByte_;

	 if (testAccDuringAdapt) {
	 float_t acc = getAccuracy(trainDataset, classes);

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
	 std::cout
	 << "The grid is becoming worse --> stop learning"
	 << std::endl;

	 break;
	 }
	 } else {
	 if ((i > 0) && (oldAcc >= acc)) {
	 if (isVerbose_)
	 std::cout
	 << "The grid is becoming worse --> stop learning"
	 << std::endl;

	 break;
	 }
	 }

	 oldAcc = acc;
	 }
	 }

	 if (isVerbose_) {
	 std::cout << "Finished Training!" << std::endl << std::endl;
	 std::cout << "Training took: " << execTime_ << " seconds" << std::endl
	 << std::endl;
	 }

	 isTrained_ = true;

	 delete myStopwatch;
	 delete myStopwatch2;
	 delete myCG;
	 delete DMSystem;

	 return result;
	 */
}

}
}
