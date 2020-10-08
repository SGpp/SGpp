// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/base/tools/PrecisionConverter.hpp>

#include <sgpp/solver/sle/BiCGStabSP.hpp>
#include <sgpp/solver/sle/ConjugateGradientsSP.hpp>

#include <sgpp/datadriven/application/LearnerBaseSP.hpp>

#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/globaldef.hpp>

#include <iostream>
#include <string>

namespace sgpp {
namespace datadriven {

LearnerBaseSP::LearnerBaseSP(const bool isRegression, const bool isVerbose)
    : alpha_(nullptr),
      grid_(nullptr),
      isVerbose_(isVerbose),
      isRegression_(isRegression),
      isTrained_(false),
      execTime_(0.0),
      GFlop_(0.0),
      GByte_(0.0) {}

LearnerBaseSP::LearnerBaseSP(const std::string tGridFilename, const std::string tAlphaFilename,
                             const bool isRegression, const bool isVerbose)
    : alpha_(nullptr),
      grid_(nullptr),
      isVerbose_(isVerbose),
      isRegression_(isRegression),
      isTrained_(false),
      execTime_(0.0),
      GFlop_(0.0),
      GByte_(0.0) {
  throw base::application_exception(
      "LearnerBaseSP::LearnerBaseSP: This construct isn't "
      "implemented, yet!");
}

LearnerBaseSP::LearnerBaseSP(const LearnerBaseSP& copyMe)
    : isVerbose_(copyMe.isVerbose_),
      isRegression_(copyMe.isRegression_),
      isTrained_(false),
      execTime_(0.0),
      GFlop_(0.0),
      GByte_(0.0) {
  this->isRegression_ = copyMe.isRegression_;
  this->isTrained_ = false;
  this->GFlop_ = 0.0;
  this->GByte_ = 0.0;
  this->execTime_ = 0.0;

  // can be solved better with a grid copy constructor
  grid_ = sgpp::base::Grid::unserialize(copyMe.grid_->serialize());
  alpha_ = new sgpp::base::DataVectorSP(*(copyMe.alpha_));
}

LearnerBaseSP::~LearnerBaseSP() {
  // if user does no cleaning
  if (alpha_ != nullptr) delete alpha_;

  if (grid_ != nullptr) delete grid_;
}

void LearnerBaseSP::InitializeGrid(const sgpp::base::RegularGridConfiguration& gridConfig) {
  if (gridConfig.type_ == sgpp::base::GridType::LinearBoundary) {
    grid_ = new sgpp::base::LinearBoundaryGrid(gridConfig.dim_);
  } else if (gridConfig.type_ == sgpp::base::GridType::ModLinear) {
    grid_ = new sgpp::base::ModLinearGrid(gridConfig.dim_);
  } else if (gridConfig.type_ == sgpp::base::GridType::Linear) {
    grid_ = new sgpp::base::LinearGrid(gridConfig.dim_);
  } else {
    grid_ = nullptr;
    throw base::application_exception(
        "LearnerBaseSP::InitializeGrid: An unsupported grid type was "
        "chosen!");
  }

  // Generate regular Grid with LEVELS Levels
  grid_->getGenerator().regular(gridConfig.level_);

  // Create alpha
  alpha_ = new sgpp::base::DataVectorSP(grid_->getSize());
  alpha_->setAll(0.0);
}

void LearnerBaseSP::preProcessing() {}

void LearnerBaseSP::postProcessing(const sgpp::base::DataMatrixSP& trainDataset,
                                   const sgpp::solver::SLESolverType& solver,
                                   const size_t numNeededIterations) {
  if (this->isVerbose_) {
    std::cout << std::endl;
    std::cout << "Current Execution Time: " << execTime_ << std::endl;
    std::cout << std::endl;
  }
}

LearnerTiming LearnerBaseSP::train(sgpp::base::DataMatrixSP& trainDataset,
                                   sgpp::base::DataVectorSP& classes,
                                   const sgpp::base::RegularGridConfiguration& gridConfig,
                                   const sgpp::solver::SLESolverSPConfiguration& SolverConfigRefine,
                                   const sgpp::solver::SLESolverSPConfiguration& SolverConfigFinal,
                                   const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                                   const bool testAccDuringAdapt,
                                   const float lambdaRegularization) {
  LearnerTiming result;

  if (trainDataset.getNrows() != classes.getSize()) {
    throw base::application_exception(
        "LearnerBaseSP::train: length of classes vector does not "
        "match to dataset!");
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
  if (alpha_ != nullptr) delete alpha_;

  if (grid_ != nullptr) delete grid_;

  if (isTrained_ == true) isTrained_ = false;

  InitializeGrid(gridConfig);

  // check if grid was created
  if (grid_ == nullptr) return result;

  // create DMSystem
  sgpp::datadriven::DMSystemMatrixBaseSP* DMSystem =
      createDMSystem(trainDataset, lambdaRegularization);

  // check if System was created
  if (DMSystem == nullptr) return result;

  sgpp::solver::SLESolverSP* myCG;

  if (SolverConfigRefine.type_ == sgpp::solver::SLESolverType::CG) {
    myCG = new sgpp::solver::ConjugateGradientsSP(SolverConfigRefine.maxIterations_,
                                                  SolverConfigRefine.eps_);
  } else if (SolverConfigRefine.type_ == sgpp::solver::SLESolverType::BiCGSTAB) {
    myCG = new sgpp::solver::BiCGStabSP(SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
  } else {
    throw base::application_exception(
        "LearnerBaseSP::train: An unsupported SLE solver type was "
        "chosen!");
  }

  // Pre-Procession
  preProcessing();

  if (isVerbose_) std::cout << "Starting Learning...." << std::endl;

  // execute adaptsteps
  sgpp::base::SGppStopwatch* myStopwatch = new sgpp::base::SGppStopwatch();

  for (size_t i = 0; i < adaptivityConfig.numRefinements_ + 1; i++) {
    if (isVerbose_) std::cout << std::endl << "Doing refinement: " << i << std::endl;

    myStopwatch->start();

    // Do Refinements
    if (i > 0) {
      sgpp::base::DataVector alphaDP(alpha_->getSize());
      sgpp::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_, alphaDP);
      sgpp::base::SurplusRefinementFunctor myRefineFunc(
          alphaDP, adaptivityConfig.numRefinementPoints_, adaptivityConfig.refinementThreshold_);
      grid_->getGenerator().refine(myRefineFunc);
      DMSystem->rebuildLevelAndIndex();

      if (isVerbose_) std::cout << "New Grid Size: " << grid_->getSize() << std::endl;

      alpha_->resizeZero(grid_->getSize());
    } else {
      if (isVerbose_) std::cout << "Grid Size: " << grid_->getSize() << std::endl;
    }

    sgpp::base::DataVectorSP b(alpha_->getSize());
    DMSystem->generateb(classes, b);

    if (i == adaptivityConfig.numRefinements_) {
      myCG->setMaxIterations(SolverConfigFinal.maxIterations_);
      myCG->setEpsilon(SolverConfigFinal.eps_);
    }

    myCG->solve(*DMSystem, *alpha_, b, true, false, 0.0);

    execTime_ += myStopwatch->stop();

    if (isVerbose_) {
      std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
      std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;
    }

    // use post-processing to determine Flops and time
    if (i < adaptivityConfig.numRefinements_) {
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
    result.timeRegularization_ = 0.0;
    result.GFlop_ = GFlop_;
    result.GByte_ = GByte_;

    if (testAccDuringAdapt) {
      double acc = getAccuracy(trainDataset, classes);

      if (isVerbose_) {
        if (isRegression_) {
          if (isVerbose_) std::cout << "MSE (train): " << acc << std::endl;
        } else {
          if (isVerbose_) std::cout << "Acc (train): " << acc << std::endl;
        }
      }

      if (isRegression_) {
        if ((i > 0) && (oldAcc <= acc)) {
          if (isVerbose_) std::cout << "The grid is becoming worse --> stop learning" << std::endl;

          break;
        }
      } else {
        if ((i > 0) && (oldAcc >= acc)) {
          if (isVerbose_) std::cout << "The grid is becoming worse --> stop learning" << std::endl;

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

LearnerTiming LearnerBaseSP::train(sgpp::base::DataMatrixSP& trainDataset,
                                   sgpp::base::DataVectorSP& classes,
                                   const sgpp::base::RegularGridConfiguration& gridConfig,
                                   const sgpp::solver::SLESolverSPConfiguration& SolverConfig,
                                   const float lambdaRegularization) {
  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 0;
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.percent_ = 0.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  return train(trainDataset, classes, gridConfig, SolverConfig, SolverConfig, adaptivityConfig,
               false, lambdaRegularization);
}

sgpp::base::DataVectorSP LearnerBaseSP::predict(sgpp::base::DataMatrixSP& testDataset) {
  sgpp::base::DataVectorSP classesComputed(testDataset.getNrows());

  sgpp::base::DataVector classesComputedDP(testDataset.getNrows());
  sgpp::base::DataVector alphaDP(grid_->getSize());
  sgpp::base::DataMatrix testDatasetDP(testDataset.getNrows(), testDataset.getNcols());

  sgpp::base::PrecisionConverter::convertDataMatrixSPToDataMatrix(testDataset, testDatasetDP);
  sgpp::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_, alphaDP);

  sgpp::base::OperationMultipleEval* MultEval =
      sgpp::op_factory::createOperationMultipleEval(*grid_, testDatasetDP);
  MultEval->mult(alphaDP, classesComputedDP);
  delete MultEval;

  sgpp::base::PrecisionConverter::convertDataVectorToDataVectorSP(classesComputedDP,
                                                                  classesComputed);

  return classesComputed;
}

void LearnerBaseSP::store(std::string tGridFilename, std::string tAlphaFilename) {
  throw base::application_exception("LearnerBaseSP::store: This method isn't implemented, yet!");
}

double LearnerBaseSP::getAccuracy(sgpp::base::DataMatrixSP& testDataset,
                                  const sgpp::base::DataVectorSP& classesReference,
                                  const float threshold) {
  // evaluate test dataset
  sgpp::base::DataVectorSP classesComputed = predict(testDataset);

  return getAccuracy(classesComputed, classesReference, threshold);
}

double LearnerBaseSP::getAccuracy(const sgpp::base::DataVectorSP& classesComputed,
                                  const sgpp::base::DataVectorSP& classesReference,
                                  const float threshold) {
  double result = -1.0;

  if (classesComputed.getSize() != classesReference.getSize()) {
    throw base::application_exception(
        "LearnerBaseSP::getAccuracy: lengths of classes vectors do "
        "not match!");
  }

  if (isRegression_) {
    sgpp::base::DataVectorSP tmp(classesComputed);
    tmp.sub(classesReference);
    tmp.sqr();
    result = static_cast<double>(tmp.sum());
    result /= static_cast<double>(tmp.getSize());
  } else {
    size_t correct = 0;

    for (size_t i = 0; i < classesComputed.getSize(); i++) {
      if ((classesComputed.get(i) >= threshold && classesReference.get(i) >= 0.0f) ||
          (classesComputed.get(i) < threshold && classesReference.get(i) < 0.0f)) {
        correct++;
      }
    }

    result = static_cast<double>(correct) / static_cast<double>(classesComputed.getSize());
  }

  return result;
}

ClassificatorQuality LearnerBaseSP::getCassificatorQuality(
    sgpp::base::DataMatrixSP& testDataset, const sgpp::base::DataVectorSP& classesReference,
    const float threshold) {
  // evaluate test dataset
  sgpp::base::DataVectorSP classesComputed = predict(testDataset);

  return getCassificatorQuality(classesComputed, classesReference, threshold);
}

ClassificatorQuality LearnerBaseSP::getCassificatorQuality(
    const sgpp::base::DataVectorSP& classesComputed,
    const sgpp::base::DataVectorSP& classesReference, const float threshold) {
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
    } else {  // ( (classesComputed.get(i) < threshold && classesReference.get(i) >= 0) )
      result.falseNegative_++;
    }
  }

  return result;
}

void LearnerBaseSP::dumpGrid(std::string tFilename) {
  if (isTrained_) {
    sgpp::base::GridPrinter myPlotter(*grid_);
    sgpp::base::DataVector tmp_alpha(alpha_->getSize());
    sgpp::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_, tmp_alpha);
    myPlotter.printSparseGrid(tmp_alpha, tFilename, false);
  }
}

void LearnerBaseSP::dumpFunction(std::string tFilename, size_t resolution) {
  if (isTrained_ && grid_->getDimension() <= 2) {
    sgpp::base::GridPrinter myPlotter(*grid_);
    sgpp::base::DataVector tmp_alpha(alpha_->getSize());
    sgpp::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_, tmp_alpha);
    myPlotter.printGrid(tmp_alpha, tFilename, resolution);
  }
}

bool LearnerBaseSP::getIsRegression() const { return isRegression_; }

bool LearnerBaseSP::getIsVerbose() const { return isVerbose_; }

void LearnerBaseSP::setIsVerbose(const bool isVerbose) { isVerbose_ = isVerbose; }

}  // namespace datadriven
}  // namespace sgpp
