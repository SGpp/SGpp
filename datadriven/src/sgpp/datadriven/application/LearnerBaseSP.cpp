// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/PrecisionConverter.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>

#include <sgpp/solver/sle/ConjugateGradientsSP.hpp>
#include <sgpp/solver/sle/BiCGStabSP.hpp>

#include <sgpp/datadriven/application/LearnerBaseSP.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>

#include <iostream>
#include <string>

namespace SGPP {

namespace datadriven {

LearnerBaseSP::LearnerBaseSP(const bool isRegression, const bool isVerbose) :
  alpha_(NULL), grid_(NULL), isVerbose_(isVerbose), isRegression_(isRegression),
  isTrained_(false), execTime_(
    0.0), GFlop_(0.0), GByte_(0.0) {
}

LearnerBaseSP::LearnerBaseSP(const std::string tGridFilename,
                             const std::string tAlphaFilename, const bool isRegression,
                             const bool isVerbose) :
  alpha_(NULL), grid_(NULL), isVerbose_(isVerbose), isRegression_(isRegression),
  isTrained_(false), execTime_(
    0.0), GFlop_(0.0), GByte_(0.0) {
  throw base::application_exception("LearnerBaseSP::LearnerBaseSP: This construct isn't "
    "implemented, yet!");
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
  isVerbose_(copyMe.isVerbose_), isRegression_(copyMe.isRegression_),
  isTrained_(false), execTime_(0.0), GFlop_(
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

  // can be solved better with a grid copy constructor
  grid_ = SGPP::base::Grid::unserialize(copyMe.grid_->serialize()).release();
  alpha_ = new SGPP::base::DataVectorSP(*(copyMe.alpha_));
}

LearnerBaseSP::~LearnerBaseSP() {
  // if user does no cleaning
  if (alpha_ != NULL)
    delete alpha_;

  if (grid_ != NULL)
    delete grid_;
}

void LearnerBaseSP::InitializeGrid(const SGPP::base::RegularGridConfiguration&
                                   GridConfig) {
  if (GridConfig.type_ == SGPP::base::GridType::LinearBoundary) {
    grid_ = new SGPP::base::LinearBoundaryGrid(GridConfig.dim_);
  } else if (GridConfig.type_ == SGPP::base::GridType::ModLinear) {
    grid_ = new SGPP::base::ModLinearGrid(GridConfig.dim_);
  } else if (GridConfig.type_ == SGPP::base::GridType::Linear) {
    grid_ = new SGPP::base::LinearGrid(GridConfig.dim_);
  } else {
    grid_ = NULL;
    throw base::application_exception("LearnerBaseSP::InitializeGrid: An unsupported grid type was "
      "chosen!");
  }

  // Generate regular Grid with LEVELS Levels
  grid_->getGenerator().regular(GridConfig.level_);

  // Create alpha
  alpha_ = new SGPP::base::DataVectorSP(grid_->getSize());
  alpha_->setAll(0.0);
}

void LearnerBaseSP::preProcessing() {
}

void LearnerBaseSP::postProcessing(const SGPP::base::DataMatrixSP& trainDataset,
                                   const SGPP::solver::SLESolverType& solver,
                                   const size_t numNeededIterations) {
  if (this->isVerbose_) {
    std::cout << std::endl;
    std::cout << "Current Execution Time: " << execTime_ << std::endl;
    std::cout << std::endl;
  }
}

LearnerTiming LearnerBaseSP::train(SGPP::base::DataMatrixSP& trainDataset,
                                   SGPP::base::DataVectorSP& classes,
                                   const SGPP::base::RegularGridConfiguration& GridConfig,
                                   const SGPP::solver::SLESolverSPConfiguration& SolverConfigRefine,
                                   const SGPP::solver::SLESolverSPConfiguration& SolverConfigFinal,
                                   const SGPP::base::AdpativityConfiguration& AdaptConfig,
                                   const bool testAccDuringAdapt,
                                   const float lambdaRegularization) {
  LearnerTiming result;

  if (trainDataset.getNrows() != classes.getSize()) {
    throw base::application_exception("LearnerBaseSP::train: length of classes vector does not "
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
  SGPP::datadriven::DMSystemMatrixBaseSP* DMSystem = createDMSystem(trainDataset,
      lambdaRegularization);

  // check if System was created
  if (DMSystem == NULL)
    return result;

  SGPP::solver::SLESolverSP* myCG;

  if (SolverConfigRefine.type_ == SGPP::solver::SLESolverType::CG) {
    myCG = new SGPP::solver::ConjugateGradientsSP(SolverConfigRefine.maxIterations_,
        SolverConfigRefine.eps_);
  } else if (SolverConfigRefine.type_ == SGPP::solver::SLESolverType::BiCGSTAB) {
    myCG = new SGPP::solver::BiCGStabSP(SolverConfigRefine.maxIterations_,
                                        SolverConfigRefine.eps_);
  } else {
    throw base::application_exception("LearnerBaseSP::train: An unsupported SLE solver type was "
      "chosen!");
  }

  // Pre-Procession
  preProcessing();

  if (isVerbose_)
    std::cout << "Starting Learning...." << std::endl;

  // execute adaptsteps
  SGPP::base::SGppStopwatch* myStopwatch = new SGPP::base::SGppStopwatch();

  for (size_t i = 0; i < AdaptConfig.numRefinements_ + 1; i++) {
    if (isVerbose_)
      std::cout << std::endl << "Doing refinement: " << i << std::endl;

    myStopwatch->start();

    // Do Refinements
    if (i > 0) {
      SGPP::base::DataVector alphaDP(alpha_->getSize());
      SGPP::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_,
          alphaDP);
      SGPP::base::SurplusRefinementFunctor myRefineFunc(
          alphaDP, AdaptConfig.noPoints_, AdaptConfig.threshold_);
      grid_->getGenerator().refine(myRefineFunc);
      DMSystem->rebuildLevelAndIndex();

      if (isVerbose_)
        std::cout << "New Grid Size: " << grid_->getSize() << std::endl;

      alpha_->resizeZero(grid_->getSize());
    } else {
      if (isVerbose_)
        std::cout << "Grid Size: " << grid_->getSize() << std::endl;
    }

    SGPP::base::DataVectorSP b(alpha_->getSize());
    DMSystem->generateb(classes, b);

    if (i == AdaptConfig.numRefinements_) {
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
    std::cout << "Training took: " << execTime_ << " seconds" << std::endl <<
              std::endl;
  }

  isTrained_ = true;

  delete myStopwatch;
  delete myCG;
  delete DMSystem;

  return result;
}

LearnerTiming LearnerBaseSP::train(SGPP::base::DataMatrixSP& trainDataset,
                                   SGPP::base::DataVectorSP& classes,
                                   const SGPP::base::RegularGridConfiguration& GridConfig,
                                   const SGPP::solver::SLESolverSPConfiguration& SolverConfig,
                                   const float lambdaRegularization) {
  SGPP::base::AdpativityConfiguration AdaptConfig;

  AdaptConfig.maxLevelType_ = false;
  AdaptConfig.noPoints_ = 0;
  AdaptConfig.numRefinements_ = 0;
  AdaptConfig.percent_ = 0.0;
  AdaptConfig.threshold_ = 0.0;

  return train(trainDataset, classes, GridConfig, SolverConfig, SolverConfig,
               AdaptConfig, false, lambdaRegularization);
}

SGPP::base::DataVectorSP LearnerBaseSP::predict(SGPP::base::DataMatrixSP&
    testDataset) {
  SGPP::base::DataVectorSP classesComputed(testDataset.getNrows());

  SGPP::base::DataVector classesComputedDP(testDataset.getNrows());
  SGPP::base::DataVector alphaDP(grid_->getSize());
  SGPP::base::DataMatrix testDatasetDP(testDataset.getNrows(),
                                       testDataset.getNcols());

  SGPP::base::PrecisionConverter::convertDataMatrixSPToDataMatrix(testDataset,
      testDatasetDP);
  SGPP::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_,
      alphaDP);

  SGPP::base::OperationMultipleEval* MultEval =
    SGPP::op_factory::createOperationMultipleEval(*grid_, testDatasetDP).release();
  MultEval->mult(alphaDP, classesComputedDP);
  delete MultEval;

  SGPP::base::PrecisionConverter::convertDataVectorToDataVectorSP(
    classesComputedDP, classesComputed);

  return classesComputed;
}

void LearnerBaseSP::store(std::string tGridFilename,
                          std::string tAlphaFilename) {
  throw base::application_exception("LearnerBaseSP::store: This method isn't implemented, yet!");
}

double LearnerBaseSP::getAccuracy(SGPP::base::DataMatrixSP& testDataset,
                                  const SGPP::base::DataVectorSP& classesReference,
                                  const float threshold) {
  // evaluate test dataset
  SGPP::base::DataVectorSP classesComputed = predict(testDataset);

  return getAccuracy(classesComputed, classesReference, threshold);
}

double LearnerBaseSP::getAccuracy(const SGPP::base::DataVectorSP&
                                  classesComputed,
                                  const SGPP::base::DataVectorSP& classesReference,
                                  const float threshold) {
  double result = -1.0;

  if (classesComputed.getSize() != classesReference.getSize()) {
    throw base::application_exception("LearnerBaseSP::getAccuracy: lengths of classes vectors do "
      "not match!");
  }

  if (isRegression_) {
    SGPP::base::DataVectorSP tmp(classesComputed);
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

    result = static_cast<double>(correct) / static_cast<double>
             (classesComputed.getSize());
  }

  return result;
}

ClassificatorQuality LearnerBaseSP::getCassificatorQuality(
  SGPP::base::DataMatrixSP& testDataset,
  const SGPP::base::DataVectorSP& classesReference, const float threshold) {
  // evaluate test dataset
  SGPP::base::DataVectorSP classesComputed = predict(testDataset);

  return getCassificatorQuality(classesComputed, classesReference, threshold);
}

ClassificatorQuality LearnerBaseSP::getCassificatorQuality(
  const SGPP::base::DataVectorSP& classesComputed,
  const SGPP::base::DataVectorSP& classesReference, const float threshold) {
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
    } else if ((classesComputed.get(i) < threshold
                && classesReference.get(i) < 0.0f)) {
      result.trueNegative_++;
    } else if ((classesComputed.get(i) >= threshold
                && classesReference.get(i) < 0.0f)) {
      result.falsePositive_++;
    } else {  // ( (classesComputed.get(i) < threshold && classesReference.get(i) >= 0) )
      result.falseNegative_++;
    }
  }

  return result;
}

void LearnerBaseSP::dumpGrid(std::string tFilename) {
  if (isTrained_) {
    SGPP::base::GridPrinter myPlotter(*grid_);
    SGPP::base::DataVector tmp_alpha(alpha_->getSize());
    SGPP::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_,
        tmp_alpha);
    myPlotter.printSparseGrid(tmp_alpha, tFilename, false);
  }
}

void LearnerBaseSP::dumpFunction(std::string tFilename, size_t resolution) {
  if (isTrained_ && grid_->getDimension() <= 2) {
    SGPP::base::GridPrinter myPlotter(*grid_);
    SGPP::base::DataVector tmp_alpha(alpha_->getSize());
    SGPP::base::PrecisionConverter::convertDataVectorSPToDataVector(*alpha_,
        tmp_alpha);
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

}  // namespace datadriven
}  // namespace SGPP

