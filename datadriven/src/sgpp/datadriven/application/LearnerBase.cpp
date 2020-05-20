// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/datadriven/application/LearnerBase.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

LearnerBase::LearnerBase(const bool isRegression, const bool isVerbose)
    : isVerbose(isVerbose),
      isRegression(isRegression),
      reuseCoefficients(true),
      solverVerbose(false),
      isTrained(false),
      execTime(0.0),
      stepExecTime(0.0),
      GFlop(0.0),
      stepGFlop(0.0),
      GByte(0.0),
      stepGByte(0.0),
      currentRefinementStep(0) {}

// LearnerBase::LearnerBase(const std::string tGridFilename, const std::string
// tAlphaFilename,
//                         const bool isRegression, const bool isVerbose)
//    : isVerbose(isVerbose),
//      isRegression(isRegression),
//      reuseCoefficients(false),
//      isTrained(false),
//      execTime(0.0),
//      stepExecTime(0.0),
//      GFlop(0.0),
//      stepGFlop(0.0),
//      GByte(0.0),
//      stepGByte(0.0),
//      currentRefinementStep(0) {
//  throw base::application_exception(
//      "LearnerBase::LearnerBase: This construct isn't implemented, yet!");
//}

LearnerBase::LearnerBase(const LearnerBase& copyMe) {
  this->isVerbose = copyMe.isVerbose;
  this->reuseCoefficients = copyMe.reuseCoefficients;
  this->solverVerbose = copyMe.solverVerbose;
  this->isTrained = false;
  this->isRegression = copyMe.isRegression;
  this->execTime = -1.0;
  this->GFlop = -1.0;
  this->GByte = -1.0;
  this->stepExecTime = -1.0;
  this->stepGFlop = -1.0;
  this->stepGByte = -1.0;
  this->currentRefinementStep = 0;

  // TODO(pfandedd): don't use grid serialization to not have to implement a
  // copy constructor!
  grid.reset(sgpp::base::Grid::unserialize(copyMe.grid->serialize()));
  alpha = std::make_unique<sgpp::base::DataVector>(*(copyMe.alpha));
}

LearnerBase::~LearnerBase() {}

void LearnerBase::InitializeGrid(const sgpp::base::RegularGridConfiguration& gridConfig) {
  if (gridConfig.type_ == sgpp::base::GridType::LinearBoundary) {
    grid = std::make_unique<sgpp::base::LinearBoundaryGrid>(gridConfig.dim_);
  } else if (gridConfig.type_ == sgpp::base::GridType::ModLinear) {
    grid = std::make_unique<sgpp::base::ModLinearGrid>(gridConfig.dim_);
  } else if (gridConfig.type_ == sgpp::base::GridType::Linear) {
    grid = std::make_unique<sgpp::base::LinearGrid>(gridConfig.dim_);
  } else {
    throw base::application_exception(
        "LearnerBase::InitializeGrid: An unsupported grid type was chosen!");
  }

  // Generate regular Grid with LEVELS Levels
  grid->getGenerator().regular(gridConfig.level_);

  // Create alpha
  alpha = std::make_unique<sgpp::base::DataVector>(grid->getSize());
  alpha->setAll(0.0);
}

void LearnerBase::preProcessing() {}

void LearnerBase::postProcessing(const sgpp::base::DataMatrix& trainDataset,
                                 const sgpp::solver::SLESolverType& solver,
                                 const size_t numNeededIterations) {
  if (this->isVerbose) {
    std::cout << std::endl;
    std::cout << "Current Execution Time: " << execTime << std::endl;
    std::cout << std::endl;
  }
}

LearnerTiming LearnerBase::train(sgpp::base::DataMatrix& trainDataset,
                                 sgpp::base::DataVector& classes,
                                 const sgpp::base::RegularGridConfiguration& gridConfig,
                                 const sgpp::solver::SLESolverConfiguration& SolverConfigRefine,
                                 const sgpp::solver::SLESolverConfiguration& SolverConfigFinal,
                                 const sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                                 const bool testAccDuringAdapt, const double lambdaRegularization,
                                 sgpp::base::DataMatrix* testDataset,
                                 sgpp::base::DataVector* testClasses) {
  LearnerTiming result;

  if (trainDataset.getNrows() != classes.getSize()) {
    throw base::application_exception(
        "LearnerBase::train: length of classes vector does not match to "
        "dataset!");
  }

  result.timeComplete_ = 0.0;
  result.timeMultComplete_ = 0.0;
  result.timeMultCompute_ = 0.0;
  result.timeMultTransComplete_ = 0.0;
  result.timeMultTransCompute_ = 0.0;
  result.timeRegularization_ = 0.0;
  result.GFlop_ = 0.0;
  result.GByte_ = 0.0;

  execTime = 0.0;
  GFlop = 0.0;
  GByte = 0.0;

  double oldAcc = 0.0;

  // Construct Grid
  //  if (alpha_ != nullptr) delete alpha_;
  //  if (grid_ != nullptr) delete grid_;

  if (isTrained == true) isTrained = false;

  InitializeGrid(gridConfig);

  // check if grid was created
  if (!grid.operator bool()) {
    throw base::application_exception("error: couldn't create grid");
  }

  // create DMSystem
  std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase> DMSystem =
      createDMSystem(trainDataset, lambdaRegularization);

  // check if System was created
  if (!DMSystem.operator bool()) {
    throw base::application_exception("error: couldn't create DMSystem");
  }

  std::unique_ptr<sgpp::solver::SLESolver> myCG;

  if (SolverConfigRefine.type_ == sgpp::solver::SLESolverType::CG) {
    myCG = std::make_unique<sgpp::solver::ConjugateGradients>(SolverConfigRefine.maxIterations_,
                                                              SolverConfigRefine.eps_);
  } else if (SolverConfigRefine.type_ == sgpp::solver::SLESolverType::BiCGSTAB) {
    myCG = std::make_unique<sgpp::solver::BiCGStab>(SolverConfigRefine.maxIterations_,
                                                    SolverConfigRefine.eps_);
  } else {
    throw base::application_exception(
        "LearnerBase::train: An unsupported SLE solver type was chosen!");
  }

  // Pre-Procession
  preProcessing();

  if (isVerbose) std::cout << "Starting Learning...." << std::endl;

  // execute adaptsteps
  std::unique_ptr<sgpp::base::SGppStopwatch> myStopwatch =
      std::make_unique<sgpp::base::SGppStopwatch>();
  std::unique_ptr<sgpp::base::SGppStopwatch> myStopwatch2 =
      std::make_unique<sgpp::base::SGppStopwatch>();

  for (size_t i = 0; i < adaptivityConfig.numRefinements_ + 1; i++) {
    if (isVerbose) std::cout << std::endl << "Doing refinement: " << i << std::endl;

    this->currentRefinementStep = i;

    myStopwatch->start();

    // Do Refinements
    if (i > 0) {
      myStopwatch2->start();

      // disable refinement here!
      if (adaptivityConfig.errorBasedRefinement_) {
        std::unique_ptr<sgpp::base::DataVector> residuals =
            std::make_unique<sgpp::base::DataVector>(alpha->getSize());
        this->predict(trainDataset, *residuals);
        residuals->sub(classes);
        residuals->sqr();
        std::unique_ptr<sgpp::base::DataVector> mseResiduals =
            std::make_unique<sgpp::base::DataVector>(grid->getSize());
        multTranspose(trainDataset, *residuals, *mseResiduals);
        mseResiduals->componentwise_mult(*alpha);
        sgpp::base::SurplusRefinementFunctor myRefineFunc(*mseResiduals,
                                                          adaptivityConfig.numRefinementPoints_,
                                                          adaptivityConfig.refinementThreshold_);
        grid->getGenerator().refine(myRefineFunc);
      } else {
        sgpp::base::SurplusRefinementFunctor myRefineFunc(
            *alpha, adaptivityConfig.numRefinementPoints_, adaptivityConfig.refinementThreshold_);
        grid->getGenerator().refine(myRefineFunc);
      }

      // tell the SLE manager that the grid changed (for interal data
      // structures)
      DMSystem->prepareGrid();

      alpha->resizeZero(grid->getSize());
      double refineTime = myStopwatch2->stop();

      if (isVerbose)
        std::cout << "New Grid Size: " << grid->getSize() << " (Refinement took " << refineTime
                  << " secs)" << std::endl;
    } else {
      if (isVerbose) std::cout << "Grid Size: " << grid->getSize() << std::endl;
    }

    sgpp::base::DataVector b(alpha->getSize());
    DMSystem->generateb(classes, b);

    if (i == adaptivityConfig.numRefinements_) {
      myCG->setMaxIterations(SolverConfigFinal.maxIterations_);
      myCG->setEpsilon(SolverConfigFinal.eps_);
    }

    myCG->solve(*DMSystem, *alpha, b, reuseCoefficients, solverVerbose, 0.0);

    double stopTime = myStopwatch->stop();
    this->execTime += stopTime;
    this->stepExecTime = stopTime;

    if (isVerbose) {
      std::cout << std::endl;
      std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
      std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;
    }

    // use post-processing to determine Flops and time
    if (i < adaptivityConfig.numRefinements_) {
      postProcessing(trainDataset, SolverConfigRefine.type_, myCG->getNumberIterations());
    } else {
      postProcessing(trainDataset, SolverConfigFinal.type_, myCG->getNumberIterations());
    }

    double timeMult, computeMult, timeMultTrans, computeMultTrans;
    DMSystem->getTimers(timeMult, computeMult, timeMultTrans, computeMultTrans);
    result.timeComplete_ = execTime;
    result.timeMultComplete_ = timeMult;
    result.timeMultCompute_ = computeMult;
    result.timeMultTransComplete_ = timeMultTrans;
    result.timeMultTransCompute_ = computeMultTrans;
    result.timeRegularization_ = 0.0;
    result.GFlop_ = GFlop;
    result.GByte_ = GByte;

    if (testAccDuringAdapt) {
      double acc = getAccuracy(trainDataset, classes);

      if (isVerbose) {
        if (isRegression) {
          if (isVerbose) std::cout << "MSE (train): " << acc << std::endl;
        } else {
          if (isVerbose) std::cout << "Acc (train): " << acc << std::endl;
        }
      }

      if (testDataset != nullptr && testClasses != nullptr) {
        double testAcc = getAccuracy(*testDataset, *testClasses);

        if (isVerbose) {
          if (isRegression) {
            if (isVerbose) std::cout << "MSE (test): " << testAcc << std::endl;
          } else {
            if (isVerbose) std::cout << "Acc (test): " << testAcc << std::endl;
          }
        }
      }

      if (isRegression) {
        if ((i > 0) && (oldAcc <= acc)) {
          if (isVerbose) std::cout << "The grid is becoming worse --> stop learning" << std::endl;

          break;
        }
      } else {
        if ((i > 0) && (oldAcc >= acc)) {
          if (isVerbose) std::cout << "The grid is becoming worse --> stop learning" << std::endl;

          break;
        }
      }
      oldAcc = acc;
    }
  }

  if (isVerbose) {
    std::cout << "Finished Training!" << std::endl << std::endl;
    std::cout << "Training took: " << execTime << " seconds" << std::endl << std::endl;
  }

  isTrained = true;

  //  delete myStopwatch;
  //  delete myStopwatch2;
  //  delete myCG;
  //  delete DMSystem;

  return result;
}

LearnerTiming LearnerBase::train(sgpp::base::DataMatrix& trainDataset,
                                 sgpp::base::DataVector& classes,
                                 const sgpp::base::RegularGridConfiguration& gridConfig,
                                 const sgpp::solver::SLESolverConfiguration& SolverConfig,
                                 const double lambdaRegularization) {
  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 0;
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.percent_ = 0.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  return train(trainDataset, classes, gridConfig, SolverConfig, SolverConfig, adaptivityConfig,
               false, lambdaRegularization);
}

void LearnerBase::predict(sgpp::base::DataMatrix& testDataset,
                          sgpp::base::DataVector& classesComputed) {
  classesComputed.resize(testDataset.getNrows());

  std::unique_ptr<sgpp::base::OperationMultipleEval> MultEval(
      sgpp::op_factory::createOperationMultipleEval(*grid, testDataset));
  MultEval->mult(*alpha, classesComputed);
}

void LearnerBase::multTranspose(sgpp::base::DataMatrix& dataset, sgpp::base::DataVector& multiplier,
                                sgpp::base::DataVector& result) {
  result.resize(grid->getSize());

  std::unique_ptr<sgpp::base::OperationMultipleEval> MultEval(
      sgpp::op_factory::createOperationMultipleEval(*grid, dataset));
  MultEval->multTranspose(multiplier, result);
}

void LearnerBase::store(std::string tGridFilename, std::string tAlphaFilename) {
  throw base::application_exception("LearnerBase::store: This method isn't implemented, yet!");
}

double LearnerBase::getAccuracy(sgpp::base::DataMatrix& testDataset,
                                const sgpp::base::DataVector& classesReference,
                                const double threshold) {
  // evaluate test dataset

  sgpp::base::DataVector classesComputed(testDataset.getNrows());
  predict(testDataset, classesComputed);

  return getAccuracy(classesComputed, classesReference, threshold);
}

double LearnerBase::getAccuracy(const sgpp::base::DataVector& classesComputed,
                                const sgpp::base::DataVector& classesReference,
                                const double threshold) {
  double result = -1.0;

  if (classesComputed.getSize() != classesReference.getSize()) {
    throw base::application_exception(
        "LearnerBase::getAccuracy: lengths of classes vectors do not match!");
  }

  if (isRegression) {
    sgpp::base::DataVector tmp(classesComputed);
    tmp.sub(classesReference);
    tmp.sqr();
    result = tmp.sum();
    result /= static_cast<double>(tmp.getSize());
  } else {
    size_t correct = 0;

    for (size_t i = 0; i < classesComputed.getSize(); i++) {
      if ((classesComputed.get(i) >= threshold && classesReference.get(i) >= 0.0) ||
          (classesComputed.get(i) < threshold && classesReference.get(i) < 0.0)) {
        correct++;
      }
    }

    result = static_cast<double>(correct) / static_cast<double>(classesComputed.getSize());
  }

  return result;
}

ClassificatorQuality LearnerBase::getCassificatorQuality(
    sgpp::base::DataMatrix& testDataset, const sgpp::base::DataVector& classesReference,
    const double threshold) {
  // evaluate test dataset
  sgpp::base::DataVector classesComputed(testDataset.getNrows());
  predict(testDataset, classesComputed);

  return getCassificatorQuality(classesComputed, classesReference, threshold);
}

ClassificatorQuality LearnerBase::getCassificatorQuality(
    const sgpp::base::DataVector& classesComputed, const sgpp::base::DataVector& classesReference,
    const double threshold) {
  ClassificatorQuality result;

  if (isRegression) {
    throw base::application_exception(
        "LearnerBase::getCassificatorQuality: this method is not valid for "
        "regression problems!");
  }

  if (classesComputed.getSize() != classesReference.getSize()) {
    throw base::application_exception(
        "LearnerBase::getCassificatorQuality: lengths of classes vectors do "
        "not match!");
  }

  result.truePositive_ = 0;
  result.trueNegative_ = 0;
  result.falsePositive_ = 0;
  result.falseNegative_ = 0;

  for (size_t i = 0; i < classesComputed.getSize(); i++) {
    if ((classesComputed.get(i) >= threshold && classesReference.get(i) >= 0.0)) {
      result.truePositive_++;
    } else if ((classesComputed.get(i) < threshold && classesReference.get(i) < 0.0)) {
      result.trueNegative_++;
    } else if ((classesComputed.get(i) >= threshold && classesReference.get(i) < 0.0)) {
      result.falsePositive_++;
    } else {  // ( (classesComputed.get(i) < threshold &&
              // classesReference.get(i) >= 0) )
      result.falseNegative_++;
    }
  }

  return result;
}

void LearnerBase::dumpGrid(std::string tFilename) {
  if (isTrained) {
    sgpp::base::GridPrinter myPlotter(*grid);
    myPlotter.printSparseGrid(*alpha, tFilename, false);
  }
}

void LearnerBase::dumpFunction(std::string tFilename, size_t resolution) {
  if (isTrained && grid->getDimension() <= 2) {
    sgpp::base::GridPrinter myPlotter(*grid);
    myPlotter.printGrid(*alpha, tFilename, resolution);
  }
}

bool LearnerBase::getIsRegression() const { return isRegression; }

bool LearnerBase::getIsVerbose() const { return isVerbose; }

void LearnerBase::setIsVerbose(const bool isVerbose) { this->isVerbose = isVerbose; }

std::vector<std::pair<size_t, double> > LearnerBase::getRefinementExecTimes() {
  return this->ExecTimeOnStep;
}

sgpp::base::Grid& LearnerBase::getGrid() {
  if (this->grid.get() == nullptr) {
    throw;
  }
  return *this->grid;
}

sgpp::base::DataVector& LearnerBase::getAlpha() {
  if (this->alpha == nullptr) {
    throw;
  }
  return *this->alpha;
}

void LearnerBase::setReuseCoefficients(bool reuseCoefficients) {
  this->reuseCoefficients = reuseCoefficients;
}

void LearnerBase::setSolverVerbose(bool solverVerbose) { this->solverVerbose = solverVerbose; }

}  // namespace datadriven
}  // namespace sgpp
