// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>

#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>

#include <sgpp/datadriven/application/LearnerBase.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>


namespace SGPP {

  namespace datadriven {

    LearnerBase::LearnerBase(const bool isRegression, const bool isVerbose) :
      alpha_(NULL), grid_(NULL), isVerbose_(isVerbose), isRegression_(isRegression), isTrained_(false), execTime_(
        0.0), stepExecTime_(0.0), GFlop_(0.0), stepGFlop_(0.0), GByte_(0.0), stepGByte_(0.0), currentRefinementStep(0) {
    }

    LearnerBase::LearnerBase(const std::string tGridFilename, const std::string tAlphaFilename, const bool isRegression,
                             const bool isVerbose) :
      alpha_(NULL), grid_(NULL), isVerbose_(isVerbose), isRegression_(isRegression), isTrained_(false), execTime_(
        0.0), stepExecTime_(0.0), GFlop_(0.0), stepGFlop_(0.0), GByte_(0.0), stepGByte_(0.0), currentRefinementStep(0) {
      throw base::application_exception("LearnerBase::LearnerBase: This construct isn't implemented, yet!");
    }

    LearnerBase::LearnerBase(const LearnerBase& copyMe) {
      this->isVerbose_ = copyMe.isVerbose_;
      this->isTrained_ = false;
      this->isRegression_ = copyMe.isRegression_;
      this->execTime_ = -1.0;
      this->GFlop_ = -1.0;
      this->GByte_ = -1.0;
      this->stepExecTime_ = -1.0;
      this->stepGFlop_ = -1.0;
      this->stepGByte_ = -1.0;
      this->currentRefinementStep = 0;

      // safety, should not happen
      if (alpha_ != NULL)
        delete alpha_;

      if (grid_ != NULL)
        delete grid_;

      grid_ = SGPP::base::Grid::unserialize(copyMe.grid_->serialize());
      alpha_ = new SGPP::base::DataVector(*(copyMe.alpha_));
    }

    LearnerBase::~LearnerBase() {
      // if user does no cleaning
      if (alpha_ != NULL)
        delete alpha_;

      if (grid_ != NULL)
        delete grid_;
    }

    void LearnerBase::InitializeGrid(const SGPP::base::RegularGridConfiguration& GridConfig) {
      if (GridConfig.type_ == SGPP::base::GridType::LinearBoundary) {
        grid_ = new SGPP::base::LinearBoundaryGrid(GridConfig.dim_);
      } else if (GridConfig.type_ == SGPP::base::GridType::ModLinear) {
        grid_ = new SGPP::base::ModLinearGrid(GridConfig.dim_);
      } else if (GridConfig.type_ == SGPP::base::GridType::Linear) {
        grid_ = new SGPP::base::LinearGrid(GridConfig.dim_);
      } else {
        grid_ = NULL;
        throw base::application_exception("LearnerBase::InitializeGrid: An unsupported grid type was chosen!");
      }

      // Generate regular Grid with LEVELS Levels
      SGPP::base::GridGenerator* myGenerator = grid_->createGridGenerator();
      myGenerator->regular(GridConfig.level_);
      delete myGenerator;

      // Create alpha
      alpha_ = new SGPP::base::DataVector(grid_->getSize());
      alpha_->setAll(0.0);
    }

    void LearnerBase::preProcessing() {
    }

    void LearnerBase::postProcessing(const SGPP::base::DataMatrix& trainDataset, const SGPP::solver::SLESolverType& solver,
                                     const size_t numNeededIterations) {
      if (this->isVerbose_) {
        std::cout << std::endl;
        std::cout << "Current Execution Time: " << execTime_ << std::endl;
        std::cout << std::endl;
      }
    }

    LearnerTiming LearnerBase::train(SGPP::base::DataMatrix& trainDataset, SGPP::base::DataVector& classes,
                                     const SGPP::base::RegularGridConfiguration& GridConfig,
                                     const SGPP::solver::SLESolverConfiguration& SolverConfigRefine,
                                     const SGPP::solver::SLESolverConfiguration& SolverConfigFinal,
                                     const SGPP::base::AdpativityConfiguration& AdaptConfig, const bool testAccDuringAdapt, const float_t lambdaRegularization) {
      LearnerTiming result;

      if (trainDataset.getNrows() != classes.getSize()) {
        throw base::application_exception("LearnerBase::train: length of classes vector does not match to dataset!");
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
      SGPP::datadriven::DMSystemMatrixBase* DMSystem = createDMSystem(trainDataset, lambdaRegularization);

      // check if System was created
      if (DMSystem == NULL)
        return result;

      SGPP::solver::SLESolver* myCG;

      if (SolverConfigRefine.type_ == SGPP::solver::SLESolverType::CG) {
        myCG = new SGPP::solver::ConjugateGradients(SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
      } else if (SolverConfigRefine.type_ == SGPP::solver::SLESolverType::BiCGSTAB) {
        myCG = new SGPP::solver::BiCGStab(SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
      } else {
        throw base::application_exception("LearnerBase::train: An unsupported SLE solver type was chosen!");
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
          SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(alpha_,
              AdaptConfig.noPoints_, AdaptConfig.threshold_);
          grid_->createGridGenerator()->refine(myRefineFunc);
          delete myRefineFunc;

          //tell the SLE manager that the grid changed (for interal data structures)
          DMSystem->prepareGrid();

          alpha_->resizeZero(grid_->getSize());
          float_t refineTime = myStopwatch2->stop();

          if (isVerbose_)
            std::cout << "New Grid Size: " << grid_->getSize() << " (Refinement took " << refineTime << " secs)"
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
          std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
          std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;
        }

        // use post-processing to determine Flops and time
        if (i < AdaptConfig.numRefinements_) {
          postProcessing(trainDataset, SolverConfigRefine.type_, myCG->getNumberIterations());
        } else {
          postProcessing(trainDataset, SolverConfigFinal.type_, myCG->getNumberIterations());
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
      delete myStopwatch2;
      delete myCG;
      delete DMSystem;

      return result;
    }

    LearnerTiming LearnerBase::train(SGPP::base::DataMatrix& trainDataset, SGPP::base::DataVector& classes,
                                     const SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverConfiguration& SolverConfig,
                                     const float_t lambdaRegularization) {
      SGPP::base::AdpativityConfiguration AdaptConfig;

      AdaptConfig.maxLevelType_ = false;
      AdaptConfig.noPoints_ = 0;
      AdaptConfig.numRefinements_ = 0;
      AdaptConfig.percent_ = 0.0;
      AdaptConfig.threshold_ = 0.0;

      return train(trainDataset, classes, GridConfig, SolverConfig, SolverConfig, AdaptConfig, false, lambdaRegularization);
    }

    void LearnerBase::predict(SGPP::base::DataMatrix& testDataset, SGPP::base::DataVector& classesComputed) {
      classesComputed.resize(testDataset.getNrows());

      SGPP::base::OperationMultipleEval* MultEval = SGPP::op_factory::createOperationMultipleEval(*grid_, testDataset);
      MultEval->mult(*alpha_, classesComputed);
      delete MultEval;
    }

    void LearnerBase::store(std::string tGridFilename, std::string tAlphaFilename) {
      throw base::application_exception("LearnerBase::store: This method isn't implemented, yet!");
    }

    float_t LearnerBase::getAccuracy(SGPP::base::DataMatrix& testDataset, const SGPP::base::DataVector& classesReference,
                                     const float_t threshold) {
      // evaluate test dataset

      SGPP::base::DataVector classesComputed(testDataset.getNrows());
      predict(testDataset, classesComputed);

      return getAccuracy(classesComputed, classesReference, threshold);
    }

    float_t LearnerBase::getAccuracy(const SGPP::base::DataVector& classesComputed,
                                     const SGPP::base::DataVector& classesReference, const float_t threshold) {
      float_t result = -1.0;

      if (classesComputed.getSize() != classesReference.getSize()) {
        throw base::application_exception("LearnerBase::getAccuracy: lengths of classes vectors do not match!");
      }

      if (isRegression_) {
        SGPP::base::DataVector tmp(classesComputed);
        tmp.sub(classesReference);
        tmp.sqr();
        result = tmp.sum();
        result /= static_cast<float_t>(tmp.getSize());
      } else {
        size_t correct = 0;

        for (size_t i = 0; i < classesComputed.getSize(); i++) {
          if ((classesComputed.get(i) >= threshold && classesReference.get(i) >= 0.0)
              || (classesComputed.get(i) < threshold && classesReference.get(i) < 0.0)) {
            correct++;
          }
        }

        result = static_cast<float_t>(correct) / static_cast<float_t>(classesComputed.getSize());
      }

      return result;
    }

    ClassificatorQuality LearnerBase::getCassificatorQuality(SGPP::base::DataMatrix& testDataset,
        const SGPP::base::DataVector& classesReference, const float_t threshold) {
      // evaluate test dataset
      SGPP::base::DataVector classesComputed(testDataset.getNrows());
      predict(testDataset, classesComputed);

      return getCassificatorQuality(classesComputed, classesReference, threshold);
    }

    ClassificatorQuality LearnerBase::getCassificatorQuality(const SGPP::base::DataVector& classesComputed,
        const SGPP::base::DataVector& classesReference, const float_t threshold) {
      ClassificatorQuality result;

      if (isRegression_) {
        throw base::application_exception(
          "LearnerBase::getCassificatorQuality: this method is not valid for regression problems!");
      }

      if (classesComputed.getSize() != classesReference.getSize()) {
        throw base::application_exception(
          "LearnerBase::getCassificatorQuality: lengths of classes vectors do not match!");
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
        } else { // ( (classesComputed.get(i) < threshold && classesReference.get(i) >= 0) )
          result.falseNegative_++;
        }
      }

      return result;
    }

    void LearnerBase::dumpGrid(std::string tFilename) {
      if (isTrained_) {
        SGPP::base::GridPrinter myPlotter(*grid_);
        myPlotter.printSparseGrid(*alpha_, tFilename, false);
      }
    }

    void LearnerBase::dumpFunction(std::string tFilename, size_t resolution) {
      if (isTrained_ && grid_->getStorage()->dim() <= 2) {
        SGPP::base::GridPrinter myPlotter(*grid_);
        myPlotter.printGrid(*alpha_, tFilename, resolution);
      }
    }

    bool LearnerBase::getIsRegression() const {
      return isRegression_;
    }

    bool LearnerBase::getIsVerbose() const {
      return isVerbose_;
    }

    void LearnerBase::setIsVerbose(const bool isVerbose) {
      isVerbose_ = isVerbose;
    }

    std::vector<std::pair<size_t, float_t> > LearnerBase::getRefinementExecTimes() {
      return this->ExecTimeOnStep;
    }

    std::shared_ptr<base::Grid> LearnerBase::getGridCopy() {
      if (this->grid_ == nullptr) {
        throw;
      }

      return std::shared_ptr<base::Grid>(base::Grid::unserialize(this->grid_->serialize()));
    }

    std::shared_ptr<base::DataVector> LearnerBase::getAlphaCopy() {
      if (this->alpha_ == nullptr) {
        throw;
      }

      return std::shared_ptr<base::DataVector>(new base::DataVector(*this->alpha_));
    }

  }

}
