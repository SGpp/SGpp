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
#include "base/exception/application_exception.hpp"
#include "base/tools/GridPrinter.hpp"

#include "solver/sle/ConjugateGradients.hpp"
#include "solver/sle/BiCGStab.hpp"

#include "datadriven/application/LearnerBase.hpp"

#include <iostream>

#include "parallel/tools/MPI/SGppMPITools.hpp"

namespace sg {

  namespace datadriven {

    LearnerBase::LearnerBase(const bool isRegression, const bool isVerbose)
      : alpha_(NULL), grid_(NULL), isVerbose_(isVerbose), isRegression_(isRegression), isTrained_(false), execTime_(0.0), GFlop_(0.0), GByte_(0.0) {
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

    LearnerBase::LearnerBase(const std::string tGridFilename, const std::string tAlphaFilename, const bool isRegression, const bool isVerbose)
      : alpha_(NULL), grid_(NULL), isVerbose_(isVerbose), isRegression_(isRegression), isTrained_(false), execTime_(0.0), GFlop_(0.0), GByte_(0.0) {
      // @TODO (heinecke) implement
      throw base::application_exception("LearnerBase::LearnerBase: This construct isn't implemented, yet!");
    }

    LearnerBase::LearnerBase(const LearnerBase& copyMe) {
      // safety, should not happen
      if (alpha_ != NULL)
        delete alpha_;

      if (grid_ != NULL)
        delete grid_;

      // @TODO (heinecke) grid copy constructor
      grid_ = sg::base::Grid::unserialize(copyMe.grid_->serialize());
      alpha_ = new sg::base::DataVector(*(copyMe.alpha_));
    }

    LearnerBase::~LearnerBase() {
      // if user does no cleaning
      if (alpha_ != NULL)
        delete alpha_;

      if (grid_ != NULL)
        delete grid_;
    }

    void LearnerBase::InitializeGrid(const sg::base::RegularGridConfiguration& GridConfig) {
      if (GridConfig.type_ == sg::base::LinearTrapezoidBoundary) {
        grid_ = new sg::base::LinearTrapezoidBoundaryGrid(GridConfig.dim_);
      } else if (GridConfig.type_ == sg::base::ModLinear) {
        grid_ = new sg::base::ModLinearGrid(GridConfig.dim_);
      } else if (GridConfig.type_ == sg::base::Linear) {
        grid_ = new sg::base::LinearGrid(GridConfig.dim_);
      } else {
        grid_ = NULL;
        throw base::application_exception("LearnerBase::InitializeGrid: An unsupported grid type was chosen!");
      }

      // Generate regular Grid with LEVELS Levels
      sg::base::GridGenerator* myGenerator = grid_->createGridGenerator();
      myGenerator->regular(GridConfig.level_);
      delete myGenerator;

      // Create alpha
      alpha_ = new sg::base::DataVector(grid_->getSize());
      alpha_->setAll(0.0);
    }

    void LearnerBase::preProcessing() {
    }

    void LearnerBase::postProcessing(const sg::base::DataMatrix& trainDataset, const sg::solver::SLESolverType& solver,
                                     const size_t numNeededIterations) {
      if (this->isVerbose_) {
        std::cout << std::endl;
        std::cout << "Current Execution Time: " << execTime_ << std::endl;
        std::cout << std::endl;
      }
    }

    LearnerTiming LearnerBase::train(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes,
                                     const sg::base::RegularGridConfiguration& GridConfig, const sg::solver::SLESolverConfiguration& SolverConfigRefine,
                                     const sg::solver::SLESolverConfiguration& SolverConfigFinal, const sg::base::AdpativityConfiguration& AdaptConfig,
                                     const bool testAccDuringAdapt, const double lambda) {
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
      sg::datadriven::DMSystemMatrixBase* DMSystem = createDMSystem(trainDataset, lambda);

      // check if System was created
      if (DMSystem == NULL)
        return result;

      sg::solver::SLESolver* myCG;

      if (SolverConfigRefine.type_ == sg::solver::CG) {
        myCG = new sg::solver::ConjugateGradients(SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
      } else if (SolverConfigRefine.type_ == sg::solver::BiCGSTAB) {
        myCG = new sg::solver::BiCGStab(SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
      } else {
        throw base::application_exception("LearnerBaseSP::train: An unsupported SLE solver type was chosen!");
      }

      // Pre-Procession
      preProcessing();

      if (isVerbose_)
        std::cout << "Starting Learning...." << std::endl;

      // execute adaptsteps
      sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();
      sg::base::SGppStopwatch* myStopwatch2 = new sg::base::SGppStopwatch();

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
          myStopwatch2->start();

#ifdef USE_MPI

          if (parallel::myGlobalMPIComm->getMyRank() == 0) {
#endif
            sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(alpha_, AdaptConfig.noPoints_, AdaptConfig.threshold_);
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

          alpha_->resizeZero(grid_->getSize());
          double refineTime = myStopwatch2->stop();

          if (isVerbose_)
            std::cout << "New Grid Size: " << grid_->getSize() << " (Refinement took " << refineTime << " secs)" << std::endl;
        } else {
          if (isVerbose_)
            std::cout << "Grid Size: " << grid_->getSize() << std::endl;
        }

        sg::base::DataVector b(alpha_->getSize());
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
          std::cout << std::endl;
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
      delete myStopwatch2;
      delete myCG;
      delete DMSystem;

      return result;
    }

    LearnerTiming LearnerBase::train(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes,
                                     const sg::base::RegularGridConfiguration& GridConfig, const sg::solver::SLESolverConfiguration& SolverConfig,
                                     const double lambda) {
      sg::base::AdpativityConfiguration AdaptConfig;

      AdaptConfig.maxLevelType_ = false;
      AdaptConfig.noPoints_ = 0;
      AdaptConfig.numRefinements_ = 0;
      AdaptConfig.percent_ = 0.0;
      AdaptConfig.threshold_ = 0.0;

      return train(trainDataset, classes, GridConfig, SolverConfig, SolverConfig, AdaptConfig, false, lambda);
    }

    sg::base::DataVector LearnerBase::predict(sg::base::DataMatrix& testDataset) {
      sg::base::DataVector classesComputed(testDataset.getNrows());

      sg::base::OperationMultipleEval* MultEval = sg::op_factory::createOperationMultipleEval(*grid_, &testDataset);
      MultEval->mult(*alpha_, classesComputed);
      delete MultEval;

      return classesComputed;
    }

    void LearnerBase::store(std::string tGridFilename, std::string tAlphaFilename) {
      // @TODO (heinecke) implement
      throw base::application_exception("LearnerBase::store: This method isn't implemented, yet!");
    }

    double LearnerBase::getAccuracy(sg::base::DataMatrix& testDataset, const sg::base::DataVector& classesReference, const double threshold) {
      // evaluate test dataset
      sg::base::DataVector classesComputed = predict(testDataset);

      return getAccuracy(classesComputed, classesReference, threshold);
    }

    double LearnerBase::getAccuracy(const sg::base::DataVector& classesComputed, const sg::base::DataVector& classesReference, const double threshold) {
      double result = -1.0;

      if (classesComputed.getSize() != classesReference.getSize()) {
        throw base::application_exception("LearnerBase::getAccuracy: lengths of classes vectors do not match!");
      }

      if (isRegression_) {
        sg::base::DataVector tmp(classesComputed);
        tmp.sub(classesReference);
        tmp.sqr();
        result = tmp.sum();
        result /= static_cast<double>(tmp.getSize());
      } else {
        size_t correct = 0;

        for (size_t i = 0; i < classesComputed.getSize(); i++) {
          if ( (classesComputed.get(i) >= threshold && classesReference.get(i) >= 0.0)
               || (classesComputed.get(i) < -threshold && classesReference.get(i) < 0.0) ) {
            correct++;
            //std::cout << classesComputed.get(i) << " " <<classesReference.get(i) << std::endl;
          }
        }

        result = static_cast<double>(correct) / static_cast<double>(classesComputed.getSize());
      }

      return result;
    }

    ClassificatorQuality LearnerBase::getCassificatorQuality(sg::base::DataMatrix& testDataset, const sg::base::DataVector& classesReference, const double threshold) {
      // evaluate test dataset
      sg::base::DataVector classesComputed = predict(testDataset);

      return getCassificatorQuality(classesComputed, classesReference, threshold);
    }

    ClassificatorQuality LearnerBase::getCassificatorQuality(const sg::base::DataVector& classesComputed, const sg::base::DataVector& classesReference, const double threshold) {
      ClassificatorQuality result;

      if (isRegression_) {
        throw base::application_exception("LearnerBase::getCassificatorQuality: this method is not valid for regression problems!");
      }

      if (classesComputed.getSize() != classesReference.getSize()) {
        throw base::application_exception("LearnerBase::getCassificatorQuality: lengths of classes vectors do not match!");
      }

      result.truePositive_ = 0;
      result.trueNegative_ = 0;
      result.falsePositive_ = 0;
      result.falseNegative_ = 0;

      for (size_t i = 0; i < classesComputed.getSize(); i++) {
        if ( (classesComputed.get(i) >= threshold && classesReference.get(i) >= 0.0) ) {
          result.truePositive_++;
        } else if ( (classesComputed.get(i) < threshold && classesReference.get(i) < 0.0) ) {
          result.trueNegative_++;
        } else if ( (classesComputed.get(i) >= threshold && classesReference.get(i) < 0.0) ) {
          result.falsePositive_++;
        } else { // ( (classesComputed.get(i) < threshold && classesReference.get(i) >= 0) )
          result.falseNegative_++;
        }
      }

      return result;
    }

    void LearnerBase::dumpGrid(std::string tFilename) {
      if (isTrained_) {
        sg::base::GridPrinter myPlotter(*grid_);
        myPlotter.printSparseGrid(*alpha_, tFilename, false);
      }
    }

    void LearnerBase::dumpFunction(std::string tFilename, size_t resolution) {
      if (isTrained_ && grid_->getStorage()->dim() <= 2) {
        sg::base::GridPrinter myPlotter(*grid_);
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

  }

}
