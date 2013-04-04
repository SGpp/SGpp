#include "LearnerDensityBased.hpp"
#include "base/exception/application_exception.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "solver/sle/BiCGStab.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "datadriven/algorithm/DensitySystemMatrix.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "pde/operation/PdeOpFactory.hpp"
#include "base/exception/data_exception.hpp"
#include <limits>

namespace sg {

  namespace datadriven {

    LearnerDensityBased::LearnerDensityBased(sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression,
        const bool isVerbose) :
      LearnerBase(isRegression, isVerbose), CMode_(regularization) {
      this->withPrior = true;
      this->nrClasses = 0;
    }

    LearnerDensityBased::~LearnerDensityBased() {
      if (!CVec_.empty()) {
        for (size_t i = 0; i < CVec_.size(); i++)
          delete CVec_[i];
      }

      CVec_.clear();

      if (!gridVec_.empty()) {
        for (size_t i = 0; i < gridVec_.size(); i++)
          delete gridVec_[i];
      }

      gridVec_.clear();
      alphaVec_.clear();
    }

    sg::datadriven::DMSystemMatrixBase* LearnerDensityBased::createDMSystem(
      sg::base::DataMatrix& trainDataset, double lambda) {
      // Is not used
      return NULL;
    }

    void LearnerDensityBased::InitializeGrid(const sg::base::RegularGridConfiguration& GridConfig) {
      if (this->nrClasses == 0) {
        throw base::application_exception("LearnerDensityBased::InitializeGrid: Number of classes not set!");
      }

      if (!gridVec_.empty()) {
        for (size_t i = 0; i < gridVec_.size(); i++)
          delete gridVec_[i];
      }

      gridVec_.clear();

      if (!alphaVec_.empty()) {
        alphaVec_.clear();
      }

      for (size_t i = 0; i < nrClasses; i++) {
        if (GridConfig.type_ == sg::base::LinearTrapezoidBoundary) {
          gridVec_.push_back(new sg::base::LinearTrapezoidBoundaryGrid(GridConfig.dim_));
        } else if (GridConfig.type_ == sg::base::ModLinear) {
          gridVec_.push_back(new sg::base::ModLinearGrid(GridConfig.dim_));
        } else if (GridConfig.type_ == sg::base::Linear) {
          gridVec_.push_back(new sg::base::LinearGrid(GridConfig.dim_));
        } else {
          gridVec_.push_back(NULL);
          throw base::application_exception("LearnerDensityBased::InitializeGrid: An unsupported grid type was chosen!");
        }

        // Generate regular Grid with LEVELS Levels
        sg::base::GridGenerator* myGenerator = gridVec_[i]->createGridGenerator();
        myGenerator->regular(GridConfig.level_);
        delete myGenerator;

        // Create alpha
        alphaVec_.push_back(sg::base::DataVector(gridVec_[i]->getSize()));
        alphaVec_[i].setAll(0.0);
      }
    }

    // LearnerTiming LearnerDensityBased::train(sg::base::DataMatrix& trainDataset,
    //    sg::base::DataVector& classes,
    //    const sg::base::RegularGridConfiguration& GridConfig,
    //    const sg::solver::SLESolverConfiguration& SolverConfigRefine,
    //    const sg::solver::SLESolverConfiguration& SolverConfigFinal,
    //    const sg::base::AdpativityConfiguration& AdaptConfig,
    //    const bool testAccDuringAdapt, const double lambda)
    // {
    //   return train(trainDataset, classes, GridConfig, SolverConfigRefine, SolverConfigFinal, AdaptConfig, testAccDuringAdapt, lambda, true);
    // }


    LearnerTiming LearnerDensityBased::train(sg::base::DataMatrix& trainDataset,
        sg::base::DataVector& classes,
        const sg::base::RegularGridConfiguration& GridConfig,
        const sg::solver::SLESolverConfiguration& SolverConfigRefine,
        const sg::solver::SLESolverConfiguration& SolverConfigFinal,
        const sg::base::AdpativityConfiguration& AdaptConfig,
        const bool testAccDuringAdapt, const double lambda) {
      LearnerTiming result;

      if (trainDataset.getNrows() != classes.getSize()) {
        throw base::application_exception(
          "LearnerBase::train: length of classes vector does not match to dataset!");
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

      sg::solver::SLESolver* myCG;

      if (SolverConfigRefine.type_ == sg::solver::CG) {
        myCG = new sg::solver::ConjugateGradients(
          SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
      } else if (SolverConfigRefine.type_ == sg::solver::BiCGSTAB) {
        myCG = new sg::solver::BiCGStab(SolverConfigRefine.maxIterations_,
                                        SolverConfigRefine.eps_);
      } else {
        throw base::application_exception(
          "LearnerBaseSP::train: An unsupported SLE solver type was chosen!");
      }

      // Pre-Procession
      //preProcessing();

      if (isVerbose_)
        std::cout << "Starting Learning...." << std::endl;

      sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();

      int dim = (int)trainDataset.getNcols();

      //Compute all occurring class labels and how many data points exist per label:
      std::map<double, int> entriesPerClass;

      for (unsigned int i = 0; i < classes.getSize(); i++) {
        double classNum = classes.get(i);

        if (entriesPerClass.find(classNum) == entriesPerClass.end()) {
          entriesPerClass.insert(std::pair<double, int>(classNum, 1));
        } else {
          entriesPerClass[classNum]++;
        }
      }

      //Create an empty matrix for every class:
      std::vector<sg::base::DataMatrix> trainDataClasses;
      std::map<double, int> class_indeces; //Maps class numbers to indices

      std::map<double, int>::iterator it;
      int index = 0;

      for (it = entriesPerClass.begin(); it != entriesPerClass.end(); it++) {
        sg::base::DataMatrix m(0/*(*it).second*/, dim);
        trainDataClasses.push_back(m);
        class_indeces[(*it).first] = index;
        index_to_class_.insert(std::pair<int, double>(index, (*it).first));

        //compute prior
        if (this->withPrior) {
          prior.push_back((*it).second / (double)classes.getSize());
        } else {
          prior.push_back(1.);
        }

        index++;
      }

      //Split the data into the different classes:
      for (unsigned int i = 0; i < trainDataset.getNrows(); i++) {
        double classLabel = classes[i];
        sg::base::DataVector vec(dim);
        trainDataset.getRow(i, vec);
        trainDataClasses[class_indeces[classLabel]].appendRow(vec);
      }

      // Construct Grid
      if (!alphaVec_.empty())
        alphaVec_.clear();

      if (grid_ != NULL)
        delete grid_;

      if (isTrained_ == true)
        isTrained_ = false;

      setNrClasses(class_indeces.size());
      InitializeGrid(GridConfig);


      // check if grids were created
      if (gridVec_.size() != class_indeces.size())
        return result;

      for (size_t i = 0; i < AdaptConfig.numRefinements_ + 1; i++) {
        if (isVerbose_)
          std::cout << std::endl << "Doing refinement: " << i << std::endl;

        myStopwatch->start();

        // Do Refinements
        if (i > 0) {
          for (size_t ii = 0; ii < alphaVec_.size(); ii++) {
            sg::base::SurplusRefinementFunctor* myRefineFunc =
              new sg::base::SurplusRefinementFunctor(&alphaVec_[ii], AdaptConfig.noPoints_);
            gridVec_[ii]->createGridGenerator()->refine(myRefineFunc);
            delete myRefineFunc;

            //DMSystem->rebuildLevelAndIndex();   not implemented

            if (isVerbose_)
              std::cout << "New Grid Size[" << ii << "]: " << gridVec_[ii]->getSize() << std::endl;
          }
        } else {
          for (size_t ii = 0; ii < gridVec_.size(); ii++) {
            if (isVerbose_)
              std::cout << "Grid Size[" << ii << "]: " << gridVec_[ii]->getSize() << std::endl;
          }
        }


        //Create regularization operator
        if (!CVec_.empty()) {
          for (size_t ii = 0; ii < CVec_.size(); ii++)
            delete CVec_[ii];

          CVec_.clear();
        }

        for (size_t ii = 0; ii < class_indeces.size(); ii++) {
          if (this->CMode_ == Laplace) {
            CVec_.push_back(sg::op_factory::createOperationLaplace(*this->gridVec_[ii]));
          } else if (this->CMode_ == Identity) {
            CVec_.push_back(sg::op_factory::createOperationIdentity(*this->gridVec_[ii]));
          } else {
            throw base::application_exception("LearnerDensityBased::train: Unknown regularization operator");
          }
        }

        //Solve the system for every class and store coefficients:
        for (size_t ii = 0; ii < trainDataClasses.size(); ii++) {
          sg::datadriven::DensitySystemMatrix DMatrix(*gridVec_[ii], trainDataClasses[ii], *(CVec_[ii]), lambda);
          sg::base::DataVector rhs(gridVec_[ii]->getStorage()->size());
          sg::base::DataVector alpha(gridVec_[ii]->getStorage()->size());
          alpha.setAll(0.0);
          DMatrix.generateb(rhs);

          if (i == AdaptConfig.numRefinements_) {
            myCG->setMaxIterations(SolverConfigFinal.maxIterations_);
            myCG->setEpsilon(SolverConfigFinal.eps_);
          }

          myCG->solve(DMatrix, alpha, rhs, false, false, -1.0);
          alphaVec_[ii].resize(alpha.getSize());
          alphaVec_[ii].copyFrom(alpha);
        }

        execTime_ += myStopwatch->stop();

        if (isVerbose_) {
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

        /*double tmp1, tmp2, tmp3, tmp4;
          DMSystem->getTimers(tmp1, tmp2, tmp3, tmp4);
          result.timeComplete_ = execTime_;
          result.timeMultComplete_ = tmp1;
          result.timeMultCompute_ = tmp2;
          result.timeMultTransComplete_ = tmp3;
          result.timeMultTransCompute_ = tmp4;*/
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
      delete myCG;
      //delete DMSystem;

      return result;
    }

    sg::base::DataVector LearnerDensityBased::predict(
      sg::base::DataMatrix& testDataset) {
      if (isTrained_) {
        sg::base::DataVector result(testDataset.getNrows());

        for (unsigned int i = 0; i < testDataset.getNrows(); i++) {
          sg::base::DataVector p(testDataset.getNcols());
          testDataset.getRow(i, p);

          //Compute maximum of all density functions:
          std::vector<sg::base::DataVector>::iterator it;
          int max_index = -1;
          double max = std::numeric_limits<double>::min();
          int class_index = 0;

          for (it = alphaVec_.begin(); it != alphaVec_.end(); it++) {
            sg::base::OperationEval* Eval = sg::op_factory::createOperationEval(*gridVec_[class_index]);
            sg::base::DataVector alpha = *it;
            //posterior = likelihood*prior
            double res = Eval->eval(alpha, p) * this->prior[class_index];

            if (res > max) {
              max = res;
              max_index = class_index;
            }

            class_index++;
            delete Eval;
          }

          result[i] = index_to_class_[max_index];
        }

        return result;
      } else {
        throw base::data_exception ("Cannot predict. Learner has to be trained first!");
      }

    }

    time_t LearnerDensityBased::getExecTime() {
      return (time_t)this->execTime_;
    }

    size_t LearnerDensityBased::getNrGridPoints() {
      size_t maxGrid = 0;

      for (size_t i = 0; i < gridVec_.size(); i++) {
        if (gridVec_[i]->getSize() > maxGrid)
          maxGrid = gridVec_[i]->getSize();
      }

      return maxGrid;
    }

  }
}

