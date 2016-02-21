// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerDensityBased.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <limits>
#include <map>
#include <utility>
#include <vector>


namespace SGPP {
namespace datadriven {

LearnerDensityBased::LearnerDensityBased(SGPP::datadriven::RegularizationType&
    regularization, const bool isRegression,
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

SGPP::datadriven::DMSystemMatrixBase* LearnerDensityBased::createDMSystem(
  SGPP::base::DataMatrix& trainDataset, float_t lambda) {
  // Is not used
  return NULL;
}

void LearnerDensityBased::InitializeGrid(const
    SGPP::base::RegularGridConfiguration& GridConfig) {
  if (this->nrClasses == 0) {
    throw base::application_exception("LearnerDensityBased::InitializeGrid: Number of classes not "
      "set!");
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
    if (GridConfig.type_ == SGPP::base::GridType::LinearBoundary) {
      gridVec_.push_back(new SGPP::base::LinearBoundaryGrid(GridConfig.dim_));
    } else if (GridConfig.type_ == SGPP::base::GridType::ModLinear) {
      gridVec_.push_back(new SGPP::base::ModLinearGrid(GridConfig.dim_));
    } else if (GridConfig.type_ == SGPP::base::GridType::Linear) {
      gridVec_.push_back(new SGPP::base::LinearGrid(GridConfig.dim_));
    } else {
      gridVec_.push_back(NULL);
      throw base::application_exception("LearnerDensityBased::InitializeGrid: An unsupported grid "
        "type was chosen!");
    }

    // Generate regular Grid with LEVELS Levels
    gridVec_[i]->getGenerator().regular(GridConfig.level_);

    // Create alpha
    alphaVec_.push_back(SGPP::base::DataVector(gridVec_[i]->getSize()));
    alphaVec_[i].setAll(0.0);
  }
}

// LearnerTiming LearnerDensityBased::train(SGPP::base::DataMatrix& trainDataset,
//    SGPP::base::DataVector& classes,
//    const SGPP::base::RegularGridConfiguration& GridConfig,
//    const SGPP::solver::SLESolverConfiguration& SolverConfigRefine,
//    const SGPP::solver::SLESolverConfiguration& SolverConfigFinal,
//    const SGPP::base::AdpativityConfiguration& AdaptConfig,
//    const bool testAccDuringAdapt, const float_t lambda)
// {
//   return train(trainDataset, classes, GridConfig, SolverConfigRefine, SolverConfigFinal,
//     AdaptConfig, testAccDuringAdapt, lambda, true);
// }


LearnerTiming LearnerDensityBased::train(SGPP::base::DataMatrix& trainDataset,
    SGPP::base::DataVector& classes,
    const SGPP::base::RegularGridConfiguration& GridConfig,
    const SGPP::solver::SLESolverConfiguration& SolverConfigRefine,
    const SGPP::solver::SLESolverConfiguration& SolverConfigFinal,
    const SGPP::base::AdpativityConfiguration& AdaptConfig,
    const bool testAccDuringAdapt, const float_t lambda) {
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

  float_t oldAcc = 0.0;

  SGPP::solver::SLESolver* myCG;

  if (SolverConfigRefine.type_ == SGPP::solver::SLESolverType::CG) {
    myCG = new SGPP::solver::ConjugateGradients(
      SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
  } else if (SolverConfigRefine.type_ == SGPP::solver::SLESolverType::BiCGSTAB) {
    myCG = new SGPP::solver::BiCGStab(SolverConfigRefine.maxIterations_,
                                      SolverConfigRefine.eps_);
  } else {
    throw base::application_exception(
      "LearnerDensityBased::train: An unsupported SLE solver type was chosen!");
  }

  // Pre-Procession
  // preProcessing();

  if (isVerbose_)
    std::cout << "Starting Learning...." << std::endl;

  SGPP::base::SGppStopwatch* myStopwatch = new SGPP::base::SGppStopwatch();

  int dim = static_cast<int>(trainDataset.getNcols());

  // Compute all occurring class labels and how many data points exist per label:
  std::map<float_t, int> entriesPerClass;

  for (unsigned int i = 0; i < classes.getSize(); i++) {
    float_t classNum = classes.get(i);

    if (entriesPerClass.find(classNum) == entriesPerClass.end()) {
      entriesPerClass.insert(std::pair<float_t, int>(classNum, 1));
    } else {
      entriesPerClass[classNum]++;
    }
  }

  // Create an empty matrix for every class:
  std::vector<SGPP::base::DataMatrix> trainDataClasses;
  std::map<float_t, int> class_indeces;  // Maps class numbers to indices

  std::map<float_t, int>::iterator it;
  int index = 0;

  for (it = entriesPerClass.begin(); it != entriesPerClass.end(); it++) {
    SGPP::base::DataMatrix m(0/*(*it).second*/, dim);
    trainDataClasses.push_back(m);
    class_indeces[(*it).first] = index;
    index_to_class_.insert(std::pair<int, float_t>(index, (*it).first));

    // compute prior
    if (this->withPrior) {
      prior.push_back(static_cast<float_t>((*it).second) /
                      static_cast<float_t>(classes.getSize()));
    } else {
      prior.push_back(1.);
    }

    index++;
  }

  // Split the data into the different classes:
  for (unsigned int i = 0; i < trainDataset.getNrows(); i++) {
    float_t classLabel = classes[i];
    SGPP::base::DataVector vec(dim);
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
        SGPP::base::SurplusRefinementFunctor myRefineFunc(alphaVec_[ii], AdaptConfig.noPoints_);
        gridVec_[ii]->getGenerator().refine(myRefineFunc);

        // DMSystem->rebuildLevelAndIndex();   not implemented

        if (isVerbose_)
          std::cout << "New Grid Size[" << ii << "]: " << gridVec_[ii]->getSize() <<
                    std::endl;
      }
    } else {
      for (size_t ii = 0; ii < gridVec_.size(); ii++) {
        if (isVerbose_)
          std::cout << "Grid Size[" << ii << "]: " << gridVec_[ii]->getSize() <<
                    std::endl;
      }
    }


    // Create regularization operator
    if (!CVec_.empty()) {
      for (size_t ii = 0; ii < CVec_.size(); ii++)
        delete CVec_[ii];

      CVec_.clear();
    }

    for (size_t ii = 0; ii < class_indeces.size(); ii++) {
      if (this->CMode_ == SGPP::datadriven::RegularizationType::Laplace) {
        CVec_.push_back(SGPP::op_factory::createOperationLaplace(*this->gridVec_[ii]).release());
      } else if (this->CMode_ == SGPP::datadriven::RegularizationType::Identity) {
        CVec_.push_back(SGPP::op_factory::createOperationIdentity(*this->gridVec_[ii]).release());
      } else {
        throw base::application_exception("LearnerDensityBased::train: Unknown regularization "
          "operator");
      }
    }

    // Solve the system for every class and store coefficients:
    for (size_t ii = 0; ii < trainDataClasses.size(); ii++) {
      SGPP::datadriven::DensitySystemMatrix DMatrix(*gridVec_[ii],
          trainDataClasses[ii], *(CVec_[ii]), lambda);
      SGPP::base::DataVector rhs(gridVec_[ii]->getStorage().size());
      SGPP::base::DataVector alpha(gridVec_[ii]->getStorage().size());
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

    /*float_t tmp1, tmp2, tmp3, tmp4;
      DMSystem->getTimers(tmp1, tmp2, tmp3, tmp4);
      result.timeComplete_ = execTime_;
      result.timeMultComplete_ = tmp1;
      result.timeMultCompute_ = tmp2;
      result.timeMultTransComplete_ = tmp3;
      result.timeMultTransCompute_ = tmp4;*/
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
  delete myCG;
  // delete DMSystem;

  return result;
}

SGPP::base::DataVector LearnerDensityBased::predict(
  SGPP::base::DataMatrix& testDataset) {
  if (isTrained_) {
    SGPP::base::DataVector result(testDataset.getNrows());

    for (unsigned int i = 0; i < testDataset.getNrows(); i++) {
      SGPP::base::DataVector p(testDataset.getNcols());
      testDataset.getRow(i, p);

      // Compute maximum of all density functions:
      std::vector<SGPP::base::DataVector>::iterator it;
      int max_index = -1;
      float_t max = std::numeric_limits<float_t>::min();
      int class_index = 0;

      for (it = alphaVec_.begin(); it != alphaVec_.end(); it++) {
        std::unique_ptr<SGPP::base::OperationEval> Eval(
            SGPP::op_factory::createOperationEval(*gridVec_[class_index]));
        SGPP::base::DataVector alpha = *it;
        // posterior = likelihood*prior
        float_t res = Eval->eval(alpha, p) * this->prior[class_index];

        if (res > max) {
          max = res;
          max_index = class_index;
        }

        class_index++;
      }

      result[i] = index_to_class_[max_index];
    }

    return result;
  } else {
    throw base::data_exception("Cannot predict. Learner has to be trained first!");
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

}  // namespace datadriven
}  // namespace SGPP

