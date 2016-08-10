// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/datadriven/application/LearnerDensityBased.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/globaldef.hpp>

#include <limits>
#include <map>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

LearnerDensityBased::LearnerDensityBased(datadriven::RegularizationType& regularization,
                                         const bool isRegression, const bool isVerbose)
    : LearnerBase(isRegression, isVerbose), CMode(regularization) {
  this->withPrior = true;
  this->nrClasses = 0;
}

LearnerDensityBased::~LearnerDensityBased() {
  CVec.clear();
  gridVec.clear();
  alphaVec.clear();
}

std::unique_ptr<datadriven::DMSystemMatrixBase> LearnerDensityBased::createDMSystem(
    base::DataMatrix& trainDataset, double lambda) {
  // Is not used
  return nullptr;
}

void LearnerDensityBased::InitializeGrid(const base::RegularGridConfiguration& GridConfig) {
  if (this->nrClasses == 0) {
    throw base::application_exception(
        "LearnerDensityBased::InitializeGrid: Number of classes not "
        "set!");
  }

  gridVec.clear();
  alphaVec.clear();

  for (size_t i = 0; i < nrClasses; i++) {
    if (GridConfig.type_ == base::GridType::LinearBoundary) {
      gridVec.push_back(std::make_unique<base::LinearBoundaryGrid>(GridConfig.dim_));
    } else if (GridConfig.type_ == base::GridType::ModLinear) {
      gridVec.push_back(std::make_unique<base::ModLinearGrid>(GridConfig.dim_));
    } else if (GridConfig.type_ == base::GridType::Linear) {
      gridVec.push_back(std::make_unique<base::LinearGrid>(GridConfig.dim_));
    } else {
      gridVec.push_back(nullptr);
      throw base::application_exception(
          "LearnerDensityBased::InitializeGrid: An unsupported grid "
          "type was chosen!");
    }

    // Generate regular Grid with LEVELS Levels
    gridVec[i]->getGenerator().regular(GridConfig.level_);

    // Create alpha
    alphaVec.push_back(base::DataVector(gridVec[i]->getSize()));
    alphaVec[i].setAll(0.0);
  }
}

LearnerTiming LearnerDensityBased::train(base::DataMatrix& trainDataset, base::DataVector& classes,
                                         const base::RegularGridConfiguration& GridConfig,
                                         const solver::SLESolverConfiguration& SolverConfigRefine,
                                         const solver::SLESolverConfiguration& SolverConfigFinal,
                                         const base::AdpativityConfiguration& AdaptConfig,
                                         const bool testAccDuringAdapt, const double lambda) {
  LearnerTiming timing = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  if (trainDataset.getNrows() != classes.getSize()) {
    throw base::application_exception(
        "LearnerBase::train: length of classes vector does not match to dataset!");
  }

  double oldAcc = 0.0;

  std::unique_ptr<solver::SLESolver> myCG;

  if (SolverConfigRefine.type_ == solver::SLESolverType::CG) {
    myCG = std::make_unique<solver::ConjugateGradients>(SolverConfigRefine.maxIterations_,
                                                        SolverConfigRefine.eps_);
  } else if (SolverConfigRefine.type_ == solver::SLESolverType::BiCGSTAB) {
    myCG = std::make_unique<solver::BiCGStab>(SolverConfigRefine.maxIterations_,
                                              SolverConfigRefine.eps_);
  } else {
    throw base::application_exception(
        "LearnerDensityBased::train: An unsupported SLE solver type was chosen!");
  }

  if (isVerbose) std::cout << "Starting Learning...." << std::endl;

  auto myStopwatch = std::make_unique<base::SGppStopwatch>();

  size_t dim = trainDataset.getNcols();

  // Compute all occurring class labels and how many data points exist per label:
  std::map<double, size_t> entriesPerClass;

  for (size_t i = 0; i < classes.getSize(); i++) {
    double classNum = classes.get(i);

    if (entriesPerClass.find(classNum) == entriesPerClass.end()) {
      entriesPerClass.insert(std::pair<double, int>(classNum, 1));
    } else {
      entriesPerClass[classNum]++;
    }
  }

  // Create an empty matrix for every class:
  std::vector<base::DataMatrix> trainDataClasses;
  std::map<double, size_t> classIndeces;  // Maps class numbers to indices

  std::map<double, size_t>::iterator it;
  int index = 0;

  for (it = entriesPerClass.begin(); it != entriesPerClass.end(); it++) {
    base::DataMatrix m(0 /*(*it).second*/, dim);
    trainDataClasses.push_back(m);
    classIndeces[(*it).first] = index;
    indexToClass.insert(std::pair<size_t, double>(index, (*it).first));

    // compute prior
    if (this->withPrior) {
      prior.push_back(static_cast<double>((*it).second) / static_cast<double>(classes.getSize()));
    } else {
      prior.push_back(1.);
    }

    index++;
  }

  // Split the data into the different classes:
  for (size_t i = 0; i < trainDataset.getNrows(); i++) {
    double classLabel = classes[i];
    base::DataVector vec(dim);
    trainDataset.getRow(i, vec);
    trainDataClasses[classIndeces[classLabel]].appendRow(vec);
  }

  alphaVec.clear();

  if (isTrained == true) isTrained = false;

  setNrClasses(classIndeces.size());
  InitializeGrid(GridConfig);

  // check if grids were created
  if (gridVec.size() != classIndeces.size()) return timing;

  for (size_t i = 0; i < AdaptConfig.numRefinements_ + 1; i++) {
    if (isVerbose) std::cout << std::endl << "Doing refinement: " << i << std::endl;

    myStopwatch->start();

    // Do Refinements
    if (i > 0) {
      for (size_t ii = 0; ii < alphaVec.size(); ii++) {
        base::SurplusRefinementFunctor myRefineFunc(alphaVec[ii], AdaptConfig.noPoints_);
        gridVec[ii]->getGenerator().refine(myRefineFunc);

        // DMSystem->rebuildLevelAndIndex();   not implemented

        if (isVerbose)
          std::cout << "New Grid Size[" << ii << "]: " << gridVec[ii]->getSize() << std::endl;
      }
    } else {
      for (size_t ii = 0; ii < gridVec.size(); ii++) {
        if (isVerbose)
          std::cout << "Grid Size[" << ii << "]: " << gridVec[ii]->getSize() << std::endl;
      }
    }

    // Create regularization operator
    CVec.clear();

    for (size_t j = 0; j < classIndeces.size(); j++) {
      if (CMode == datadriven::RegularizationType::Laplace) {
        CVec.push_back(std::unique_ptr<base::OperationMatrix>(
            op_factory::createOperationLaplace(*this->gridVec[j])));
      } else if (CMode == datadriven::RegularizationType::Identity) {
        CVec.push_back(std::unique_ptr<base::OperationMatrix>(
            op_factory::createOperationIdentity(*this->gridVec[j])));
      } else {
        throw base::application_exception(
            "LearnerDensityBased::train: Unknown regularization "
            "operator");
      }
    }

    // Solve the system for every class and store coefficients:
    for (size_t j = 0; j < trainDataClasses.size(); j++) {
      datadriven::DensitySystemMatrix DMatrix(*gridVec[j], trainDataClasses[j], CVec[j].release(),
                                              lambda);
      base::DataVector rhs(gridVec[j]->getSize());
      base::DataVector alpha(gridVec[j]->getSize());
      alpha.setAll(0.0);
      DMatrix.generateb(rhs);

      if (i == AdaptConfig.numRefinements_) {
        myCG->setMaxIterations(SolverConfigFinal.maxIterations_);
        myCG->setEpsilon(SolverConfigFinal.eps_);
      }

      myCG->solve(DMatrix, alpha, rhs, false, false, -1.0);
      alphaVec[j].resize(alpha.getSize());
      alphaVec[j].copyFrom(alpha);
    }

    execTime += myStopwatch->stop();

    if (isVerbose) {
      std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
      std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;
    }

    // use post-processing to determine Flops and time
    if (i < AdaptConfig.numRefinements_) {
      postProcessing(trainDataset, SolverConfigRefine.type_, myCG->getNumberIterations());
    } else {
      postProcessing(trainDataset, SolverConfigFinal.type_, myCG->getNumberIterations());
    }

    if (testAccDuringAdapt) {
      double acc = getAccuracy(trainDataset, classes);

      if (isVerbose) {
        if (isRegression) {
          if (isVerbose) std::cout << "MSE (train): " << acc << std::endl;
        } else {
          if (isVerbose) std::cout << "Acc (train): " << acc << std::endl;
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

  timing.timeComplete_ = execTime;
  timing.GFlop_ = GFlop;
  timing.GByte_ = GByte;

  return timing;
}

base::DataVector LearnerDensityBased::predict(base::DataMatrix& testDataset) {
  base::DataVector result(testDataset.getNrows());
  if (isTrained) {
    for (size_t i = 0; i < testDataset.getNrows(); i++) {
      base::DataVector p(testDataset.getNcols());
      testDataset.getRow(i, p);

      // Compute maximum of all density functions:
      std::vector<base::DataVector>::iterator it;
      int maxIndex = -1;
      double max = std::numeric_limits<double>::min();
      int class_index = 0;

      for (auto alpha : alphaVec) {
        std::unique_ptr<base::OperationEval> Eval(
            op_factory::createOperationEval(*gridVec[class_index]));
        // posterior = likelihood*prior
        double res = Eval->eval(alpha, p) * this->prior[class_index];

        if (res > max) {
          max = res;
          maxIndex = class_index;
        }

        class_index++;
      }

      result[i] = indexToClass[maxIndex];
    }

    return result;
  } else {
    throw base::data_exception("Cannot predict. Learner has to be trained first!");
    return result;
  }
}

size_t LearnerDensityBased::getNrGridPoints() {
  size_t maxGrid = 0;

  for (size_t i = 0; i < gridVec.size(); i++) {
    if (gridVec[i]->getSize() > maxGrid) maxGrid = gridVec[i]->getSize();
  }

  return maxGrid;
}

}  // namespace datadriven
}  // namespace sgpp
