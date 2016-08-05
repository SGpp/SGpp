// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/datadriven/application/LearnerDensityBasedReg.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

LearnerDensityBasedReg::LearnerDensityBasedReg(datadriven::RegularizationType& regularization,
                                               double border)
    : LearnerBase(true), CMode(regularization), maxValue(0.), minValue(0.), border(border) {}

LearnerDensityBasedReg::~LearnerDensityBasedReg() {}

std::unique_ptr<datadriven::DMSystemMatrixBase> LearnerDensityBasedReg::createDMSystem(
    base::DataMatrix& trainDataset, double lambda) {
  // Is not used
  return nullptr;
}

LearnerTiming LearnerDensityBasedReg::train(
    base::DataMatrix& trainDataset, base::DataVector& classes,
    const base::RegularGridConfiguration& GridConfig,
    const solver::SLESolverConfiguration& SolverConfigRefine,
    const solver::SLESolverConfiguration& SolverConfigFinal,
    const base::AdpativityConfiguration& AdaptConfig, bool testAccDuringAdapt,
    const double lambda) {
  LearnerTiming timing = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  size_t dim = trainDataset.getNcols();
  size_t m = trainDataset.getNrows();

  if (m != classes.getSize()) {
    throw base::application_exception(
        "LearnerBase::train: length of classes vector does not match to dataset!");
  }

  execTime = 0.0;
  GFlop = 0.0;
  GByte = 0.0;

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
        "LearnerDensityBasedReg::train: An unsupported SLE solver type was chosen!");
  }

  std::unique_ptr<base::SGppStopwatch> myStopwatch = std::make_unique<base::SGppStopwatch>();

  if (isVerbose) std::cout << "Starting Learning...." << std::endl;

  // Normalize "class" vector and append it to the train data matrix to create a m x (n+1) matrix
  // for estimation of Pr[x, t]:
  base::DataVector classesCopy(classes);
  classesCopy.minmax(&minValue, &maxValue);
  classesCopy.normalize(border);
  base::DataMatrix densityMatrix(m, dim + 1);
  double* densityData = densityMatrix.getPointer();
  double* trainData = trainDataset.getPointer();
  double* classesData = classesCopy.getPointer();
  size_t acc = 0;

  for (size_t i = 0; i < m; i++) {
    size_t j = 0;

    for (; j < dim; j++) {
      densityData[i * (dim + 1) + j] = trainData[acc + j];
    }

    densityData[i * (dim + 1) + j] = classesData[i];
    acc += dim;
  }

  // Adjust grid configuration to the higher dimension:
  base::RegularGridConfiguration adaptivityConfig;
  adaptivityConfig.dim_ = GridConfig.dim_ + 1;
  adaptivityConfig.level_ = GridConfig.level_;
  adaptivityConfig.type_ = GridConfig.type_;

  if (isTrained == true) isTrained = false;

  InitializeGrid(adaptivityConfig);

  // check if grid was created
  if (grid == nullptr) return timing;

  for (size_t i = 0; i < AdaptConfig.numRefinements_ + 1; i++) {
    if (isVerbose) std::cout << std::endl << "Doing refinement: " << i << std::endl;

    myStopwatch->start();

    // Do Refinements
    if (i > 0) {
      base::SurplusRefinementFunctor myRefineFunc(*alpha, AdaptConfig.noPoints_);
      grid->getGenerator().refine(myRefineFunc);

      alpha->resize(grid->getSize());

      // DMSystem->rebuildLevelAndIndex();   not implemented

      if (isVerbose) std::cout << "New Grid Size: " << grid->getSize() << std::endl;
    } else {
      if (isVerbose) std::cout << "Grid Size: " << grid->getSize() << std::endl;
    }

    // regularization term:

    if (this->CMode == datadriven::RegularizationType::Laplace) {
      C = std::unique_ptr<base::OperationMatrix>(op_factory::createOperationLaplace(*grid));
    } else if (this->CMode == datadriven::RegularizationType::Identity) {
      C = std::unique_ptr<base::OperationMatrix>(op_factory::createOperationIdentity(*grid));
    } else {
      throw base::application_exception(
          "LearnerDensityBased::train: Unknown regularization "
          "operator");
    }

    datadriven::DensitySystemMatrix DMatrix(*grid, densityMatrix, *C, lambda);
    base::DataVector rhs(grid->getSize());
    DMatrix.generateb(rhs);

    if (i == AdaptConfig.numRefinements_) {
      myCG->setMaxIterations(SolverConfigFinal.maxIterations_);
      myCG->setEpsilon(SolverConfigFinal.eps_);
    }

    myCG->solve(DMatrix, *alpha, rhs, false, false, -1.0);

    if (isVerbose) {
      std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
      std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;
    }

    if (testAccDuringAdapt) {
      double acc = getAccuracy(trainDataset, classes);

      if (isVerbose) std::cout << "MSE (train): " << acc << std::endl;

      if ((i > 0) && (oldAcc <= acc)) {
        if (isVerbose) std::cout << "The grid is becoming worse --> stop learning" << std::endl;

        break;
      }

      oldAcc = acc;
    }

    execTime += myStopwatch->stop();

    // use post-processing to determine Flops and time
    if (i < AdaptConfig.numRefinements_) {
      postProcessing(trainDataset, SolverConfigRefine.type_, myCG->getNumberIterations());
    } else {
      postProcessing(trainDataset, SolverConfigFinal.type_, myCG->getNumberIterations());
    }
  }

  timing.timeComplete_ = execTime;

  isTrained = true;

  if (isVerbose) {
    std::cout << "Finished Training!" << std::endl << std::endl;
    std::cout << "Training took: " << execTime << " seconds" << std::endl << std::endl;
  }

  return timing;
}

base::DataVector LearnerDensityBasedReg::predict(base::DataMatrix& testDataset) {
  base::DataVector res(testDataset.getNrows());

  double delta = (maxValue - minValue) / (1 - 2 * border);

  size_t dim = testDataset.getNcols();
  size_t m = testDataset.getNrows();
  base::DataVector point(dim);

  // TODO(franzfn): might be buggy thanks to smart pointer transition (by David)
  for (size_t i = 0; i < m; i++) {
    testDataset.getRow(i, point);

    base::Grid* tempGrid = grid.get();
    base::Grid* lastGrid = NULL;
    base::DataVector* lastAlpha = alpha.get();

    // Conditionalize for all dimensions, but the last one:
    for (size_t j = 0; j < dim; j++) {
      base::DataVector* tempAlpha = new base::DataVector(1);
      op_factory::createOperationDensityConditional(*grid)->doConditional(
          *lastAlpha, tempGrid, *tempAlpha, 0, point.get(j));

      if (j > 0) {
        delete lastAlpha;
        delete lastGrid;
      }

      // TODO(franzfn): please check, all these assignments seem to without effect (by David)
      lastGrid = tempGrid;
      lastAlpha = tempAlpha;
    }

    // Compute conditional expectation:
    double valueNormalized =
        op_factory::createOperationFirstMoment(*lastGrid)->doQuadrature(*lastAlpha);

    res.set(i, ((valueNormalized - border) * delta) + minValue);

    delete lastAlpha;
    delete lastGrid;
  }

  return res;
}

void LearnerDensityBasedReg::dumpDensityAtPoint(base::DataVector& point, std::string fileName,
                                                unsigned int resolution) {
  size_t dim = point.getSize();
  base::Grid* tempGrid = grid.get();
  base::Grid* lastGrid = NULL;
  base::DataVector* lastAlpha = alpha.get();

  // Conditionalize for all dimensions, but the last one:
  for (size_t j = 0; j < dim; j++) {
    base::DataVector* tempAlpha = new base::DataVector(1);
    op_factory::createOperationDensityConditional(*tempGrid)->doConditional(
        *lastAlpha, tempGrid, *tempAlpha, 0, point.get(j));

    if (j > 0) {
      delete lastAlpha;
      delete lastGrid;
    }

    lastGrid = tempGrid;
    lastAlpha = tempAlpha;
  }

  base::GridPrinter myPlotter(*lastGrid);
  myPlotter.printGrid(*lastAlpha, fileName, resolution);
}

}  // namespace datadriven
}  // namespace sgpp
