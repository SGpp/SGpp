// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSGDE.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

namespace sgpp {
namespace datadriven {

LearnerSGDE::LearnerSGDE(sgpp::base::RegularGridConfiguration& gridConfig,
                         sgpp::base::AdpativityConfiguration& adaptivityConfig,
                         sgpp::solver::SLESolverConfiguration& solverConfig,
                         sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                         LearnerSGDEConfiguration& learnerSGDEConfig)
    : grid(nullptr),
      alpha(nullptr),
      samples(nullptr),
      gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      solverConfig(solverConfig),
      regularizationConfig(regularizationConfig),
      learnerSGDEConfig(learnerSGDEConfig) {}

LearnerSGDE::~LearnerSGDE() {}

void LearnerSGDE::initialize(base::DataMatrix& psamples) {
  samples = std::make_shared<base::DataMatrix>(psamples);
  size_t ndim = psamples.getNcols();
  grid = createRegularGrid(ndim);
  alpha = std::make_shared<base::DataVector>(grid->getSize());

  // optimize the regularization parameter
  double lambdaReg = 0.0;

  if (learnerSGDEConfig.doCrossValidation_) {
    lambdaReg = optimizeLambdaCV();
  } else {
    lambdaReg = learnerSGDEConfig.lambda_;
  }

  // learn the data -> do the density estimation
  train(*grid, *alpha, *samples, lambdaReg);
}

// ---------------------------------------------------------------------------

double LearnerSGDE::pdf(base::DataVector& x) {
  return op_factory::createOperationEval(*grid)->eval(*alpha, x);
}

void LearnerSGDE::pdf(base::DataMatrix& points, base::DataVector& res) {
  op_factory::createOperationMultipleEval(*grid, points)->eval(*alpha, res);
}

double LearnerSGDE::mean() {
  return op_factory::createOperationFirstMoment(*grid)->doQuadrature(*alpha);
}

double LearnerSGDE::variance() {
  double secondMoment = op_factory::createOperationSecondMoment(*grid)->doQuadrature(*alpha);

  // use Steiners translation theorem to compute the variance
  double firstMoment = mean();
  double res = secondMoment - firstMoment * firstMoment;
  return res;
}

void LearnerSGDE::cov(base::DataMatrix& cov) { return; }

std::shared_ptr<base::DataVector> LearnerSGDE::getSamples(size_t dim) {
  std::shared_ptr<base::DataVector> isamples = std::make_shared<base::DataVector>(getNsamples());
  samples->getColumn(dim, *isamples);
  return isamples;
}

std::shared_ptr<base::DataMatrix> LearnerSGDE::getSamples() { return samples; }

size_t LearnerSGDE::getDim() { return samples->getNcols(); }

size_t LearnerSGDE::getNsamples() { return samples->getNrows(); }

std::shared_ptr<base::DataVector> LearnerSGDE::getSurpluses() { return alpha; }

std::shared_ptr<base::Grid> LearnerSGDE::getGrid() { return grid; }

// ---------------------------------------------------------------------------

std::shared_ptr<base::Grid> LearnerSGDE::createRegularGrid(size_t ndim) {
  // load grid
  std::unique_ptr<base::Grid> uGrid;
  if (gridConfig.type_ == base::GridType::Linear) {
    uGrid = base::Grid::createLinearGrid(ndim);
  } else if (gridConfig.type_ == base::GridType::LinearL0Boundary) {
    uGrid = base::Grid::createLinearBoundaryGrid(ndim, 0);
  } else if (gridConfig.type_ == base::GridType::LinearBoundary) {
    uGrid = base::Grid::createLinearBoundaryGrid(ndim);
  } else {
    throw base::application_exception("LeanerSGDE::initialize : grid type is not supported");
  }

  uGrid->getGenerator().regular(gridConfig.level_);

  // move the grid to be shared
  std::shared_ptr<base::Grid> sGrid{std::move(uGrid)};

  return sGrid;
}

double LearnerSGDE::optimizeLambdaCV() {
  double curLambda;
  double bestLambda = 0;
  double curMean = 0;
  double curMeanAcc = 0;
  double bestMeanAcc = 0;

  size_t kfold = learnerSGDEConfig.kfold_;

  std::vector<std::shared_ptr<base::DataMatrix> > kfold_train(kfold);
  std::vector<std::shared_ptr<base::DataMatrix> > kfold_test(kfold);
  splitset(kfold_train, kfold_test);

  double lambdaStart = learnerSGDEConfig.lambdaStart_;
  double lambdaEnd = learnerSGDEConfig.lambdaEnd_;

  if (learnerSGDEConfig.logScale_) {
    lambdaStart = std::log(lambdaStart);
    lambdaEnd = std::log(lambdaEnd);
  }

  for (size_t i = 0; i < learnerSGDEConfig.lambdaSteps_; i++) {
    // compute current lambda
    curLambda = lambdaStart +
                static_cast<double>(i) * (lambdaEnd - lambdaStart) /
                    static_cast<double>(learnerSGDEConfig.lambdaSteps_ - 1);

    if (learnerSGDEConfig.logScale_) curLambda = exp(curLambda);

    if (i % static_cast<size_t>(
                std::max(static_cast<double>(learnerSGDEConfig.lambdaSteps_) / 10.0f,
                         static_cast<double>(1.0f))) ==
        0) {
      if (!learnerSGDEConfig.silent_) {
        std::cout << i + 1 << "/" << learnerSGDEConfig.lambdaSteps_ << " (lambda = " << curLambda
                  << ") " << std::endl;
        std::cout.flush();
      }
    }

    // cross-validation
    curMeanAcc = 0.0;
    curMean = 0.0;
    std::shared_ptr<base::Grid> grid;
    base::DataVector alpha(getNsamples());

    for (size_t j = 0; j < kfold; j++) {
      // initialize standard grid and alpha vector
      grid = createRegularGrid(getDim());
      alpha.setAll(0.0);

      // compute density
      train(*grid, alpha, *(kfold_train[j]), curLambda);
      // get L2 norm of residual for test set
      curMean = computeResidual(*grid, alpha, *(kfold_test[j]), 0.0);
      curMeanAcc += curMean;

      if (!learnerSGDEConfig.silent_) {
        std::cout << "# " << curLambda << " " << i << " " << j << " " << curMeanAcc << " "
                  << curMean << std::endl;
      }
    }

    curMeanAcc /= static_cast<double>(kfold);

    if (i == 0 || curMeanAcc < bestMeanAcc) {
      bestMeanAcc = curMeanAcc;
      bestLambda = curLambda;
    }

    if (!learnerSGDEConfig.silent_) {
      std::cout << "# " << curLambda << " " << bestLambda << " " << i << " " << curMeanAcc
                << std::endl;
    }
  }

  if (!learnerSGDEConfig.silent_) {
    std::cout << "# -> best lambda = " << bestLambda << std::endl;
  }

  return bestLambda;
}

void LearnerSGDE::train(base::Grid& grid, base::DataVector& alpha, base::DataMatrix& train,
                        double lambdaReg) {
  size_t dim = train.getNcols();

  base::GridStorage& gridStorage = grid.getStorage();
  base::GridGenerator& gridGen = grid.getGenerator();
  base::DataVector rhs(grid.getSize());
  alpha.resize(grid.getSize());
  alpha.setAll(0.0);

  if (!learnerSGDEConfig.silent_) {
    std::cout << "# LearnerSGDE: grid points " << grid.getSize() << std::endl;
  }

  for (size_t ref = 0; ref <= adaptivityConfig.numRefinements_; ref++) {
    std::unique_ptr<base::OperationMatrix> C = computeRegularizationMatrix(grid);

    datadriven::DensitySystemMatrix SMatrix(grid, train, *C, lambdaReg);
    SMatrix.generateb(rhs);

    if (!learnerSGDEConfig.silent_) {
      std::cout << "# LearnerSGDE: Solving " << std::endl;
    }

    solver::ConjugateGradients myCG(solverConfig.maxIterations_, solverConfig.eps_);
    myCG.solve(SMatrix, alpha, rhs, false, false, solverConfig.threshold_);

    if (ref < adaptivityConfig.numRefinements_) {
      if (!learnerSGDEConfig.silent_) {
        std::cout << "# LearnerSGDE: Refine grid ... ";
      }

      // Weight surplus with function evaluation at grid points
      std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(grid));
      base::GridIndex* gp;
      base::DataVector p(dim);
      base::DataVector alphaWeight(alpha.getSize());

      for (size_t i = 0; i < grid.getSize(); i++) {
        gp = gridStorage.get(i);
        gp->getCoords(p);
        alphaWeight[i] = alpha[i] * opEval->eval(alpha, p);
      }

      base::SurplusRefinementFunctor srf(alphaWeight, adaptivityConfig.noPoints_,
                                         adaptivityConfig.threshold_);
      gridGen.refine(srf);

      if (!learnerSGDEConfig.silent_) {
        std::cout << "# LearnerSGDE: ref " << ref << "/" << adaptivityConfig.numRefinements_ - 1
                  << ": " << grid.getSize() << std::endl;
      }

      alpha.resize(grid.getSize());
      rhs.resize(grid.getSize());
      alpha.setAll(0.0);
      rhs.setAll(0.0);
    }
  }

  return;
}

double LearnerSGDE::computeResidual(base::Grid& grid, base::DataVector& alpha,
                                     base::DataMatrix& test, double lambdaReg) {
  std::unique_ptr<base::OperationMatrix> C = computeRegularizationMatrix(grid);

  base::DataVector rhs(grid.getSize());
  base::DataVector res(grid.getSize());
  datadriven::DensitySystemMatrix SMatrix(grid, test, *C, lambdaReg);
  SMatrix.generateb(rhs);

  SMatrix.mult(alpha, res);

  for (size_t i = 0; i < res.getSize(); i++) {
    res[i] = res[i] - rhs[i];
  }
  return res.l2Norm();
}

std::unique_ptr<base::OperationMatrix> LearnerSGDE::computeRegularizationMatrix(base::Grid& grid) {
  std::unique_ptr<base::OperationMatrix> C;

  if (regularizationConfig.regType_ == datadriven::RegularizationType::Identity) {
    C = op_factory::createOperationIdentity(grid);
  } else if (regularizationConfig.regType_ == datadriven::RegularizationType::Laplace) {
    C = op_factory::createOperationLaplace(grid);
  } else {
    throw base::application_exception("LearnerSGDE::train : unknown regularization type");
  }

  return C;
}

void LearnerSGDE::splitset(std::vector<std::shared_ptr<base::DataMatrix> >& strain,
                           std::vector<std::shared_ptr<base::DataMatrix> >& stest) {
  std::shared_ptr<base::DataMatrix> mydata = std::make_shared<base::DataMatrix>(*samples);
  base::DataVector p(samples->getNcols());
  base::DataVector tmp(samples->getNcols());

  size_t kfold = learnerSGDEConfig.kfold_;

  std::vector<size_t> s(kfold);        // size of partition
  std::vector<size_t> ind(kfold + 1);  // index of partition
  size_t n = mydata->getNrows();       // size of data

  if (learnerSGDEConfig.shuffle_) {
    if (learnerSGDEConfig.seed_ == -1)
      srand(static_cast<unsigned int>(time(0)));
    else
      srand(learnerSGDEConfig.seed_);

    for (size_t i = 0; i < mydata->getNrows(); i++) {
      size_t r = i + (static_cast<size_t>(rand()) % (mydata->getNrows() - i));
      mydata->getRow(i, p);
      mydata->getRow(r, tmp);
      mydata->setRow(r, p);
      mydata->setRow(i, tmp);
    }
  }

  // set size of partitions
  if (!learnerSGDEConfig.silent_) std::cout << "# kfold: ";

  ind[0] = 0;

  for (size_t i = 0; i < kfold - 1; i++) {
    s[i] = n / kfold;
    ind[i + 1] = ind[i] + s[i];

    if (!learnerSGDEConfig.silent_) std::cout << s[i] << " ";
  }

  ind[kfold] = n;
  s[kfold - 1] = n - (kfold - 1) * (n / kfold);

  if (!learnerSGDEConfig.silent_) std::cout << s[kfold - 1] << std::endl;

  if (!learnerSGDEConfig.silent_) {
    std::cout << "# kfold ind: ";

    for (size_t i = 0; i <= kfold; i++) std::cout << ind[i] << " ";

    std::cout << std::endl;
  }

  // fill data
  for (size_t i = 0; i < kfold; i++) {
    // allocate memory
    strain[i] = std::make_shared<base::DataMatrix>(mydata->getNrows() - s[i], mydata->getNcols());
    stest[i] = std::make_shared<base::DataMatrix>(s[i], mydata->getNcols());

    size_t local_test = 0;
    size_t local_train = 0;

    for (size_t j = 0; j < mydata->getNrows(); j++) {
      mydata->getRow(j, p);

      if (ind[i] <= j && j < ind[i + 1]) {
        stest[i]->setRow(local_test, p);
        local_test++;
      } else {
        strain[i]->setRow(local_train, p);
        local_train++;
      }
    }
  }
}

}  // namespace datadriven
}  // namespace sgpp
