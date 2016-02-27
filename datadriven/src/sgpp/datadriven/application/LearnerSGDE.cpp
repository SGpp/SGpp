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

#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

// using namespace std;
// using namespace sgpp::base;

using sgpp::base::DataVector;
using sgpp::base::DataMatrix;
using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::OperationMatrix;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::GridIndex;
using sgpp::base::OperationEval;
using sgpp::base::GridGenerator;
using sgpp::base::GridType;
using sgpp::base::OperationSecondMoment;
using sgpp::base::OperationFirstMoment;
using sgpp::base::OperationMultipleEval;

using std::vector;
using std::cout;
using std::endl;

namespace sgpp {
namespace datadriven {

LearnerSGDE::LearnerSGDE(sgpp::base::RegularGridConfiguration& gridConfig,
                         sgpp::base::AdpativityConfiguration& adaptivityConfig,
                         sgpp::solver::SLESolverConfiguration& solverConfig,
                         sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                         LearnerSGDEConfiguration& learnerSGDEConfig) :
  grid(NULL), alpha(1), samples(NULL), gridConfig(gridConfig), adaptivityConfig(
    adaptivityConfig), solverConfig(solverConfig), regularizationConfig(
      regularizationConfig), learnerSGDEConfig(learnerSGDEConfig) {
}

LearnerSGDE::~LearnerSGDE() {
  if (samples != NULL) {
    delete samples;
  }

  if (grid != NULL) {
    delete grid;
  }
}

void LearnerSGDE::initialize(sgpp::base::DataMatrix& psamples) {
  samples = new DataMatrix(psamples);
  size_t ndim = psamples.getNcols();
  createRegularGrid(grid, ndim);
  alpha.resize(grid->getSize());

  // optimize the regularization parameter
  double lambdaReg = 0.0;

  if (learnerSGDEConfig.doCrossValidation_) {
    lambdaReg = optimizeLambdaCV();
  } else {
    lambdaReg = learnerSGDEConfig.lambda_;
  }

  // learn the data -> do the density estimation
  train(*grid, alpha, *samples, lambdaReg);
}

// ---------------------------------------------------------------------------

double LearnerSGDE::pdf(DataVector& x) {
  return sgpp::op_factory::createOperationEval(*grid)->eval(alpha, x);
}

void LearnerSGDE::pdf(DataMatrix& points, DataVector& res) {
  sgpp::op_factory::createOperationMultipleEval(*grid, points)->eval(alpha, res);
}

double LearnerSGDE::mean() {
  return op_factory::createOperationFirstMoment(*grid)->doQuadrature(alpha);
}

double LearnerSGDE::variance() {
  double secondMoment = op_factory::createOperationSecondMoment(*grid)->doQuadrature(alpha);

  // use Steiners translation theorem to compute the variance
  double firstMoment = mean();
  double res = secondMoment - firstMoment * firstMoment;
  return res;
}

void LearnerSGDE::cov(DataMatrix& cov) {
  return;
}

DataVector* LearnerSGDE::getSamples(size_t dim) {
  DataVector* isamples = new DataVector(getNsamples());
  samples->getColumn(dim, *isamples);
  return isamples;
}

DataMatrix* LearnerSGDE::getSamples() {
  return new DataMatrix(*samples);
}

size_t LearnerSGDE::getDim() {
  return samples->getNcols();
}

size_t LearnerSGDE::getNsamples() {
  return samples->getNrows();
}

DataVector LearnerSGDE::getSurpluses() {
  return alpha;
}

GridStorage* LearnerSGDE::getGridStorage() {
  return &grid->getStorage();
}
// ---------------------------------------------------------------------------

void LearnerSGDE::createRegularGrid(Grid*& grid, size_t ndim) {
  // load grid
  if (gridConfig.type_ == GridType::Linear) {
    grid = Grid::createLinearGrid(ndim).release();
  } else if (gridConfig.type_ == GridType::LinearL0Boundary) {
    grid = Grid::createLinearBoundaryGrid(ndim, 0).release();
  } else if (gridConfig.type_ == GridType::LinearBoundary) {
    grid = Grid::createLinearBoundaryGrid(ndim).release();
  } else {
    throw base::application_exception("LeanerSGDE::initialize : grid type is not supported");
  }

  grid->getGenerator().regular(gridConfig.level_);
}

double LearnerSGDE::optimizeLambdaCV() {
  Grid* grid = NULL;
  DataVector* alpha = NULL;

  double curLambda;
  double bestLambda = 0;
  double curMean = 0;
  double curMeanAcc = 0;
  double bestMeanAcc = 0;

  size_t kfold = learnerSGDEConfig.kfold_;

  vector<DataMatrix*> kfold_train(kfold);
  vector<DataMatrix*> kfold_test(kfold);
  splitset(kfold_train, kfold_test);

  double lambdaStart = learnerSGDEConfig.lambdaStart_;
  double lambdaEnd = learnerSGDEConfig.lambdaEnd_;

  if (learnerSGDEConfig.logScale_) {
    lambdaStart = log(lambdaStart);
    lambdaEnd = log(lambdaEnd);
  }

  for (size_t i = 0; i < learnerSGDEConfig.lambdaSteps_; i++) {
    // compute current lambda
    curLambda = lambdaStart
                + static_cast<double>(i) * (lambdaEnd - lambdaStart)
                / static_cast<double>(learnerSGDEConfig.lambdaSteps_
                                       - 1);

    if (learnerSGDEConfig.logScale_)
      curLambda = exp(curLambda);

    if (i
        % static_cast<size_t>(std::max(
                                static_cast<double>(learnerSGDEConfig.lambdaSteps_)
                                / 10.0f, static_cast<double>(1.0f))) == 0) {
      if (!learnerSGDEConfig.silent_) {
        cout << i + 1 << "/" << learnerSGDEConfig.lambdaSteps_
             << " (lambda = " << curLambda << ") " << endl;
        cout.flush();
      }
    }

    // cross-validation
    curMeanAcc = 0.0;
    curMean = 0.0;

    for (size_t j = 0; j < kfold; j++) {
      // initialize standard grid and alpha vector
      createRegularGrid(grid, getDim());
      alpha = new DataVector(getNsamples());
      // compute density
      train(*grid, *alpha, *(kfold_train[j]), curLambda);
      // get L2 norm of residual for test set
      curMean = computeResidual(*grid, *alpha, *(kfold_test[j]), 0.0);
      curMeanAcc += curMean;

      if (!learnerSGDEConfig.silent_) {
        cout << "# " << curLambda << " " << i << " " << j << " "
             << curMeanAcc << " " << curMean << endl;
      }

      // free space
      delete grid;
      delete alpha;
    }

    curMeanAcc /= static_cast<double>(kfold);

    if (i == 0 || curMeanAcc < bestMeanAcc) {
      bestMeanAcc = curMeanAcc;
      bestLambda = curLambda;
    }

    if (!learnerSGDEConfig.silent_) {
      cout << "# " << curLambda << " " << bestLambda << " " << i << " "
           << curMeanAcc << endl;
    }
  }

  if (!learnerSGDEConfig.silent_) {
    cout << "# -> best lambda = " << bestLambda << endl;
  }

  // free splitted sets
  for (size_t i = 0; i < kfold; i++) {
    delete kfold_train[i];
    delete kfold_test[i];
  }

  return bestLambda;
}

void LearnerSGDE::train(Grid& grid, DataVector& alpha, DataMatrix& train,
                        double lambdaReg) {
  size_t dim = train.getNcols();

  GridStorage* gridStorage = &grid.getStorage();
  GridGenerator& gridGen = grid.getGenerator();
  DataVector rhs(grid.getSize());
  alpha.resize(grid.getSize());
  alpha.setAll(0.0);

  if (!learnerSGDEConfig.silent_) {
    cout << "# LearnerSGDE: grid points " << grid.getSize() << endl;
  }

  for (size_t ref = 0; ref <= adaptivityConfig.numRefinements_; ref++) {
    OperationMatrix* C = computeRegularizationMatrix(grid);

    sgpp::datadriven::DensitySystemMatrix SMatrix(grid, train, *C, lambdaReg);
    SMatrix.generateb(rhs);

    if (!learnerSGDEConfig.silent_) {
      cout << "# LearnerSGDE: Solving " << endl;
    }

    sgpp::solver::ConjugateGradients myCG(solverConfig.maxIterations_,
                                          solverConfig.eps_);
    myCG.solve(SMatrix, alpha, rhs, false, false, solverConfig.threshold_);

    if (ref < adaptivityConfig.numRefinements_) {
      if (!learnerSGDEConfig.silent_) {
        cout << "# LearnerSGDE: Refine grid ... ";
      }

      // Weight surplus with function evaluation at grid points
      std::unique_ptr<OperationEval> opEval(sgpp::op_factory::createOperationEval(grid));
      GridIndex* gp;
      DataVector p(dim);
      DataVector alphaWeight(alpha.getSize());

      for (size_t i = 0; i < gridStorage->getSize(); i++) {
        gp = gridStorage->get(i);
        gp->getCoords(p);
        alphaWeight[i] = alpha[i] * opEval->eval(alpha, p);
      }

      SurplusRefinementFunctor srf(alphaWeight,
                                   adaptivityConfig.noPoints_, adaptivityConfig.threshold_);
      gridGen.refine(srf);

      if (!learnerSGDEConfig.silent_) {
        cout << "# LearnerSGDE: ref " << ref << "/"
             << adaptivityConfig.numRefinements_ - 1 << ": "
             << grid.getSize() << endl;
      }

      alpha.resize(grid.getSize());
      rhs.resize(grid.getSize());
      alpha.setAll(0.0);
      rhs.setAll(0.0);
    }

    delete C;
  }

  return;
}

double LearnerSGDE::computeResidual(Grid& grid, DataVector& alpha,
                                     DataMatrix& test, double lambdaReg) {
  OperationMatrix* C = computeRegularizationMatrix(grid);

  DataVector rhs(grid.getSize());
  DataVector res(grid.getSize());
  sgpp::datadriven::DensitySystemMatrix SMatrix(grid, test, *C, lambdaReg);
  SMatrix.generateb(rhs);

  SMatrix.mult(alpha, res);

  for (size_t i = 0; i < res.getSize(); i++)
    res[i] = res[i] - rhs[i];

  delete C;
  return res.l2Norm();
}

OperationMatrix* LearnerSGDE::computeRegularizationMatrix(
  sgpp::base::Grid& grid) {
  OperationMatrix* C = NULL;

  if (regularizationConfig.regType_
      == sgpp::datadriven::RegularizationType::Identity) {
    C = sgpp::op_factory::createOperationIdentity(grid).release();
  } else if (regularizationConfig.regType_
             == sgpp::datadriven::RegularizationType::Laplace) {
    C = sgpp::op_factory::createOperationLaplace(grid).release();
  } else {
    throw base::application_exception("LearnerSGDE::train : unknown regularization type");
  }

  return C;
}

void LearnerSGDE::splitset(vector<DataMatrix*>& strain,
                           vector<DataMatrix*>& stest) {
  DataMatrix* mydata = new DataMatrix(*samples);
  DataVector p(samples->getNcols());
  DataVector tmp(samples->getNcols());

  size_t kfold = learnerSGDEConfig.kfold_;

  vector<size_t> s(kfold);  // size of partition
  vector<size_t> ind(kfold + 1);  // index of partition
  size_t n = mydata->getNrows();  // size of data

  if (learnerSGDEConfig.shuffle_) {
    if (learnerSGDEConfig.seed_ == -1)
      srand(static_cast<unsigned int>(time(0)));
    else
      srand(learnerSGDEConfig.seed_);

    for (size_t i = 0; i < mydata->getNrows(); i++) {
      size_t r = i
                 + (static_cast<size_t>(rand()) % (mydata->getNrows() - i));
      mydata->getRow(i, p);
      mydata->getRow(r, tmp);
      mydata->setRow(r, p);
      mydata->setRow(i, tmp);
    }
  }

  // set size of partitions
  if (!learnerSGDEConfig.silent_)
    cout << "# kfold: ";

  ind[0] = 0;

  for (size_t i = 0; i < kfold - 1; i++) {
    s[i] = n / kfold;
    ind[i + 1] = ind[i] + s[i];

    if (!learnerSGDEConfig.silent_)
      cout << s[i] << " ";
  }

  ind[kfold] = n;
  s[kfold - 1] = n - (kfold - 1) * (n / kfold);

  if (!learnerSGDEConfig.silent_)
    cout << s[kfold - 1] << endl;

  if (!learnerSGDEConfig.silent_) {
    cout << "# kfold ind: ";

    for (size_t i = 0; i <= kfold; i++)
      cout << ind[i] << " ";

    cout << endl;
  }

  // fill data
  for (size_t i = 0; i < kfold; i++) {
    // allocate memory
    strain[i] = new DataMatrix(mydata->getNrows() - s[i],
                               mydata->getNcols());
    stest[i] = new DataMatrix(s[i], mydata->getNcols());

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

  delete mydata;
}

base::Grid* LearnerSGDE::getGrid() {
  return grid;
}

base::DataVector* LearnerSGDE::getAlpha() {
  return &alpha;
}

}  // namespace datadriven
}  // namespace sgpp

