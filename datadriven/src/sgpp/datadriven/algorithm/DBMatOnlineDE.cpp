// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinear.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalInterModLinear.hpp>
#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>

#ifdef USE_GSL
#include <gsl/gsl_blas.h>
#endif /* USE_GSL */

#include <list>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::algorithm_exception;

DBMatOnlineDE::DBMatOnlineDE(DBMatOffline& offline, Grid& grid, double lambda, double beta)
    : DBMatOnline{offline},
      beta(beta),
      totalPoints(0),
      canCV(false),
      lambdaStep(0),
      lambdaStart(0),
      lambdaEnd(0),
      testMat(nullptr),
      testMatRes(nullptr),
      cvLogscale(false),
      normFactor(1.),
      lambda(lambda) {
  functionComputed = false;
  bSave = DataVector(offlineObject.getDecomposedMatrix().getNcols());
  bTotalPoints = DataVector(offlineObject.getDecomposedMatrix().getNcols(), 0.0);
  oDim = grid.getDimension();
}

void DBMatOnlineDE::computeDensityFunction(DataVector& alpha, DataMatrix& m, Grid& grid,
    DensityEstimationConfiguration& densityEstimationConfig, bool save_b, bool do_cv,
    std::list<size_t>* deletedPoints,
    size_t newPoints) {
  std::cout << "Computing density function..." << std::endl;
  if (m.getNrows() > 0) {
    DataMatrix& lhsMatrix = offlineObject.getDecomposedMatrix();

    // in case OrthoAdapt, the current size is not lhs size, but B size
    bool use_B_size = false;
    sgpp::datadriven::DBMatOnlineDEOrthoAdapt* thisOrthoAdaptPtr;
    if (densityEstimationConfig.decomposition_ ==
        sgpp::datadriven::MatrixDecompositionType::OrthoAdapt) {
      thisOrthoAdaptPtr = static_cast<sgpp::datadriven::DBMatOnlineDEOrthoAdapt*>(&*this);
      if (thisOrthoAdaptPtr->getB().getNcols() > 1) {
        use_B_size = true;
      }
    }

    // Compute right hand side of the equation:
    size_t numberOfPoints = m.getNrows();
    totalPoints++;
    DataVector b(use_B_size ? thisOrthoAdaptPtr->getB().getNcols() : lhsMatrix.getNcols());
    b.setAll(0);
    if (b.getSize() != offlineObject.getGridSize()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatOnlineDE::computeDensityFunction: b doesn't match size of system matrix");
    }

    std::unique_ptr<sgpp::base::OperationMultipleEval> B(
        (offlineObject.interactions.size() == 0)
            ? sgpp::op_factory::createOperationMultipleEval(grid, m)
            : sgpp::op_factory::createOperationMultipleEvalInter(grid, m,
                                                                 offlineObject.interactions));

    DataVector y(numberOfPoints);
    y.setAll(1.0);
    // Bt * 1
    B->multTranspose(y, b);

    // Perform permutation because of decomposition (LU)
    if (densityEstimationConfig.decomposition_ == MatrixDecompositionType::LU) {
#ifdef USE_GSL
      static_cast<DBMatOfflineLU&>(offlineObject).permuteVector(b);
#else
      throw algorithm_exception("built withot GSL");
#endif /*USE_GSL*/
    }

    // std::cout << b.getSize() << std::endl;

    if (save_b) {
      if (functionComputed) {
        // double tmpBeta = std::max(beta, (1./(double)totalPoints));
        // b.mult(tmpBeta);
        // bSave->mult(1.-tmpBeta);

        // Delete indices when grid got coarsend-> reduce 'bSave'
        if (deletedPoints != nullptr && !deletedPoints->empty()) {
          std::vector<size_t> v{std::begin(*deletedPoints), std::end(*deletedPoints)};
          std::vector<size_t> v1(bSave.getSize() - deletedPoints->size());
          size_t old_size = bSave.getSize();

          size_t index_coarse = 0;
          size_t index_remain = 0;
          size_t temp;
          for (size_t j = 0; j < old_size; j++) {
            temp = v[index_coarse];
            if (temp == j) {
              index_coarse++;
              continue;
            } else {
              v1[index_remain] = j;
              index_remain++;
            }
          }
          bSave.restructure(v1);
          bTotalPoints.restructure(v1);
        }

        // Expand 'b_save' when grid got refined
        if (newPoints > 0) {
          bSave.resizeZero(b.getSize());
          bTotalPoints.resizeZero(b.getSize());
        }

        b.add(bSave);
        // b.mult(beta);
      }

      // Update weighting based on processed data points
      for (size_t i = 0; i < b.getSize(); i++) {
        bSave.set(i, b.get(i));
        bTotalPoints.set(i, static_cast<double>(numberOfPoints) + bTotalPoints.get(i));
        b.set(i, bSave.get(i) * (1. / bTotalPoints.get(i)));
      }
    } else {
      // 1 / M * Bt * 1
      b.mult(1. / static_cast<double>(numberOfPoints));
    }

    solveSLE(alpha, b, grid, densityEstimationConfig, do_cv);

    functionComputed = true;
  }
}

double DBMatOnlineDE::resDensity(DataVector& alpha, Grid& grid) {
  auto C = sgpp::op_factory::createOperationIdentity(grid);
  DataVector rhs(grid.getSize());
  DataVector res(grid.getSize());
  sgpp::datadriven::DensitySystemMatrix SMatrix(grid, *testMat, C, 0.0);

  SMatrix.generateb(rhs);

  SMatrix.mult(alpha, res);

  for (size_t i = 0; i < res.getSize(); i++) res[i] -= rhs[i];
  return res.l2Norm();
}

double DBMatOnlineDE::computeL2Error(DataVector& alpha, Grid& grid) {
  size_t nRows = testMatRes->getNrows();
  DataVector r(nRows);
  DataVector tmp(testMat->getNcols());
  for (size_t i = 0; i < nRows; i++) {
    testMat->getRow(i, tmp);
    r[i] = this->eval(alpha, tmp, grid, true);
  }
  double l2err = 0;
  for (size_t i = 0; i < nRows; i++) {
    l2err += (testMatRes->get(i, 0) - r[i]) * (testMatRes->get(i, 0) - r[i]);
  }
  return sqrt(l2err) / static_cast<double>(nRows);
}

double DBMatOnlineDE::eval(DataVector& alpha, const DataVector& p, Grid& grid, bool force) {
  if (functionComputed || force == true) {
    double res;
    std::unique_ptr<sgpp::base::OperationEval> opEval(
        sgpp::op_factory::createOperationEval(grid));
    res = opEval->eval(alpha, p);
    return res * normFactor;
  } else {
    throw algorithm_exception("Density function not computed, yet!");
  }
}

void DBMatOnlineDE::eval(DataVector& alpha, DataMatrix& values, DataVector& results, Grid& grid,
    bool force) {
  if (functionComputed || force == true) {
    std::unique_ptr<sgpp::base::OperationMultipleEval> opEval(
        (offlineObject.interactions.size() == 0)
            ? sgpp::op_factory::createOperationMultipleEval(grid, values)
            : sgpp::op_factory::createOperationMultipleEvalInter(grid, values,
                                                                 offlineObject.interactions));
    opEval->eval(alpha, results);
    results.mult(normFactor);
  } else {
    throw algorithm_exception("Density function not computed, yet!");
  }
}


bool DBMatOnlineDE::isComputed() { return functionComputed; }

void DBMatOnlineDE::setCrossValidationParameters(int lambda_step, double lambda_start,
                                                 double lambda_end, DataMatrix* test,
                                                 DataMatrix* test_cc, bool logscale) {
  lambdaStep = lambda_step;
  cvLogscale = logscale;
  if (cvLogscale) {
    lambdaStart = std::log(lambda_start);
    lambdaEnd = std::log(lambda_end);
  } else {
    lambdaStart = lambda_start;
    lambdaEnd = lambda_end;
  }
  if (test != nullptr) testMat = test;
  if (test_cc != nullptr) testMatRes = test_cc;
  canCV = true;
}

double DBMatOnlineDE::getBestLambda() { return lambda; }

void DBMatOnlineDE::setBeta(double newBeta) { beta = newBeta; }

double DBMatOnlineDE::getBeta() { return beta; }

double DBMatOnlineDE::normalize(DataVector& alpha, Grid& grid, size_t samples) {
  this->normFactor = 1.;
  double sum = 0.;
  DataVector p(this->oDim);
  srand(static_cast<unsigned int>(time(nullptr)));
  for (size_t i = 0; i < samples; i++) {
    for (size_t j = 0; j < this->oDim; j++) p[j] = (static_cast<double>(rand()) / RAND_MAX);
    sum += this->eval(alpha, p, grid);
  }
  return this->normFactor = static_cast<double>(samples) / sum;
}

double DBMatOnlineDE::normalizeQuadrature(DataVector& alpha, Grid& grid) {
  this->normFactor = 1.;
  double quadrature =
      sgpp::op_factory::createOperationQuadrature(grid)->doQuadrature(alpha);

  return this->normFactor /= quadrature;
}

}  // namespace datadriven
}  // namespace sgpp
