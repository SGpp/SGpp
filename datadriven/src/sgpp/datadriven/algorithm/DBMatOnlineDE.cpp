// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalInterModLinear.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinear.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineLU.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE_SMW.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>

#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalDistributed.hpp>

#ifdef USE_GSL
#include <gsl/gsl_blas.h>
#endif /* USE_GSL */

#include <algorithm>
#include <iostream>
#include <list>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::algorithm_exception;

DBMatOnlineDE::DBMatOnlineDE(DBMatOffline& offline, Grid& grid, double lambda,
                             double beta)
    : DBMatOnline{offline},
      localVectorsInitialized(false),
      distributedVectorsInitialized(false),
      beta(beta),
      totalPoints(0),
      testMat(nullptr),
      testMatRes(nullptr),
      normFactor(1.),
      lambda(lambda) {
  functionComputed = false;
  oDim = grid.getDimension();
}

void DBMatOnlineDE::updateRhs(size_t gridSize,
                              std::vector<size_t>& deletedPoints) {
  if (functionComputed) {
    // Coarsening -> remove all idx in deletedPoints
    if (deletedPoints.size() > 0) {
      if (localVectorsInitialized) {
        bSave.remove(deletedPoints);
        bTotalPoints.remove(deletedPoints);
      }

      if (distributedVectorsInitialized) {
        DataVector tmpbSave = bSaveDistributed->toLocalDataVector();
        DataVector tmpbTotalPoints =
            bTotalPointsDistributed->toLocalDataVector();

        auto processGrid = bSaveDistributed->getProcessGrid();
        if (processGrid->getCurrentRow() == 0 &&
            processGrid->getCurrentColumn() == 0) {
          tmpbSave.remove(deletedPoints);
          tmpbTotalPoints.remove(deletedPoints);
        }
        bSaveDistributed->distribute(tmpbSave.data());
        bTotalPointsDistributed->distribute(tmpbTotalPoints.data());
      }
    }
    // Refinement -> append newPoints zeros to b
    if (localVectorsInitialized && gridSize > bSave.size()) {
      bSave.resizeZero(gridSize);
      bTotalPoints.resizeZero(gridSize);
    }

    if (distributedVectorsInitialized &&
        gridSize > bSaveDistributed->getGlobalRows()) {
      bSaveDistributed->resize(gridSize);
      bTotalPointsDistributed->resize(gridSize);
    }
  }
}

void DBMatOnlineDE::computeDensityFunction(
    DataVector& alpha, Grid& grid,
    DensityEstimationConfiguration& densityEstimationConfig, bool do_cv) {
  if (functionComputed) {
    if (alpha.size() == bSave.size() && alpha.size() == bTotalPoints.size()) {
      // Assemble the rhs
      DataVector b(alpha.size());
      for (size_t i = 0; i < b.size(); i++) {
        b.set(i, bSave.get(i) * (1. / bTotalPoints.get(i)));
      }
      // Resolve the SLE
      solveSLE(alpha, b, grid, densityEstimationConfig, do_cv);
    } else {
      throw sgpp::base::algorithm_exception(
          "Recomputation of density function with mismatching alpha size and b "
          "size");
    }
  } else {
    throw sgpp::base::algorithm_exception(
        "Density function can not be recomputed without any b stored in "
        "DBMatOnlineDE");
  }
}

void DBMatOnlineDE::computeDensityFunction(
    DataVector& alpha, DataMatrix& m, Grid& grid,
    DensityEstimationConfiguration& densityEstimationConfig, bool save_b,
    bool do_cv) {
  if (!localVectorsInitialized) {
    // init bsave and bTotalPoints only here, as they are not needed in the
    // parallel version
    bSave = DataVector(offlineObject.getDecomposedMatrix().getNcols(), 0.0);
    bTotalPoints =
        DataVector(offlineObject.getDecomposedMatrix().getNcols(), 0.0);

    localVectorsInitialized = true;
  }

  if (m.getNrows() > 0) {
    DataVector b = computeBFromBatch(m, grid, densityEstimationConfig);
    size_t numberOfPoints = m.getNrows();
    totalPoints++;

    if (save_b) {
      // Old rhs is weighted by beta
      bSave.mult(beta);
      b.add(bSave);

      // Update weighting based on processed data points
      for (size_t i = 0; i < b.getSize(); i++) {
        bSave.set(i, b.get(i));
        bTotalPoints.set(
            i, static_cast<double>(numberOfPoints) + bTotalPoints.get(i));
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

void DBMatOnlineDE::computeDensityFunctionParallel(
    DataVectorDistributed& alpha, Grid& grid,
    DensityEstimationConfiguration& densityEstimationConfig,
    const ParallelConfiguration& parallelConfig,
    std::shared_ptr<BlacsProcessGrid> processGrid, bool do_cv) {
  if (functionComputed) {
    if (alpha.getGlobalRows() == bSaveDistributed->getGlobalRows() &&
        alpha.getGlobalRows() == bTotalPointsDistributed->getGlobalRows()) {
      // Assemble the rhs (on each process)
      DataVectorDistributed b(processGrid, alpha.getGlobalRows(),
                              parallelConfig.rowBlockSize_);
      for (size_t i = 0; i < b.getLocalRows(); i++) {
        b.getLocalPointer()[i] =
            bSaveDistributed->getLocalPointer()[i] *
            (1. / bTotalPointsDistributed->getLocalPointer()[i]);
      }
      // Resolve the SLE
      solveSLEParallel(alpha, b, grid, densityEstimationConfig, do_cv);
    } else {
      throw sgpp::base::algorithm_exception(
          "Recomputation of density function with mismatching alpha size and b "
          "size");
    }
  } else {
    throw sgpp::base::algorithm_exception(
        "Density function can not be recomputed without any b stored in "
        "DBMatOnlineDE");
  }
}

void DBMatOnlineDE::computeDensityFunctionParallel(
    DataVectorDistributed& alpha, DataMatrix& m, Grid& grid,
    DensityEstimationConfiguration& densityEstimationConfig,
    const ParallelConfiguration& parallelConfig,
    std::shared_ptr<BlacsProcessGrid> processGrid, bool save_b, bool do_cv,
    std::list<size_t>* deletedPoints, size_t newPoints) {
  if (save_b && !distributedVectorsInitialized) {
    // init bSaveDistributed and bTotalPointsDistributed only here, as they are
    // not needed in the
    // local version
    bSaveDistributed = std::make_unique<DataVectorDistributed>(
        processGrid, offlineObject.getDecomposedMatrix().getNcols(),
        parallelConfig.rowBlockSize_);
    bTotalPointsDistributed = std::make_unique<DataVectorDistributed>(
        processGrid, offlineObject.getDecomposedMatrix().getNcols(),
        parallelConfig.rowBlockSize_);

    distributedVectorsInitialized = true;
  }

  DataVectorDistributed b = computeBFromBatchParallel(
      m, grid, densityEstimationConfig, parallelConfig, processGrid);
  size_t numberOfPoints = m.getNrows();
  totalPoints++;

  if (save_b) {
    // Old rhs is weighted by beta
    bSaveDistributed->scale(beta);
    b.add(*bSaveDistributed);

    // Update weighting based on processed data points
    for (size_t i = 0; i < b.getGlobalRows(); i++) {
      bSaveDistributed->set(i, b.get(i));
      bTotalPointsDistributed->set(i, static_cast<double>(numberOfPoints) +
                                          bTotalPointsDistributed->get(i));
      b.set(i,
            bSaveDistributed->get(i) * (1. / bTotalPointsDistributed->get(i)));
    }
  } else {
    // 1/M * (Bt * 1)
    b.scale(1. / static_cast<double>(numberOfPoints));
  }

  solveSLEParallel(alpha, b, grid, densityEstimationConfig, do_cv);

  functionComputed = true;
}

DataVector DBMatOnlineDE::computeBFromBatch(
    DataMatrix& m, Grid& grid,
    DensityEstimationConfiguration& densityEstimationConfig) {
  if (m.getNrows() > 0) {
    DataMatrix& lhsMatrix = offlineObject.getDecomposedMatrix();

    // in case OrthoAdapt or both SMW_, the current size is not lhs size, but B
    // size
    bool use_B_size = false;
    size_t B_size = 0;
    sgpp::datadriven::DBMatOnlineDEOrthoAdapt* this_OrthoAdapt_pointer;
    if (densityEstimationConfig.decomposition_ ==
        sgpp::datadriven::MatrixDecompositionType::OrthoAdapt) {
      this_OrthoAdapt_pointer =
          static_cast<sgpp::datadriven::DBMatOnlineDEOrthoAdapt*>(&*this);
      if (this_OrthoAdapt_pointer->getB().getNcols() > 1) {
        use_B_size = true;
        B_size = this_OrthoAdapt_pointer->getB().getNcols();
      }
    }

    sgpp::datadriven::DBMatOnlineDE_SMW* this_SMW_pointer;
    if (densityEstimationConfig.decomposition_ ==
            sgpp::datadriven::MatrixDecompositionType::SMW_ortho ||
        densityEstimationConfig.decomposition_ ==
            sgpp::datadriven::MatrixDecompositionType::SMW_chol) {
      this_SMW_pointer =
          static_cast<sgpp::datadriven::DBMatOnlineDE_SMW*>(&*this);
      if (this_SMW_pointer->getB().getNcols() > 1) {
        use_B_size = true;
        B_size = this_SMW_pointer->getB().getNcols();
      }
    }

    // Compute right hand side of the equation:
    size_t numberOfPoints = m.getNrows();
    DataVector b(use_B_size ? B_size : lhsMatrix.getNcols());

    if (b.getSize() != grid.getSize()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatOnlineDE::computeBFromBatch: b doesn't match size of system "
          "matrix");
    }

    std::unique_ptr<sgpp::base::OperationMultipleEval> B(
        (offlineObject.interactions.size() == 0)
            ? sgpp::op_factory::createOperationMultipleEval(grid, m)
            : sgpp::op_factory::createOperationMultipleEvalInter(
                  grid, m, offlineObject.interactions));

    // Bt * 1
    DataVector y(numberOfPoints, 1.0);
    B->multTranspose(y, b);

    // Perform permutation because of decomposition (LU)
    if (densityEstimationConfig.decomposition_ == MatrixDecompositionType::LU) {
#ifdef USE_GSL
      static_cast<DBMatOfflineLU&>(offlineObject).permuteVector(b);
#else
      throw algorithm_exception("built without GSL");
#endif /*USE_GSL*/
    }
    return b;
  }
  return DataVector();
}

DataVectorDistributed DBMatOnlineDE::computeBFromBatchParallel(
    DataMatrix& m, Grid& grid,
    const DensityEstimationConfiguration& densityEstimationConfig,
    const ParallelConfiguration& parallelConfig,
    std::shared_ptr<BlacsProcessGrid> processGrid) {
  if (m.getNrows() > 0) {
    DataMatrix& lhsMatrix = offlineObject.getDecomposedMatrix();

    // in case OrthoAdapt, the current size is not lhs size, but B size
    bool use_B_size = false;
    size_t B_size = 0;
    sgpp::datadriven::DBMatOnlineDEOrthoAdapt* this_OrthoAdapt_pointer;
    if (densityEstimationConfig.decomposition_ ==
        sgpp::datadriven::MatrixDecompositionType::OrthoAdapt) {
      this_OrthoAdapt_pointer =
          static_cast<sgpp::datadriven::DBMatOnlineDEOrthoAdapt*>(&*this);
      if (this_OrthoAdapt_pointer->getB().getNcols() > 1) {
        use_B_size = true;
        B_size = this_OrthoAdapt_pointer->getB().getNcols();
      }
    }

    // in case of SMW_ortho or SMW_chol, current size is B size and
    // also B is not used, bus B_distributed_
    sgpp::datadriven::DBMatOnlineDE_SMW* this_SMW_pointer;
    if (densityEstimationConfig.decomposition_ ==
            sgpp::datadriven::MatrixDecompositionType::SMW_ortho ||
        densityEstimationConfig.decomposition_ ==
            sgpp::datadriven::MatrixDecompositionType::SMW_chol) {
      this_SMW_pointer =
          static_cast<sgpp::datadriven::DBMatOnlineDE_SMW*>(&*this);
      if (this_SMW_pointer->getBDistributed().getGlobalCols() > 1) {
        use_B_size = true;
        B_size = this_SMW_pointer->getBDistributed().getGlobalCols();
      }
    }

    // Compute right hand side of the equation:
    size_t numberOfPoints = m.getNrows();

    size_t bSize = use_B_size ? B_size : lhsMatrix.getNcols();

    DataVectorDistributed b(processGrid, bSize, parallelConfig.rowBlockSize_);

    if (b.getGlobalRows() != grid.getSize()) {
      throw sgpp::base::algorithm_exception(
          "In DBMatOnlineDE::computeDensityFunctionParallel: b doesn't match "
          "size of system "
          "matrix");
    }

    OperationMultipleEvalConfiguration opConfig(
        OperationMultipleEvalType::SCALAPACK);

    if (offlineObject.interactions.size() != 0) {
      throw sgpp::base::not_implemented_exception(
          "Parallel evaluation operation not yet implemented for offline "
          "objects with "
          "interations");
    }
    std::unique_ptr<OperationMultipleEvalDistributed> B(
        static_cast<OperationMultipleEvalDistributed*>(
            sgpp::op_factory::createOperationMultipleEval(grid, m, opConfig)));

    // Bt * 1
    DataVector y(numberOfPoints, 1.0);
    B->multTransposeDistributed(y, b);

    return b;
  }
  return DataVectorDistributed(processGrid, 0, 1);
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

double DBMatOnlineDE::eval(DataVector& alpha, const DataVector& p, Grid& grid,
                           bool force) {
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

void DBMatOnlineDE::eval(DataVector& alpha, DataMatrix& values,
                         DataVector& results, Grid& grid, bool force) {
  if (functionComputed || force == true) {
    std::unique_ptr<sgpp::base::OperationMultipleEval> opEval(
        (offlineObject.interactions.size() == 0)
            ? sgpp::op_factory::createOperationMultipleEval(grid, values)
            : sgpp::op_factory::createOperationMultipleEvalInter(
                  grid, values, offlineObject.interactions));
    opEval->eval(alpha, results);
    results.mult(normFactor);
  } else {
    throw algorithm_exception("Density function not computed, yet!");
  }
}

void DBMatOnlineDE::evalParallel(DataVector& alpha, DataMatrix& values,
                                 DataVectorDistributed& results, Grid& grid,
                                 bool force) {
  if (functionComputed || force == true) {
    OperationMultipleEvalConfiguration opConfig(
        OperationMultipleEvalType::SCALAPACK);
    if (offlineObject.interactions.size() != 0) {
      throw sgpp::base::not_implemented_exception(
          "Parallel evaluation operation not yet implemented for offline "
          "objects with "
          "interations");
    }
    std::unique_ptr<OperationMultipleEvalDistributed> opEval(
        static_cast<OperationMultipleEvalDistributed*>(
            sgpp::op_factory::createOperationMultipleEval(grid, values,
                                                          opConfig)));

    opEval->evalParallel(alpha, results);
    results.scale(normFactor);
  } else {
    throw algorithm_exception("Density function not yet computed!");
  }
}

bool DBMatOnlineDE::isComputed() { return functionComputed; }

void DBMatOnlineDE::setBeta(double newBeta) { beta = newBeta; }

double DBMatOnlineDE::getBeta() { return beta; }

double DBMatOnlineDE::normalize(DataVector& alpha, Grid& grid, size_t samples) {
  this->normFactor = 1.;
  double sum = 0.;
  DataVector p(this->oDim);
  srand(static_cast<unsigned int>(time(nullptr)));
  for (size_t i = 0; i < samples; i++) {
    for (size_t j = 0; j < this->oDim; j++)
      p[j] = (static_cast<double>(rand()) / RAND_MAX);
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

void DBMatOnlineDE::syncDistributedDecomposition(
    std::shared_ptr<BlacsProcessGrid> processGrid,
    const ParallelConfiguration& parallelConfig) {
  offlineObject.syncDistributedDecomposition(processGrid, parallelConfig);
}

void DBMatOnlineDE::resetTraining() {
  functionComputed = false;
  localVectorsInitialized = false;
  distributedVectorsInitialized = false;
}

}  // namespace datadriven
}  // namespace sgpp
