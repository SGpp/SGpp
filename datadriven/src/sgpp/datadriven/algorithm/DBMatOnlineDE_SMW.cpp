// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#endif /* USE_GSL */

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMS_SMW.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE_SMW.hpp>

#include <algorithm>
#include <functional>
#include <list>
#include <vector>

namespace sgpp {
namespace datadriven {

DBMatOnlineDE_SMW::DBMatOnlineDE_SMW(DBMatOffline& offline, Grid& grid, double lambda, double beta)
    : sgpp::datadriven::DBMatOnlineDE(offline, grid, lambda, beta) {
  if (offline.getDecompositionType() != sgpp::datadriven::MatrixDecompositionType::OrthoAdapt &&
      offline.getDecompositionType() != sgpp::datadriven::MatrixDecompositionType::Chol) {
    throw sgpp::base::algorithm_exception(
        "In DBMatOnlineDE_SMW::DBMatOnlineDE_SMW: !!offline!! object has wrong "
        "decomposition type: DecompositionType::OrthoAdapt or ::Chol needed!");
  } else {
    // initial dummy space for information of refinement/coarsening
    this->b_adapt_matrix_ = sgpp::base::DataMatrix(1, 1);
    // note: if the size of b_adapt_matrix_ here is not initialized with 1, other functions
    // which rely on first initialization with size 1 may not function correctly
    this->b_is_refined = false;
    this->refined_points_ = {};
    this->current_refine_index = 0;
  }
}

std::vector<size_t> DBMatOnlineDE_SMW::updateSystemMatrixDecomposition(
    DensityEstimationConfiguration& densityEstimationConfig, Grid& grid, size_t numAddedGridPoints,
    std::list<size_t> deletedGridPointIndices, double lambda) {
  // points not possible to coarsen
  std::vector<size_t> return_vector = {};

  // coarsening:
  // split the valid coarsen indices and the non valid ones and do sherman-morrison coarsening
  if (!deletedGridPointIndices.empty()) {
    // indices of coarsened points and their corresponding slot
    std::vector<size_t> coarsen_points = {};
    size_t offMatrixSize = offlineObject.getGridSize();
    while (!deletedGridPointIndices.empty()) {
      size_t cur = deletedGridPointIndices.back();
      // check, which points can/cannot be coarsened
      if (cur < offMatrixSize) {  // cannot be coarsened, because in offline's initial lhs matrix
        return_vector.push_back(cur);
      } else {
        coarsen_points.push_back(cur);
      }
      deletedGridPointIndices.pop_back();  // throw away already considered indices
    }

    if (!coarsen_points.empty()) {
      // todo(): optimize size of matrix X when coarsening
      // break-even point, k times single point VS. k points at once
      DataMatrix X(grid.getSize(), this->current_refine_index, 0.0);
      compute_L2_coarsen_matrix(X, grid, coarsen_points);
      this->smw_adapt(X, 0, false, coarsen_points);
    }
  }

  // refinement:
  if (numAddedGridPoints > 0) {
    // refine/coarsen matrix X for smw formula
    DataMatrix X(grid.getSize(), numAddedGridPoints);
    compute_L2_refine_matrix(X, grid, numAddedGridPoints, lambda);
    this->smw_adapt(X, numAddedGridPoints, true);
  }

  return return_vector;
}

void DBMatOnlineDE_SMW::solveSLE(DataVector& alpha, DataVector& b, Grid& grid,
                                 DensityEstimationConfiguration& densityEstimationConfig,
                                 bool do_cv) {
  // create solver
  sgpp::datadriven::DBMatDMS_SMW* solver = new sgpp::datadriven::DBMatDMS_SMW();
  // solve the created system
  alpha.resizeZero(b.getSize());
  solver->solve(this->offlineObject.getInverseMatrix(), this->getB(), b, alpha);

  free(solver);
}

void DBMatOnlineDE_SMW::solveSLEParallel(DataVectorDistributed& alpha, DataVectorDistributed& b,
                                         Grid& grid,
                                         DensityEstimationConfiguration& densityEstimationConfig,
                                         bool do_cv) {
  sgpp::datadriven::DBMatOfflineOrthoAdapt* offline =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&this->offlineObject);
  DataMatrixDistributed TinvDistributed = offline->getTinvDistributed();

  DataMatrixDistributed QDistributed = offline->getQDistributed();

  DataMatrixDistributed BDistributed = this->getBDistributed();

  // create solver
  std::unique_ptr<sgpp::datadriven::DBMatDMSOrthoAdapt> solver =
      std::make_unique<sgpp::datadriven::DBMatDMSOrthoAdapt>();

  alpha.resize(b.getGlobalRows());
  solver->solveParallel(TinvDistributed, QDistributed, BDistributed, b, alpha);
}

void DBMatOnlineDE_SMW::smw_adapt(DataMatrix& X, size_t newPoints, bool refine,
                                  std::vector<size_t> coarsenIndices) {
#ifdef USE_GSL
  // dimension of offline's lhs matrix and its inverse
  size_t offMatrixSize = this->offlineObject.getGridSize();

  // check, if offline object has been decomposed yet
  if (offMatrixSize == 0) {
    throw sgpp::base::algorithm_exception(
        "In DBMatOnlineDE_SMW::smw_adapt:\noffline object wasn't decomposed "
        "yet, so can't perform refinement/coarsening.");
  }

  // determine the final size of the B matrix (additive component in smw formula)
  size_t oldSize = this->b_is_refined ? this->b_adapt_matrix_.getNcols() : offMatrixSize;
  size_t newSize = refine ? (oldSize + newPoints) : (oldSize - coarsenIndices.size());
  size_t adaptSteps = (newSize > oldSize) ? (newSize - oldSize) : (oldSize - newSize);

  // allocate additional space for b_adapt_matrix_ and fill new diagonal entries with ones
  // note: only done in refining! Coarsening will resize at the end of function
  if (refine) {
    this->b_adapt_matrix_.resizeQuadratic(newSize);
    for (size_t i = oldSize; i < newSize; i++) {
      this->b_adapt_matrix_.set(i, i, 1.0);
    }
  }

  /************************************************************
   * BEGIN OF SHERMAN-MORRISON-WOODBURRY
   *
   * Phase 1: B~ = B - (A^-1 + B) X (I + E^t (A^-1 + B) X)^-1 E^t (A^-1 + B)
   *
   ************************************************************/

  // create view of B, will store values of B~
  gsl_matrix_view b_adapt_full_view =
      gsl_matrix_view_array(this->b_adapt_matrix_.getPointer(), this->b_adapt_matrix_.getNrows(),
                            this->b_adapt_matrix_.getNcols());

  // adapt X's diagonal, lambda already added before
  for (size_t k = X.getNrows() - X.getNcols(); k < X.getNrows(); k++) {
    // <x, x> = 0 => x=0, therefore:
    // diag_k = 0 -> ignoring k-th row/column while coarsening
    double val = (X.get(k, k - X.getNrows() + X.getNcols()) != 0)
                     ? (X.get(k, k - X.getNrows() + X.getNcols()) - 1) * 0.5
                     : 0;
    X.set(k, k - X.getNrows() + X.getNcols(), val);
  }

  // create E
  // todo(): replace operations done with E/E^t by matrix_subviews
  //         could give performance boost
  DataMatrix E(X.getNrows(), X.getNcols(), 0.0);
  for (size_t k = E.getNrows() - E.getNcols(); k < E.getNrows(); k++) {
    E.set(k, k - E.getNrows() + E.getNcols(), 1.0);
  }

  // create view of A^-1
  gsl_matrix_view A_inv_view = gsl_matrix_view_array(
      this->offlineObject.getInverseMatrix().getPointer(), offMatrixSize, offMatrixSize);

  // A^-1 + B (calculate in size of A^-1, store in size of B)
  DataMatrix AB(this->b_adapt_matrix_);
  gsl_matrix_view AB_view = gsl_matrix_view_array(AB.getPointer(), this->b_adapt_matrix_.getNrows(),
                                                  this->b_adapt_matrix_.getNcols());
  gsl_matrix_view AB_sub_view =
      gsl_matrix_submatrix(&AB_view.matrix, 0, 0, offMatrixSize, offMatrixSize);
  gsl_matrix_add(&AB_sub_view.matrix, &A_inv_view.matrix);

  // (A^-1 + B) * X
  gsl_matrix* AX = gsl_matrix_alloc(X.getNrows(), X.getNcols());
  gsl_matrix_view X_view = gsl_matrix_view_array(X.getPointer(), X.getNrows(), X.getNcols());
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &AB_view.matrix, &X_view.matrix, 0.0, AX);

  // E^t * (A^-1 + B)
  gsl_matrix* EA = gsl_matrix_alloc(X.getNcols(), X.getNrows());
  gsl_matrix_view E_view = gsl_matrix_view_array(E.getPointer(), X.getNrows(), X.getNcols());
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &E_view.matrix, &AB_view.matrix, 0.0, EA);

  // I + E^t * (A^-1 + B) * X
  gsl_matrix* TO_INV = gsl_matrix_alloc(X.getNcols(), X.getNcols());
  gsl_matrix_set_identity(TO_INV);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &E_view.matrix, AX, 1.0, TO_INV);

  // (I + E^t * (A^-1 + B) * X)^-1
  gsl_permutation* P = gsl_permutation_alloc(X.getNcols());
  int* signum = new int[X.getNcols()];
  gsl_linalg_LU_decomp(TO_INV, P, signum);
  gsl_matrix* INV = gsl_matrix_alloc(X.getNcols(), X.getNcols());
  gsl_linalg_LU_invert(TO_INV, P, INV);

  // (A^-1 + B) X (I + E^t (A^-1 + B) X)^-1
  gsl_matrix* AXINV = gsl_matrix_alloc(X.getNrows(), X.getNcols());
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AX, INV, 0.0, AXINV);

  // (A^-1 + B) X (I + E^t (A^-1 + B) X)^-1 E^t (A^-1 + B)
  gsl_matrix* AXINVEA = gsl_matrix_alloc(X.getNrows(), X.getNrows());
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AXINV, EA, 0.0, AXINVEA);

  // B - (A^-1 + B) X (I + E^t (A^-1 + B) X)^-1 E^t (A^-1 + B)
  // effectively stores final values of Phase 1 in b_adapt_matrix_
  gsl_matrix_sub(&b_adapt_full_view.matrix, AXINVEA);

  /************************************************************
   *
   * Phase 2: B = B~ - (A^-1 + B~) E (I + X^t (A^-1 + B~) E)^-1 X^t (A^-1 + B~)
   *
   ************************************************************/

  // A^-1 + B~ (calculate in size of A^-1, store in size of B~)
  DataMatrix ABtilde(this->b_adapt_matrix_);
  gsl_matrix_view ABtilde_view = gsl_matrix_view_array(
      ABtilde.getPointer(), this->b_adapt_matrix_.getNrows(), this->b_adapt_matrix_.getNcols());
  gsl_matrix_view ABtilde_sub_view =
      gsl_matrix_submatrix(&ABtilde_view.matrix, 0, 0, offMatrixSize, offMatrixSize);
  gsl_matrix_add(&ABtilde_sub_view.matrix, &A_inv_view.matrix);

  // (A^-1 + B~) * E    und    X^t*(A+B)
  gsl_matrix* AE = gsl_matrix_alloc(X.getNrows(), X.getNcols());
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &ABtilde_view.matrix, &E_view.matrix, 0.0, AE);

  // X^t * (A^-1 + B~)
  gsl_matrix* XA = gsl_matrix_alloc(X.getNcols(), X.getNrows());
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &X_view.matrix, &ABtilde_view.matrix, 0.0, XA);

  // I + X^t * (A^-1 + B) * E
  gsl_matrix_set_identity(TO_INV);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &X_view.matrix, AE, 1.0, TO_INV);

  // (I + X^t * (A^-1 + B) * E)^-1
  gsl_permutation* Ptilde = gsl_permutation_alloc(X.getNcols());
  int* signumtilde = new int[X.getNcols()];
  gsl_linalg_LU_decomp(TO_INV, Ptilde, signumtilde);
  gsl_linalg_LU_invert(TO_INV, Ptilde, INV);

  // (A^-1 + B) E (I + X^t * (A^-1 + B) * E)^-1
  gsl_matrix* AEINV = gsl_matrix_alloc(X.getNrows(), X.getNcols());
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AE, INV, 0.0, AEINV);

  // (A^-1 + B) E (I + X^t (A^-1 + B) E)^-1 X^t (A^-1 + B)
  gsl_matrix* AEINVXA = gsl_matrix_alloc(X.getNrows(), X.getNrows());
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AEINV, XA, 0.0, AEINVXA);

  // B~ - (A^-1 + B~) E (I + X^t (A^-1 + B~) E)^-1 X^t (A^-1 + B~)
  // effectively stores final values of Phase 2 in b_adapt_matrix_
  gsl_matrix_sub(&b_adapt_full_view.matrix, AEINVXA);

  /*****
   *
   * END OF SHERMAN-MORRISON-WOODBURRY
   *
   *****/

  // If points were coarsened the b_adapt_matrix will now have empty rows and columns
  // on the indices of the coarsened points. In the following algorithm, the symmetry
  // of b_adapt_matrix_ is used to glue the blocks together again.
  if (!refine) {
    std::sort(coarsenIndices.begin(), coarsenIndices.end());
    size_t vi = 1;  // skips points
    sgpp::base::DataVector copied_row(oldSize);

    for (size_t i = coarsenIndices[0]; i < oldSize - coarsenIndices.size(); i++) {
      // find amount of successive indices, to skip a block of zeroes
      while (i + vi == coarsenIndices[vi]) {
        vi++;
      }

      // set current row and column to copied one
      this->b_adapt_matrix_.getRow(i + vi, copied_row);
      this->b_adapt_matrix_.setRow(i, copied_row);
      this->b_adapt_matrix_.setColumn(i, copied_row);  // symmetry

      // copy diagonal elements too
      double value = this->b_adapt_matrix_.get(i + vi, i + vi);
      this->b_adapt_matrix_.set(i, i, value);
    }

    this->b_adapt_matrix_.resizeQuadratic(newSize);
    // size fitting ends here

    // remove coarsened points from online's internal storage of refined_vectors_
    for (size_t k = adaptSteps; k > 0; k--) {
      this->refined_points_.erase(this->refined_points_.begin() + coarsenIndices[k - 1] -
                                  offMatrixSize);

      // adjust current_refine_index of online object
      this->current_refine_index--;
    }
  }

  // move the index of the container to the end of it
  if (refine) {
    this->current_refine_index += adaptSteps;
  }

  // determine, if any refined information now is contained in matrix b_adapt
  this->b_is_refined = this->b_adapt_matrix_.getNcols() > offMatrixSize;

  return;
#endif /* USE_GSL */
}

void DBMatOnlineDE_SMW::compute_L2_refine_matrix(DataMatrix& X, Grid& grid, size_t newPoints,
                                                 double newLambda) {
  size_t gridSize = grid.getSize();

  if (X.getNcols() != newPoints || X.getNrows() != gridSize) {
    throw sgpp::base::algorithm_exception(
        "In DBMatOnlineDE_SMW::compute_L2_refine_matrix:\n"
        "The passed matrix container doesn't have the correct size.\n");
  }

  this->offlineObject.compute_L2_refine_vectors(&X, &grid, newPoints);

  // add lambda to diagonal elements
  for (size_t i = gridSize - newPoints; i < gridSize; i++) {
    double res = X.get(i, i - gridSize + newPoints);
    X.set(i, i - gridSize + newPoints, res + newLambda);
  }

  // put the new points inside the internal storage
  // this is needed for coarsening points later, without recalculating
  for (size_t i = 0; i < newPoints; i++) {
    sgpp::base::DataVector vec(gridSize);
    X.getColumn(i, vec);
    this->refined_points_.push_back(vec);
  }

  // fill in the new L2 products in the older points, and resize them
  for (size_t i = 0; i < this->current_refine_index; i++) {
    this->refined_points_[i].resize(gridSize);
    for (size_t j = this->current_refine_index; j < this->current_refine_index + newPoints; j++) {
      double xxx = this->refined_points_[j].get(i);
      this->refined_points_[i].set(j, xxx);
    }
  }
}

void DBMatOnlineDE_SMW::compute_L2_coarsen_matrix(DataMatrix& X, Grid& grid,
                                                  std::vector<size_t> coarsen_indices) {
  if (X.getNrows() != this->b_adapt_matrix_.getNrows()) {
    throw sgpp::base::algorithm_exception(
        "in DBMatOnlineDE_SMW::compute_L2_coarsen_matrix:\n matrix X doesn't match B");
  }
  for (size_t i : coarsen_indices) {
    X.setColumn(i - this->offlineObject.getGridSize(),
                this->refined_points_[i - this->offlineObject.getGridSize()]);
  }
}

void DBMatOnlineDE_SMW::syncDistributedDecomposition(std::shared_ptr<BlacsProcessGrid> processGrid,
                                                     const ParallelConfiguration& parallelConfig) {
#ifdef USE_SCALAPACK
  offlineObject.syncDistributedDecomposition(processGrid, parallelConfig);
  b_adapt_matrix_distributed_ = DataMatrixDistributed::fromSharedData(
      b_adapt_matrix_.data(), processGrid, b_adapt_matrix_.getNrows(), b_adapt_matrix_.getNcols(),
      parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);
#endif
  // no action needed without scalapack
}

}  // namespace datadriven
}  // namespace sgpp