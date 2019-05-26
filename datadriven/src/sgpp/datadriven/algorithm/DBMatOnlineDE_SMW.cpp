// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE_SMW.hpp>

#ifdef USE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#endif /* USE_GSL */

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>

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

  // refine/coarsen matrix X for smw formula
  DataMatrix X(grid.getSize(), numAddedGridPoints);
  // coarsening:
  // split the valid coarsen indices and the non valid ones and do sherman-morrison coarsening
  if (!deletedGridPointIndices.empty()) {
    std::cout << "\n\nCOARSEN CASE\n\n\n";
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
      compute_L2_coarsen_matrix(X, grid, coarsen_points);
      this->smw_adapt(X, 0, false, coarsen_points);
    }
  }

  // refinement:
  if (numAddedGridPoints > 0) {
    compute_L2_refine_matrix(X, grid, numAddedGridPoints, lambda);
    this->smw_adapt(X, numAddedGridPoints, true);
  }

  return return_vector;
}

void DBMatOnlineDE_SMW::solveSLE(DataVector& alpha, DataVector& b, Grid& grid,
                                 DensityEstimationConfiguration& densityEstimationConfig,
                                 bool do_cv) {
  sgpp::datadriven::DBMatOfflineOrthoAdapt* offline =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&this->offlineObject);
  // create solver
  sgpp::datadriven::DBMatDMSOrthoAdapt* solver = new sgpp::datadriven::DBMatDMSOrthoAdapt();
  // solve the created system
  alpha.resizeZero(b.getSize());
  solver->solve(offline->getTinv(), offline->getQ(), this->getB(), b, alpha);

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
  std::cout << "\nentered smw_adapt\n";
#ifdef USE_GSL

  sgpp::datadriven::DBMatOfflineOrthoAdapt* offlinePtr =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&this->offlineObject);

  // dimension of offline's lhs matrix A^{-1} = Q * T^{-1} * Q^t
  size_t offMatrixSize = offlinePtr->getGridSize();

  // check, if offline object has been decomposed yet
  if (offMatrixSize == 0) {
    throw sgpp::base::algorithm_exception(
        "In DBMatOnlineDE_SMW::smw_adapt:\noffline object wasn't decomposed "
        "yet, so can't perform adapting");
  }

  // determine the final size of the B matrix (additive component in smw formula)
  size_t oldSize = this->b_is_refined ? this->b_adapt_matrix_.getNcols() : offMatrixSize;
  size_t newSize = refine ? (oldSize + newPoints) : (oldSize - coarsenIndices.size());
  size_t adaptSteps = (newSize > oldSize) ? (newSize - oldSize) : (oldSize - newSize);

  // allocate space for b_adapt_matrix_ and fill new diagonal entries with ones
  // note: only done in refining! Coarsening will resize at the end of function
  if (refine) {
    this->b_adapt_matrix_.resizeQuadratic(newSize);
    for (size_t i = oldSize; i < newSize; i++) {
      this->b_adapt_matrix_.set(i, i, 1.0);
    }
  }

  // create view of the whole matrix b_adapt, needed to later create submatrix_view
  gsl_matrix_view b_adapt_full_view =
      gsl_matrix_view_array(this->b_adapt_matrix_.getPointer(), this->b_adapt_matrix_.getNrows(),
                            this->b_adapt_matrix_.getNcols());

  // apply sherman-morrison-woodburry formula, which computes k refinements at once
  // TODO
  //

  // compute A^-1 explicitly
  // TODO(dima) make function in offline objects for supported types (ortho + chol)

  // create X (already done)
  // adapt X's diagonal, differ between refine/coarsen
  for (size_t k = X.getNrows() - X.getNcols(); k < X.getNrows(); k++) {
    double val = (X.get(k, k - X.getNrows() + X.getNcols()) - 1) * 0.5;
    X.set(k, k - X.getNrows() + X.getNcols(), val);
  }

  // create E, (*+1 if refine, *-1 if coarsen)
  DataMatrix E(X.getNrows(), X.getNcols(), 0.0);
  for (size_t k = E.getNrows() - E.getNcols(); k < E.getNrows(); k++) {
    double val = refine ? 1 : -1;
    E.set(k, k - E.getNrows() + E.getNcols(), val);
  }

  std::cout << "\nstarting gsl calculations\n";
  if (X.getNrows() != this->b_adapt_matrix_.getNcols()) {
    throw sgpp::base::algorithm_exception("in adapt:\nX.getNrows != B.size\n");
  }

  std::cout << "\ncreate views for Q and T_inv\n";
  if (offMatrixSize != offlinePtr->getTinv().getNcols()) {
    throw sgpp::base::algorithm_exception("in adapt:\noffMatrixSize != T und Q\n");
  }
  // create matrix views and buffers for interim values
  // view of T^{-1} of the offline object
  gsl_matrix_view t_inv_view =
      gsl_matrix_view_array(offlinePtr->getTinv().getPointer(), offMatrixSize, offMatrixSize);

  // view of Q of the offline object
  gsl_matrix_view q_view =
      gsl_matrix_view_array(offlinePtr->getQ().getPointer(), offMatrixSize, offMatrixSize);
  ////////////////////////////////////
  //
  // smw phase 1, calculate B_tilde
  std::cout << "\nSMW PHASE 1\n\n";
  //
  ////////////////////////////////////
  // calculate A^-1 explicitly
  gsl_matrix* interim = gsl_matrix_alloc(offMatrixSize, offMatrixSize);
  gsl_matrix* A_inv = gsl_matrix_alloc(offMatrixSize, offMatrixSize);

  std::cout << "\nberechne A^-1 der offline matrix\n";
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &q_view.matrix, &t_inv_view.matrix, 0.0, interim);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, interim, &q_view.matrix, 0.0, A_inv);

  std::cout << "\nA^-1 + B\n";
  // A^-1 + B (in size of B)
  DataMatrix AB(this->b_adapt_matrix_);
  gsl_matrix_view AB_view = gsl_matrix_view_array(AB.getPointer(), this->b_adapt_matrix_.getNrows(),
                                                  this->b_adapt_matrix_.getNcols());
  gsl_matrix_view AB_sub_view =
      gsl_matrix_submatrix(&AB_view.matrix, 0, 0, offMatrixSize, offMatrixSize);
  gsl_matrix_add(&AB_sub_view.matrix, A_inv);

  std::cout << "\n.*X und E^t*.\n";
  // (A+B)*X    und    E^t*(A+B)
  gsl_matrix* AX = gsl_matrix_alloc(X.getNrows(), X.getNcols());
  gsl_matrix* EA = gsl_matrix_alloc(X.getNcols(), X.getNrows());
  gsl_matrix_view X_view = gsl_matrix_view_array(X.getPointer(), X.getNrows(), X.getNcols());
  gsl_matrix_view E_view = gsl_matrix_view_array(E.getPointer(), X.getNrows(), X.getNcols());

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &AB_view.matrix, &X_view.matrix, 0.0, AX);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &E_view.matrix, &AB_view.matrix, 0.0, EA);

  std::cout << "\n (I + E^t*.)\n";
  // I + E^t * (A+B) * X
  gsl_matrix* TO_INV = gsl_matrix_alloc(X.getNcols(), X.getNcols());
  gsl_matrix_set_identity(TO_INV);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &E_view.matrix, AX, 1.0, TO_INV);

  std::cout << "\ninvertieren mit gsl_LU\n";
  // invert (I + E^t * (A + B) * X)
  gsl_permutation* P = gsl_permutation_alloc(X.getNcols());
  int* signum = new int[X.getNcols()];
  gsl_linalg_LU_decomp(TO_INV, P, signum);
  gsl_matrix* INV = gsl_matrix_alloc(X.getNcols(), X.getNcols());
  gsl_linalg_LU_invert(TO_INV, P, INV);

  std::cout << "\n AX*INV\n";
  gsl_matrix* AXINV = gsl_matrix_alloc(X.getNrows(), X.getNcols());
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AX, INV, 0.0, AXINV);

  std::cout << "\n AX*INV*EA\n";
  gsl_matrix* AXINVEA = gsl_matrix_alloc(X.getNrows(), X.getNrows());
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AXINV, EA, 0.0, AXINVEA);
  gsl_matrix_sub(&b_adapt_full_view.matrix, AXINVEA);
  ////////////////////////////////////
  //
  // smw phase 2, calculate B
  std::cout << "\nSMW PHASE 2\n\n";
  //
  ////////////////////////////////////
  std::cout << "\nA^-1 + B_tilde\n";
  // A^-1 + B (in size of B)
  DataMatrix ABtilde(this->b_adapt_matrix_);
  gsl_matrix_view ABtilde_view = gsl_matrix_view_array(
      ABtilde.getPointer(), this->b_adapt_matrix_.getNrows(), this->b_adapt_matrix_.getNcols());
  gsl_matrix_view ABtilde_sub_view =
      gsl_matrix_submatrix(&ABtilde_view.matrix, 0, 0, offMatrixSize, offMatrixSize);
  gsl_matrix_add(&ABtilde_sub_view.matrix, A_inv);

  std::cout << "\n.*E und X^t*.\n";
  // (A+B)*E    und    X^t*(A+B)
  gsl_matrix* AE = gsl_matrix_alloc(X.getNrows(), X.getNcols());
  gsl_matrix* XA = gsl_matrix_alloc(X.getNcols(), X.getNrows());

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &ABtilde_view.matrix, &E_view.matrix, 0.0, AE);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &X_view.matrix, &ABtilde_view.matrix, 0.0, XA);

  std::cout << "\n (I + X^t*.)\n";
  // I + X^t * (A+B) * E
  gsl_matrix_set_identity(TO_INV);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &X_view.matrix, AE, 1.0, TO_INV);

  std::cout << "\ninvertieren mit gsl_LU\n";
  // invert (I + X^t * (A + B) * E)
  gsl_permutation* Ptilde = gsl_permutation_alloc(X.getNcols());
  int* signumtilde = new int[X.getNcols()];
  gsl_linalg_LU_decomp(TO_INV, Ptilde, signumtilde);
  gsl_linalg_LU_invert(TO_INV, Ptilde, INV);

  std::cout << "\n AE*INV\n";
  gsl_matrix* AEINV = gsl_matrix_alloc(X.getNrows(), X.getNcols());
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AE, INV, 0.0, AEINV);

  std::cout << "\n AE*INV*XA\n";
  gsl_matrix* AEINVXA = gsl_matrix_alloc(X.getNrows(), X.getNrows());
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AEINV, XA, 0.0, AEINVXA);
  gsl_matrix_sub(&b_adapt_full_view.matrix, AEINVXA);
  // TODO END

  std::cout << "\nending gsl calculations\n";

  // If points were coarsened the b_adapt_matrix will now have empty rows and columns
  // on the indices of the coarsened points. In the following algorithm, the symmetry
  // of b_adapt_matrix is used to fit the blocks together again.
  if (!refine) {
    std::sort(coarsenIndices.begin(), coarsenIndices.end());
    size_t vi = 1;  // skipps points
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

  // std::cout << "\n\n\n B = \n\n" << this->b_adapt_matrix_.toString() << "\n\n";

  std::cout << "\nexited smw_adapt\n";
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
  // TODO(dima):
  // baue matrix X für den coarsen fall. Matrix hat full size, und an nicht angegebenen indices
  // gleicht sie der einheitsmatrix.
  // note: kein vorzeichenhandling! einfach nur die matrix X für ausgewählte indices (siehe
  // refine_matrix)
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
