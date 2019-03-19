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

#include <sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <algorithm>
#include <functional>
#include <list>
#include <vector>

namespace sgpp {
namespace datadriven {

DBMatOnlineDEOrthoAdapt::DBMatOnlineDEOrthoAdapt(DBMatOffline& offline, Grid& grid, double lambda,
                                                 double beta)
    : sgpp::datadriven::DBMatOnlineDE(offline, grid, lambda, beta) {
  if (offline.getDecompositionType() != sgpp::datadriven::MatrixDecompositionType::OrthoAdapt) {
    throw sgpp::base::algorithm_exception(
        "In DBMatOnlineDEOrthoAdapt::DBMatOnlineDEOrthoAdapt: offline object has wrong "
        "decomposition type: DecompositionType::OrthoAdapt needed!");
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

std::vector<size_t> DBMatOnlineDEOrthoAdapt::updateSystemMatrixDecomposition(
    DensityEstimationConfiguration& densityEstimationConfig, Grid& grid, size_t numAddedGridPoints,
    std::list<size_t> deletedGridPointIndices, double lambda) {
  // points not possible to coarsen
  std::vector<size_t> return_vector = {};

  // coarsening:
  // split the valid coarsen indices and the non valid ones and do sherman-morrison coarsening
  if (!deletedGridPointIndices.empty()) {
    // indices of coarsened points and their corresponding slot
    std::vector<size_t> coarsen_points = {};
    size_t dima = offlineObject.getGridSize();
    while (!deletedGridPointIndices.empty()) {
      size_t cur = deletedGridPointIndices.back();
      // check, which points can/cannot be coarsened
      if (cur < dima) {  // cannot be coarsened, because part of the offline's initial lhs matrix
        return_vector.push_back(cur);
      } else {
        coarsen_points.push_back(cur);
      }
      deletedGridPointIndices.pop_back();  // throw away already considered indices
    }

    if (!coarsen_points.empty()) {
      this->sherman_morrison_adapt(0, false, coarsen_points);
    }
  }

  // refinement:
  if (numAddedGridPoints > 0) {
    // get the L2_gridvectors from the new refined grid and do Sherman-Morrison refinement
    this->compute_L2_gridvectors(grid, numAddedGridPoints, lambda);
    this->sherman_morrison_adapt(numAddedGridPoints, true);
  }

  return return_vector;
}

void DBMatOnlineDEOrthoAdapt::solveSLE(DataVector& alpha, DataVector& b, Grid& grid,
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

void DBMatOnlineDEOrthoAdapt::solveSLEParallel(
    DataVectorDistributed& alpha, DataVectorDistributed& b, Grid& grid,
    DensityEstimationConfiguration& densityEstimationConfig,
    const ParallelConfiguration& parallelConfig, std::shared_ptr<BlacsProcessGrid> processGrid,
    bool do_cv) {
  // TODO(jan) parallelize offline version at this point?
  sgpp::datadriven::DBMatOfflineOrthoAdapt* offline =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&this->offlineObject);
  DataMatrix Tinv = offline->getTinv();
  DataMatrixDistributed TinvDistributed =
      DataMatrixDistributed(Tinv.data(), processGrid, Tinv.getNrows(), Tinv.getNcols(),
                            parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);

  DataMatrix Q = offline->getQ();
  DataMatrixDistributed QDistributed =
      DataMatrixDistributed(Q.data(), processGrid, Q.getNrows(), Q.getNcols(),
                            parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);

  DataMatrix B = this->getB();
  DataMatrixDistributed BDistributed =
      DataMatrixDistributed(B.data(), processGrid, B.getNrows(), B.getNcols(),
                            parallelConfig.rowBlockSize_, parallelConfig.columnBlockSize_);

  // create solver
  sgpp::datadriven::DBMatDMSOrthoAdapt* solver = new sgpp::datadriven::DBMatDMSOrthoAdapt();

  alpha.resize(b.getGlobalRows());
  solver->solveParallel(TinvDistributed, QDistributed, BDistributed, b, alpha, processGrid,
                        parallelConfig);
}

void DBMatOnlineDEOrthoAdapt::sherman_morrison_adapt(size_t newPoints, bool refine,
                                                     std::vector<size_t> coarsenIndices) {
#ifdef USE_GSL
  sgpp::datadriven::DBMatOfflineOrthoAdapt* offlinePtr =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&this->offlineObject);

  // dimension of offline's lhs matrix A^{-1} = Q * T^{-1} * Q^t
  size_t dima = offlinePtr->getGridSize();

  // check, if offline object has been decomposed yet
  if (dima == 0) {
    throw sgpp::base::algorithm_exception(
        "In DBMatOnlineDEOrthoAdapt::sherman_morrison_adapt:\noffline object wasn't decomposed "
        "yet, so can't perform sherman_morrison_adapt");
  }

  // determine the final size of the matrix system
  size_t oldSize = this->b_is_refined ? this->b_adapt_matrix_.getNcols() : dima;
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

  // apply sherman morrison formula k times (one for each point refining/coarsening)
  for (size_t k = 0; k < adaptSteps; k++) {
    // fetch point with according size
    // refine -> size gets bigger by 1
    // coarse -> size stays the same
    size_t current_size = oldSize + (refine ? k + 1 : 0);

    // copies grid point vector from container
    sgpp::base::DataVector x =
        this->refined_points_[refine ? oldSize - dima + k : coarsenIndices[k] - dima];

    // configure x according to Sherman-Morrison formula
    if (refine) {
      // clipping off not needed entries
      x.resize(current_size);

      // note: lambda is added already on all of refinePts
      double config_x_value = (x.get(current_size - 1) - 1) * 0.5;
      x.set(current_size - 1, config_x_value);
    } else {
      // no clipping off when coarsening
      double config_x_value = (x.get(coarsenIndices[k]) - 1) * 0.5;

      // lambda is added already on all of the points in container
      x.set(coarsenIndices[k], config_x_value);
    }

    // configure unit vector depending on refine/coarsen
    // e[unit_index] = 1 when refining, -1 when coarsening
    size_t unit_index = refine ? current_size - 1 : coarsenIndices[k];

    // view of T^{-1} of the offline object
    gsl_matrix_view t_inv_view =
        gsl_matrix_view_array(offlinePtr->getTinv().getPointer(), dima, dima);

    // view of Q of the offline object
    gsl_matrix_view q_view = gsl_matrix_view_array(offlinePtr->getQ().getPointer(), dima, dima);

    // view of B of the online object, which holds all information of refinement/coarsening
    gsl_matrix_view b_adapt_view =
        gsl_matrix_submatrix(&b_adapt_full_view.matrix, 0, 0, current_size, current_size);

    // view of current point to refine/coarsen
    gsl_vector_view x_view = gsl_vector_view_array(x.getPointer(), current_size);

    // view of current point to refine/coarsen cut to offline objects matrix size
    gsl_vector_view x_cut_view = gsl_vector_view_array(x.getPointer(), dima);

    // allocating space for buffer
    gsl_vector* buffer = gsl_vector_alloc(dima);

    // allocating space for the term: Q*T^{-1}*Q^t * x_cut
    gsl_vector* x_term = gsl_vector_alloc(dima);

    // allocating space for the term: B * x
    gsl_vector* bx_term = gsl_vector_alloc(current_size);

    //##########################################################################
    //
    // first Sherman-Morrison update/downdate, calculating B_tilde
    //
    //##########################################################################

    // calculating bx_term = B * x, (B is symmetric)
    gsl_blas_dgemv(CblasNoTrans, 1.0, &b_adapt_view.matrix, &x_view.vector, 0.0, bx_term);

    // calculating x_term = Q*T^{-1}*Q^t * x_cut (the whole matrix term is symmetric)
    gsl_blas_dgemv(CblasTrans, 1.0, &q_view.matrix, &x_cut_view.vector, 0.0,
                   x_term);  // x_term = Q^t * x_cut

    gsl_blas_dgemv(CblasNoTrans, 1.0, &t_inv_view.matrix, x_term, 0.0,
                   buffer);  // buffer = T^{-1} * Q^t * x_cut

    gsl_blas_dgemv(CblasNoTrans, 1.0, &q_view.matrix, buffer, 0.0,
                   x_term);  // x_term = Q*T^{-1}*Q^t * x_cut

    // calculating the divisor of the sherman-morrison-formula: 1 + e^t * x_term + e^t * bx_term
    // note: the term e^t * Q * T^{-1} * Q^t * x is always zero,
    // because initial gridpoints cannot be refined
    double eBx = (refine ? 1.0 : -1.0) * gsl_vector_get(bx_term, unit_index);
    double divisor = 1 + eBx;
    divisor = 1.0 / divisor;  // optional optimization

    // calculating: bx_term + x_term, where x_term is "filled" with zero to fit dimensions
    // the result is stored in bx_term, as x_term is later needed
    for (size_t i = 0; i < dima; i++) {
      bx_term->data[i] += x_term->data[i];
    }

    // subtracting the matrices B - ( x_term * e_term ) / divisor
    // where the vector product is an outer product yielding a matrix

    // saving the row of B, because it gets overwritten after unit_index-th step
    sgpp::base::DataVector b_row(refine ? newSize : oldSize);
    this->b_adapt_matrix_.getRow(unit_index, b_row);

    for (size_t i = 0; i < current_size; i++) {
      for (size_t j = 0; j < current_size; j++) {
        // calculate matrix entry new value according to Sherman-Morrison-formula
        double final_value = this->b_adapt_matrix_.get(i, j) -
                             (bx_term->data[i] * (refine ? 1.0 : -1.0) * b_row.get(j) * divisor);

        // by algorithm, the column of B at the index of coarsening
        // should be zero, except for the diagonal entry
        if (!refine) {
          if (j == unit_index && i != unit_index) {
            final_value = 0.0;
          }
        }
        this->b_adapt_matrix_.set(i, j, final_value);
      }
    }

    //##########################################################################
    //
    // second Sherman-Morrison update/downdate, calculating B
    //
    //##########################################################################

    // calculating new bx_term = x^t * B_tilde = (B_tilde^t * x)^t
    // note: B is symmetric, but B_tilde is never symmetric
    gsl_blas_dgemv(CblasTrans, 1.0, &b_adapt_view.matrix, &x_view.vector, 0.0, bx_term);

    // calculating the divisor of sherman-morrison-formula: 1 + e^t * b_tilde_x_term
    eBx = (refine ? 1.0 : -1.0) * gsl_vector_get(bx_term, unit_index);
    divisor = 1 + eBx;
    divisor = 1.0 / divisor;

    // calculating: bx_tilde_term + x_term
    for (size_t i = 0; i < dima; i++) {
      bx_term->data[i] += x_term->data[i];
    }

    // subtracting the matrices B_final - ( x_term * b_tilde_term ) / divisor
    // where the vector product is an outer product yielding a matrix
    // note: as the order is flipped, i and j have to be changed AND
    // the b_term now is taken out of b_adapt_matrix_'s column instead of row

    // saving the column of B, because it gets overwritten after unit_index-th step
    sgpp::base::DataVector b_col(refine ? newSize : oldSize);
    this->b_adapt_matrix_.getColumn(unit_index, b_col);

    for (size_t i = 0; i < current_size; i++) {
      for (size_t j = 0; j < current_size; j++) {
        double final_value = this->b_adapt_matrix_.get(i, j) -
                             (bx_term->data[j] * (refine ? 1.0 : -1.0) * b_col.get(i)) * divisor;

        // By algorithm, the row of B at the index of coarsening
        // should be now the unit vector, and can be discarded after coarsening.
        // Since it will be discarded, the unit entry will be set to zero, to not
        // meddle with other dimensions of the matrix
        if (!refine) {
          if (i == unit_index) {
            final_value = 0.0;
          }
        }
        this->b_adapt_matrix_.set(i, j, final_value);
      }
    }

    gsl_vector_free(buffer);
    gsl_vector_free(x_term);
    gsl_vector_free(bx_term);
  }

  //##########################################################################
  //
  // ending Sherman-Morrison updates/downdates
  //
  //##########################################################################

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

    // remove coarsened points from online's refined_vectors
    for (size_t k = adaptSteps; k > 0; k--) {
      this->refined_points_.erase(this->refined_points_.begin() + coarsenIndices[k - 1] - dima);

      // adjust current_refine_index of online object
      this->current_refine_index--;
    }
  }

  // move the index of the container to the end of it
  if (refine) {
    this->current_refine_index += adaptSteps;
  }

  // determine, if any refined information now is contained in matrix b_adapt
  this->b_is_refined = this->b_adapt_matrix_.getNcols() > dima;
  return;
#endif /* USE_GSL */
}

/**
 * todo: Kilian! this function can probably be refactored into some parent class,
 * because it calculates similar stuff as corresponding part in in
 * DBMatOfflineChol::choleskyModification.
 * But note, that this class is an online class, whereas the class containing the
 * choleskyModification function is an offline class.
 * If you want to refactor, also look at the adjusted part of this function by
 * searching for "//### begin adjusted part"
 *
 * greets,
 * Dima
 */
void DBMatOnlineDEOrthoAdapt::compute_L2_gridvectors(Grid& grid, size_t newPoints,
                                                     double newLambda) {
  if (newPoints > 0) {
    size_t gridSize = grid.getStorage().getSize();
    size_t gridDim = grid.getStorage().getDimension();

    // allocate the new points
    for (size_t i = 0; i < newPoints; i++) {
      sgpp::base::DataVector vec(gridSize);
      this->refined_points_.push_back(vec);
    }

    DataMatrix level(gridSize, gridDim);
    DataMatrix index(gridSize, gridDim);

    grid.getStorage().getLevelIndexArraysForEval(level, index);
    double lambda_conf = newLambda;

    // loop to calculate all L2-products of added points based on the
    // hat-function as basis function
    for (size_t i = 0; i < gridSize; i++) {
      for (size_t j = gridSize - newPoints; j < gridSize; j++) {
        double res = 1;
        for (size_t k = 0; k < gridDim; k++) {
          double lik = level.get(i, k);
          double ljk = level.get(j, k);
          double iik = index.get(i, k);
          double ijk = index.get(j, k);

          if (lik == ljk) {
            if (iik == ijk) {
              // use formula for identical ansatz functions:
              res *= 2 / lik / 3;
            } else {
              // different index, but same level => ansatz functions do not overlap
              res = 0.;
              break;
            }
          } else {
            if (std::max((iik - 1) / lik, (ijk - 1) / ljk) >=
                std::min((iik + 1) / lik, (ijk + 1) / ljk)) {
              // ansatz functions do not not overlap:
              res = 0.;
              break;
            } else {
              // use formula for different overlapping ansatz functions
              if (lik > ljk) {  // phi_i_k is the "smaller" ansatz function
                double diff = (iik / lik) - (ijk / ljk);  // x_i_k - x_j_k
                double temp_res = fabs(diff - (1 / lik)) + fabs(diff + (1 / lik)) - fabs(diff);
                temp_res *= ljk;
                temp_res = (1 - temp_res) / lik;
                res *= temp_res;
              } else {  // phi_j_k is the "smaller" ansatz function
                double diff = (ijk / ljk) - (iik / lik);  // x_j_k - x_i_k
                double temp_res = fabs(diff - (1 / ljk)) + fabs(diff + (1 / ljk)) - fabs(diff);
                temp_res *= lik;
                temp_res = (1 - temp_res) / ljk;
                res *= temp_res;
              }
            }
          }
        }

        //### begin adjusted part
        // add current lambda to lower diagonal elements of the refine points
        if (i == j) {
          this->refined_points_[j - gridSize + newPoints + this->current_refine_index].set(
              i, res + lambda_conf);
        } else {
          this->refined_points_[j - gridSize + newPoints + this->current_refine_index].set(i, res);
        }
      }
    }

    // fill in the new L2 products in the older points, and resize them
    // loop all the old points of this->refined_points_
    for (size_t i = 0; i < this->current_refine_index; i++) {
      this->refined_points_[i].resize(gridSize);
      for (size_t j = this->current_refine_index; j < this->current_refine_index + newPoints; j++) {
        double xxx = this->refined_points_[j].get(i);
        this->refined_points_[i].set(j, xxx);
      }
    }
    //### end adjusted part
  }
}

}  // namespace datadriven
}  // namespace sgpp
