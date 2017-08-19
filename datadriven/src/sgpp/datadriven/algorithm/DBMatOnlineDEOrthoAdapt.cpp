// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <algorithm>
#include <iomanip>
#include <vector>

namespace sgpp {
namespace datadriven {

// print datamatrices for debugging
static void printMatrix(sgpp::base::DataMatrix a) {
  for (size_t i = 0; i < a.getNrows(); i++) {
    for (size_t j = 0; j < a.getNcols(); j++) {
      std::cout << std::setprecision(5) << std::fixed << a.get(i, j) << "  ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

static void printMatrixView(gsl_matrix_view a) {
  for (size_t i = 0; i < (&a.matrix)->size1; i++) {
    for (size_t j = 0; j < (&a.matrix)->size2; j++) {
      std::cout << std::setprecision(5) << std::fixed
                << (&a.matrix)->data[i * (&a.matrix)->size1 + j] << "  ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

DBMatOnlineDEOrthoAdapt::DBMatOnlineDEOrthoAdapt(DBMatOffline& offline, double beta)
    : sgpp::datadriven::DBMatOnlineDE(offline, beta) {
  if (offline.getConfig().decomp_type_ != sgpp::datadriven::DBMatDecompostionType::OrthoAdapt) {
    throw sgpp::base::algorithm_exception(
        "offline object has wrong decomposition type: OrthoAdapt needed!");
  } else {
    // dummy space for information of refinement/coarsening
    this->b_adapt_matrix_ = sgpp::base::DataMatrix(1, 1);
    this->b_is_refined = false;
  }
}

// not tested yet
void DBMatOnlineDEOrthoAdapt::getAdaptPointsCartesian(bool refine,
                                                      sgpp::base::DataMatrix& adapt_matrix,
                                                      std::vector<size_t> coarsenIndices) {
  // at this point, the old grid size equals to the dimension of the b_adapt_matrix_
  sgpp::datadriven::DBMatOfflineOrthoAdapt* childPtr =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&this->offlineObject);
  size_t oldGridSize = this->b_is_refined ? this->b_adapt_matrix_.getNcols() : childPtr->getDimA();
  size_t newGridSize = this->getOfflineObject().getGrid().getStorage().getSize();
  size_t gridDim;

  if (refine) {
    // when refining, the dimension of the point is equal to the size of new grid
    gridDim = this->getOfflineObject().getGrid().getStorage().getDimension();
  } else {
    // when coarsening, the dimension of the points is equal to the prior grid
    // which is equal to b_adapt_matrix
    gridDim = this->b_adapt_matrix_.getNcols();
  }

  size_t numberAdaptPoints = refine ? newGridSize - oldGridSize : oldGridSize - newGridSize;
  if (numberAdaptPoints <= 0) {
    std::cout << "ain't nuthin' to adapt, yo\n";
    return;
  }
  sgpp::base::DataMatrix level(newGridSize, gridDim);
  sgpp::base::DataMatrix index(newGridSize, gridDim);

  this->getOfflineObject().getGrid().getStorage().getLevelIndexArraysForEval(level, index);
  double lambda_conf = this->lambda;
  // Loop to calculate all L2-products of added points based on the
  // hat-function as basis function
  for (size_t i = 0; i < newGridSize; i++) {
    // Loops all new added points
    for (size_t j = oldGridSize; j < newGridSize; j++) {
      double res = 1;
      for (size_t k = 0; k < gridDim; k++) {
        double lik = level.get(i, k);
        double ljk = level.get(j, k);
        double iik = index.get(i, k);
        double ijk = index.get(j, k);

        if (lik == ljk) {
          if (iik == ijk) {
            // Use formula for identical ansatz functions:
            res *= 2 / lik / 3;
          } else {
            // Different index, but same level => ansatz functions do not
            // overlap:
            res = 0.;
            break;
          }
        } else {
          if (std::max((iik - 1) / lik, (ijk - 1) / ljk) >=
              std::min((iik + 1) / lik, (ijk + 1) / ljk)) {
            // Ansatz functions do not not overlap:
            res = 0.;
            break;
          } else {
            // Use formula for different overlapping ansatz functions:
            if (lik > ljk) {                            // Phi_i_k is the "smaller" ansatz function
              double diff = (iik / lik) - (ijk / ljk);  // x_i_k - x_j_k
              double temp_res = fabs(diff - (1 / lik)) + fabs(diff + (1 / lik)) - fabs(diff);
              temp_res *= ljk;
              temp_res = (1 - temp_res) / lik;
              res *= temp_res;
            } else {                                    // Phi_j_k is the "smaller" ansatz function
              double diff = (ijk / ljk) - (iik / lik);  // x_j_k - x_i_k
              double temp_res = fabs(diff - (1 / ljk)) + fabs(diff + (1 / ljk)) - fabs(diff);
              temp_res *= lik;
              temp_res = (1 - temp_res) / ljk;
              res *= temp_res;
            }
          }
        }
      }
      // The new Rows/Cols are stored in mat_refine

      // add current lambda to lower diagonal elements of mat_refine
      if (i == j) {
        adapt_matrix.set(i, j - newGridSize + numberAdaptPoints, res + lambda_conf);
      } else {
        adapt_matrix.set(i, j - newGridSize + numberAdaptPoints, res);
      }
    }
  }
}

void DBMatOnlineDEOrthoAdapt::adapt(sgpp::base::DataMatrix& adapt_points, bool refine,
                                    std::vector<size_t> coarsenIndices) {
  // todo: add parameter "list of coarsen points"
  // and gather the points to refine/coarsen in this very function
  if (refine) {
    this->sherman_morrison_adapt(adapt_points, refine);
  } else {
    if (adapt_points.getNcols() != coarsenIndices.size()) {
      throw sgpp::base::algorithm_exception("coarsen point indices and amount don't match");
    }
    this->sherman_morrison_adapt(adapt_points, false, coarsenIndices);
  }
}

void DBMatOnlineDEOrthoAdapt::solveSLE(DataVector& b, bool do_cv) {
  sgpp::datadriven::DBMatOfflineOrthoAdapt* childPtr =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&this->offlineObject);
  // create solver
  sgpp::datadriven::DBMatDMSOrthoAdapt* solver = new sgpp::datadriven::DBMatDMSOrthoAdapt();
  // solve the created system
  this->alpha = sgpp::base::DataVector(b.getSize());
  solver->solve(childPtr->getTinv(), childPtr->getQ(), this->getB(), b, this->alpha);
}

void DBMatOnlineDEOrthoAdapt::sherman_morrison_adapt(sgpp::base::DataMatrix& adapt_points,
                                                     bool refine,
                                                     std::vector<size_t> coarsenIndices) {
  // check, if offline object has been decomposed
  sgpp::datadriven::DBMatOfflineOrthoAdapt* childPtr =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&this->offlineObject);
  if (childPtr->getDimA() == 0) {
    throw sgpp::base::algorithm_exception(
        "offline object wasn't decomposed yet ...\ncan't apply sherman_morrison this way");
  }

  // determine the final size of the matrix system
  size_t oldSize = this->b_is_refined ? this->b_adapt_matrix_.getNcols() : childPtr->getDimA();
  size_t newSize = refine ? (oldSize + adapt_points.getNcols()) : (oldSize - coarsenIndices.size());
  size_t adaptSteps = (newSize > oldSize) ? (newSize - oldSize) : (oldSize - newSize);

  // todo: If instead these points should just be ignored, change this to kill some coarsePts
  // If coarsening points in base grid or out of grid
  for (size_t i = 0; i < coarsenIndices.size(); i++) {
    if (coarsenIndices[i] < childPtr->getDimA() ||
        coarsenIndices[i] > this->b_adapt_matrix_.getNcols()) {
      // todo: split condition and kill points here with coarsePts.delete
      throw sgpp::base::algorithm_exception(
          "OrthoAdapt can't coarsen points, which were not refined before");
    }
  }

  // allocate space for b_adapt_matrix_ and fill new diagonal entries with ones
  // note: only done in refining! Coarsening will resize at the end of function
  if (refine) {
    this->b_adapt_matrix_.resizeQuadratic(newSize);
    for (size_t i = oldSize; i < newSize; i++) {
      this->b_adapt_matrix_.set(i, i, 1.0);
    }
  }
  // create view of the whole matrix b_adapt
  gsl_matrix_view b_adapt_full_view =
      gsl_matrix_view_array(this->b_adapt_matrix_.getPointer(), this->b_adapt_matrix_.getNrows(),
                            this->b_adapt_matrix_.getNcols());

  // apply sherman morrison formula k times
  // ( one for each point
  // refining/coarsening )
  for (size_t k = 0; k < adaptSteps; k++) {
    // fetch point with according size
    // refine -> size gets bigger by 1
    // coarse -> size stays the same
    size_t current_size = oldSize + (refine ? k + 1 : 0);
    sgpp::base::DataVector x((refine ? newSize : oldSize));

    // configure x according to sherman morrison formula
    if (refine) {
      adapt_points.getColumn(k, x);
      x.resize(current_size);  // clipping off not needed entries
      double config_x_value =
          (x.get(current_size - 1) - 1) * 0.5;  // lambda is added already on all of refinePts
      x.set(current_size - 1, config_x_value);
    } else {
      adapt_points.getColumn(k, x);
      // no clipping off when coarsening
      double config_x_value =
          (x.get(coarsenIndices[k]) - 1) * 0.5;  // lambda is added already on all of refinePts
      x.set(coarsenIndices[k], config_x_value);
    }

    // configure unit vector depending on
    // refine/coarsen
    // e[unit_index] = unit_value, others are zero
    size_t unit_index = refine ? current_size - 1 : coarsenIndices[k];

    // view of T^{-1} of the offline object
    gsl_matrix_view t_inv_view = gsl_matrix_view_array(childPtr->getTinv().getPointer(),
                                                       childPtr->getDimA(), childPtr->getDimA());
    // view of Q of the offline object
    gsl_matrix_view q_view = gsl_matrix_view_array(childPtr->getQ().getPointer(),
                                                   childPtr->getDimA(), childPtr->getDimA());
    // view of B of the online object, which holds all information of refinement/coarsening
    gsl_matrix_view b_adapt_view =
        gsl_matrix_submatrix(&b_adapt_full_view.matrix, 0, 0, current_size, current_size);

    // view of current point to refine/coarsen
    gsl_matrix_view x_view = gsl_matrix_view_array(x.getPointer(), current_size, 1);

    // view of current point to refine/coarsen cut to offline objects matrix size
    gsl_matrix_view x_cut_view = gsl_matrix_view_array(x.getPointer(), childPtr->getDimA(), 1);

    // allocating space for buffer
    gsl_matrix* buffer = gsl_matrix_alloc(childPtr->getDimA(), 1);

    // allocating space for the term: Q*T^{-1}*Q^t * x_cut
    gsl_matrix* x_term = gsl_matrix_alloc(childPtr->getDimA(), 1);

    // allocating space for the term: B * x
    gsl_matrix* bx_term = gsl_matrix_alloc(current_size, 1);

    /*
     * starting calculations: fist sherman morrison update/downdate
     */
    // calculating bx_term = B * x
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &b_adapt_view.matrix, &x_view.matrix, 0.0,
                   bx_term);

    // calculating x_term = Q*T^{-1}*Q^t * x_cut
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &q_view.matrix, &x_cut_view.matrix, 0.0,
                   x_term);  // x_term = Q^t * x_cut

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &t_inv_view.matrix, x_term, 0.0,
                   buffer);  // buffer = T^{-1} * Q^t * x_cut

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &q_view.matrix, buffer, 0.0,
                   x_term);  // x_term = Q*T^{-1}*Q^t * x_cut

    // calculating the divisor of sherman morrison formula: 1 - e^t * x_term - e^t * bx_term
    double eQTQx = 0.0;  // is always zero, because initial gridpoints cannot be refined
    double eBx = (refine ? 1.0 : -1.0) * gsl_matrix_get(bx_term, unit_index, 0);
    double divisor = 1 + eQTQx + eBx;

    // calculating: bx_term + x_term, where x_term is filled with zero to fit dimension of bx
    // the result is stored in bx_term, as x_term is later needed
    for (size_t i = 0; i < childPtr->getDimA(); i++) {
      bx_term->data[i] += x_term->data[i];
    }

    // subtracting the matrices B - ( x_term * e_term ) / divisor
    // where the vector product is an outer product yielding a matrix
    for (size_t i = 0; i < current_size; i++) {
      for (size_t j = 0; j < current_size; j++) {
        double final_value =
            this->b_adapt_matrix_.get(i, j) -
            (bx_term->data[i] * (refine ? 1.0 : -1.0) * this->b_adapt_matrix_.get(unit_index, j)) /
                divisor;
        this->b_adapt_matrix_.set(i, j, final_value);
      }
    }
    /*
     * starting calculations: second sherman morrison update/downdate
     * reusing value of first step
     */
    // calculating new bx_term = x^t * B_tilde = B_tilde^t * x
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &b_adapt_view.matrix, &x_view.matrix, 0.0,
                   bx_term);

    // calculating the divisor of sherman morrison formula: 1 - e^t * x_term - e^t * b_tilde_x_term
    // note: the eQTQx term is equal to xQTQe, so it is ok taking the old value
    eBx = (refine ? 1.0 : -1.0) * gsl_matrix_get(bx_term, unit_index, 0);
    divisor = 1 + eQTQx + eBx;

    // calculating: bx_tilde_term + x_term
    for (size_t i = 0; i < childPtr->getDimA(); i++) {
      bx_term->data[i] += x_term->data[i];
    }

    // subtracting the matrices B_tilde - ( x_term * b_tilde_term ) / divisor
    // where the vector product is an outer product yielding a matrix
    // note: as the order is flipped, i and j have to be changed AND b column now
    for (size_t i = 0; i < current_size; i++) {
      for (size_t j = 0; j < current_size; j++) {
        double final_value =
            this->b_adapt_matrix_.get(i, j) -
            (bx_term->data[j] * (refine ? 1.0 : -1.0) * this->b_adapt_matrix_.get(i, unit_index)) /
                divisor;
        this->b_adapt_matrix_.set(i, j, final_value);
      }
    }
  } /*** ending sherman morrison update/downdate calculations ***/

  // If points were coarsened the b_adapt_matrix will now have empty rows and columns
  // on the indices of the coarsened points. In the following algorithm, the symmetry
  // of b_adapt_matrix is used to fit the blocks together again.
  if (!refine) {
    std::sort(coarsenIndices.begin(), coarsenIndices.end());
    size_t vj = 0;  // number of indices already considered in rows
    size_t vi = 0;  // number of indices already considered in columns
    for (size_t i = 0; i < oldSize - coarsenIndices.size(); i++) {
      if (i == coarsenIndices[vi] - vi) {  // then skip this row
        vi++;
      }
      for (size_t j = i; j < oldSize - coarsenIndices.size(); j++) {
        if (j == coarsenIndices[vj] - vj) {  // then skip this column
          vj++;
        }
        double value = this->b_adapt_matrix_.get(i + vi, j + vj);
        this->b_adapt_matrix_.set(i, j, value);
        this->b_adapt_matrix_.set(j, i, value);
      }

      vj = 0;  // new row, new count
    }
    this->b_adapt_matrix_.resizeQuadratic(newSize);
  }

  // determine, if any refined information is in matrix b_adapt
  this->b_is_refined = this->b_adapt_matrix_.getNcols() <= childPtr->getDimA() ? false : true;
}
}  // namespace datadriven
}  // namespace sgpp
