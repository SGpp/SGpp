/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOnlineDEOrthoAdapt.cpp
 */

// #ifdef USE_GSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>

#include <algorithm>
#include <vector>

namespace sgpp {
namespace datadriven {

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

void DBMatOnlineDEOrthoAdapt::adapt(sgpp::base::DataMatrix& refinePts,
                                    std::vector<size_t>& coarsePts) {
  if (refinePts.getNcols() > 0) {
    this->sherman_morrison_adapt(refinePts, coarsePts, true);
  }
  if (coarsePts.size() > 0) {
    this->sherman_morrison_adapt(refinePts, coarsePts, false);
  }
}
void DBMatOnlineDEOrthoAdapt::solveSLE(DataVector& b, bool do_cv) {}

void DBMatOnlineDEOrthoAdapt::sherman_morrison_adapt(sgpp::base::DataMatrix& refinePts,
                                                     std::vector<size_t>& coarsePts, bool refine) {
  // check, if offline object has been decomposed
  sgpp::datadriven::DBMatOfflineOrthoAdapt* childPtr =
      static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt*>(&this->offlineObject);
  if (childPtr->getDimA() == 0) {
    throw sgpp::base::algorithm_exception(
        "offline object wasn't decomposed yet ...\ncan't apply sherman_morrison this way");
  }

  // determine the final size of the matrix system
  size_t oldSize = this->b_is_refined ? this->b_adapt_matrix_.getNcols() : childPtr->getDimA();
  size_t newSize = refine ? (oldSize + refinePts.getNcols()) : (oldSize - coarsePts.size());
  size_t adaptSteps = newSize > oldSize ? (newSize - oldSize) : (oldSize - newSize);

  // todo: If instead these points should just be ignored, change this to kill some coarsePts
  // If coarsening points in base grid or out of grid
  for (size_t i = 0; i < coarsePts.size(); i++) {
    if (coarsePts[i] < childPtr->getDimA() || coarsePts[i] > this->b_adapt_matrix_.getNcols()) {
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

  // apply sherman morrison formula k times ( one for each point refining/coarsening )
  for (size_t k = 0; k < adaptSteps; k++) {
    // fetch point with according size
    // refine -> size gets bigger by 1
    // coarse -> size stays the same
    size_t current_size = oldSize + (refine ? k + 1 : 0);
    sgpp::base::DataVector x(current_size);
    if (refine) {
      // configure x according to sherman morrison formula
      refinePts.getColumn(k, x);
      double config_x_value =
          (x.get(current_size - 1) - 1) * 0.5;  // lambda is added already on all of refinePts
      x.set(current_size - 1, config_x_value);
    } else {
      // determine former value of vector x, given only its index in the grid
      // todo: retrieve x vector out of index
    }

    // configure unit vector depending on refine/coarsen
    // e[unit_index] = unit_value, others are zero
    size_t unit_index = refine ? current_size - 1 : coarsePts[k];

    // view of T^{-1} of the offline object
    gsl_matrix_view t_inv_view = gsl_matrix_view_array(childPtr->getTinv().getPointer(),
                                                       childPtr->getDimA(), childPtr->getDimA());
    // view of Q of the offline object
    gsl_matrix_view q_view = gsl_matrix_view_array(childPtr->getQ().getPointer(),
                                                   childPtr->getDimA(), childPtr->getDimA());
    // view of B of the online object, which holds all information of refinement/coarsening
    gsl_matrix_view b_adapt_view =
        gsl_matrix_view_array(this->b_adapt_matrix_.getPointer(), current_size, current_size);
    // view of current point to refine/coarsen
    gsl_matrix_view x_view = gsl_matrix_view_array(x.getPointer(), current_size, 1);

    // view of current point to refine/coarsen cut to offline objects matrix size
    gsl_matrix_view x_cut_view = gsl_matrix_view_array(x.getPointer(), childPtr->getDimA(), 1);

    // allocating space for buffer
    gsl_matrix* buffer = gsl_matrix_alloc(childPtr->getDimA(), 1);

    // allocating space for the term: Q*T^{-1}*Q^t * x_cut
    gsl_matrix* x_term = gsl_matrix_alloc(current_size, 1);

    // allocating space for the term: B * x
    gsl_matrix* bx_term = gsl_matrix_alloc(current_size, 1);

    /*
     * starting calculations: fist sherman morrison update/downdate
     * note: paying attention to sign of unitvector !!!
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
    double eQTQx = (refine ? 0.0 : gsl_matrix_get(x_term, unit_index - 1, 1));
    double eBx = (refine ? 1.0 : -1.0) * gsl_matrix_get(bx_term, unit_index - 1, 1);
    double divisor = (refine ? 1.0 : -1.0) * (1 - eQTQx - eBx);

    // calculating: bx_term + x_term, where x_term is filled with zero to fit dimension of bx
    // the result is stored in bx_term, as x_term is later needed
    for (size_t i = 0; i < childPtr->getDimA(); i++) {
      bx_term->data[i] += x_term->data[i];
    }

    // subtracting the matrices B - ( x_term * e_term ) / divisor
    // where the vector product is an outer product yielding a matrix
    for (size_t i = 0; i < current_size - 1; i++) {
      for (size_t j = 0; j < current_size - 1; j++) {
        double final_value =
            this->b_adapt_matrix_.get(i, j) -
            (bx_term->data[i] * this->b_adapt_matrix_.get(unit_index - 1, j)) / divisor;
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
    eBx = (refine ? 1.0 : -1.0) * gsl_matrix_get(bx_term, unit_index - 1, 1);
    divisor = (refine ? 1.0 : -1.0) * (1 - eQTQx - eBx);

    // calculating: bx_tilde_term + x_term
    for (size_t i = 0; i < childPtr->getDimA(); i++) {
      bx_term->data[i] += x_term->data[i];
    }

    // subtracting the matrices B_tilde - ( x_term * b_tilde_term ) / divisor
    // where the vector product is an outer product yielding a matrix
    // note: as the order is flipped, i and j have to be changed AND b column now
    for (size_t i = 0; i < current_size - 1; i++) {
      for (size_t j = 0; j < current_size - 1; j++) {
        double final_value =
            this->b_adapt_matrix_.get(i, j) -
            (bx_term->data[j] * this->b_adapt_matrix_.get(i, unit_index - 1)) / divisor;
        this->b_adapt_matrix_.set(i, j, final_value);
      }
    }
  } /*** ending sherman morrison update/downdate calculations ***/

  // If points were coarsened the b_adapt_matrix will now have empty rows and columns
  // on the indices of the coarsened points. In the following algorithm, the symmetry
  // of b_adapt_matrix is used to fit the blocks together again.
  if (!refine) {
    std::sort(coarsePts.begin(), coarsePts.end());
    size_t vj = 0;  // number of indices already considered in rows
    size_t vi = 0;  // number of indices already considered in columns
    for (size_t i = 0; i < newSize - coarsePts.size(); i++) {
      for (size_t j = i; j < newSize - coarsePts.size(); j++) {
        double value = this->b_adapt_matrix_.get(i + vi, j + vj);
        this->b_adapt_matrix_.set(i, j, value);
        this->b_adapt_matrix_.set(j, i, value);
        if (j == coarsePts[vj] - vj) {  // then skip this column
          vj++;
        }
      }
      if (i == coarsePts[vi] - vi) {  // then skip this row
        vi++;
      }
      vj = 0;  // new row, new count
    }
    this->b_adapt_matrix_.resizeQuadratic(newSize);
  }
}
}  // namespace datadriven
}  // namespace sgpp
