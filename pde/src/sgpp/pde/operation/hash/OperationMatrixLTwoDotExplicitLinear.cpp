// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitLinear.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

#include <string.h>
#include <cmath>
#include <vector>
#include <algorithm>

namespace SGPP {
namespace pde {

OperationMatrixLTwoDotExplicitLinear::OperationMatrixLTwoDotExplicitLinear(
    SGPP::base::DataMatrix* m, SGPP::base::Grid* grid)
    : ownsMatrix_(false) {
  m_ = m;
  buildMatrix(grid);
}

OperationMatrixLTwoDotExplicitLinear::OperationMatrixLTwoDotExplicitLinear(SGPP::base::Grid* grid)
    : ownsMatrix_(true) {
  m_ = new SGPP::base::DataMatrix(grid->getStorage().size(), grid->getStorage().size());
  buildMatrix(grid);
}

void OperationMatrixLTwoDotExplicitLinear::buildMatrix(SGPP::base::Grid* grid) {
  size_t gridSize = grid->getStorage().size();
  size_t gridDim = grid->getStorage().dim();

  SGPP::base::DataMatrix level(gridSize, gridDim);
  SGPP::base::DataMatrix index(gridSize, gridDim);

  grid->getStorage().getLevelIndexArraysForEval(level, index);

  for (size_t i = 0; i < gridSize; i++) {
    for (size_t j = i; j < gridSize; j++) {
      float_t res = 1;

      for (size_t k = 0; k < gridDim; k++) {
        float_t lik = level.get(i, k);
        float_t ljk = level.get(j, k);
        float_t iik = index.get(i, k);
        float_t ijk = index.get(j, k);

        if (lik == ljk) {
          if (iik == ijk) {
            // Use formula for identical ansatz functions:
            res *= 2 / lik / 3;
          } else {
            // Different index, but same level => ansatz functions do not overlap:
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
            if (lik > ljk) {                             // Phi_i_k is the "smaller" ansatz function
              float_t diff = (iik / lik) - (ijk / ljk);  // x_i_k - x_j_k
              float_t temp_res = fabs(diff - (1 / lik)) + fabs(diff + (1 / lik)) - fabs(diff);
              temp_res *= ljk;
              temp_res = (1 - temp_res) / lik;
              res *= temp_res;
            } else {                                     // Phi_j_k is the "smaller" ansatz function
              float_t diff = (ijk / ljk) - (iik / lik);  // x_j_k - x_i_k
              float_t temp_res = fabs(diff - (1 / ljk)) + fabs(diff + (1 / ljk)) - fabs(diff);
              temp_res *= lik;
              temp_res = (1 - temp_res) / ljk;
              res *= temp_res;
            }
          }
        }
      }

      m_->set(i, j, res);
      m_->set(j, i, res);
    }
  }
}

OperationMatrixLTwoDotExplicitLinear::~OperationMatrixLTwoDotExplicitLinear() {
  if (ownsMatrix_) delete m_;
}

void OperationMatrixLTwoDotExplicitLinear::mult(SGPP::base::DataVector& alpha,
                                                SGPP::base::DataVector& result) {
  size_t nrows = m_->getNrows();
  size_t ncols = m_->getNcols();

  if (alpha.getSize() != ncols || result.getSize() != nrows) {
    throw SGPP::base::data_exception("Dimensions do not match!");
  }

  float_t* data = m_->getPointer();

  // Standard matrix multiplication:
  float_t temp = 0.;
  size_t acc = 0;

  for (size_t i = 0; i < nrows; i++) {
    for (size_t j = 0; j < ncols; j++) {
      temp += data[j + acc] * alpha[j];
    }

    result[i] = temp;
    temp = 0.;
    acc += ncols;
  }
}

}  // namespace pde
}  // namespace SGPP
