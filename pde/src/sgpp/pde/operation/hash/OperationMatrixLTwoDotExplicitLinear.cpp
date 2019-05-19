// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitLinear.hpp>

#include <sgpp/globaldef.hpp>

#include <string.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace sgpp {
namespace pde {

OperationMatrixLTwoDotExplicitLinear::OperationMatrixLTwoDotExplicitLinear(
    sgpp::base::DataMatrix* m, sgpp::base::Grid* grid)
    : ownsMatrix_(false) {
  m_ = m;
  buildMatrix(grid);
}

OperationMatrixLTwoDotExplicitLinear::OperationMatrixLTwoDotExplicitLinear(sgpp::base::Grid* grid)
    : ownsMatrix_(true) {
  m_ = new sgpp::base::DataMatrix(grid->getSize(), grid->getSize());
  buildMatrix(grid);
}

inline void OperationMatrixLTwoDotExplicitLinear::buildMatrixWithBounds(
    sgpp::base::DataMatrix* mat, sgpp::base::Grid* grid, size_t i_start, size_t i_end,
    size_t j_start, size_t j_end) {
  // todo: dima, bounds einbauen
  size_t gridSize = grid->getSize();
  size_t gridDim = grid->getDimension();

  sgpp::base::DataMatrix level(gridSize, gridDim);
  sgpp::base::DataMatrix index(gridSize, gridDim);

  grid->getStorage().getLevelIndexArraysForEval(level, index);

  // init standard values
  i_start = i_start == 0 ? 0 : i_start;
  i_end = i_end == 0 ? gridSize : i_end;

  for (size_t i = i_start; i < i_end; i++) {
    j_start = j_start == 0 ? i : j_start;
    j_end = j_end == 0 ? gridSize : j_end;
#pragma omp parallel for schedule(guided)
    for (size_t j = j_start; j < j_end; j++) {
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

      mat->set(i, j, res);
      mat->set(j, i, res);
    }
  }
}

void OperationMatrixLTwoDotExplicitLinear::buildMatrix(sgpp::base::Grid* grid) {
  this->buildMatrixWithBounds(this->m_, grid);
}

OperationMatrixLTwoDotExplicitLinear::~OperationMatrixLTwoDotExplicitLinear() {
  if (ownsMatrix_) delete m_;
}

void OperationMatrixLTwoDotExplicitLinear::mult(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataVector& result) {
  size_t nrows = m_->getNrows();
  size_t ncols = m_->getNcols();

  if (alpha.getSize() != ncols || result.getSize() != nrows) {
    throw sgpp::base::data_exception("Dimensions do not match!");
  }

  double* data = m_->getPointer();

  // Standard matrix multiplication:
  double temp = 0.;
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
}  // namespace sgpp
