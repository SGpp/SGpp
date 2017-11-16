// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitModifiedLinear.hpp>

#include <sgpp/globaldef.hpp>

#include <string.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace sgpp {
namespace pde {

OperationMatrixLTwoDotExplicitModifiedLinear::OperationMatrixLTwoDotExplicitModifiedLinear(
    sgpp::base::DataMatrix* m, sgpp::base::Grid* grid)
    : ownsMatrix_(false) {
  m_ = m;
  buildMatrix(grid);
}

OperationMatrixLTwoDotExplicitModifiedLinear::OperationMatrixLTwoDotExplicitModifiedLinear
    (sgpp::base::Grid* grid)
    : ownsMatrix_(true) {
  m_ = new sgpp::base::DataMatrix(grid->getSize(), grid->getSize());
  buildMatrix(grid);
}

void OperationMatrixLTwoDotExplicitModifiedLinear::buildMatrix(sgpp::base::Grid* grid) {
  size_t gridSize = grid->getSize();
  size_t gridDim = grid->getDimension();

  sgpp::base::DataMatrix level(gridSize, gridDim);
  sgpp::base::DataMatrix index(gridSize, gridDim);

  grid->getStorage().getLevelIndexArraysForEval(level, index);

  for (size_t i = 0; i < gridSize; i++) {
// #pragma omp parallel for schedule(guided)
    for (size_t j = i; j < gridSize; j++) {
      double res = 1;

      for (size_t k = 0; k < gridDim; k++) {
        double lik = level.get(i, k);
        double ljk = level.get(j, k);
        double iik = index.get(i, k);
        double ijk = index.get(j, k);
        if (lik == 2.) {
          if (ljk == 2.) {
            // do nothing (multiply with 1)
          } else if (ljk -1 == static_cast<int>(ijk) || ijk == 1.) {
            // right or left modified basis
            res *= 2./ ljk;
          } else {
            // regular ansatz function
            res *= 1./ljk;
          }
        } else if (ljk == 2.) {
          if (lik -1 == iik || iik == 1.) {
            res *= 2./ lik;
          } else {
            res *= 1./lik;
          }
        } else if (lik == ljk) {
          if (iik == ijk) {
            if (lik -1 == iik || iik == 1.) {
              // Identical modified basis
              res *= 16.0/(3.0* lik);
            } else {
              // Use formula for identical ansatz functions:
              res *= 2 / lik / 3;
            }
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
              if (ljk -1 == static_cast<int>(ijk) || ijk == 1.) {  // "larger" function is modified
                double ldiff = lik/ljk;
                // "smaller" function is modified too
                if (lik -1 == static_cast<int>(iik) || iik == 1.) {
                  res *= 2/3. * (2+ 2*(2-1./ldiff)) /lik;
                } else {  // smaller function is a inner function
                  auto phi = [] (double x) -> double{return 2-2*x;};
                  double tmpi = (iik > lik/2.) ? lik - iik : iik;
                  double xi = tmpi*0.5/ldiff;
                  double delta = 0.25/ldiff;
                  res *= 4./6. * (phi(xi -delta) + phi(xi) + phi(xi+delta)) /lik;
                }
              } else {  // two inner ansatz functions
                double diff = (iik / lik) - (ijk / ljk);  // x_i_k - x_j_k
                double temp_res = fabs(diff - (1 / lik)) + fabs(diff + (1 / lik)) - fabs(diff);
                temp_res *= ljk;
                temp_res = (1 - temp_res) / lik;
                res *= temp_res;
              }
            } else {                                    // Phi_j_k is the "smaller" ansatz function
              if (lik -1 == static_cast<int>(iik) || iik == 1.) {  // "larger" function is modified
                double ldiff = ljk/lik;
                // "smaller" function is modified too
                if (ljk -1 == static_cast<int>(ijk) || ijk == 1.) {
                  res *= 2/3. * (2+ 2*(2-1./ldiff)) /ljk;
                } else {  // smaller function is a inner function
                  auto phi = [] (double x) -> double{return 2-2*x;};
                  double tmpi = (ijk > ljk/2.)? ljk - ijk: ijk;
                  double xi = tmpi*0.5/ldiff;
                  double delta = 0.25/ldiff;
                  res *= 4./6. * (phi(xi -delta) + phi(xi) + phi(xi+delta)) /ljk;
                }
              } else {
                double diff = (ijk / ljk) - (iik / lik);  // x_j_k - x_i_k
                double temp_res = fabs(diff - (1 / ljk)) + fabs(diff + (1 / ljk)) - fabs(diff);
                temp_res *= lik;
                temp_res = (1 - temp_res) / ljk;
                res *= temp_res;
              }
            }
          }
        }
      }
      m_->set(i, j, res);
      m_->set(j, i, res);
    }
  }
}

OperationMatrixLTwoDotExplicitModifiedLinear::~OperationMatrixLTwoDotExplicitModifiedLinear() {
  if (ownsMatrix_) delete m_;
}

void OperationMatrixLTwoDotExplicitModifiedLinear::mult(sgpp::base::DataVector& alpha,
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
