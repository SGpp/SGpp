// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OperationMatrixLTwoDotExplicitModifiedLinear_HPP_
#define OperationMatrixLTwoDotExplicitModifiedLinear_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Explicit representation of the matrix \f$(\Phi_i,\Phi_j)_{L2}\f$ for a sparse grid
 */
class OperationMatrixLTwoDotExplicitModifiedLinear : public sgpp::base::OperationMatrix {
 public:
  /**
   * Constructor that only builds the object, without initialization
   */
  OperationMatrixLTwoDotExplicitModifiedLinear();
  /**
   * Constructor that uses a external matrix pointer to construct the matrix,
   * i.e. matrix is NOT destroyed by the destructor of
   * OperationMatrixLTwoDotExplicitModifiedLinearFullGrid
   *
   * @param m pointer to datamatrix of size (number of grid point) x (number of grid points)
   * @param grid the sparse grid
   */
  OperationMatrixLTwoDotExplicitModifiedLinear(sgpp::base::DataMatrix* m, sgpp::base::Grid* grid);
  /**
   * Constructor that creates an own matrix
   * i.e. matrix is destroyed by the destructor of
   * OperationMatrixLTwoDotExplicitModifiedLinearFullGrid
   *
   * @param grid the sparse grid
   */
  explicit OperationMatrixLTwoDotExplicitModifiedLinear(sgpp::base::Grid* grid);

  /**
   * Destructor
   */
  virtual ~OperationMatrixLTwoDotExplicitModifiedLinear();

  /**
   * Implementation of standard matrix multiplication
   *
   * @param alpha DataVector that is multiplied to the matrix
   * @param result DataVector into which the result of multiplication is stored
   */
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  /**
   * generalization of "buildMatrix" function, creates L2-dot-product matrix for specified bounds
   * @param mat matrix for storage of L2 products
   * @param grid the underlying grid
   * @param i_start start index for row iteration
   * @param i_end end index for row iteration
   * @param j_start start index for column iteration
   * @param j_end end index for column iteration
   */
  inline void buildMatrixWithBounds(sgpp::base::DataMatrix* mat, sgpp::base::Grid* grid,
                                    size_t i_start = 0, size_t i_end = 0, size_t j_start = 0,
                                    size_t j_end = 0) {
    size_t gridSize = grid->getSize();
    size_t gridDim = grid->getDimension();

    sgpp::base::DataMatrix level(gridSize, gridDim);
    sgpp::base::DataMatrix index(gridSize, gridDim);

    grid->getStorage().getLevelIndexArraysForEval(level, index);

    // needed for non-quadratic matrix cases
    bool mat_quadratic = (i_start + i_end + j_start + j_end == 0);

    // init standard values
    i_end = i_end == 0 ? gridSize : i_end;

    for (size_t i = 0; i < gridSize; i++) {
      j_start = j_start == 0 ? i : j_start;
      j_end = j_end == 0 ? gridSize : j_end;
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
            } else if (ljk - 1 == static_cast<int>(ijk) || ijk == 1.) {
              // right or left modified basis
              res *= 2. / ljk;
            } else {
              // regular ansatz function
              res *= 1. / ljk;
            }
          } else if (ljk == 2.) {
            if (lik - 1 == iik || iik == 1.) {
              res *= 2. / lik;
            } else {
              res *= 1. / lik;
            }
          } else if (lik == ljk) {
            if (iik == ijk) {
              if (lik - 1 == iik || iik == 1.) {
                // Identical modified basis
                res *= 16.0 / (3.0 * lik);
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
              if (lik > ljk) {  // Phi_i_k is the "smaller" ansatz function
                if (ljk - 1 == static_cast<int>(ijk) ||
                    ijk == 1.) {  // "larger" function is modified
                  double ldiff = lik / ljk;
                  // "smaller" function is modified too
                  if (lik - 1 == static_cast<int>(iik) || iik == 1.) {
                    res *= 2 / 3. * (2 + 2 * (2 - 1. / ldiff)) / lik;
                  } else {  // smaller function is a inner function
                    auto phi = [](double x) -> double { return 2 - 2 * x; };
                    double tmpi = (iik > lik / 2.) ? lik - iik : iik;
                    double xi = tmpi * 0.5 / ldiff;
                    double delta = 0.25 / ldiff;
                    res *= 4. / 6. * (phi(xi - delta) + phi(xi) + phi(xi + delta)) / lik;
                  }
                } else {                                    // two inner ansatz functions
                  double diff = (iik / lik) - (ijk / ljk);  // x_i_k - x_j_k
                  double temp_res = fabs(diff - (1 / lik)) + fabs(diff + (1 / lik)) - fabs(diff);
                  temp_res *= ljk;
                  temp_res = (1 - temp_res) / lik;
                  res *= temp_res;
                }
              } else {  // Phi_j_k is the "smaller" ansatz function
                if (lik - 1 == static_cast<int>(iik) ||
                    iik == 1.) {  // "larger" function is modified
                  double ldiff = ljk / lik;
                  // "smaller" function is modified too
                  if (ljk - 1 == static_cast<int>(ijk) || ijk == 1.) {
                    res *= 2 / 3. * (2 + 2 * (2 - 1. / ldiff)) / ljk;
                  } else {  // smaller function is a inner function
                    auto phi = [](double x) -> double { return 2 - 2 * x; };
                    double tmpi = (ijk > ljk / 2.) ? ljk - ijk : ijk;
                    double xi = tmpi * 0.5 / ldiff;
                    double delta = 0.25 / ldiff;
                    res *= 4. / 6. * (phi(xi - delta) + phi(xi) + phi(xi + delta)) / ljk;
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
        if (mat_quadratic) {
          mat->set(i, j, res);
          mat->set(j, i, res);
        } else {
          mat->set(i, j - j_start, res);
        }
      }
    }
  };

 private:
  /**
   * This method is used by both constructors to build the matrix
   */
  void buildMatrix(sgpp::base::Grid* grid);

  sgpp::base::DataMatrix* m_;
  bool ownsMatrix_;
};

}  // namespace pde
}  // namespace sgpp
#endif /* OperationMatrixLTwoDotExplicitModifiedLinear_HPP_ */
