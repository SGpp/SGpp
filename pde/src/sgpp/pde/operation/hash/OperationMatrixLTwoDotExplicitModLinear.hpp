// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OperationMatrixLTwoDotExplicitModLinear_HPP_
#define OperationMatrixLTwoDotExplicitModLinear_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Explicit representation of the matrix \f$(\Phi_i,\Phi_j)_{L2}\f$ for a sparse grid
 */
class OperationMatrixLTwoDotExplicitModLinear : public sgpp::base::OperationMatrix {
 public:
  /**
   * Constructor that only builds the object, without initialization
   */
  OperationMatrixLTwoDotExplicitModLinear();

  /**
   * Constructor that uses a external matrix pointer to construct the matrix,
   * i.e. matrix is NOT destroyed by the destructor of
   * OperationMatrixLTwoDotExplicitModLinearFullGrid
   *
   * @param m pointer to datamatrix of size (number of grid point) x (number of grid points)
   * @param grid the sparse grid
   */
  OperationMatrixLTwoDotExplicitModLinear(sgpp::base::DataMatrix* m, sgpp::base::Grid* grid);
  /**
   * Constructor that creates an own matrix
   * i.e. matrix is destroyed by the destructor of OperationMatrixLTwoDotExplicitModLinearFullGrid
   *
   * @param grid the sparse grid
   */
  explicit OperationMatrixLTwoDotExplicitModLinear(sgpp::base::Grid* grid);

  /**
   * Destructor
   */
  virtual ~OperationMatrixLTwoDotExplicitModLinear();

  /**
   * Implementation of standard matrix multiplication
   *
   * @param alpha DataVector that is multiplied to the matrix
   * @param result DataVector into which the result of multiplication is stored
   */
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  /**
   * generalization of "buildMatrix" function, creates L2-dot-product matrix for specified bounds
   * @param mat matrix for storage of L2 producs
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
    base::GridStorage& storage = grid->getStorage();
    base::SLinearModifiedBase& basis = const_cast<base::SLinearModifiedBase&>(
        dynamic_cast<const base::SLinearModifiedBase&>(grid->getBasis()));

    // needed for non-quadratic matrix cases
    bool mat_quadratic = (i_start + i_end + j_start + j_end == 0);

    // init standard values
    i_end = i_end == 0 ? gridSize : i_end;

    for (size_t i = i_start; i < i_end; i++) {
      j_start = j_start == 0 ? i : j_start;
      j_end = j_end == 0 ? gridSize : j_end;
#pragma omp parallel for schedule(guided)
      for (size_t j = j_start; j < j_end; j++) {
        double res = 1;

        for (size_t k = 0; k < gridDim; k++) {
          const base::level_t lik = storage[i].getLevel(k);
          const base::level_t ljk = storage[j].getLevel(k);
          const base::index_t iik = storage[i].getIndex(k);
          const base::index_t ijk = storage[j].getIndex(k);
          base::index_t hInvi = (1 << lik);
          base::index_t hInvj = (1 << ljk);
          double hInviDbl = static_cast<double>(hInvi);
          double hInvjDbl = static_cast<double>(hInvj);
          double temp_res;

          if (lik == ljk) {
            if (lik == 1) {
              continue;
            } else if (iik == ijk) {
              if (iik == 1 || iik == hInvi - 1) {
                // Use formula for identical modified ansatz functions:
                temp_res = 8 / (hInviDbl * 3);
              } else {
                // Use formula for identical ansatz functions:
                temp_res = 2 / (hInviDbl * 3);
              }
            } else {
              // Different index, but same level => ansatz functions do not overlap:
              res = 0.;
              break;
            }
          } else {
            // if one of the basis functions is from level 1 it's easy
            if (lik == 1) {
              temp_res = basis.getIntegral(ljk, ijk);
            } else if (ljk == 1) {
              temp_res = basis.getIntegral(lik, iik);
            } else if ((iik - 1) / hInviDbl >= (ijk + 1) / hInvjDbl ||
                       (iik + 1) / hInviDbl <= (ijk - 1) / hInvjDbl) {
              // Ansatz functions do not not overlap:
              res = 0.;
              break;
            } else {
              // use formula for different overlapping ansatz functions:
              if (lik > ljk) {  // Phi_i_k is the "smaller" ansatz function
                if ((iik == 1 && ijk == 1) || (iik == hInvi - 1 && ijk == hInvj - 1)) {
                  // integrate modified basis prdouct from 0 to 2^(-lik + 1)
                  temp_res = 4 * ((1 / hInviDbl) - (hInvjDbl / 3 / (hInviDbl * hInviDbl)));
                } else if (ijk == 1) {
                  // integrate product of modified Phi_i_k with
                  // regular Phi_j_k from (ijk-1)/2^(ljk) to  (ijk+1)/2^(ljk)
                  temp_res = (1 / hInviDbl) * (1 / hInviDbl) * (2 * hInviDbl - iik * hInvjDbl);
                } else if (ijk == hInvj - 1) {
                  // symmetric to ijk == 1
                  temp_res =
                      (1 / hInviDbl) * (1 / hInviDbl) * (2 * hInviDbl - (hInvi - iik) * hInvjDbl);
                } else {
                  double diff = (iik / hInviDbl) - (ijk / hInvjDbl);  // x_i_k - x_j_k
                  temp_res = fabs(diff - (1 / hInviDbl)) + fabs(diff + (1 / hInviDbl)) - fabs(diff);
                  temp_res *= hInvjDbl;
                  temp_res = (1 - temp_res) / hInviDbl;
                }
              } else {  // Phi_j_k is the "smaller" ansatz function
                // symmetric to case above
                if ((iik == 1 && ijk == 1) || (iik == hInvi - 1 && ijk == hInvj - 1)) {
                  // both basis functions are modified
                  // integrate modified basis prdouct from 0 to 2^(-ljk + 1)
                  temp_res = 4 * ((1 / hInvjDbl) - (hInviDbl / (3 * hInvjDbl * hInvjDbl)));
                } else if (iik == 1) {
                  // integrate product of modified Phi_i_k with
                  // regular Phi_j_k from (ijk-1)/2^(ljk) to  (ijk+1)/2^(ljk)
                  temp_res = (1 / hInvjDbl) * (1 / hInvjDbl) * (2 * hInvjDbl - ijk * hInviDbl);
                } else if (iik == hInvi - 1) {
                  // symmetric to iik == 1
                  temp_res =
                      (1 / hInvjDbl) * (1 / hInvjDbl) * (2 * hInvjDbl - (hInvj - ijk) * hInviDbl);
                } else {
                  double diff = (ijk / hInvjDbl) - (iik / hInviDbl);  // x_j_k - x_i_k
                  temp_res = fabs(diff - (1 / hInvjDbl)) + fabs(diff + (1 / hInvjDbl)) - fabs(diff);
                  temp_res *= hInviDbl;
                  temp_res = (1 - temp_res) / hInvjDbl;
                }
              }
            }
          }
          res *= temp_res;
        }
        if (mat_quadratic) {
          mat->set(i, j, res);
          mat->set(j, i, res);
        } else {
          mat->set(i, j - j_start, res);
        }
      }
    }
  }

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
#endif /* OperationMatrixLTwoDotExplicitModLinear_HPP_ */
