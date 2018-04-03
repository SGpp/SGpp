// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLaplaceModPolyClenshawCurtis.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/type/ModPolyClenshawCurtisGrid.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedClenshawCurtisBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <string.h>
#include <cmath>
#include <vector>
#include <algorithm>

namespace sgpp {
namespace pde {

OperationLaplaceModPolyClenshawCurtis::OperationLaplaceModPolyClenshawCurtis(sgpp::base::Grid* grid)
  : grid(grid),
    clenshawCurtisTable(base::ClenshawCurtisTable::getInstance()) {}

OperationLaplaceModPolyClenshawCurtis::~OperationLaplaceModPolyClenshawCurtis() {}

void OperationLaplaceModPolyClenshawCurtis::mult(sgpp::base::DataVector& alpha,
                                sgpp::base::DataVector& result) {
  size_t gridSize = grid->getSize();
  size_t gridDim = grid->getDimension();
  const size_t p = dynamic_cast<sgpp::base::ModPolyClenshawCurtisGrid*>(grid)->getDegree();
  const size_t quadOrder = p + 1;
  sgpp::base::SPolyModifiedClenshawCurtisBase& basis
    = dynamic_cast<sgpp::base::SPolyModifiedClenshawCurtisBase&>(grid->getBasis());
  sgpp::base::GridStorage& storage = grid->getStorage();

  sgpp::base::DataVector integrals1D(gridDim);
  sgpp::base::DataVector integralsDeriv1D(gridDim);

  sgpp::base::DataVector coordinates;
  sgpp::base::DataVector weights;
  sgpp::base::GaussLegendreQuadRule1D gauss;
  gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);

  size_t nrows = storage.getSize();
  size_t ncols = storage.getSize();

  if (alpha.getSize() != ncols || result.getSize() != nrows) {
    throw sgpp::base::data_exception("Dimensions do not match!");
  }

  for (size_t i = 0; i < gridSize; i++) {
    result[i] = 0;
  }

  for (size_t i = 0; i < gridSize; i++) {
    for (size_t j = i; j < gridSize; j++) {
      for (size_t k = 0; k < gridDim; k++) {
        const base::level_t lik = storage[i].getLevel(k);
        const base::level_t ljk = storage[j].getLevel(k);
        const base::index_t iik = storage[i].getIndex(k);
        const base::index_t ijk = storage[j].getIndex(k);
        const double left_i = clenshawCurtisTable.getPoint(lik, iik - 1);
        const double right_i = clenshawCurtisTable.getPoint(lik, iik + 1);
        const double left_j = clenshawCurtisTable.getPoint(ljk, ijk - 1);
        const double right_j = clenshawCurtisTable.getPoint(ljk, ijk + 1);

        if (left_j >= right_i || left_i >= right_j) {
          // Ansatz functions do not not overlap:
          integrals1D[k] = 0.0;
          integralsDeriv1D[k] = 0.0;
          break;
        } else {
          // Use formula for different overlapping ansatz functions:

          double temp_res = 0.0;
          double temp_res_deriv = 0.0;

          const double left = std::max(left_i, left_j);
          const double right = std::min(right_i, right_j);
          double scaling = right - left;

          for (size_t c = 0; c < quadOrder; c++) {
            const double x = left + scaling * (coordinates[c]);
            temp_res += weights[c] * basis.eval(lik, iik, x) * basis.eval(ljk, ijk, x);
            temp_res_deriv += weights[c] * basis.evalDx(lik, iik, x) *
              basis.evalDx(ljk, ijk, x);
          }

          integrals1D[k] = scaling * temp_res;
          integralsDeriv1D[k] = scaling * temp_res_deriv;
        }
      }

      /**
       * int nabla phi_i(x) * nabla phi_j(x) dx
       * = sum_k int (dx_k phi_i(x)) * (dx_k phi_j(x)) dx
       * = sum_k int (phi'_{i_k}(x_k) * phi'_{i_k}(x_k) *
       *              prod_{l!=k} phi_{i_l}(x_l) * phi_{j_l}(x_l)) dx
       * = sum_k int (phi'_{i_k}(x_k) * phi'_{j_k}(x_k)) dx_k *
       *         prod_{l!=k} int (phi_{i_l}(x_l) * phi_{j_l}(x_l)) dx_l
       */
      double temp_ij = 0.0;

      // sum due to scalar product (of gradients)
      for (size_t k = 0; k < gridDim; k++) {
        double temp_res = 1.0;

        // integral of product of partial derivatives w.r.t. dimension k
        // = product of all 1D integrals except dimension k, times 1D integral of derivatives
        for (size_t l = 0; l < gridDim; l++) {
          if (l == k) {
            temp_res *= integralsDeriv1D[l];
          } else {
            temp_res *= integrals1D[l];
          }

          if (temp_res == 0.0) {
            break;
          }
        }

        temp_ij += temp_res;
      }

      // multiplication and summation of results
      result[i] += temp_ij * alpha[j];
      if (i != j)
        result[j] += temp_ij * alpha[i];
    }
  }
}

}  // namespace pde
}  // namespace sgpp
