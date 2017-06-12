// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceBsplineClenshawCurtis.hpp>

#include <sgpp/globaldef.hpp>

#include <string.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace sgpp {
namespace pde {

OperationLaplaceBsplineClenshawCurtis::OperationLaplaceBsplineClenshawCurtis(sgpp::base::Grid* grid)
    : grid(grid), clenshawCurtisTable(base::ClenshawCurtisTable::getInstance()) {}

OperationLaplaceBsplineClenshawCurtis::~OperationLaplaceBsplineClenshawCurtis() {}

void OperationLaplaceBsplineClenshawCurtis::mult(sgpp::base::DataVector& alpha,
                                                 sgpp::base::DataVector& result) {
  size_t gridSize = grid->getSize();
  size_t gridDim = grid->getDimension();
  const size_t p = dynamic_cast<sgpp::base::BsplineClenshawCurtisGrid*>(grid)->getDegree();
  const size_t pp1h = (p + 1) / 2;
  const size_t quadOrder = p + 1;
  sgpp::base::SBsplineClenshawCurtisBase& basis =
      dynamic_cast<sgpp::base::SBsplineClenshawCurtisBase&>(grid->getBasis());
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
        const sgpp::base::level_t lik = storage[i].getLevel(k);
        const sgpp::base::level_t ljk = storage[j].getLevel(k);
        const sgpp::base::index_t iik = storage[i].getIndex(k);
        const sgpp::base::index_t ijk = storage[j].getIndex(k);
        const int left_iik = static_cast<int>(iik) - static_cast<int>(pp1h);
        const int left_ijk = static_cast<int>(ijk) - static_cast<int>(pp1h);
        // clenshawCurtisPoint returns 0.0 if point is right of 1.0
        const double right_iik_point =
            basis.clenshawCurtisPoint(lik, iik + static_cast<base::index_t>(pp1h));
        const double right_ijk_point =
            basis.clenshawCurtisPoint(ljk, ijk + static_cast<base::index_t>(pp1h));
        // points are not uniformly distributed thus we need to find the left and right boundarys
        const double left_i = ((left_iik > 0) ? clenshawCurtisTable.getPoint(lik, left_iik) : 0.0);
        const double right_i =
            (right_iik_point == 0.0 || (right_iik_point >= 1.0)) ? 1.0 : right_iik_point;
        const double left_j = ((left_ijk > 0) ? clenshawCurtisTable.getPoint(ljk, left_ijk) : 0.0);
        const double right_j =
            (right_ijk_point == 0.0 || (right_ijk_point >= 1.0)) ? 1.0 : right_ijk_point;

        if (left_j > right_i && left_i > right_j) {
          // Ansatz functions do not not overlap:
          integrals1D[k] = 0.0;
          integralsDeriv1D[k] = 0.0;
          break;
        } else {
          // Use formula for different overlapping ansatz functions:
          double scaling;
          size_t start;
          size_t stop;

          // find the finer one of the two levels and calculate the first and last intervall
          const base::level_t finest_l = std::max(lik, ljk);
          // start and stop are the *absolute* index values of the interval we want to sum up
          if (lik >= ljk) {
            start = ((iik < pp1h) ? 0 : (iik - pp1h));
            stop = std::min(iik + pp1h - 1, static_cast<size_t>((1 << lik) - 1));
          } else {
            start = ((ijk < pp1h) ? 0 : (ijk - pp1h));
            stop = std::min(ijk + pp1h - 1, static_cast<size_t>((1 << ljk) - 1));
          }

          double temp_ij = 0.0;
          double temp_ij_deriv = 0.0;

          for (size_t n = start; n <= stop; n++) {
            double left = std::max(
                clenshawCurtisTable.getPoint(finest_l, static_cast<base::index_t>(n)), 0.0);
            double right = std::min(
                clenshawCurtisTable.getPoint(finest_l, static_cast<base::index_t>(n + 1)), 1.0);
            scaling = right - left;
            for (size_t c = 0; c < quadOrder; c++) {
              const double x = left + scaling * coordinates[c];
              temp_ij += scaling * weights[c] * basis.eval(lik, iik, x) * basis.eval(ljk, ijk, x);
              temp_ij_deriv +=
                  scaling * weights[c] * basis.evalDx(lik, iik, x) * basis.evalDx(ljk, ijk, x);
            }
          }
          integrals1D[k] = temp_ij;
          integralsDeriv1D[k] = temp_ij_deriv;
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
      double res = 0.0;

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
        res += temp_res;
      }
      // multiplication and summation of results
      result[i] += res * alpha[j];
      if (i != j) result[j] += res * alpha[i];
    }
  }
}

}  // namespace pde
}  // namespace sgpp
