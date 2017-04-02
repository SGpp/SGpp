// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLaplaceExplicitBspline.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <string.h>
#include <cmath>
#include <vector>
#include <algorithm>

namespace sgpp {
namespace pde {

OperationLaplaceExplicitBspline::OperationLaplaceExplicitBspline(
    sgpp::base::DataMatrix* m, sgpp::base::Grid* grid)
    : ownsMatrix_(false) {
  m_ = m;
  buildMatrix(grid);
}

OperationLaplaceExplicitBspline::OperationLaplaceExplicitBspline(sgpp::base::Grid* grid)
    : ownsMatrix_(true) {
  m_ = new sgpp::base::DataMatrix(grid->getSize(), grid->getSize());
  buildMatrix(grid);
}

void OperationLaplaceExplicitBspline::buildMatrix(sgpp::base::Grid* grid) {
  size_t gridSize = grid->getSize();
  size_t gridDim = grid->getDimension();
  const size_t p = dynamic_cast<sgpp::base::BsplineGrid*>(grid)->getDegree();
  const size_t pp1h = (p + 1) / 2;
  const double pp1hDbl = static_cast<double>(pp1h);
  const size_t quadOrder = p + 1;
  sgpp::base::SBsplineBase& basis = dynamic_cast<sgpp::base::SBsplineBase&>(grid->getBasis());
  sgpp::base::GridStorage& storage = grid->getStorage();

  sgpp::base::DataVector integrals1D(gridDim);
  sgpp::base::DataVector integralsDeriv1D(gridDim);

  sgpp::base::DataVector coordinates;
  sgpp::base::DataVector weights;
  sgpp::base::GaussLegendreQuadRule1D gauss;
  gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);

  for (size_t i = 0; i < gridSize; i++) {
    for (size_t j = i; j < gridSize; j++) {
      for (size_t k = 0; k < gridDim; k++) {
        const sgpp::base::level_t lik = storage[i].getLevel(k);
        const sgpp::base::level_t ljk = storage[j].getLevel(k);
        const sgpp::base::index_t iik = storage[i].getIndex(k);
        const sgpp::base::index_t ijk = storage[j].getIndex(k);
        const sgpp::base::index_t hInvik = 1 << lik;
        const sgpp::base::index_t hInvjk = 1 << ljk;
        const double hik = 1.0 / static_cast<double>(hInvik);
        const double hjk = 1.0 / static_cast<double>(hInvjk);

        if (std::max((static_cast<double>(iik) - pp1hDbl) * hik,
                     (static_cast<double>(ijk) - pp1hDbl) * hjk) >=
            std::min((static_cast<double>(iik) + pp1hDbl) * hik,
                     (static_cast<double>(ijk) + pp1hDbl) * hjk)) {
          // Ansatz functions do not not overlap:
          integrals1D[k] = 0.0;
          integralsDeriv1D[k] = 0.0;
          break;
        } else {
          // Use formula for different overlapping ansatz functions:
          double offset;
          double scaling;
          size_t start;
          size_t stop;

          if (lik >= ljk) {
            offset = (static_cast<double>(iik) - pp1hDbl) * hik;
            scaling = hik;
            start = ((iik > pp1h) ? 0 : (pp1h - iik));
            stop = std::min(p, hInvik + pp1h - iik - 1);
          } else {
            offset = (static_cast<double>(ijk) - pp1hDbl) * hjk;
            scaling = hjk;
            start = ((ijk > pp1h) ? 0 : (pp1h - ijk));
            stop = std::min(p, hInvjk + pp1h - ijk - 1);
          }

          double temp_res = 0.0;
          double temp_res_deriv = 0.0;

          for (size_t n = start; n <= stop; n++) {
            for (size_t c = 0; c < quadOrder; c++) {
              const double x = offset + scaling * (coordinates[c] + static_cast<double>(n));
              temp_res += weights[c] * basis.eval(lik, iik, x) * basis.eval(ljk, ijk, x);
              temp_res_deriv += weights[c] * basis.evalDx(lik, iik, x) *
                  basis.evalDx(ljk, ijk, x);
            }
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
      m_->set(i, j, res);
      m_->set(j, i, res);
    }
  }
}

OperationLaplaceExplicitBspline::~OperationLaplaceExplicitBspline() {
  if (ownsMatrix_) delete m_;
}

void OperationLaplaceExplicitBspline::mult(sgpp::base::DataVector& alpha,
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
