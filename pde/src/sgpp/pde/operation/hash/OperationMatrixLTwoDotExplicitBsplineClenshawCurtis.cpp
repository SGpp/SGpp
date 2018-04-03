// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitBsplineClenshawCurtis.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

#include <string.h>
#include <cmath>
#include <vector>
#include <algorithm>

namespace sgpp {
namespace pde {

OperationMatrixLTwoDotExplicitBsplineClenshawCurtis::
  OperationMatrixLTwoDotExplicitBsplineClenshawCurtis(
    sgpp::base::DataMatrix* m, sgpp::base::Grid* grid)
    : ownsMatrix_(false) {
  m_ = m;
  buildMatrix(grid);
}

OperationMatrixLTwoDotExplicitBsplineClenshawCurtis::
  OperationMatrixLTwoDotExplicitBsplineClenshawCurtis(
    sgpp::base::Grid* grid)
    : ownsMatrix_(true) {
  m_ = new sgpp::base::DataMatrix(grid->getSize(), grid->getSize());
  buildMatrix(grid);
}

void OperationMatrixLTwoDotExplicitBsplineClenshawCurtis::buildMatrix(sgpp::base::Grid* grid) {
  size_t gridSize = grid->getSize();
  size_t gridDim = grid->getDimension();
  const size_t p = dynamic_cast<sgpp::base::BsplineClenshawCurtisGrid*>(grid)->getDegree();
  const size_t pp1h = (p + 1) >> 1;  // (p + 1) / 2
  // const double pp1hDbl = static_cast<double>(pp1h);
  const size_t quadOrder = p + 1;
  // clenshawCurtisPoint Method modifies base so we need to cast constness away
  base::SBsplineClenshawCurtisBase& basis =
    const_cast<base::SBsplineClenshawCurtisBase&> (
      dynamic_cast<const base::SBsplineClenshawCurtisBase&>(grid->getBasis()));
  base::GridStorage& storage = grid->getStorage();
  base::DataVector coordinates;
  base::DataVector weights;
  base::GaussLegendreQuadRule1D gauss;
  gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);
  for (size_t i = 0; i < gridSize; i++) {
    for (size_t j = i; j < gridSize; j++) {
      double res = 1.0;

      for (size_t k = 0; k < gridDim; k++) {
        const base::level_t lik = storage[i].getLevel(k);
        const base::level_t ljk = storage[j].getLevel(k);
        const base::index_t iik = storage[i].getIndex(k);
        const base::index_t ijk = storage[j].getIndex(k);
        const int left_iik = static_cast<int>(iik) - static_cast<int>(pp1h);
        const int left_ijk = static_cast<int>(ijk) - static_cast<int>(pp1h);
        // clenshawCurtisPoint returns 0.0 if point is right of 1.0
        const double right_iik_point =
          basis.clenshawCurtisPoint(lik, iik + static_cast<base::index_t>(pp1h));
        const double right_ijk_point =
          basis.clenshawCurtisPoint(ljk, ijk + static_cast<base::index_t>(pp1h));
        // points are not uniformly distributed thus we need to find the left and right boundarys
        const double left_i = ((left_iik > 0)? basis.clenshawCurtisPoint(lik, left_iik) : 0.0);
        const double right_i = (right_iik_point == 0.0 || (right_iik_point >= 1.0))
                               ? 1.0 : right_iik_point;
        const double left_j = ((left_ijk > 0)? basis.clenshawCurtisPoint(ljk, left_ijk) : 0.0);
        const double right_j = (right_ijk_point == 0.0 || (right_ijk_point >= 1.0))
                               ? 1.0 : right_ijk_point;

        // Check if ansatz functions overlap. We need to use the actual position of the
        // boundaries because the index values iik and ijk might be for different levels.
        if (left_j > right_i && left_i > right_j) {
          // Ansatz functions do not not overlap:
          res = 0.0;
          break;
        } else {
          size_t start;
          size_t stop;
          double scaling;
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
          // std::cout << "start: " << start << std::endl;
          // std::cout << "stop: " << stop << std::endl;
          double temp_res = 0.0;
          for (size_t n = start; n <= stop; n++) {
            double left = std::max(basis.clenshawCurtisPoint(
                                                      finest_l,
                                                      static_cast<base::index_t>(n)), 0.0);
            double right = std::min(basis.clenshawCurtisPoint(
                                                       finest_l,
                                                       static_cast<base::index_t>(n+1)), 1.0);
            scaling = right - left;
            for (size_t c = 0; c < quadOrder; c++) {
              const double x = left + scaling * coordinates[c];
              temp_res += scaling *
                         (weights[c] * basis.eval(lik, iik, x) * basis.eval(ljk, ijk, x));
            }
          }
          res *= temp_res;
        }
      }
      // std::cout << "res:" << res << std::endl;
      m_->set(i, j, res);
      m_->set(j, i, res);
    }
  }
}

OperationMatrixLTwoDotExplicitBsplineClenshawCurtis::
  ~OperationMatrixLTwoDotExplicitBsplineClenshawCurtis() {
  if (ownsMatrix_) delete m_;
}

void OperationMatrixLTwoDotExplicitBsplineClenshawCurtis::mult(sgpp::base::DataVector& alpha,
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
