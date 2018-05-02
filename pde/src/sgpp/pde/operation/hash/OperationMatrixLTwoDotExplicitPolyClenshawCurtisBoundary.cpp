// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary.hpp>
#include <sgpp/base/grid/type/PolyClenshawCurtisBoundaryGrid.hpp>

#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <string.h>
#include <cmath>
#include <vector>
#include <algorithm>

namespace sgpp {
namespace pde {

OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary::
    OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary(
    sgpp::base::DataMatrix* m, sgpp::base::Grid* grid)
    : ownsMatrix_(false),
    clenshawCurtisTable(base::ClenshawCurtisTable::getInstance()) {
  m_ = m;
  buildMatrix(grid);
}

OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary::
    OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary(sgpp::base::Grid* grid)
    : ownsMatrix_(true),
    clenshawCurtisTable(base::ClenshawCurtisTable::getInstance()) {
  m_ = new sgpp::base::DataMatrix(grid->getSize(), grid->getSize());
  buildMatrix(grid);
}

void OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary::buildMatrix(sgpp::base::Grid* grid) {
  size_t gridSize = grid->getSize();
  size_t gridDim = grid->getDimension();
  const size_t p = dynamic_cast<sgpp::base::PolyClenshawCurtisBoundaryGrid*>(grid)->getDegree();
  const size_t quadOrder = p + 1;
  base::SBasis& basis = const_cast<base::SBasis&>(grid->getBasis());
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
        // points are not uniformly distributed thus we need to find the left and right boundarys
        const double left_i = (iik == 0) ? 0 : clenshawCurtisTable.getPoint(lik, iik - 1);
        const double left_j = (ijk == 0) ? 0 : clenshawCurtisTable.getPoint(ljk, ijk - 1);
        const double right_i = (iik == static_cast<base::index_t>(1 << lik))
          ? 1.0
          : clenshawCurtisTable.getPoint(lik, iik + 1);
        const double right_j = (ijk == static_cast<base::index_t>(1 << ljk))
          ? 1.0
          : clenshawCurtisTable.getPoint(ljk, ijk + 1);

        // Check if ansatz functions overlap. We need to use the actual position of the
        // boundaries because the index values iik and ijk might be for different levels.
        if (left_j >= right_i || left_i >= right_j) {
          // Ansatz functions do not not overlap:
          res = 0.0;
          break;
        } else {
          const double left = std::max(left_i, left_j);
          const double right = std::min(right_i, right_j);
          const double scaling = right - left;
          double temp_res = 0.0;
          for (size_t c = 0; c < quadOrder; c++) {
            const double x = left + scaling * coordinates[c];
            temp_res += scaling * (weights[c] * basis.eval(lik, iik, x) * basis.eval(ljk, ijk, x));
          }
          res *= temp_res;
        }
      }
      m_->set(i, j, res);
      m_->set(j, i, res);
    }
  }
}

OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary::
    ~OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary() {
  if (ownsMatrix_) delete m_;
}

void OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary::mult(
     sgpp::base::DataVector& alpha,
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
