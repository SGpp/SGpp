// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitPolyBoundary.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

#include <string.h>
#include <cmath>
#include <vector>
#include <algorithm>

namespace sgpp {
namespace pde {

OperationMatrixLTwoDotExplicitPolyBoundary::OperationMatrixLTwoDotExplicitPolyBoundary(
    sgpp::base::DataMatrix* m, sgpp::base::Grid* grid)
    : ownsMatrix_(false) {
  m_ = m;
  buildMatrix(grid);
}

OperationMatrixLTwoDotExplicitPolyBoundary::
    OperationMatrixLTwoDotExplicitPolyBoundary(sgpp::base::Grid* grid)
    : ownsMatrix_(true) {
  m_ = new sgpp::base::DataMatrix(grid->getSize(), grid->getSize());
  buildMatrix(grid);
}

void OperationMatrixLTwoDotExplicitPolyBoundary::buildMatrix(sgpp::base::Grid* grid) {
  size_t gridSize = grid->getSize();
  size_t gridDim = grid->getDimension();
  const size_t p = dynamic_cast<sgpp::base::PolyBoundaryGrid*>(grid)->getDegree();
  // const double pp1hDbl = static_cast<double>(pp1h);
  const size_t quadOrder = p + 1;
  // clenshawCurtisPoint Method does not leave base const so we cast const'ness away
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
        // double left_i;
        // double right_i;
        // double left_j;
        // double right_j;
        // // correct boundary cases
        // // left i
        // if (iik == 0)
        //   left_i = 0;
        // else
        //   left_i = 1.0/(1 << lik) * (iik - 1);
        // // left j
        // if (ijk == 0)
        //   left_j = 0;
        // else
        //   left_j = 1.0/(1 << ljk) * (ijk - 1);
        // // right i
        // if (iik == static_cast<base::index_t>(1 << lik))
        //   right_i = 1.0/(1 << lik) * iik;
        // else
        //   right_i = 1.0/(1 << lik) * (iik + 1);
        // // right j
        // if (ijk == static_cast<base::index_t>(1 << ljk))
        //   right_j = 1.0/(1 << ljk) * ijk;
        // else
        //   right_j = 1.0/(1 << ljk) * (ijk + 1);

        const double left_i = (iik == 0) ? 0 : 1.0/(1 << lik) * (iik - 1);
        const double right_i =
          (iik ==  static_cast<base::index_t>(1 << lik))
          ? 1.0/(1 << lik) * iik
          : 1.0/(1 << lik) * (iik + 1);
        const double left_j = (ijk == 0) ? 0 : 1.0/(1 << ljk) * (ijk - 1);
        const double right_j =
          (ijk ==  static_cast<base::index_t>(1 << ljk))
          ? 1.0/(1 << ljk) * ijk
          : 1.0/(1 << ljk) * (ijk + 1);


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
            temp_res += weights[c] * basis.eval(lik, iik, x) * basis.eval(ljk, ijk, x);
          }
          res *= scaling*temp_res;
        }
      }
      // std::cout << "res:" << res << std::endl;
      m_->set(i, j, res);
      m_->set(j, i, res);
    }
  }
}

OperationMatrixLTwoDotExplicitPolyBoundary::~OperationMatrixLTwoDotExplicitPolyBoundary() {
  if (ownsMatrix_) delete m_;
}

void OperationMatrixLTwoDotExplicitPolyBoundary::mult(sgpp::base::DataVector& alpha,
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
