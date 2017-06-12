// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotModPoly.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>

namespace sgpp {
namespace pde {

OperationMatrixLTwoDotModPoly::OperationMatrixLTwoDotModPoly(
    sgpp::base::Grid* grid) {
  this->grid = grid;
}

OperationMatrixLTwoDotModPoly::~OperationMatrixLTwoDotModPoly() {}

void OperationMatrixLTwoDotModPoly::mult(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result) {
  const size_t p = dynamic_cast<sgpp::base::ModPolyGrid*>(grid)->getDegree();
  const size_t quadOrder = p + 1;
  // clenshawCurtisPoint Method does not leave base const so we cast const'ness away
  base::SBasis& basis = const_cast<base::SBasis&>(grid->getBasis());
  base::GridStorage& storage = grid->getStorage();
  base::DataVector coordinates;
  base::DataVector weights;
  base::GaussLegendreQuadRule1D gauss;
  gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);

  size_t nrows = storage.getSize();
  size_t ncols = storage.getSize();

  if (alpha.getSize() != ncols || result.getSize() != nrows) {
    throw sgpp::base::data_exception("Dimensions do not match!");
  }

  size_t gridSize = storage.getSize();
  size_t gridDim = storage.getDimension();

  for (size_t i = 0; i < gridSize; i++) {
    result[i] = 0;
  }

  for (size_t i = 0; i < gridSize; i++) {
    for (size_t j = i; j < gridSize; j++) {
      double temp_ij = 1;

      for (size_t k = 0; k < gridDim; k++) {
        const base::level_t lik = storage[i].getLevel(k);
        const base::level_t ljk = storage[j].getLevel(k);
        const base::index_t iik = storage[i].getIndex(k);
        const base::index_t ijk = storage[j].getIndex(k);
        const double left_i = 1.0/(1 << lik) * (iik - 1);
        const double left_j = 1.0/(1 << ljk) * (ijk - 1);
        const double right_i = 1.0/(1 << lik) * (iik + 1);
        const double right_j = 1.0/(1 << ljk) * (ijk + 1);

        if (left_j >= right_i || left_i >= right_j) {
          // Ansatz functions do not not overlap:
          temp_ij = 0.0;
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
          temp_ij *= scaling*temp_res;
        }
      }
      result[i] += temp_ij * alpha[j];
      if (i != j)
        result[j] += temp_ij * alpha[i];
    }
  }
}
}  // namespace pde
}  // namespace sgpp
