// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotNakBsplineBoundaryCombigrid.hpp>

#include <algorithm>

namespace sgpp {
namespace pde {

OperationMatrixLTwoDotNakBsplineBoundaryCombigrid::
    OperationMatrixLTwoDotNakBsplineBoundaryCombigrid(sgpp::base::Grid* grid) {
  this->grid = grid;
}

OperationMatrixLTwoDotNakBsplineBoundaryCombigrid::
    ~OperationMatrixLTwoDotNakBsplineBoundaryCombigrid() {}

void OperationMatrixLTwoDotNakBsplineBoundaryCombigrid::mult(sgpp::base::DataVector& alpha,
                                                             sgpp::base::DataVector& result) {
  const size_t p = dynamic_cast<sgpp::base::NakBsplineBoundaryCombigridGrid*>(grid)->getDegree();
  if ((p != 1) && (p != 3) && (p != 5)) {
    std::cerr << "OperationMatrixLTwoDotNakBsplineBoundary: only B spline degrees 1, 3 and 5 are "
                 "supported."
              << std::endl;
  }

  const size_t pp1h = (p + 1) >> 1;  //  =|_p/2_|
  const double pp1hDbl = static_cast<double>(pp1h);
  const size_t quadOrder = p + 1;
  //  base::SBasis& basis = const_cast<base::SBasis&>(grid->getBasis());
  sgpp::base::SNakBsplineBoundaryCombigridBase basis(p);
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

  result.setAll(0.0);

// use either schedule static or set #pragma omp atomic before the two result[j/i] =
// tempij*alpha[i/j] lines
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < gridSize; i++) {
    for (size_t j = i; j < gridSize; j++) {
      double temp_ij = 1;

      for (size_t k = 0; k < gridDim; k++) {
        const base::level_t lik = storage[i].getLevel(k);
        const base::level_t ljk = storage[j].getLevel(k);
        const base::index_t iik = storage[i].getIndex(k);
        const base::index_t ijk = storage[j].getIndex(k);
        const sgpp::base::index_t hInvik = 1 << lik;  // = 2^lik
        const sgpp::base::index_t hInvjk = 1 << ljk;
        const double hik = 1.0 / static_cast<double>(hInvik);
        const double hjk = 1.0 / static_cast<double>(hInvjk);
        double offseti_left = (static_cast<double>(iik) - pp1hDbl) * hik;
        double offseti_right = (static_cast<double>(iik) + pp1hDbl) * hik;
        double offsetj_left = (static_cast<double>(ijk) - pp1hDbl) * hjk;
        double offsetj_right = (static_cast<double>(ijk) + pp1hDbl) * hjk;

        //        std::cout << k << "| (" << lik << " " << iik << ") (" << ljk << " " << ijk << ") "
        //                  << std::endl;

        double temp_res = 0.0;
        double offset = 0.0, scaling = 0.0;
        size_t start = 0, stop = 0;
        if (p == 3) {
          if (iik == 3) offseti_left -= hik;
          if (iik == hInvik - 3) offseti_right += hik;
          if (ijk == 3) offsetj_left -= hjk;
          if (ijk == hInvjk - 3) offsetj_right += hjk;
        } else if (p == 5) {
          if ((iik == 3) || (iik == 5)) offseti_left -= 2 * hik;
          if ((iik == hInvik - 3) || (iik == hInvik - 5)) offseti_right += 2 * hik;
          if ((ijk == 3) || (ijk == 5)) offsetj_left -= 2 * hjk;
          if ((ijk == hInvjk - 3) || (ijk == hInvjk - 5)) offsetj_right += 2 * hjk;
        }

        if (std::max(offseti_left, offsetj_left) >= std::min(offseti_right, offsetj_right)) {
          // B spline supports do not not overlap:
          temp_ij = 0.0;

          break;
        } else {
          if (lik >= ljk) {
            offset = offseti_left;
            scaling = hik;
            start = ((iik > pp1h) ? 0 : (pp1h - iik));
            stop = std::min(p, hInvik + pp1h - iik - 1);
            if (p == 3) {
              if ((iik == 3) || (iik == hInvik - 3)) stop += 1;
            } else if (p == 5) {
              if ((iik == 3) || (iik == 5) || (iik == hInvik - 3) || (iik == hInvik - 5)) stop += 2;
            }
            if (lik == 2) {
              start = 1;
              stop = 4;
              offset = -0.25;
              scaling = 0.25;
            }
            if ((p == 5) && (lik == 3)) {
              start = 1;
              stop = 8;
              offset = -0.125;
              scaling = 0.125;
            }
          } else {
            // if lik <= ljk
            offset = offsetj_left;
            scaling = hjk;
            start = ((ijk > pp1h) ? 0 : (pp1h - ijk));
            stop = std::min(p, hInvjk + pp1h - ijk - 1);
            if (p == 3) {
              if ((ijk == 3) || (ijk == hInvjk - 3)) stop += 1;
            } else if (p == 5) {
              if ((ijk == 3) || (ijk == 5) || (ijk == hInvjk - 3) || (ijk == hInvjk - 5)) stop += 2;
            }
            if (ljk == 2) {
              start = 1;
              stop = 4;
              offset = -0.25;
              scaling = 0.25;
            }
            if ((p == 5) && (ljk == 3)) {
              start = 1;
              stop = 8;
              offset = -0.125;
              scaling = 0.125;
            }
          }
        }

        for (size_t n = start; n <= stop; n++) {
          for (size_t c = 0; c < quadOrder; c++) {
            const double x = offset + scaling * (coordinates[c] + static_cast<double>(n));
            temp_res += weights[c] * basis.eval(lik, iik, x) * basis.eval(ljk, ijk, x);
          }
        }
        temp_ij *= scaling * temp_res;
      }

      //#pragma omp atomic
      result[i] += temp_ij * alpha[j];
      if (i != j) {
        //#pragma omp atomic
        result[j] += temp_ij * alpha[i];
      }
    }
  }
}
}  // namespace pde
}  // namespace sgpp
