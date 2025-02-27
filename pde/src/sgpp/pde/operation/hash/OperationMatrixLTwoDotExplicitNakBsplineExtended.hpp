// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>

namespace sgpp {
namespace pde {

/**
 * Explicit representation of the matrix \f$(\Phi_i,\Phi_j)_{L2}\f$ for a sparse grid
 */
class OperationMatrixLTwoDotExplicitNakBsplineExtended : public sgpp::base::OperationMatrix {
 public:
  /**
   * Constructor that only builds the object, without initialization
   */
  OperationMatrixLTwoDotExplicitNakBsplineExtended();
  /**
   * Constructor that uses a external matrix pointer to construct the matrix,
   * i.e. matrix is NOT destroyed by the destructor of
   * OperationMatrixLTwoDotExplicitModNakBsplineFullGrid
   *
   * @param m pointer to datamatrix of size (number of grid point) x (number of grid points)
   * @param grid the sparse grid
   */
  OperationMatrixLTwoDotExplicitNakBsplineExtended(sgpp::base::DataMatrix* m,
                                                   sgpp::base::Grid* grid);
  /**
   * Constructor that creates an own matrix
   * i.e. matrix is destroyed by the destructor of OperationMatrixLTwoDotExplicitModBsplineFullGrid
   *
   * @param grid the sparse grid
   */
  explicit OperationMatrixLTwoDotExplicitNakBsplineExtended(sgpp::base::Grid* grid);

  /**
   * Destructor
   */
  virtual ~OperationMatrixLTwoDotExplicitNakBsplineExtended();

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

    const size_t p = dynamic_cast<sgpp::base::NakBsplineExtendedGrid*>(grid)->getDegree();
    const size_t pp1h = (p + 1) / 2;
    const double pp1hDbl = static_cast<double>(pp1h);
    const size_t quadOrder = p + 1;

    base::SNakBsplineExtendedBase& basis = const_cast<base::SNakBsplineExtendedBase&>(
        dynamic_cast<const base::SNakBsplineExtendedBase&>(grid->getBasis()));
    sgpp::base::GridStorage& storage = grid->getStorage();

    sgpp::base::DataVector coordinates;
    sgpp::base::DataVector weights;
    sgpp::base::GaussLegendreQuadRule1D gauss;
    gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);

    // needed for non-quadratic matrix cases
    bool mat_quadratic = (i_start + i_end + j_start + j_end == 0);

    // init standard values
    i_end = i_end == 0 ? gridSize : i_end;

    for (size_t i = i_start; i < i_end; i++) {
      j_start = j_start == 0 ? i : j_start;
      j_end = j_end == 0 ? gridSize : j_end;
#pragma omp parallel for schedule(guided)
      for (size_t j = j_start; j < j_end; j++) {
        double res = 1.;

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
            res = 0.;
            break;
          } else {
            double temp_res = 0.0;

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

            for (size_t n = start; n <= stop; n++) {
              for (size_t c = 0; c < quadOrder; c++) {
                const double x = offset + scaling * (coordinates[c] + static_cast<double>(n));
                temp_res += weights[c] * basis.eval(lik, iik, x) * basis.eval(ljk, ijk, x);
              }
            }

            res *= scaling * temp_res;
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
