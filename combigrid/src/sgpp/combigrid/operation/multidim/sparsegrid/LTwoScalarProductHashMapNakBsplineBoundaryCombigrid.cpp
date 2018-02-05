// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/sparsegrid/LTwoScalarProductHashMapNakBsplineBoundaryCombigrid.hpp>

#include <algorithm>
#include <map>
#include <vector>

namespace sgpp {
namespace combigrid {

LTwoScalarProductHashMapNakBsplineBoundaryCombigrid::
    LTwoScalarProductHashMapNakBsplineBoundaryCombigrid(sgpp::base::Grid* grid)
    : grid(grid) {
  // initilaize collection
  sgpp::combigrid::SingleFunction constant_weight_function =
      sgpp::combigrid::SingleFunction(sgpp::combigrid::constantFunction<double>(1.0));
  weightFunctionsCollection =
      sgpp::combigrid::WeightFunctionsCollection(grid->getDimension(), constant_weight_function);
  isCustomWeightFunction = false;
  // dummy bounds. in evaluation [0,1] is used
  bounds = sgpp::base::DataVector(std::vector<double>(1, 0));
  numAdditionalPoints = 0;
  degree = dynamic_cast<sgpp::base::NakBsplineBoundaryCombigridGrid*>(grid)->getDegree();
}

LTwoScalarProductHashMapNakBsplineBoundaryCombigrid::
    LTwoScalarProductHashMapNakBsplineBoundaryCombigrid(
        sgpp::base::Grid* grid,
        sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection)
    : grid(grid),
      weightFunctionsCollection(weightFunctionsCollection),
      isCustomWeightFunction(true) {
  // dummy bounds. in evaluation [0,1] is used
  bounds = sgpp::base::DataVector(std::vector<double>(1, 0));
  numAdditionalPoints = 0;
  degree = dynamic_cast<sgpp::base::NakBsplineBoundaryCombigridGrid*>(grid)->getDegree();
}

LTwoScalarProductHashMapNakBsplineBoundaryCombigrid::
    LTwoScalarProductHashMapNakBsplineBoundaryCombigrid(
        sgpp::base::Grid* grid,
        sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection,
        sgpp::base::DataVector bounds)
    : grid(grid),
      weightFunctionsCollection(weightFunctionsCollection),
      isCustomWeightFunction(true),
      bounds(bounds),
      numAdditionalPoints(0) {
  degree = dynamic_cast<sgpp::base::NakBsplineBoundaryCombigridGrid*>(grid)->getDegree();
}

LTwoScalarProductHashMapNakBsplineBoundaryCombigrid::
    LTwoScalarProductHashMapNakBsplineBoundaryCombigrid(
        sgpp::base::Grid* grid,
        sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection,
        sgpp::base::DataVector bounds, size_t numAdditionalPoints)
    : grid(grid),
      weightFunctionsCollection(weightFunctionsCollection),
      isCustomWeightFunction(true),
      bounds(bounds),
      numAdditionalPoints(numAdditionalPoints) {
  degree = dynamic_cast<sgpp::base::NakBsplineBoundaryCombigridGrid*>(grid)->getDegree();
}

LTwoScalarProductHashMapNakBsplineBoundaryCombigrid::
    ~LTwoScalarProductHashMapNakBsplineBoundaryCombigrid() {}

void LTwoScalarProductHashMapNakBsplineBoundaryCombigrid::updateGrid(sgpp::base::Grid* grid) {
  this->grid = grid;
}

void LTwoScalarProductHashMapNakBsplineBoundaryCombigrid::hashLevelIndex(base::level_t li,
                                                                         base::index_t ii,
                                                                         base::level_t lj,
                                                                         base::index_t ij, size_t d,
                                                                         MultiIndex& hashMI) {
  if (li > lj) {
    hashMI[0] = li;
    hashMI[1] = ii;
    hashMI[2] = lj;
    hashMI[3] = ij;
  } else {
    hashMI[0] = lj;
    hashMI[1] = ij;
    hashMI[2] = li;
    hashMI[3] = ii;
  }
  hashMI[4] = d;
}
double LTwoScalarProductHashMapNakBsplineBoundaryCombigrid::calculateScalarProduct(
    base::level_t lid, base::index_t iid, base::level_t ljd, base::index_t ijd, size_t quadOrder,
    size_t d, double offseti_left, double offseti_right, double offsetj_left, double offsetj_right,
    sgpp::base::index_t hInvik, sgpp::base::index_t hInvjk, double hik, double hjk, size_t pp1h) {
  sgpp::base::SNakBsplineBoundaryCombigridBase basis(degree);

  base::DataVector coordinates;
  base::DataVector weights;
  base::GaussLegendreQuadRule1D gauss;
  gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);

  double temp_res = 0.0;
  size_t start = 0, stop = 0;
  double scaling = 0.0, offset = 0.0;

  if (lid >= ljd) {
    offset = offseti_left;
    scaling = hik;
    start = ((iid > pp1h) ? 0 : (pp1h - iid));
    stop = std::min(degree, hInvik + pp1h - iid - 1);
    if (degree == 3) {
      if ((iid == 3) || (iid == hInvik - 3)) stop += 1;
    } else if (degree == 5) {
      if ((iid == 3) || (iid == 5) || (iid == hInvik - 3) || (iid == hInvik - 5)) stop += 2;
    }
    if (lid == 2) {
      start = 1;
      stop = 4;
      offset = -0.25;
      scaling = 0.25;
    }
    if ((degree == 5) && (lid == 3)) {
      start = 1;
      stop = 8;
      offset = -0.125;
      scaling = 0.125;
    }
  } else {
    // if lid <= ljd
    offset = offsetj_left;
    scaling = hjk;
    start = ((ijd > pp1h) ? 0 : (pp1h - ijd));
    stop = std::min(degree, hInvjk + pp1h - ijd - 1);
    if (degree == 3) {
      if ((ijd == 3) || (ijd == hInvjk - 3)) stop += 1;
    } else if (degree == 5) {
      if ((ijd == 3) || (ijd == 5) || (ijd == hInvjk - 3) || (ijd == hInvjk - 5)) stop += 2;
    }
    if (ljd == 2) {
      start = 1;
      stop = 4;
      offset = -0.25;
      scaling = 0.25;
    }
    if ((degree == 5) && (ljd == 3)) {
      start = 1;
      stop = 8;
      offset = -0.125;
      scaling = 0.125;
    }
  }

  for (size_t n = start; n <= stop; n++) {
    for (size_t c = 0; c < quadOrder; c++) {
      // transform  the quadrature points to the segment on which the Bspline is
      // evaluated and from there to the interval[a,b] on which the weight function is defined
      const double x = offset + scaling * (coordinates[c] + static_cast<double>(n));
      double transX = x;
      if (bounds.size() > 1) {
        transX = bounds[2 * d] + (bounds[2 * d + 1] - bounds[2 * d]) * x;
      }
      temp_res += weights[c] * basis.eval(lid, iid, x) * basis.eval(ljd, ijd, x) *
                  weightFunctionsCollection[d](transX);
    }
  }

  return temp_res * scaling;
}

void LTwoScalarProductHashMapNakBsplineBoundaryCombigrid::mult(sgpp::base::DataVector& alpha,
                                                               sgpp::base::DataVector& result) {
  if ((degree != 1) && (degree != 3) && (degree != 5)) {
    std::cerr << "OperationMatrixLTwoDotNakBsplineBoundary: only B spline degrees 1, 3 and 5 are "
                 "supported."
              << std::endl;
  }

  base::GridStorage& storage = grid->getStorage();

  size_t lastNumAdditionalPoints = 0;
  size_t quadOrder = degree + 1 + numAdditionalPoints;

  size_t nrows = storage.getSize();
  size_t ncols = storage.getSize();

  if (alpha.getSize() != ncols || result.getSize() != nrows) {
    throw sgpp::base::data_exception("Dimensions do not match!");
  }

  size_t gridSize = storage.getSize();
  size_t gridDim = storage.getDimension();

  result.setAll(0.0);

  MultiIndex hashMI(5);

  //#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < gridSize; i++) {
    for (size_t j = i; j < gridSize; j++) {
      double temp_ij = 1;

      for (size_t d = 0; d < gridDim; d++) {
        const base::level_t lid = storage[i].getLevel(d);
        const base::level_t ljd = storage[j].getLevel(d);
        const base::index_t iid = storage[i].getIndex(d);
        const base::index_t ijd = storage[j].getIndex(d);

        const size_t pp1h = (degree + 1) >> 1;  //  =|_p/2_|
        const double pp1hDbl = static_cast<double>(pp1h);
        const sgpp::base::index_t hInvik = 1 << lid;  // = 2^lid
        const sgpp::base::index_t hInvjk = 1 << ljd;
        const double hik = 1.0 / static_cast<double>(hInvik);
        const double hjk = 1.0 / static_cast<double>(hInvjk);
        double offseti_left = (static_cast<double>(iid) - pp1hDbl) * hik;
        double offseti_right = (static_cast<double>(iid) + pp1hDbl) * hik;
        double offsetj_left = (static_cast<double>(ijd) - pp1hDbl) * hjk;
        double offsetj_right = (static_cast<double>(ijd) + pp1hDbl) * hjk;

        if (degree == 3) {
          if (iid == 3) offseti_left -= hik;
          if (iid == hInvik - 3) offseti_right += hik;
          if (ijd == 3) offsetj_left -= hjk;
          if (ijd == hInvjk - 3) offsetj_right += hjk;
        } else if (degree == 5) {
          if ((iid == 3) || (iid == 5)) offseti_left -= 2 * hik;
          if ((iid == hInvik - 3) || (iid == hInvik - 5)) offseti_right += 2 * hik;
          if ((ijd == 3) || (ijd == 5)) offsetj_left -= 2 * hjk;
          if ((ijd == hInvjk - 3) || (ijd == hInvjk - 5)) offsetj_right += 2 * hjk;
        }
        if (std::max(offseti_left, offsetj_left) >= std::min(offseti_right, offsetj_right)) {
          // B spline supports do not not overlap:
          temp_ij = 0.0;
          break;
        } else {
          // if (level-Index) is already in hash map use it
          // else calculate 1D integral int bi bj f dx
          double temp_res = 0.0;

          hashLevelIndex(lid, iid, ljd, ijd, d, hashMI);

          auto it = innerProducts.find(hashMI);

          if (it != innerProducts.end()) {
            //            std::cout << hashMI[0] << " " << hashMI[1] << " " << hashMI[2] << " " <<
            //            hashMI[3]
            //                      << " " << hashMI[4] << std::endl;
            temp_res = it->second;
          } else {
            temp_res = calculateScalarProduct(lid, iid, ljd, ijd, quadOrder, d, offseti_left,
                                              offseti_right, offsetj_left, offsetj_right, hInvik,
                                              hInvjk, hik, hjk, pp1h);
            numAdditionalPoints = lastNumAdditionalPoints;
            size_t incrementQuadraturePoints = 3;  // this might be variable
                                                   //            if (isCustomWeightFunction) {
                                                   //              double tol = 1e-14;
                                                   //              double err = 1e14;
                                                   //
                                                   //              while (err > tol) {
            //                lastNumAdditionalPoints = numAdditionalPoints;
            //                numAdditionalPoints += incrementQuadraturePoints;
            //                quadOrder = degree + 1 + numAdditionalPoints;
            //
            //                double finer_temp_res = calculateScalarProduct(
            //                    lid, iid, ljd, ijd, quadOrder, d, offseti_left, offseti_right,
            //                    offsetj_left,
            //                    offsetj_right, hInvik, hInvjk, hik, hjk, pp1h);
            //                err = fabs(temp_res - finer_temp_res);
            //                temp_res = finer_temp_res;
            //
            //                std::cout << numAdditionalPoints << " " << err << std::endl;
            //              }
            //            }

            innerProducts[hashMI] = temp_res;
          }
          temp_ij *= temp_res;
        }
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
}  // namespace combigrid
}  // namespace sgpp
