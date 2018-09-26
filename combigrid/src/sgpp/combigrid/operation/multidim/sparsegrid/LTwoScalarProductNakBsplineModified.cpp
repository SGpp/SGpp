// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "LTwoScalarProductNakBsplineModified.hpp"

#include <sgpp/combigrid/threading/ThreadPool.hpp>

#include <omp.h>
#include <algorithm>
#include <map>
#include <vector>

namespace sgpp {
namespace combigrid {

LTwoScalarProductNakBsplineModified::LTwoScalarProductNakBsplineModified(
    sgpp::base::Grid* grid)
    : grid(grid) {
  // initilaize collection
  sgpp::combigrid::SingleFunction constant_weight_function =
      sgpp::combigrid::SingleFunction(sgpp::combigrid::constantFunction<double>(1.0));
  weightFunctionsCollection =
      sgpp::combigrid::WeightFunctionsCollection(grid->getDimension(), constant_weight_function);
  isCustomWeightFunction = false;
  bounds = sgpp::base::DataVector(0);
  for (size_t d = 0; d < grid->getDimension(); d++) {
    bounds.push_back(0);
    bounds.push_back(1);
  }
  numAdditionalPoints = 0;
  incrementQuadraturePoints = 1;
  degree = dynamic_cast<sgpp::base::NakBsplineModifiedGrid*>(grid)->getDegree();
}

LTwoScalarProductNakBsplineModified::LTwoScalarProductNakBsplineModified(
    sgpp::base::Grid* grid, sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection)
    : grid(grid),
      weightFunctionsCollection(weightFunctionsCollection),
      isCustomWeightFunction(true) {
  bounds = sgpp::base::DataVector(0);
  for (size_t d = 0; d < grid->getDimension(); d++) {
    bounds.push_back(0);
    bounds.push_back(1);
  }
  numAdditionalPoints = 0;
  incrementQuadraturePoints = 1;
  degree = dynamic_cast<sgpp::base::NakBsplineModifiedGrid*>(grid)->getDegree();
}

LTwoScalarProductNakBsplineModified::LTwoScalarProductNakBsplineModified(
    sgpp::base::Grid* grid, sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection,
    sgpp::base::DataVector bounds)
    : grid(grid),
      weightFunctionsCollection(weightFunctionsCollection),
      isCustomWeightFunction(true),
      bounds(bounds),
      numAdditionalPoints(0),
      incrementQuadraturePoints(1) {
  degree = dynamic_cast<sgpp::base::NakBsplineModifiedGrid*>(grid)->getDegree();
}

LTwoScalarProductNakBsplineModified::LTwoScalarProductNakBsplineModified(
    sgpp::base::Grid* grid, sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection,
    sgpp::base::DataVector bounds, size_t numAdditionalPoints, size_t incrementQuadraturePoints)
    : grid(grid),
      weightFunctionsCollection(weightFunctionsCollection),
      isCustomWeightFunction(true),
      bounds(bounds),
      numAdditionalPoints(numAdditionalPoints),
      incrementQuadraturePoints(incrementQuadraturePoints) {
  degree = dynamic_cast<sgpp::base::NakBsplineModifiedGrid*>(grid)->getDegree();
}

LTwoScalarProductNakBsplineModified::~LTwoScalarProductNakBsplineModified() {}

MultiIndex LTwoScalarProductNakBsplineModified::hashLevelIndex(base::level_t li,
                                                                    base::index_t ii,
                                                                    base::level_t lj,
                                                                    base::index_t ij, size_t d) {
  MultiIndex hashMI(5);

  // use symmetry <(li, ii) , (lj,ij)> = <(lj,ij) , (li,ii)>
  if ((li > lj) || ((li == lj) && (ii > ij)) || ((li == lj) && (ii == ij))) {
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
  return hashMI;
}

double LTwoScalarProductNakBsplineModified::calculateScalarProduct(
    base::level_t lid, base::index_t iid, base::level_t ljd, base::index_t ijd,
    base::DataVector coordinates, base::DataVector weights,
    sgpp::base::SNakBsplineModifiedBase basis, size_t d, double offseti_left,
    double offsetj_left, sgpp::base::index_t hInvik, sgpp::base::index_t hInvjk, double hik,
    double hjk, size_t pp1h) {
  double temp_res = 0.0, scaling = 0.0, offset = 0.0;
  size_t start = 0, stop = 0;

  if (lid >= ljd) {
    offset = offseti_left;
    scaling = hik;
    start = ((iid > pp1h) ? 0 : (pp1h - iid));
    stop = std::min(degree, hInvik + pp1h - iid - 1);
    if (degree == 3) {
      if ((iid == 3) || (iid == hInvik - 3)) stop += 1;
    } else if (degree == 5) {
      //      if ((iid == 3) || (iid == 5) || (iid == hInvik - 3) || (iid == hInvik - 5)) stop += 2;
      if ((iid == 5) || (iid == hInvik - 5)) stop += 2;
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
    // if lid < ljd
    offset = offsetj_left;
    scaling = hjk;
    start = ((ijd > pp1h) ? 0 : (pp1h - ijd));
    stop = std::min(degree, hInvjk + pp1h - ijd - 1);
    if (degree == 3) {
      if ((ijd == 3) || (ijd == hInvjk - 3)) stop += 1;
    } else if (degree == 5) {
      //      if ((ijd == 3) || (ijd == 5) || (ijd == hInvjk - 3) || (ijd == hInvjk - 5)) stop += 2;
      if ((ijd == 5) || (ijd == hInvjk - 5)) stop += 2;
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
    for (size_t c = 0; c < coordinates.size(); c++) {
      // transform  the quadrature points to the segment on which the Bspline is evaluated
      // the weight functions must be defined on [0,1] (otherwise they have to be transformed to
      // [0,1] and the original interval must be handed over in 'bounds'
      const double x = offset + scaling * (coordinates[c] + static_cast<double>(n));

      double transX = x;
      temp_res += weights[c] * basis.eval(lid, iid, transX) * basis.eval(ljd, ijd, transX) *
                  weightFunctionsCollection[d](transX);
    }
  }

  return temp_res * scaling;
}

void LTwoScalarProductNakBsplineModified::mult(sgpp::base::DataVector& alpha,
                                                    sgpp::base::DataVector& result) {
  size_t count = 0;

  if ((degree != 1) && (degree != 3) && (degree != 5)) {
    std::cerr << "OperationMatrixLTwoDotNakBsplineBoundary: only B spline degrees 1, 3 and 5 are "
                 "supported."
              << std::endl;
  }

  sgpp::base::SNakBsplineModifiedBase basis(degree);
  base::GridStorage& storage = grid->getStorage();
  size_t nrows = storage.getSize();
  size_t ncols = storage.getSize();
  size_t gridSize = storage.getSize();
  size_t gridDim = storage.getDimension();
  if (alpha.getSize() != ncols || result.getSize() != nrows) {
    throw sgpp::base::data_exception("Dimensions do not match!");
  }

  //  size_t lastNumAdditionalPoints = 0;
  size_t quadOrder = degree + 1 + numAdditionalPoints;

  base::DataVector coordinates, weights;
  base::GaussLegendreQuadRule1D gauss;

  result.setAll(0.0);

  // ToDo (rehmemk) This parallelization is slower than serial execution. The criticals below are
  // way too big, find out what really causes parallelization problems. While the complete
  // calculateScalarProduct routine is marked critical this is not parallel at all
  // #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < gridSize; i++) {
    //    std::cout << "OMP uses " << omp_get_num_threads() << " thread(s)" << std::endl;
    for (size_t j = i; j < gridSize; j++) {
      double temp_ij = 1;

      for (size_t d = 0; d < gridDim; d++) {
        count++;

        const base::level_t lid = storage[i].getLevel(d);
        const base::level_t ljd = storage[j].getLevel(d);
        const base::index_t iid = storage[i].getIndex(d);
        const base::index_t ijd = storage[j].getIndex(d);

        // setting the segments of the support according to the degree, level and index of the nak
        // B-splines
        const size_t pp1h = (degree + 1) >> 1;  //  =|_(p+1)/2_|
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
          if (iid == 5) offseti_left -= 2 * hik;
          if ((iid == hInvik - 3) || (iid == hInvik - 5)) offseti_right += 2 * hik;
          if (ijd == 5) offsetj_left -= 2 * hjk;
          if ((ijd == hInvjk - 3) || (ijd == hInvjk - 5)) offsetj_right += 2 * hjk;
        }
        if (std::max(offseti_left, offsetj_left) >= std::min(offseti_right, offsetj_right)) {
          // B spline supports do not not overlap:
          temp_ij = 0.0;
          break;
        } else {
          // if (level-Index) is already in the hash map use it
          // else calculate 1D integral int bi bj f dx and save it in the hash map
          double temp_res = 0.0;

          std::map<MultiIndex, double>::iterator it;
          MultiIndex hashMI = hashLevelIndex(lid, iid, ljd, ijd, d);

          // #pragma omp critical
          //          {
          it = innerProducts.find(hashMI);
          //          }
          if (it != innerProducts.end()) {
            temp_res = it->second;
          } else {
            // #pragma omp critical
            //            {

            gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);
            temp_res =
                calculateScalarProduct(lid, iid, ljd, ijd, coordinates, weights, basis, d,
                                       offseti_left, offsetj_left, hInvik, hInvjk, hik, hjk, pp1h);
            double width = bounds[2 * d + 1] - bounds[2 * d];
            temp_res *= width;
            //            }
            // ToDo (rehmemk) now for every B spline the numAdditionalPoints is reseted. Use
            // previous numAdditionalPoints (but not simply numAdditionalPoints =
            // lastNumAdditionalPoints, that didn't work)
            numAdditionalPoints = 10;  // lastNumAdditionalPoints;

            // currently this is always true
            if (isCustomWeightFunction) {
              double tol = 1e-14;
              double err = 1e14;

              while (err > tol) {
                numAdditionalPoints += incrementQuadraturePoints;
                quadOrder = degree + 1 + numAdditionalPoints;
                // This leads to problems with the simple OMP parallelization if the abort
                // criterion is too close to 500
                if (quadOrder > 480) {
                  break;
                }
                double finer_temp_res = 1e+14;
                // #pragma omp critical
                //                {
                gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);
                finer_temp_res = calculateScalarProduct(lid, iid, ljd, ijd, coordinates, weights,
                                                        basis, d, offseti_left, offsetj_left,
                                                        hInvik, hInvjk, hik, hjk, pp1h);
                finer_temp_res *= width;
                //                }
                err = fabs(temp_res - finer_temp_res);
                //                std::cout << temp_res << " " << finer_temp_res << " " << err <<
                //                std::endl;
                temp_res = finer_temp_res;
              }
            }
            // must this be synchronized for OMP?
            innerProducts[hashMI] = temp_res;
          }

          // #pragma omp atomic
          temp_ij *= temp_res;
        }
      }

      // #pragma omp atomic
      result[i] += temp_ij * alpha[j];
      if (i != j) {
        // #pragma omp atomic
        result[j] += temp_ij * alpha[i];
      }
    }
  }

  // std::cout << "map size: " << innerProducts.size() << std::endl;
  //  std::cout << "LTwoScalarProduct numAddP: " << numAdditionalPoints << std::endl;
}
}  // namespace combigrid
}  // namespace sgpp
