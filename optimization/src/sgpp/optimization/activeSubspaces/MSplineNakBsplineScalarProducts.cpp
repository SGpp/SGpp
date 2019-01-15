// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/MSplineNakBsplineScalarProducts.hpp>

namespace sgpp {
namespace optimization {

std::unique_ptr<sgpp::base::SBasis> MSplineNakBsplineScalarProducts::initializeBasis(
    sgpp::base::GridType gridType, size_t degree) {
  if (gridType == sgpp::base::GridType::NakBspline) {
    return std::make_unique<sgpp::base::SNakBsplineBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    return std::make_unique<sgpp::base::SNakBsplineBoundaryBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
    return std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
    return std::make_unique<sgpp::base::SNakBsplineExtendedBase>(degree);
  } else {
    throw sgpp::base::generation_exception("ASMatrixNakBspline: gridType not supported.");
  }
}

sgpp::base::DataVector MSplineNakBsplineScalarProducts::getCommonSupport(unsigned int level,
                                                                         unsigned int index) {
  // knots in the common support of the B-spline and the M-spline
  sgpp::base::DataVector supp = nakBSplineSupport(level, index, degree);

  // find first and last  element of xi inside supp
  sgpp::base::DataVector::iterator first;
  first = std::find_if(xi.begin(), xi.end(),
                       bind(std::greater<double>(), std::placeholders::_1, supp[0]));
  if (*first < supp.back()) {
    sgpp::base::DataVector::iterator last;
    last = std::find_if(xi.begin(), xi.end(),
                        bind(std::greater<double>(), std::placeholders::_1, supp.back()));
    --last;
    supp.append(first, last);

    // sort
    std::sort(supp.begin(), supp.end());
    //    sgpp::base::DataVector::iterator it = std::unique(supp.begin(), supp.end());
    //    supp.resize(std::distance(supp.begin(), it));
  }
  return supp;
}

double MSplineNakBsplineScalarProducts::basisScalarProduct(unsigned int level, unsigned int index) {
  sgpp::base::DataVector commonSupport = getCommonSupport(level, index);

  // loop over supp segments
  // in each segment perform quadrature of quadOrder
  double res = 0.0;
  for (size_t i = 0; i < commonSupport.getSize() - 1; i++) {
    // scale coordinates from [0,1] to [commonSupport[i],commonSupport[i+1]]
    sgpp::base::DataVector segmentCoordinates = coordinates;
    segmentCoordinates.mult(commonSupport[i + 1] - commonSupport[i]);
    base::DataVector leftVector(segmentCoordinates.getSize(), commonSupport[i]);
    segmentCoordinates.add(leftVector);

    double segmentIntegral = 0;
    for (size_t j = 0; j < segmentCoordinates.getSize(); j++) {
      double x = segmentCoordinates[j];
      segmentIntegral +=
          weights[j] * basis->eval(level, index, x) * mSplineBasis.eval(xi.getSize() - 1, 0, x);
    }
    res += segmentIntegral * (commonSupport[i + 1] - commonSupport[i]);
  }
  return res;
}

double MSplineNakBsplineScalarProducts::calculateScalarProduct(
    std::shared_ptr<sgpp::base::Grid> grid, sgpp::base::DataVector coeff) {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  double sp = 0.0;
  for (size_t k = 0; k < coeff.getSize(); k++) {
    sp += coeff[k] *
          basisScalarProduct(gridStorage.getPointLevel(k, 0), gridStorage.getPointIndex(k, 0));
  }
  return sp;
}

sgpp::base::DataVector MSplineNakBsplineScalarProducts::nakBSplineSupport(size_t level,
                                                                          size_t index,
                                                                          size_t degree) {
  double indexD = static_cast<double>(index);
  double levelD = static_cast<double>(level);
  double pp1h = floor((static_cast<double>(degree) + 1.0) / 2.0);
  double width = 1.0 / (std::pow(2.0, levelD));
  double lindex = indexD - pp1h;
  double rindex = indexD + pp1h;

  if ((indexD == 1) || (indexD == 3) || levelD <= 2) {
    lindex = 0;
  }
  if ((indexD == std::pow(2.0, levelD) - 3) || (indexD == std::pow(2.0, levelD) - 1) ||
      levelD <= 2) {
    rindex = std::pow(2.0, level);
  }

  if (degree == 5) {  // everything above is for degree 3 and 5
    if ((indexD == 5) || (levelD == 3)) {
      lindex = 0;
    }
    if ((indexD == std::pow(2.0, levelD) - 5) || (levelD == 3)) {
      rindex = std::pow(2.0, levelD);
    }
  }

  sgpp::base::DataVector support(static_cast<int>(rindex - lindex) + 1);
  for (size_t j = 0; j < support.getSize(); j++) {
    support[j] = (lindex + static_cast<double>(j)) * width;
  }
  return support;
}

}  // namespace optimization
}  // namespace sgpp

// #endif /* USE_EIGEN */
