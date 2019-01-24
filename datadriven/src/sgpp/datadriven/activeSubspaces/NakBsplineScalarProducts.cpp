// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_EIGEN

#include "../../../../../datadriven/src/sgpp/datadriven/activeSubspaces/NakBsplineScalarProducts.hpp"

namespace sgpp {
namespace datadriven {

std::unique_ptr<sgpp::base::SBasis> NakBsplineScalarProducts::initializeBasis(
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

double NakBsplineScalarProducts::basisScalarProduct(unsigned int level1, unsigned int index1,
                                                    bool dx1, unsigned int level2,
                                                    unsigned int index2, bool dx2) {
  std::map<asMatrixHashType, double>::iterator it;
  asMatrixHashType hashKey = std::make_tuple(level1, index1, dx1, level2, index2, dx2);
  // Check if this scalar products has already been caluclated and is stored
  it = innerProducts.find(hashKey);
  if (it != innerProducts.end()) {
    return it->second;
  } else {
    sgpp::base::DataVector supp1 = nakBSplineSupport(level1, index1, degree1);
    sgpp::base::DataVector supp2 = nakBSplineSupport(level2, index2, degree2);
    sgpp::base::DataVector commonSupport;

    // if supports do not overlap
    if ((supp1[0] >= supp2.back()) || (supp1.back() <= supp2[0])) {
      return 0;
    } else {
      if (level1 >= level2)
        commonSupport = supp1;
      else if (level1 < level2)
        commonSupport = supp2;
      else
        commonSupport = supp1;
    }

    std::function<double(double)> func1;
    if (dx1) {
      func1 = [&](double x) { return basis1->evalDx(level1, index1, x); };
    } else {
      func1 = [&](double x) { return basis1->eval(level1, index1, x); };
    }
    std::function<double(double)> func2;
    if (dx2) {
      func2 = [&](double x) { return basis2->evalDx(level2, index2, x); };
    } else {
      func2 = [&](double x) { return basis2->eval(level2, index2, x); };
    }

    std::function<double(double)> func = [&](double x) { return func1(x) * func2(x); };

    double result = 0.0;
    for (size_t i = 0; i < commonSupport.getSize() - 1; i++) {
      // scale coordinates from [0,1] to [commonSupport[i],commonSupport[i+1]]
      sgpp::base::DataVector segmentCoordinates = coordinates;
      segmentCoordinates.mult(commonSupport[i + 1] - commonSupport[i]);
      base::DataVector leftVector(segmentCoordinates.getSize(), commonSupport[i]);
      segmentCoordinates.add(leftVector);

      double segmentIntegral = 0;
      for (size_t j = 0; j < segmentCoordinates.getSize(); j++) {
        segmentIntegral += weights[j] * func(segmentCoordinates[j]);
      }
      segmentIntegral *= (commonSupport[i + 1] - commonSupport[i]);
      result += segmentIntegral;
    }
    innerProducts.insert(std::pair<asMatrixHashType, double>(hashKey, result));
    return result;
  }
}

double NakBsplineScalarProducts::calculateScalarProduct(std::shared_ptr<sgpp::base::Grid> grid1,
                                                        sgpp::base::DataVector coeff1,
                                                        std::shared_ptr<sgpp::base::Grid> grid2,
                                                        sgpp::base::DataVector coeff2) {
  sgpp::base::GridStorage& gridStorage1 = grid1->getStorage();
  sgpp::base::GridStorage& gridStorage2 = grid2->getStorage();
  double sp = 0.0;
  for (size_t k = 0; k < coeff1.getSize(); k++) {
    for (size_t l = 0; l < coeff2.getSize(); l++) {
      sp += coeff1[k] * coeff2[l] *
            basisScalarProduct(gridStorage1.getPointLevel(k, 0), gridStorage1.getPointIndex(k, 0),
                               false, gridStorage2.getPointLevel(l, 0),
                               gridStorage2.getPointIndex(l, 0), false);
    }
  }
  return sp;
}

sgpp::base::DataVector NakBsplineScalarProducts::nakBSplineSupport(size_t level, size_t index,
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

}  // namespace datadriven
}  // namespace sgpp

// #endif /* USE_EIGEN */
