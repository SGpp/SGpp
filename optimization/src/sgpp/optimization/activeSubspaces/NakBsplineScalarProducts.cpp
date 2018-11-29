// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ASBsplineScalarProducts.hpp>

namespace sgpp {
namespace optimization {

double ASBsplineScalarProducts::univariateScalarProduct(size_t level1, size_t index1, bool dx1,
                                                        size_t level2, size_t index2, bool dx2) {
  std::map<asMatrixHashType, double>::iterator it;
  asMatrixHashType hashKey = std::make_tuple(level1, index1, dx1, level2, index2, dx2);
  it = innerProducts.find(hashKey);
  if (it != innerProducts.end()) {
    return it->second;
  } else {
    sgpp::base::DataVector supp1 = nakBSplineSupport(level1, index1);
    sgpp::base::DataVector supp2 = nakBSplineSupport(level2, index2);
    sgpp::base::DataVector commonSupport;

    // if supports do not overlap
    if ((supp1[0] >= supp2.back()) || (supp1.back() <= supp2[0])) {
      return 0;
    } else {
      if (level1 > level2)
        commonSupport = supp1;
      else if (level1 < level2)
        commonSupport = supp2;
      else
        commonSupport = supp1;
    }

    unsigned int l1 = static_cast<unsigned int>(level1);
    unsigned int i1 = static_cast<unsigned int>(index1);
    unsigned int l2 = static_cast<unsigned int>(level2);
    unsigned int i2 = static_cast<unsigned int>(index2);
    std::function<double(double)> func1;
    if (dx1) {
      func1 = [&](double x) { return basis->evalDx(l1, i1, x); };
    } else {
      func1 = [&](double x) { return basis->eval(l1, i1, x); };
    }
    std::function<double(double)> func2;
    if (dx2) {
      func2 = [&](double x) { return basis->evalDx(l2, i2, x); };
    } else {
      func2 = [&](double x) { return basis->eval(l2, i2, x); };
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

sgpp::base::DataVector ASBsplineScalarProducts::nakBSplineSupport(size_t level, size_t index) {
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
