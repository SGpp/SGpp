// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <map>
#include <sgpp/base/operation/hash/OperationWeightedSecondMomentNak.hpp>
#include <utility>

namespace sgpp {
namespace base {

std::shared_ptr<sgpp::base::SBasis> OperationWeightedSecondMomentNak::initializeBasis(
    sgpp::base::GridType gridType, size_t degree) {
  if (gridType == sgpp::base::GridType::NakBspline) {
    return std::make_shared<sgpp::base::SNakBsplineBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    return std::make_shared<sgpp::base::SNakBsplineBoundaryBase>(degree);
  } else if (gridType == sgpp::base::GridType::ModNakBspline) {
    return std::make_shared<sgpp::base::SNakBsplineModifiedBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
    return std::make_shared<sgpp::base::SNakBsplineExtendedBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakPBspline) {
    return std::make_shared<sgpp::base::SNakPBsplineBase>(degree);
  } else {
    throw sgpp::base::generation_exception(
        "OperationWeightedSecondMomentNak: gridType not supported.");
  }
}

double OperationWeightedSecondMomentNak::doWeightedQuadrature(
    DataVector& alpha, sgpp::base::DistributionsVector pdfs) {
  //  std::cout << "alpha: " << alpha.toString() << "\n";
  double weightedSP = 0.0;
  size_t dim = storage.getDimension();
  size_t numPoints = alpha.getSize();

  for (size_t k = 0; k < numPoints; k++) {
    for (size_t l = 0; l < numPoints; l++) {
      double sp = 1;
      for (size_t d = 0; d < dim; d++) {
        double temp = weightedBasisScalarProduct(
            storage.getPointLevel(k, d), storage.getPointIndex(k, d), storage.getPointLevel(l, d),
            storage.getPointIndex(l, d), pdfs.get(d));
        sp *= temp;
      }
      weightedSP += alpha[k] * alpha[l] * sp;
    }
  }
  return weightedSP;
}

double OperationWeightedSecondMomentNak::weightedBasisScalarProduct(
    unsigned int level1, unsigned int index1, unsigned int level2, unsigned int index2,
    std::shared_ptr<sgpp::base::Distribution> pdf) {
  sgpp::base::DataVector bounds = pdf->getBounds();
  double leftPDF = bounds[0];
  double rightPDF = bounds[1];
  sgpp::base::DataVector pdfCharacteristics = pdf->getCharacteristics();
  std::map<hashType, double>::iterator it;
  hashType hashKey = std::make_tuple(level1, index1, level2, index2, pdf->getType(),
                                     pdfCharacteristics[0], pdfCharacteristics[1]);
  // Check if this scalar products has already been caluclated and is stored
  it = innerProducts.find(hashKey);
  if (it != innerProducts.end()) {
    return it->second;
  } else {
    sgpp::base::DataVector supp1 = nakBSplineSupport(level1, index1, degree);
    sgpp::base::DataVector supp2 = nakBSplineSupport(level2, index2, degree);
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
    func1 = [&](double x) { return basis->eval(level1, index1, x); };
    std::function<double(double)> func2;
    func2 = [&](double x) { return basis->eval(level2, index2, x); };

    std::function<double(double)> func = [&](double x) { return func1(x) * func2(x); };

    double result = 0.0;
    for (size_t i = 0; i < commonSupport.getSize() - 1; i++) {
      // scale coordinates from [0,1] to [commonSupport[i],commonSupport[i+1]]
      sgpp::base::DataVector segmentCoordinates = coordinates;
      segmentCoordinates.mult(commonSupport[i + 1] - commonSupport[i]);
      base::DataVector leftVector(segmentCoordinates.getSize(), commonSupport[i]);
      segmentCoordinates.add(leftVector);

      // scale again to supp(pdf)
      sgpp::base::DataVector transformedSegmentCoordinates = segmentCoordinates;
      base::DataVector leftVectorPDF(segmentCoordinates.getSize(), leftPDF);
      transformedSegmentCoordinates.mult(rightPDF - leftPDF);
      transformedSegmentCoordinates.add(leftVectorPDF);

      // width of the transformed segment, i.e. width of [commonSupport[i],commonSupport[i+1]]
      // transformed to supp(pdf)
      double transformedSegmentWidth =
          (commonSupport[i + 1] - commonSupport[i]) * (rightPDF - leftPDF);

      double segmentIntegral = 0;
      for (size_t j = 0; j < segmentCoordinates.getSize(); j++) {
        segmentIntegral +=
            weights[j] * func(segmentCoordinates[j]) * pdf->eval(transformedSegmentCoordinates[j]);
      }
      segmentIntegral *= transformedSegmentWidth;
      result += segmentIntegral;
    }
    innerProducts.insert(std::pair<hashType, double>(hashKey, result));
    return result;
  }
}  // namespace base

sgpp::base::DataVector OperationWeightedSecondMomentNak::nakBSplineSupport(size_t level,
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

}  // namespace base
}  // namespace sgpp
