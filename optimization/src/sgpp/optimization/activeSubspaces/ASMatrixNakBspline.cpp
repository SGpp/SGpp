// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ASMatrixNakBspline.hpp>

namespace sgpp {
namespace optimization {
void ASMatrixNakBspline::buildRegularInterpolant(size_t level) {
  grid->getGenerator().regular(level);
  this->calculateInterpolationCoefficients();
}

void ASMatrixNakBspline::buildAdaptiveInterpolant(size_t maxNumGridPoints, size_t initialLevel,
                                                  size_t refinementsNum) {
  // number of points to be refined in each step
  grid->getGenerator().regular(initialLevel);
  this->calculateInterpolationCoefficients();
  while (grid->getSize() < maxNumGridPoints) {
    this->refineSurplusAdaptive(refinementsNum);
  }
}

void ASMatrixNakBspline::createMatrix(size_t numPoints) { this->createMatrixMonteCarlo(numPoints); }

void ASMatrixNakBspline::createMatrixMonteCarlo(size_t numMCPoints) {
  sgpp::optimization::InterpolantScalarFunctionGradient interpolantGradient(*grid, coefficients);

  RandomNumberGenerator::getInstance().setSeed();
  C.resize(numDim, numDim);
  C.setZero();
  for (size_t i = 0; i < numMCPoints; ++i) {
    sgpp::base::DataVector randomVector(numDim, 1);
    RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    sgpp::base::DataVector gradient(numDim);
    interpolantGradient.eval(randomVector, gradient);
    Eigen::VectorXd e = DataVectorToEigen(gradient);
    this->C += e * e.transpose();
  }
  this->C /= static_cast<double>(numMCPoints);
}

void ASMatrixNakBspline::createMatrixGauss() {
  // prepare gauss quadrature rule
  base::DataVector coordinates, weights;
  base::GaussLegendreQuadRule1D gauss;
  size_t quadOrder = static_cast<size_t>(std::ceil(static_cast<double>(degree) + 1.0 / 2.0));
  gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);
  auto pCoordinates = std::make_shared<sgpp::base::DataVector>(coordinates);
  auto pWeights = std::make_shared<sgpp::base::DataVector>(weights);
  C.resize(numDim, numDim);
  for (unsigned int i = 0; i <= C.cols(); i++) {
    for (unsigned int j = i; j < C.rows(); j++) {
      double entry = matrixEntryGauss(i, j, pCoordinates, pWeights);
      // todo (rehmemk) Can the Eigen library take advantage of C being symmetrical? If so initilize
      // it as such
      C(i, j) = entry;
      C(j, i) = entry;
    }
  }
}

// ----------------- auxiliary routines -----------

void ASMatrixNakBspline::refineSurplusAdaptive(size_t refinementsNum) {
  sgpp::base::SurplusRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  this->calculateInterpolationCoefficients();
}

void ASMatrixNakBspline::calculateInterpolationCoefficients() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  evaluationPoints.resizeZero(gridStorage.getSize(), numDim);
  // ToDo (rehmemk) when refining iteratively use the function values from the last steps!
  functionValues.resizeZero(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector gridPointVector(gridStorage.getDimension());
    gp.getStandardCoordinates(gridPointVector);
    evaluationPoints.setRow(i, gridPointVector);
    functionValues[i] = objectiveFunc->eval(gridPointVector);
  }

  // solve linear system
  sgpp::base::DataVector alpha(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, alpha)) {
    std::cout << "ASMatrixNakBspline: Solving failed.\n";
    return;
  }
  coefficients = alpha;
}

double ASMatrixNakBspline::matrixEntryGauss(size_t i, size_t j,
                                            std::shared_ptr<sgpp::base::DataVector> pCoordinates,
                                            std::shared_ptr<sgpp::base::DataVector> pWeights) {
  double entry = 0.0;
  for (size_t k = 0; k < coefficients.getSize(); k++) {
    for (size_t l = 0; l < coefficients.getSize(); l++) {
      entry += coefficients[k] * coefficients[l] *
               scalarProductDxbiDxbj(i, j, k, l, pCoordinates, pWeights);
    }
  }
  return entry;
}

// Todo (rehmemk) Check already here if (multidimensional) B-spline supports overlap and return 0 if
// so
double ASMatrixNakBspline::scalarProductDxbiDxbj(
    size_t i, size_t j, size_t k, size_t l, std::shared_ptr<sgpp::base::DataVector> pCoordinates,
    std::shared_ptr<sgpp::base::DataVector> pWeights) {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  double integral = 1.0;
  for (size_t d = 0; d < numDim; d++) {
    sgpp::base::GridPoint& gpk = gridStorage.getPoint(k);
    sgpp::base::GridPoint& gpl = gridStorage.getPoint(l);

    size_t indexk = gpk.getIndex(d);
    size_t indexl = gpl.getIndex(d);
    size_t levelk = gpk.getLevel(d);
    size_t levell = gpl.getLevel(d);

    double integral1D = 0.0;
    if ((d == i && i != j)) {
      // int d/dxi b_{k_i} (x_i) b_{l_i} (x_i) dxi
      integral1D = univariateScalarProduct(levelk, indexk, true, levell, indexl, false,
                                           pCoordinates, pWeights);
    } else if (d == j && i != j) {
      // int b_{k_j} (x_j) d/dxj b_{l_j} (x_j) dxj
      integral1D = univariateScalarProduct(levelk, indexk, false, levell, indexl, true,
                                           pCoordinates, pWeights);
    } else if (d == i && i == j) {
      // int d/dxi b_{k_i} (xi) d/dxi b_{l_i} (xi) dxi
      integral1D = univariateScalarProduct(levelk, indexk, true, levell, indexl, true, pCoordinates,
                                           pWeights);
    } else {
      // int b_{k_d} (x_d) b_{l_d} (x_d) dxd
      integral1D = univariateScalarProduct(levelk, indexk, false, levell, indexl, false,
                                           pCoordinates, pWeights);
    }
    if (integral1D == 0) return 0.0;
    integral *= integral1D;
  }
  return integral;
}

double ASMatrixNakBspline::univariateScalarProduct(
    size_t level1, size_t index1, bool dx1, size_t level2, size_t index2, bool dx2,
    std::shared_ptr<sgpp::base::DataVector> pCoordinates,
    std::shared_ptr<sgpp::base::DataVector> pWeights) {
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
    if ((supp1[0] >= supp2.back()) || (supp1.back() <= supp2[0]))
      return 0;
    else {
      if (level1 > level2)
        commonSupport = supp1;
      else if (level1 < level2)
        commonSupport = supp2;
      else
        commonSupport = supp1;
    }

    sgpp::base::SBasis* basis;
    basis = &(grid->getBasis());

    unsigned int l1 = static_cast<unsigned int>(level1);
    unsigned int i1 = static_cast<unsigned int>(index1);
    unsigned int l2 = static_cast<unsigned int>(level2);
    unsigned int i2 = static_cast<unsigned int>(index2);
    std::function<double(double)> func1;
    if (dx1) {
      if (gridType == sgpp::base::GridType::NakBspline) {
        func1 = [&](double x) {
          return static_cast<sgpp::base::SNakBsplineBase*>(basis)->evalDx(l1, i1, x);
        };
      } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
        func1 = [&](double x) {
          return static_cast<sgpp::base::SNakBsplineModifiedBase*>(basis)->evalDx(l1, i1, x);
        };
      } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
        func1 = [&](double x) {
          return static_cast<sgpp::base::SNakBsplineBoundaryBase*>(basis)->evalDx(l1, i1, x);
        };
      } else {
        throw sgpp::base::generation_exception(
            "ASMatrixNakBspline: gridType does not support evalDx.");
      }
    } else {
      func1 = [&](double x) { return basis->eval(l1, i1, x); };
    }
    std::function<double(double)> func2;
    if (dx2) {
      if (gridType == sgpp::base::GridType::NakBspline) {
        func2 = [&](double x) {
          return static_cast<sgpp::base::SNakBsplineBase*>(basis)->evalDx(l2, i2, x);
        };
      } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
        func2 = [&](double x) {
          return static_cast<sgpp::base::SNakBsplineModifiedBase*>(basis)->evalDx(l2, i2, x);
        };
      } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
        func2 = [&](double x) {
          return static_cast<sgpp::base::SNakBsplineBoundaryBase*>(basis)->evalDx(l2, i2, x);
        };
      } else {
        throw sgpp::base::generation_exception(
            "ASMatrixNakBspline: gridType does not support evalDx.");
      }
    } else {
      func2 = [&](double x) { return basis->eval(l2, i2, x); };
    }
    std::function<double(double)> func = [&](double x) { return func1(x) * func2(x); };

    double result = 0.0;
    for (size_t i = 0; i < commonSupport.getSize() - 1; i++) {
      // scale coordinates from [0,1] to [commonSupport[i],commonSupport[i+1]]
      sgpp::base::DataVector segmentCoordinates = *pCoordinates;
      sgpp::base::DataVector segmentWeights = *pWeights;
      segmentCoordinates.mult(commonSupport[i + 1] - commonSupport[i]);
      base::DataVector leftVector(segmentCoordinates.getSize(), commonSupport[i]);
      segmentCoordinates.add(leftVector);
      double segmentIntegral = 0;
      for (size_t i = 0; i < segmentCoordinates.getSize(); i++) {
        segmentIntegral += segmentWeights[i] * func(segmentCoordinates[i]);
      }
      segmentIntegral *= (commonSupport[i + 1] - commonSupport[i]);
      result += segmentIntegral;
    }

    innerProducts.insert(std::pair<asMatrixHashType, double>(hashKey, result));
    it->second = result;
    return result;
  }
}

sgpp::base::DataVector ASMatrixNakBspline::nakBSplineSupport(size_t level, size_t index) {
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

//#endif /* USE_EIGEN */
