// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ASMatrixNakBspline.hpp>

class Rand_double {
 public:
  Rand_double(double low, double high)
      : r(std::bind(std::uniform_real_distribution<>(low, high), std::default_random_engine())) {}

  double operator()() { return r(); }

 private:
  std::function<double()> r;
};

namespace sgpp {
namespace optimization {
void ASMatrixNakBspline::buildRegularInterpolant(size_t level) {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);
  evaluationPoints.resizeZero(gridStorage.getSize(), numDim);
  functionValues.resizeZero(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector gridPointVector(gridStorage.getDimension());
    gp.getStandardCoordinates(gridPointVector);
    evaluationPoints.setRow(i, gridPointVector);
    functionValues[i] = objectiveFunc.eval(gridPointVector);
  }

  // solve linear system
  sgpp::base::DataVector alpha(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, alpha)) {
    std::cout << "ASMatrixNakBspline: Solving failed.\n";
    return;
  }
  this->coefficients = alpha;
  this->interpolantFlag = 1;
}

void ASMatrixNakBspline::createMatrix(size_t numPoints) { this->createMatrixMonteCarlo(numPoints); }

void ASMatrixNakBspline::createMatrixMonteCarlo(size_t numPoints) {
  if (interpolantFlag == 0) {
    std::cout << "ASMatrixNakBspline: cannot create Matrix without interpolant!\n";
    return;
  }
  sgpp::optimization::InterpolantScalarFunctionGradient interpolantGradient(*(this->grid),
                                                                            this->coefficients);

  //  RandomNumberGenerator::getInstance().setSeed();
  C.resize(numDim, numDim);
  C.setZero();
  Rand_double rd{0, 1};
  for (size_t i = 0; i < numPoints; ++i) {
    sgpp::base::DataVector randomVector(numDim, 1);
    // todo (rehmemk) somehow the randomnumbergenerator results are much worse than the Rand_double
    // result
    RandomNumberGenerator::getInstance().getUniformRV(randomVector, 0.0, 1.0);
    //    for (size_t d = 0; d < numDim; ++d) {
    //      randomVector[d] = rd();
    //    }
    //    std::cout << randomVector.toString() << std::endl;
    sgpp::base::DataVector gradient(numDim);
    interpolantGradient.eval(randomVector, gradient);
    Eigen::VectorXd e = DataVectorToEigen(gradient);
    this->C += e * e.transpose();
  }
  this->C /= static_cast<double>(numPoints);
}

double ASMatrixNakBspline::matrixEntryGauss(size_t i, size_t j) {
  double entry = 0.0;
  for (size_t k = 0; k < coefficients.getSize(); k++) {
    for (size_t l = 0; l < coefficients.getSize(); l++) {
      entry += coefficients[k] * coefficients[l] * scalarProductDxbiDxbj(i, j, k, l);
    }
  }
  return entry;
}

double ASMatrixNakBspline::scalarProductDxbiDxbj(size_t i, size_t j, size_t k, size_t l) {
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
      integral1D = univariateScalarProduct(levelk, indexk, true, levell, indexl, false);
    } else if (d == j && i != j) {
      // int b_{k_j} (x_j) d/dxj b_{l_j} (x_j) dxj
      integral1D = univariateScalarProduct(levelk, indexk, false, levell, indexl, true);
    } else if (d == i && i == j) {
      // int d/dxi b_{k_i} (xi) d/dxi b_{l_i} (xi) dxi
      integral1D = univariateScalarProduct(levelk, indexk, true, levell, indexl, true);
    } else {
      // int b_{k_d} (x_d) b_{l_d} (x_d) dxd
      integral1D = univariateScalarProduct(levelk, indexk, false, levell, indexl, false);
    }
    if (integral1D == 0) return 0.0;
    integral *= integral1D;
  }
  return integral;
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

double ASMatrixNakBspline::univariateScalarProduct(size_t level1, size_t index1, bool dx1,
                                                   size_t level2, size_t index2, bool dx2) {
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

  size_t quadOrder = static_cast<size_t>(std::ceil(static_cast<double>(degree) + 1.0 / 2.0));
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
    double segmentIntegral =
        sgpp::optimization::gaussQuad(func, commonSupport[i], commonSupport[i + 1], quadOrder);
    result += segmentIntegral;
  }

  return result;
}

void ASMatrixNakBspline::createMatrixGauss() {
  C.resize(numDim, numDim);
  for (unsigned int i = 0; i <= C.cols(); i++) {
    for (unsigned int j = i; j < C.rows(); j++) {
      double entry = matrixEntryGauss(i, j);
      // todo (rehmemk) Can the Eigen library take advantage of C being symmetrical? If so initilize
      // it as such
      C(i, j) = entry;
      C(j, i) = entry;
    }
  }
}

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
