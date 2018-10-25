// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/activeSubspaces/ASResponseSurfaceNakBspline.hpp>

namespace sgpp {
namespace optimization {

void ASResponseSurfaceNakBspline::initialize() {
  numDim = W1.cols();
  if (gridType == sgpp::base::GridType::NakBspline) {
    grid = std::make_unique<sgpp::base::NakBsplineGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    grid = std::make_unique<sgpp::base::NakBsplineBoundaryGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineBoundaryBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
    grid = std::make_unique<sgpp::base::NakBsplineModifiedGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
  } else {
    throw sgpp::base::generation_exception("ASMatrixNakBspline: gridType not supported.");
  }
}

void ASResponseSurfaceNakBspline::createRegularReducedSurfaceFromDetectionPoints(
    sgpp::base::DataMatrix evaluationPoints, sgpp::base::DataVector functionValues, size_t level) {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

  // interpolationMatrix(i,j) = b_j (y_i) = b_j (W1T * evaluationPoint_i)
  Eigen::MatrixXd interpolationMatrix(evaluationPoints.getNrows(), gridStorage.getSize());
  for (size_t i = 0; i < evaluationPoints.getNrows(); i++) {
    sgpp::base::DataVector point(evaluationPoints.getNcols());
    evaluationPoints.getRow(i, point);
    Eigen::VectorXd y_i = W1.transpose() * DataVectorToEigen(point);
    for (size_t j = 0; j < gridStorage.getSize(); j++) {
      sgpp::base::GridPoint& gpBasis = gridStorage.getPoint(j);
      double basisEval = 1;
      for (size_t t = 0; t < gridStorage.getDimension(); t++) {
        double basisEval1D = basis->eval(gpBasis.getLevel(t), gpBasis.getIndex(t), y_i(t));
        if (basisEval1D == 0) {
          basisEval = 0;
          break;
        } else {
          basisEval *= basisEval1D;
        }
      }
      interpolationMatrix(i, j) = basisEval;
    }
  }

  Eigen::VectorXd functionValues_Eigen = DataVectorToEigen(functionValues);
  // Least Squares Fit
  Eigen::VectorXd alpha_Eigen =
      interpolationMatrix.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
          .solve(functionValues_Eigen);

  sgpp::base::DataVector alpha = EigenToDataVector(alpha_Eigen);
  interpolant = std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, alpha);
  interpolantGradient =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(*grid, alpha);
}

void ASResponseSurfaceNakBspline::createRegularReducedSurfaceWithPseudoInverse(
    size_t level, sgpp::optimization::WrapperScalarFunction objectiveFunc) {
  grid->getGenerator().regular(level);
  sgpp::base::DataVector alpha = calculateInterpolationCoefficientsWithPseudoInverse(objectiveFunc);
  interpolant = std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, alpha);
  interpolantGradient =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(*grid, alpha);

  // Check interpolation property g(p) = f(pinvW1' * p) for all grid points p
  //  Eigen::MatrixXd pinvW1 = W1.transpose().completeOrthogonalDecomposition().pseudoInverse();
  //  double diff = 0;
  //  for (size_t i = 0; i < grid->getSize(); i++) {
  //    auto gp = grid->getStorage().getPoint(i);
  //    Eigen::VectorXd p(grid->getDimension());
  //    for (size_t d = 0; d < grid->getDimension(); d++) {
  //      p[d] = gp.getStandardCoordinate(d);
  //    }
  //    sgpp::base::DataVector pinv = EigenToDataVector(pinvW1 * p);
  //    diff += fabs(objectiveFunc.eval(pinv) - interpolant->eval(EigenToDataVector(p)));
  //  }
  //  std::cout << "total diff: " << diff << std::endl;
}

void ASResponseSurfaceNakBspline::createAdaptiveReducedSurfaceWithPseudoInverse(
    size_t maxNumGridPoints, sgpp::optimization::WrapperScalarFunction objectiveFunc,
    size_t initialLevel) {
  // number of points to be refined in each step
  size_t refinementsNum = 3;
  grid->getGenerator().regular(initialLevel);
  sgpp::base::DataVector alpha = calculateInterpolationCoefficientsWithPseudoInverse(objectiveFunc);
  interpolant = std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, alpha);
  interpolantGradient =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(*grid, alpha);
  while (grid->getSize() < maxNumGridPoints) {
    this->refineSurplusAdaptive(refinementsNum, objectiveFunc, alpha);
    interpolant = std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, alpha);
    interpolantGradient =
        std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(*grid, alpha);
  }
}

double ASResponseSurfaceNakBspline::eval(sgpp::base::DataVector v) {
  Eigen::VectorXd v_Eigen = sgpp::optimization::DataVectorToEigen(v);
  Eigen::VectorXd trans_v_Eigen = W1.transpose() * v_Eigen;
  sgpp::base::DataVector trans_v_DataVector = EigenToDataVector(trans_v_Eigen);
  return interpolant->eval(trans_v_DataVector);
}

double ASResponseSurfaceNakBspline::evalGradient(sgpp::base::DataVector v,
                                                 sgpp::base::DataVector& gradient) {
  Eigen::VectorXd v_Eigen = sgpp::optimization::DataVectorToEigen(v);
  Eigen::VectorXd trans_v_Eigen = W1.transpose() * v_Eigen;
  sgpp::base::DataVector trans_v_DataVector = EigenToDataVector(trans_v_Eigen);
  return interpolantGradient->eval(trans_v_DataVector, gradient);
}

// ----------------- auxiliary routines -----------

void ASResponseSurfaceNakBspline::refineSurplusAdaptive(
    size_t refinementsNum, sgpp::optimization::WrapperScalarFunction objectiveFunc,
    sgpp::base::DataVector& alpha) {
  sgpp::base::SurplusRefinementFunctor functor(alpha, refinementsNum);
  grid->getGenerator().refine(functor);
  alpha = calculateInterpolationCoefficientsWithPseudoInverse(objectiveFunc);
}

sgpp::base::DataVector
ASResponseSurfaceNakBspline::calculateInterpolationCoefficientsWithPseudoInverse(
    sgpp::optimization::WrapperScalarFunction objectiveFunc) {
  Eigen::MatrixXd pinvW1 = W1.transpose().completeOrthogonalDecomposition().pseudoInverse();
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  Eigen::MatrixXd interpolationMatrix(gridStorage.getSize(), gridStorage.getSize());
  Eigen::VectorXd functionValues(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    Eigen::VectorXd p(static_cast<int>(gridStorage.getDimension()));
    for (size_t d = 0; d < gridStorage.getDimension(); d++) {
      p[d] = gp.getStandardCoordinate(d);
    }
    for (size_t j = 0; j < gridStorage.getSize(); j++) {
      sgpp::base::GridPoint& gpBasis = gridStorage.getPoint(j);
      double basisEval = 1;
      for (size_t t = 0; t < gridStorage.getDimension(); t++) {
        double basisEval1D = basis->eval(gpBasis.getLevel(t), gpBasis.getIndex(t), p(t));
        if (basisEval1D == 0) {
          basisEval = 0;
          break;
        } else {
          basisEval *= basisEval1D;
        }
      }
      interpolationMatrix(i, j) = basisEval;
    }

    sgpp::base::DataVector pinv =
        EigenToDataVector(pinvW1 * p);  // introduce a wrapper for eigen functinos so
                                        // we don't have to transform here every time?
    functionValues(i) = objectiveFunc.eval(pinv);
  }
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(interpolationMatrix);
  Eigen::VectorXd alpha_Eigen = dec.solve(functionValues);
  //  Eigen::VectorXf alpha_Eigen =
  //  interpolationMatrix.colPivHouseholderQr().solve(functionValues);
  sgpp::base::DataVector alpha = EigenToDataVector(alpha_Eigen);

  //  std::cout << "f\n"
  //            << functionValues << " M*a\n"
  //            << interpolationMatrix * alpha_Eigen << std::endl;

  return alpha;
}

}  // namespace optimization
}  // namespace sgpp
