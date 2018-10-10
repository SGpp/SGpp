// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/activeSubspaces/ASResponseSurfaceNakBspline.hpp>

namespace sgpp {
namespace optimization {

void ASResponseSurfaceNakBspline::createRegularSurfaceFromDetectionPoints(
    sgpp::base::DataMatrix evaluationPoints, sgpp::base::DataVector functionValues, size_t level) {
  size_t numDim = W1.cols();
  std::unique_ptr<sgpp::base::SBasis> basis;
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
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

  // interpolationMatrix(i,j) = b_j (y_i) = b_j (W1T * evaluationPoint_i)
  Eigen::MatrixXd interpolationMatrix(evaluationPoints.getNrows(), gridStorage.getSize());
  for (size_t i = 0; i < evaluationPoints.getNrows(); i++) {
    sgpp::base::DataVector point(evaluationPoints.getNcols());
    evaluationPoints.getRow(i, point);
    //    Eigen::VectorXd y_i = DataVectorToEigen(point);
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
  sgpp::base::DataVector p(numDim);
  interpolantGradient =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(*grid, alpha);
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

}  // namespace optimization
}  // namespace sgpp
